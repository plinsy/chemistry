import re
import networkx as nx
from networkx.algorithms import isomorphism
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
    TimeRemainingColumn,
)

# Creer notre dictionnaire
ELEMENTS_VALENCE = {
    # Élément : {valence, color (CPK standard)}
    "H": {"valence": 1, "color": "#FFFFFF"},  # Hydrogène : Blanc
    "C": {"valence": 4, "color": "#FF0D0D"},  # Carbone : Rouge
}

VISUALISER_EN_3D = False


def calculer_insaturation(carbone, hydrogene, azote=0, halogene=0):
    # Formule standard du Degree of Unsaturation (DoU)
    dou = (2 * carbone + 2 + azote - hydrogene - halogene) / 2
    return int(dou)


def parse_formule(formule):
    # Regex : cherche un Symbole (Majuscule + minuscule optionnelle) suivi de chiffres
    # Exemple : 'C12', 'Cl2', 'H'
    pattern = r"([A-Z][a-z]?)(\d*)"
    elements = re.findall(pattern, formule)

    dictionnaire_atomes = {}

    for symbole, quantite in elements:
        # Si aucun chiffre n'est précisé (ex: 'H'), la quantité est 1
        if quantite == "":
            valeur = 1
        else:
            valeur = int(quantite)

        # On additionne au cas où l'utilisateur écrit 'CH3CH2OH'
        dictionnaire_atomes[symbole] = dictionnaire_atomes.get(symbole, 0) + valeur

    return dictionnaire_atomes


def creer_molecule_exemple():
    # 1. Initialiser le graphe
    mol = nx.Graph()

    # 2. Ajouter les atomes (nœuds) avec leurs attributs
    # On ajoute 3 carbones
    mol.add_node(0, element="C", valence_max=4)
    mol.add_node(1, element="C", valence_max=4)
    mol.add_node(2, element="C", valence_max=4)

    # 3. Ajouter les liaisons (arêtes)
    # On lie C0-C1 et C1-C2 (chaîne linéaire)
    mol.add_edge(0, 1, type_liaison=1)  # Liaison simple
    mol.add_edge(1, 2, type_liaison=1)

    return mol


def verifier_valences(graphe):
    for node_id in graphe.nodes:
        element = graphe.nodes[node_id]["element"]
        valence_max = graphe.nodes[node_id]["valence_max"]

        # Le degré réel (nombre de liaisons connectées)
        # Note : pour les liaisons doubles/triples, il faudra sommer les 'type_liaison'
        liaisons_actuelles = sum(
            attr["type_liaison"] for _, _, attr in graphe.edges(node_id, data=True)
        )

        if liaisons_actuelles > valence_max:
            print(
                f"Erreur : L'atome {node_id} ({element}) a {liaisons_actuelles} liaisons (Max: {valence_max})"
            )
            return False

    print("Structure valide selon les règles de valence.")
    return True


def obtenir_smiles_canonique(graphe_networkx):
    # 1. Créer une molécule vide RDKit (molécule modifiable)
    rw_mol = Chem.RWMol()

    # 2. Ajouter les atomes
    node_to_idx = {}
    for node in graphe_networkx.nodes:
        idx = rw_mol.AddAtom(Chem.Atom("C"))  # On commence par les squelettes de C
        node_to_idx[node] = idx

    # 3. Ajouter les liaisons
    for u, v, data in graphe_networkx.edges(data=True):
        rw_mol.AddBond(node_to_idx[u], node_to_idx[v], Chem.BondType.SINGLE)

    # 4. Convertir en molécule finale et "Saturer" d'hydrogènes
    mol = rw_mol.GetMol()
    # Sanitize the molecule to calculate implicit valences
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)  # RDKit ajoute les H automatiquement selon la valence !

    # 5. Générer le SMILES canonique
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def generer_tous_squelettes(n_carbones, n_hydrogenes=None):
    """
    Génère tous les isomères possibles pour une formule CnHm.
    Si n_hydrogenes est None, génère uniquement les alcanes saturés.
    Sinon, génère toutes les structures (alkènes, alkines, cycles) correspondant à CnHm.
    """
    # On utilise un set pour stocker les chaînes SMILES (élimination auto des doublons)
    smiles_uniques = set()
    progress_tracker = {"count": 0, "task": None}

    def convertir_en_smiles_avec_liaisons(graphe):
        """Transforme le graphe NetworkX (avec types de liaisons) en SMILES canonique"""
        try:
            rw_mol = Chem.RWMol()
            node_to_idx = {}

            # Ajouter les atomes de carbone
            for node in graphe.nodes:
                idx = rw_mol.AddAtom(Chem.Atom("C"))
                node_to_idx[node] = idx

            # Ajouter les liaisons avec leur type
            for u, v, data in graphe.edges(data=True):
                bond_type = data.get("type_liaison", 1)
                if bond_type == 1:
                    rw_mol.AddBond(node_to_idx[u], node_to_idx[v], Chem.BondType.SINGLE)
                elif bond_type == 2:
                    rw_mol.AddBond(node_to_idx[u], node_to_idx[v], Chem.BondType.DOUBLE)
                elif bond_type == 3:
                    rw_mol.AddBond(node_to_idx[u], node_to_idx[v], Chem.BondType.TRIPLE)

            mol = rw_mol.GetMol()
            Chem.SanitizeMol(mol)

            # Vérifier le nombre d'hydrogènes si spécifié
            if n_hydrogenes is not None:
                mol_with_h = Chem.AddHs(mol)
                h_count = sum(
                    1 for atom in mol_with_h.GetAtoms() if atom.GetSymbol() == "H"
                )
                if h_count != n_hydrogenes:
                    return None

            # Generate SMILES without explicit hydrogens (they will be implicit)
            return Chem.MolToSmiles(mol)
        except:
            return None

    def compter_valences_utilisees(graphe):
        """Compte les valences utilisées pour chaque nœud"""
        valences = {}
        for node in graphe.nodes:
            valences[node] = sum(
                data["type_liaison"] for _, _, data in graphe.edges(node, data=True)
            )
        return valences

    def generer_structures_acycliques():
        """Génère les structures sans cycles (arbres)"""

        def explorer_arbre(graphe):
            # Si on a n-1 liaisons et que tout est connecté, on a un arbre valide
            if graphe.number_of_edges() == n_carbones - 1:
                if nx.is_connected(graphe):
                    # Essayer différentes combinaisons de types de liaisons
                    explorer_types_liaisons(graphe, list(graphe.edges()))
                return

            # Construire l'arbre
            for i in range(n_carbones):
                for j in range(i + 1, n_carbones):
                    if not graphe.has_edge(i, j):
                        valences = compter_valences_utilisees(graphe)
                        if valences.get(i, 0) < 4 and valences.get(j, 0) < 4:
                            graphe.add_edge(i, j, type_liaison=1)
                            explorer_arbre(graphe)
                            graphe.remove_edge(i, j)

        g = nx.Graph()
        for i in range(n_carbones):
            g.add_node(i, element="C")
        explorer_arbre(g)

    def generer_structures_cycliques():
        """Génère les structures avec cycles"""
        # Calculer le degré d'insaturation pour déterminer la plage d'exploration
        if n_hydrogenes is not None:
            dou = calculer_insaturation(n_carbones, n_hydrogenes)
            # Pour les molécules fortement insaturées (DoU >= 3), on a besoin d'explorer
            # beaucoup plus de configurations incluant des systèmes polycycliques
            # et des combinaisons complexes de liaisons multiples
            if dou >= 3:
                # Pour les molécules très insaturées, explorer jusqu'au maximum théorique
                # Le max pour n nœuds avec degré max 4 est (n * 4) / 2 = 2n
                max_edges = min(2 * n_carbones, n_carbones + dou + 5)
            else:
                max_edges = min(n_carbones + dou + 3, int(n_carbones * 1.5) + 2)
        else:
            max_edges = n_carbones + 3

        # Explorer les graphes avec différents nombres d'arêtes
        for n_edges in range(n_carbones, max_edges + 1):
            explorer_graphes_cycliques(n_edges)

    def explorer_graphes_cycliques(n_edges):
        """Explore les graphes avec un nombre donné d'arêtes"""
        # Utiliser une approche exhaustive pour générer tous les graphes possibles
        # avec n_edges arêtes, en explorant toutes les combinaisons possibles

        def construire_graphe(graphe, node_pairs, pair_idx):
            """
            Construit des graphes en essayant toutes les combinaisons d'arêtes
            node_pairs: liste de toutes les paires possibles (i, j) avec i < j
            pair_idx: index actuel dans node_pairs
            """
            n_current_edges = graphe.number_of_edges()

            # Si on a le bon nombre d'arêtes, tester ce graphe
            if n_current_edges == n_edges:
                if nx.is_connected(graphe):
                    # Essayer différentes combinaisons de types de liaisons
                    explorer_types_liaisons(graphe, list(graphe.edges()))
                return

            # Si on a trop d'arêtes ou qu'on a parcouru toutes les paires, abandonner
            if n_current_edges > n_edges or pair_idx >= len(node_pairs):
                return

            # Nombre d'arêtes restantes à ajouter
            edges_needed = n_edges - n_current_edges
            # Paires restantes à considérer
            pairs_left = len(node_pairs) - pair_idx

            # Si on n'a plus assez de paires pour atteindre n_edges, abandonner
            if pairs_left < edges_needed:
                return

            i, j = node_pairs[pair_idx]

            # Option 1: Ajouter cette arête si les valences le permettent
            valences = compter_valences_utilisees(graphe)
            if valences.get(i, 0) < 4 and valences.get(j, 0) < 4:
                graphe.add_edge(i, j, type_liaison=1)
                construire_graphe(graphe, node_pairs, pair_idx + 1)
                graphe.remove_edge(i, j)

            # Option 2: Ne pas ajouter cette arête
            construire_graphe(graphe, node_pairs, pair_idx + 1)

        # Générer toutes les paires possibles de nœuds
        node_pairs = [
            (i, j) for i in range(n_carbones) for j in range(i + 1, n_carbones)
        ]

        g = nx.Graph()
        for i in range(n_carbones):
            g.add_node(i, element="C")

        construire_graphe(g, node_pairs, 0)

    def explorer_types_liaisons(graphe, edges, edge_idx=0):
        """Essaye différentes combinaisons de types de liaisons (simple, double, triple)"""
        if edge_idx == len(edges):
            # Toutes les liaisons ont été assignées, vérifier la validité
            valences = compter_valences_utilisees(graphe)
            if all(v <= 4 for v in valences.values()):
                smiles = convertir_en_smiles_avec_liaisons(graphe)
                if smiles:
                    smiles_uniques.add(smiles)
                    progress_tracker["count"] += 1
                    if progress_tracker["task"] is not None:
                        progress_tracker["progress"].update(
                            progress_tracker["task"],
                            completed=progress_tracker["count"],
                        )
            return

        u, v = edges[edge_idx]

        # Essayer liaison simple, double, triple
        for bond_type in [1, 2, 3]:
            valences = compter_valences_utilisees(graphe)
            current_u = valences.get(u, 0)
            current_v = valences.get(v, 0)

            # Retirer la valence actuelle de cette liaison
            current_bond = graphe[u][v]["type_liaison"]
            current_u -= current_bond
            current_v -= current_bond

            # Vérifier si on peut ajouter ce type de liaison
            if current_u + bond_type <= 4 and current_v + bond_type <= 4:
                graphe[u][v]["type_liaison"] = bond_type
                explorer_types_liaisons(graphe, edges, edge_idx + 1)

        # Restaurer la liaison simple par défaut
        graphe[u][v]["type_liaison"] = 1

    # Utiliser rich.progress pour les molécules complexes
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeRemainingColumn(),
    ) as progress:
        # Estimer grossièrement le nombre de structures possibles
        estimated_total = 10 ** min(n_carbones, 4)  # Estimation très approximative
        task = progress.add_task(
            "[cyan]Génération des isomères...", total=estimated_total
        )
        progress_tracker["progress"] = progress
        progress_tracker["task"] = task

        # Générer les structures acycliques
        generer_structures_acycliques()

        # Générer les structures cycliques si on cherche des structures insaturées
        if n_hydrogenes is not None and n_carbones >= 3:
            dou = calculer_insaturation(n_carbones, n_hydrogenes)
            if dou > 0:  # S'il y a de l'insaturation, explorer les cycles
                generer_structures_cycliques()

        # Finaliser la barre de progression
        progress.update(task, completed=estimated_total)

    return list(smiles_uniques)


def smiles_vers_graphe(smiles):
    """Convertit un SMILES en graphe NetworkX pour visualisation"""
    mol = Chem.MolFromSmiles(smiles)
    # Ajouter explicitement les hydrogènes pour la visualisation
    mol = Chem.AddHs(mol)
    g = nx.Graph()

    # Ajouter les atomes
    for atom in mol.GetAtoms():
        g.add_node(
            atom.GetIdx(),
            element=atom.GetSymbol(),
            valence_max=ELEMENTS_VALENCE[atom.GetSymbol()]["valence"],
        )

    # Ajouter les liaisons
    for bond in mol.GetBonds():
        g.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            type_liaison=int(bond.GetBondTypeAsDouble()),
        )

    return g


def visualiser_molecules(liste_smiles, formule_brute, max_isomers=20):
    """Affiche toutes les structures possibles dans une fenêtre scrollable avec 2 colonnes (interactive)"""
    n_molecules = len(liste_smiles)

    if n_molecules == 0:
        print("Aucune structure trouvée")
        return

    # Limiter le nombre d'isomères affichés
    if n_molecules > max_isomers:
        print(
            f"{n_molecules} isomères trouvés, affichage limité aux {max_isomers} premiers"
        )
        liste_smiles = liste_smiles[:max_isomers]
        n_molecules = max_isomers

    # Calculer le nombre de lignes et colonnes pour la grille (2 colonnes)
    n_cols = 2
    n_rows = (n_molecules + n_cols - 1) // n_cols  # Arrondi supérieur

    # Créer une fenêtre tkinter avec scrollbar
    root = tk.Tk()
    root.title(f"{formule_brute} - {n_molecules} Isomères (Interactive)")

    # Créer un frame principal avec scrollbar
    main_frame = ttk.Frame(root)
    main_frame.pack(fill=tk.BOTH, expand=1)

    # Créer un canvas
    canvas = tk.Canvas(main_frame, width=900, height=700)
    canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

    # Ajouter scrollbar verticale
    scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=canvas.yview)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    # Configurer le canvas
    canvas.configure(yscrollcommand=scrollbar.set)
    canvas.bind(
        "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
    )

    # Créer un frame dans le canvas
    second_frame = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=second_frame, anchor="nw")

    # Créer la figure matplotlib
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 5 * n_rows))
    fig.suptitle(
        f"{formule_brute} - {n_molecules} Isomères (Drag & Drop)",
        fontsize=16,
        fontweight="bold",
    )

    # S'assurer que axes est toujours un tableau 2D
    if n_molecules == 1:
        axes = np.array([[axes, axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)

    # Stocker les données pour l'interactivité
    graph_data = []

    # Utiliser rich.progress pour montrer l'avancement
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeRemainingColumn(),
    ) as progress:
        task = progress.add_task(
            "[green]Création des images interactives...", total=n_molecules
        )

        for idx, smiles in enumerate(liste_smiles):
            row = idx // n_cols
            col = idx % n_cols
            ax = axes[row, col]

            graphe = smiles_vers_graphe(smiles)

            # Extraire les couleurs basées sur les éléments
            couleurs = [
                ELEMENTS_VALENCE[graphe.nodes[node]["element"]]["color"]
                for node in graphe.nodes
            ]

            # Créer les labels avec symboles atomiques
            labels = {node: graphe.nodes[node]["element"] for node in graphe.nodes}

            # Dessiner le graphe dans le subplot
            pos = nx.spring_layout(graphe, seed=42)

            # Stocker les données pour l'interactivité
            graph_data.append(
                {
                    "ax": ax,
                    "graphe": graphe,
                    "pos": pos,
                    "couleurs": couleurs,
                    "labels": labels,
                    "selected_node": None,
                    "smiles": smiles,
                    "idx": idx,
                }
            )

            nx.draw(
                graphe,
                pos,
                labels=labels,
                node_color=couleurs,
                node_size=800,
                font_color="white",
                font_weight="bold",
                with_labels=True,
                edgecolors="black",
                linewidths=2,
                ax=ax,
            )
            ax.set_title(f"Isomère {idx+1}\n{smiles}", fontsize=10, fontweight="bold")
            ax.axis("off")

            progress.update(task, advance=1)

    # Cacher les subplots vides
    for idx in range(n_molecules, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis("off")

    plt.tight_layout()

    # Variables pour le drag and drop
    drag_state = {"dragging": False, "graph_idx": None, "node": None}

    def find_nearest_node(event, gdata):
        """Trouve le nœud le plus proche du clic"""
        if event.inaxes != gdata["ax"]:
            return None

        # Convertir les coordonnées du clic en coordonnées du graphe
        min_dist = float("inf")
        nearest_node = None

        for node, (x, y) in gdata["pos"].items():
            # Distance entre le clic et le nœud
            dist = np.sqrt((x - event.xdata) ** 2 + (y - event.ydata) ** 2)
            if dist < min_dist and dist < 0.1:  # Seuil de distance
                min_dist = dist
                nearest_node = node

        return nearest_node

    def on_press(event):
        """Gère le clic de souris"""
        if event.inaxes is None:
            return

        # Trouver quel graphe est cliqué
        for idx, gdata in enumerate(graph_data):
            node = find_nearest_node(event, gdata)
            if node is not None:
                drag_state["dragging"] = True
                drag_state["graph_idx"] = idx
                drag_state["node"] = node
                break

    def on_motion(event):
        """Gère le mouvement de la souris pendant le drag"""
        if not drag_state["dragging"] or event.inaxes is None:
            return

        gdata = graph_data[drag_state["graph_idx"]]

        if event.inaxes != gdata["ax"]:
            return

        node = drag_state["node"]

        # Mettre à jour la position du nœud
        old_pos = gdata["pos"][node].copy()
        gdata["pos"][node] = np.array([event.xdata, event.ydata])

        # Déplacer les nœuds connectés (hydrogènes) de manière fluide
        for neighbor in gdata["graphe"].neighbors(node):
            # Si c'est un hydrogène, le déplacer avec le carbone
            if gdata["graphe"].nodes[neighbor]["element"] == "H":
                # Calculer le déplacement relatif
                delta = gdata["pos"][node] - old_pos
                gdata["pos"][neighbor] = gdata["pos"][neighbor] + delta

        # Redessiner le graphe
        gdata["ax"].clear()
        nx.draw(
            gdata["graphe"],
            gdata["pos"],
            labels=gdata["labels"],
            node_color=gdata["couleurs"],
            node_size=800,
            font_color="white",
            font_weight="bold",
            with_labels=True,
            edgecolors="black",
            linewidths=2,
            ax=gdata["ax"],
        )
        gdata["ax"].set_title(
            f"Isomère {gdata['idx']+1}\n{gdata['smiles']}",
            fontsize=10,
            fontweight="bold",
        )
        gdata["ax"].axis("off")

        canvas_widget.draw_idle()

    def on_release(event):
        """Gère le relâchement de la souris"""
        drag_state["dragging"] = False
        drag_state["graph_idx"] = None
        drag_state["node"] = None

    # Connecter les événements
    fig.canvas.mpl_connect("button_press_event", on_press)
    fig.canvas.mpl_connect("motion_notify_event", on_motion)
    fig.canvas.mpl_connect("button_release_event", on_release)

    # Intégrer matplotlib dans tkinter
    canvas_widget = FigureCanvasTkAgg(fig, master=second_frame)
    canvas_widget.draw()
    canvas_widget.get_tk_widget().pack()

    # Mettre à jour le canvas pour le scroll
    second_frame.update_idletasks()
    canvas.configure(scrollregion=canvas.bbox("all"))

    # Activer le scroll avec la molette de souris
    def _on_mousewheel(event):
        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

    canvas.bind_all("<MouseWheel>", _on_mousewheel)

    print(f"\nVisualisation créée! Glissez-déposez les atomes pour les repositionner.")
    root.mainloop()


def visualiser_molecules_3d(
    liste_smiles, formule_brute, max_isomers=20, optimize=False
):
    """Affiche toutes les structures possibles en 3D dans une fenêtre scrollable avec 2 colonnes"""
    n_molecules = len(liste_smiles)

    if n_molecules == 0:
        print("Aucune structure trouvée")
        return

    # Limiter le nombre d'isomères affichés
    if n_molecules > max_isomers:
        print(
            f"{n_molecules} isomères trouvés, affichage limité aux {max_isomers} premiers"
        )
        liste_smiles = liste_smiles[:max_isomers]
        n_molecules = max_isomers

    # Calculer le nombre de lignes et colonnes pour la grille (2 colonnes)
    n_cols = 2
    n_rows = (n_molecules + n_cols - 1) // n_cols  # Arrondi supérieur

    # Créer une fenêtre tkinter avec scrollbar
    root = tk.Tk()
    root.title(f"{formule_brute} - {n_molecules} Isomères (3D)")

    # Créer un frame principal avec scrollbar
    main_frame = ttk.Frame(root)
    main_frame.pack(fill=tk.BOTH, expand=1)

    # Créer un canvas
    canvas = tk.Canvas(main_frame, width=1000, height=700)
    canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

    # Ajouter scrollbar verticale
    scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=canvas.yview)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    # Configurer le canvas
    canvas.configure(yscrollcommand=scrollbar.set)
    canvas.bind(
        "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
    )

    # Créer un frame dans le canvas
    second_frame = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=second_frame, anchor="nw")

    # Créer la figure matplotlib 3D
    fig = plt.figure(figsize=(14, 6 * n_rows))
    fig.suptitle(
        f"{formule_brute} - {n_molecules} Isomères (3D)", fontsize=16, fontweight="bold"
    )

    # Utiliser rich.progress pour montrer l'avancement
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeRemainingColumn(),
    ) as progress:
        task = progress.add_task(
            "[magenta]Création des structures 3D...", total=n_molecules
        )

        for idx, smiles in enumerate(liste_smiles):
            # Créer un subplot 3D
            ax = fig.add_subplot(n_rows, n_cols, idx + 1, projection="3d")

            # Créer la molécule et ajouter les hydrogènes
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)

            # Générer les coordonnées 3D
            AllChem.EmbedMolecule(mol, randomSeed=42)  # type: ignore

            # Optimisation optionnelle (plus lent mais plus précis)
            if optimize:
                AllChem.MMFFOptimizeMolecule(mol)  # type: ignore

            # Obtenir les coordonnées 3D
            conf = mol.GetConformer()

            # Extraire les positions et couleurs des atomes
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                couleur = ELEMENTS_VALENCE[atom.GetSymbol()]["color"]

                ax.scatter(
                    pos.x,
                    pos.y,
                    pos.z,
                    c=couleur,
                    s=500,
                    edgecolors="black",
                    linewidths=2,
                    alpha=0.9,
                )

                # Ajouter le label de l'atome
                ax.text(
                    pos.x,
                    pos.y,
                    pos.z,
                    atom.GetSymbol(),
                    fontsize=10,
                    fontweight="bold",
                    ha="center",
                    va="center",
                )

            # Dessiner les liaisons
            for bond in mol.GetBonds():
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                begin_pos = conf.GetAtomPosition(begin_idx)
                end_pos = conf.GetAtomPosition(end_idx)

                ax.plot(
                    [begin_pos.x, end_pos.x],
                    [begin_pos.y, end_pos.y],
                    [begin_pos.z, end_pos.z],
                    "k-",
                    linewidth=2,
                )

            ax.set_xlabel("X", fontsize=8)
            ax.set_ylabel("Y", fontsize=8)
            ax.set_zlabel("Z", fontsize=8)
            ax.set_title(f"Isomère {idx+1}\n{smiles}", fontsize=10, fontweight="bold")

            # Désactiver la grille pour un meilleur rendu
            ax.grid(True, alpha=0.3)

            progress.update(task, advance=1)

    plt.tight_layout()

    # Intégrer matplotlib dans tkinter
    canvas_widget = FigureCanvasTkAgg(fig, master=second_frame)
    canvas_widget.draw()
    canvas_widget.get_tk_widget().pack()

    # Mettre à jour le canvas pour le scroll
    second_frame.update_idletasks()
    canvas.configure(scrollregion=canvas.bbox("all"))

    # Activer le scroll avec la molette de souris
    def _on_mousewheel(event):
        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

    canvas.bind_all("<MouseWheel>", _on_mousewheel)

    print(f"\nVisualisation 3D créée!")
    root.mainloop()


def analyser_entree_utilisateur():
    formule_brute = input("Entrez une formule chimique (ex: C4H10) : ").strip().upper()

    try:
        atomes = parse_formule(formule_brute)

        # Extraction des valeurs avec des valeurs par défaut à 0
        c = atomes.get("C", 0)
        h = atomes.get("H", 0)
        n = atomes.get("N", 0)
        # On regroupe les halogènes courants pour le calcul
        x = (
            atomes.get("Cl", 0)
            + atomes.get("Br", 0)
            + atomes.get("F", 0)
            + atomes.get("I", 0)
        )

        # Sécurité : Pas de carbone, pas de molécule organique !
        if c == 0:
            print("Ceci n'est pas une molécule organique standard.")
            return

        dou = calculer_insaturation(c, h, n, x)

        print(f"\n--- Résultats pour {formule_brute} ---")
        print(f"Composition : {atomes}")
        print(f"Degré d'insaturation (DoU) : {dou}")

        # Générer tous les squelettes carbonés possibles
        print(f"\nGénération des isomères pour {formule_brute}...")
        resultats = generer_tous_squelettes(c, h)
        print(f"\nNombre d'isomères trouvés : {len(resultats)}")

        # Afficher les SMILES
        for idx, smiles in enumerate(resultats, 1):
            print(f"  {idx}. {smiles}")

        # Visualiser toutes les structures en 3D (rapide par défaut)
        # Pour une optimisation géométrique complète (plus lent), utilisez: optimize=True
        if VISUALISER_EN_3D:
            visualiser_molecules_3d(
                resultats, formule_brute, max_isomers=20, optimize=False
            )
        else:
            visualiser_molecules(
                resultats,
                formule_brute,
            )

    except Exception as e:
        print(f"Erreur lors de l'analyse : {e}")


def main():
    analyser_entree_utilisateur()


if __name__ == "__main__":
    main()
