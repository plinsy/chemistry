import re
import networkx as nx
from networkx.algorithms import isomorphism
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Creer notre dictionnaire
ELEMENTS_VALENCE = {
    # Élément : {valence, color (CPK standard)}
    "H": {"valence": 1, "color": "#FFFFFF"},  # Hydrogène : Blanc
    "C": {"valence": 4, "color": "#FF0D0D"},  # Carbone : Rouge
}


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
                f"❌ Erreur : L'atome {node_id} ({element}) a {liaisons_actuelles} liaisons (Max: {valence_max})"
            )
            return False

    print("✅ Structure valide selon les règles de valence.")
    return True


# def sont_identiques(mol1, mol2):
#     # On définit ce qui rend deux nœuds identiques (ici, l'élément chimique)
#     def node_match(n1, n2):
#         return n1["element"] == n2["element"]

#     # On définit ce qui rend deux liaisons identiques (le type : simple, double...)
#     def edge_match(e1, e2):
#         return e1["type_liaison"] == e2["type_liaison"]

#     return nx.is_isomorphic(mol1, mol2, node_match=node_match, edge_match=edge_match)


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


def generer_tous_squelettes(n_carbones):
    # On utilise un set pour stocker les chaînes SMILES (élimination auto des doublons)
    smiles_uniques = set()

    def convertir_en_smiles(graphe):
        """Transforme le graphe NetworkX en SMILES canonique via RDKit"""
        rw_mol = Chem.RWMol()
        node_to_idx = {}
        for node in graphe.nodes:
            idx = rw_mol.AddAtom(Chem.Atom("C"))
            node_to_idx[node] = idx
        for u, v in graphe.edges():
            rw_mol.AddBond(node_to_idx[u], node_to_idx[v], Chem.BondType.SINGLE)

        mol = rw_mol.GetMol()
        # Sanitize the molecule to calculate implicit valences
        Chem.SanitizeMol(mol)
        # Chem.AddHs remplit les valences vides avec des H automatiquement
        mol = Chem.AddHs(mol)
        return Chem.MolToSmiles(mol)

    def explorer(graphe):
        # CONDITION DE SORTIE :
        # Si on a n-1 liaisons et que tout est connecté, on a un squelette valide
        if graphe.number_of_edges() == n_carbones - 1:
            if nx.is_connected(graphe):
                # --- C'EST ICI QU'ON AJOUTE LE SMILES ---
                smiles = convertir_en_smiles(graphe)
                smiles_uniques.add(smiles)
            return

        # Logique de construction
        for i in range(n_carbones):
            for j in range(i + 1, n_carbones):
                if (
                    not graphe.has_edge(i, j)
                    and graphe.degree(i) < 4
                    and graphe.degree(j) < 4
                ):
                    graphe.add_edge(i, j)
                    explorer(graphe)
                    graphe.remove_edge(i, j)  # Backtrack

    # Initialisation du graphe
    g = nx.Graph()
    for i in range(n_carbones):
        g.add_node(i, element="C")

    explorer(g)
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
    """Affiche toutes les structures possibles pour une formule (chaque isomère dans une figure séparée)"""
    n_molecules = len(liste_smiles)

    if n_molecules == 0:
        print("Aucune structure trouvée")
        return

    # Limiter le nombre d'isomères affichés
    if n_molecules > max_isomers:
        print(
            f"⚠️  {n_molecules} isomères trouvés, affichage limité aux {max_isomers} premiers"
        )
        liste_smiles = liste_smiles[:max_isomers]
        n_molecules = max_isomers

    print(f"\n🖼️  Génération de {n_molecules} images...")

    for idx, smiles in enumerate(liste_smiles):
        print(f"  [{idx+1}/{n_molecules}] Création de l'isomère {idx+1}...", end="\r")
        graphe = smiles_vers_graphe(smiles)

        # Extraire les couleurs basées sur les éléments
        couleurs = [
            ELEMENTS_VALENCE[graphe.nodes[node]["element"]]["color"]
            for node in graphe.nodes
        ]

        # Créer les labels avec symboles atomiques
        labels = {node: graphe.nodes[node]["element"] for node in graphe.nodes}

        # Créer une nouvelle figure pour chaque isomère
        fig, ax = plt.subplots(figsize=(8, 8))
        pos = nx.spring_layout(graphe)
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
        ax.set_title(
            f"{formule_brute} - Isomère {idx+1}\n{smiles}",
            fontsize=12,
            fontweight="bold",
        )
        plt.tight_layout()
        plt.show(
            block=False
        )  # Non-bloquant: toutes les fenêtres apparaissent en même temps

    print(f"\n✅ {n_molecules} images créées!")
    print("Fermez toutes les fenêtres pour continuer...")
    plt.show()  # Bloque jusqu'à ce que toutes les fenêtres soient fermées


def visualiser_molecules_3d(
    liste_smiles, formule_brute, max_isomers=20, optimize=False
):
    """Affiche toutes les structures possibles en 3D (chaque isomère dans une figure séparée)"""
    n_molecules = len(liste_smiles)

    if n_molecules == 0:
        print("Aucune structure trouvée")
        return

    # Limiter le nombre d'isomères affichés
    if n_molecules > max_isomers:
        print(
            f"⚠️  {n_molecules} isomères trouvés, affichage limité aux {max_isomers} premiers"
        )
        liste_smiles = liste_smiles[:max_isomers]
        n_molecules = max_isomers

    if optimize:
        print(
            f"\n🔬 Génération de {n_molecules} structures 3D avec optimisation (plus lent)..."
        )
    else:
        print(
            f"\n🔬 Génération de {n_molecules} structures 3D rapide (sans optimisation)..."
        )

    for idx, smiles in enumerate(liste_smiles):
        print(f"  [{idx+1}/{n_molecules}] Création de l'isomère {idx+1}...", end="\r")

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

        # Créer une nouvelle figure pour chaque isomère
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection="3d")

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

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(
            f"{formule_brute} - Isomère {idx+1}\n{smiles}",
            fontsize=12,
            fontweight="bold",
        )

        # Désactiver la grille pour un meilleur rendu
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show(
            block=False
        )  # Non-bloquant: toutes les fenêtres apparaissent en même temps

    print(f"\n✅ {n_molecules} structures 3D créées!")
    print("Fermez toutes les fenêtres pour continuer...")
    plt.show()  # Bloque jusqu'à ce que toutes les fenêtres soient fermées


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
        print(f"\nGénération des isomères pour {c}...")
        resultats = generer_tous_squelettes(c)
        print(f"Nombre d'isomères trouvés : {len(resultats)}")

        # Afficher les SMILES
        for idx, smiles in enumerate(resultats, 1):
            print(f"  {idx}. {smiles}")

        # Visualiser toutes les structures en 3D (rapide par défaut)
        # Pour une optimisation géométrique complète (plus lent), utilisez: optimize=True
        visualiser_molecules_3d(
            resultats, formule_brute, max_isomers=20, optimize=False
        )

    except Exception as e:
        print(f"Erreur lors de l'analyse : {e}")


def main():
    analyser_entree_utilisateur()


if __name__ == "__main__":
    main()
