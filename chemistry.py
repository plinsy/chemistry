# Creer notre dictionnaire
ELEMENTS_VALENCE = {
    # Élément : [Valence standard, (Valences alternatives possibles)]
    "H": 1,  # Hydrogène : Toujours 1 lien
    "C": 4,  # Carbone : La base de la chimie organique
    "N": 3,  # Azote : 3 liens (possède 1 paire d'électrons libre)
    "O": 2,  # Oxygène : 2 liens (possède 2 paires libres)
    "P": 3,  # Phosphore : Souvent 3 (mais peut monter à 5)
    "S": 2,  # Soufre : Souvent 2 (mais peut monter à 4 ou 6)
    "F": 1,  # Fluor : Halogène (1 lien)
    "Cl": 1,  # Chlore : Halogène (1 lien)
    "Br": 1,  # Brome : Halogène (1 lien)
    "I": 1,  # Iode : Halogène (1 lien)
}


def calculer_insaturation(carbone, hydrogene, azote=0, halogene=0):
    # Formule standard du Degree of Unsaturation (DoU)
    dou = (2 * carbone + 2 + azote - hydrogene - halogene) / 2
    return int(dou)


import re


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

    except Exception as e:
        print(f"Erreur lors de l'analyse : {e}")
