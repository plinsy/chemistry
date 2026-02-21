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
