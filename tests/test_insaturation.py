from chemistry import calculer_insaturation


def test_alcanes_satures():
    """Vérifie les molécules simples sans doubles liaisons ni cycles."""
    assert calculer_insaturation(1, 4) == 0  # Méthane
    assert calculer_insaturation(4, 10) == 0  # Butane
    assert calculer_insaturation(8, 18) == 0  # Octane
    print("✅ Test Alcanes : OK")


def test_insaturations_simples():
    """Vérifie les molécules avec une double liaison ou un cycle."""
    assert calculer_insaturation(2, 4) == 1  # Éthène (1 double)
    assert calculer_insaturation(4, 8) == 1  # Cyclobutane (1 cycle)
    assert calculer_insaturation(2, 2) == 2  # Acétylène (1 triple = 2 DoU)
    print("✅ Test Insaturations : OK")


def test_avec_heteroatomes_o_n():
    """Vérifie l'influence de l'Oxygène (neutre) et de l'Azote (+1)."""
    # L'Oxygène ne change rien au calcul
    assert calculer_insaturation(2, 6) == 0  # Éthanol (C2H6O)
    assert calculer_insaturation(3, 6) == 1  # Acétone (C3H6O, 1 double C=O)

    # L'Azote ajoute une unité au calcul
    assert calculer_insaturation(1, 5, azote=1) == 0  # Méthylamine (CH5N)
    assert (
        calculer_insaturation(5, 5, azote=1) == 4
    )  # Pyridine (C5H5N : 1 cycle + 3 doubles)
    print("✅ Test Oxygène & Azote : OK")


def test_avec_halogenès():
    """Vérifie que Cl, Br, F, I sont gérés comme des Hydrogènes."""
    # CH3Cl (un H remplacé par un Cl)
    assert calculer_insaturation(1, 3, halogene=1) == 0
    # C2H2Cl2 (Dichlorure d'éthylène, DoU doit être 1)
    assert calculer_insaturation(2, 2, halogene=2) == 1
    print("✅ Test Halogènes : OK")
