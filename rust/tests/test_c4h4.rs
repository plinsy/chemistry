/// Integration tests for C4H4 (butadiyne/diacetylene)
/// C4H4 has high unsaturation (DoU = 3) with two triple bonds
use chemistry::{calculer_insaturation, molecule_to_smiles, parse_formule, Generator};

#[test]
fn test_c4h4_parse_formula() {
    let (c, h) = parse_formule("C4H4");

    assert_eq!(c, 4, "Should have 4 carbons");
    assert_eq!(h, 4, "Should have 4 hydrogens");
}

#[test]
fn test_c4h4_degree_of_unsaturation() {
    // C4H4: DoU = (2*4 + 2 - 4) / 2 = (8 + 2 - 4) / 2 = 6 / 2 = 3
    let dou = calculer_insaturation(4, 4);
    assert_eq!(dou, 3, "C4H4 should have 3 degrees of unsaturation (e.g., two triple bonds)");
}

#[test]
fn test_c4h4_generate_isomers() {
    let mut generator = Generator::new(4, 4);
    generator.run();

    // C4H4 should have multiple isomers due to high unsaturation
    assert!(!generator.results.is_empty(), "C4H4 should generate at least one isomer");
}

#[test]
fn test_c4h4_can_form_linear_conjugated_system() {
    let mut generator = Generator::new(4, 4);
    generator.run();

    // At least one isomer should contain a triple bond (bond value 3 in the adjacency matrix)
    let has_triple = generator
        .results
        .iter()
        .any(|mol| mol.adj_matrix.iter().any(|row: &Vec<u8>| row.iter().any(|&b| b == 3)));

    assert!(has_triple, "C4H4 should generate isomers with triple bonds");
}

#[test]
fn test_c4h4_smiles_not_empty() {
    let mut generator = Generator::new(4, 4);
    generator.run();

    for mol in &generator.results {
        let smiles = molecule_to_smiles(mol);
        assert!(!smiles.is_empty(), "SMILES representation should not be empty");
        assert!(smiles.contains('C'), "SMILES should contain carbon atoms");
    }
}
