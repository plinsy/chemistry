use eframe::egui;
use petgraph::graph::Graph;
use petgraph::algo::connected_components;
use regex::Regex;
use std::collections::HashSet;

// --- LOGIQUE CHIMIQUE ---

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Molecule {
    pub adj_matrix: Vec<Vec<u8>>, // 0: pas de lien, 1: simple, 2: double, 3: triple
    pub n_atoms: usize,
}

pub fn calculer_insaturation(c: i32, h: i32) -> i32 {
    (2 * c + 2 - h) / 2
}

pub fn parse_formule(formule: &str) -> (i32, i32) {
    let re = Regex::new(r"([A-Z][a-z]?)(\d*)").unwrap();
    let mut c = 0;
    let mut h = 0;

    for cap in re.captures_iter(formule) {
        let sym = &cap[1];
        let qty = cap[2].parse::<i32>().unwrap_or(1);
        if sym == "C" { c = qty; }
        else if sym == "H" { h = qty; }
    }
    (c, h)
}

// --- GÉNÉRATEUR D'ISOMÈRES (BACKTRACKING OPTIMISÉ) ---

pub struct Generator {
    pub n_carbons: usize,
    pub n_hydrogens: usize,
    pub results: HashSet<Molecule>,
}

impl Generator {
    pub fn new(c: usize, h: usize) -> Self {
        Self { n_carbons: c, n_hydrogens: h, results: HashSet::new() }
    }

    pub fn run(&mut self) {
        let mut matrix = vec![vec![0u8; self.n_carbons]; self.n_carbons];
        self.backtrack(&mut matrix, 0, 1);
    }

    fn backtrack(&mut self, matrix: &mut Vec<Vec<u8>>, row: usize, col: usize) {
        if row == self.n_carbons - 1 {
            if self.is_valid(matrix) {
                self.results.insert(Molecule { adj_matrix: matrix.clone(), n_atoms: self.n_carbons });
            }
            return;
        }

        let (next_row, next_col) = if col == self.n_carbons - 1 {
            (row + 1, row + 2)
        } else {
            (row, col + 1)
        };

        // Essayer 0 (pas de liaison), 1 (simple), 2 (double), 3 (triple)
        for bond in 0..=3 {
            if self.can_add_bond(matrix, row, col, bond) {
                matrix[row][col] = bond;
                matrix[col][row] = bond;
                self.backtrack(matrix, next_row, next_col);
                matrix[row][col] = 0;
                matrix[col][row] = 0;
            }
        }
    }

    fn can_add_bond(&self, matrix: &[Vec<u8>], r: usize, c: usize, bond: u8) -> bool {
        let sum_r: u8 = matrix[r].iter().sum::<u8>() + bond;
        let sum_c: u8 = matrix[c].iter().sum::<u8>() + bond;
        sum_r <= 4 && sum_c <= 4
    }

    fn is_valid(&self, matrix: &[Vec<u8>]) -> bool {
        // Vérifier la connectivité
        let mut g = Graph::<usize, u8, petgraph::Undirected>::new_undirected();
        let nodes: Vec<_> = (0..self.n_carbons).map(|i| g.add_node(i)).collect();
        for i in 0..self.n_carbons {
            for j in i + 1..self.n_carbons {
                if matrix[i][j] > 0 { g.add_edge(nodes[i], nodes[j], matrix[i][j]); }
            }
        }
        
        if connected_components(&g) != 1 { return false; }

        // Vérifier le nombre d'hydrogènes : H = sum(4 - valence_carbone)
        let mut total_h = 0;
        for i in 0..self.n_carbons {
            let valence: u8 = matrix[i].iter().sum();
            total_h += 4 - valence as i32;
        }
        total_h == self.n_hydrogens as i32
    }
}

// --- SMILES GENERATOR ---

/// Convert an adjacency matrix molecule to a SMILES string via DFS.
pub fn molecule_to_smiles(mol: &Molecule) -> String {
    let n = mol.n_atoms;
    if n == 0 {
        return String::new();
    }
    let mut visited = vec![false; n];
    let mut out = String::new();
    smiles_dfs(&mol.adj_matrix, 0, n, &mut visited, &mut out);
    out
}

fn smiles_dfs(
    matrix: &[Vec<u8>],
    node: usize,
    n: usize,
    visited: &mut Vec<bool>,
    out: &mut String,
) {
    visited[node] = true;
    out.push('C');

    let neighbors: Vec<(usize, u8)> = (0..n)
        .filter(|&j| matrix[node][j] > 0 && !visited[j])
        .map(|j| (j, matrix[node][j]))
        .collect();

    for (i, &(next, bond)) in neighbors.iter().enumerate() {
        let is_last = i == neighbors.len() - 1;
        if !is_last {
            out.push('(');
        }
        match bond {
            2 => out.push('='),
            3 => out.push('#'),
            _ => {}
        }
        smiles_dfs(matrix, next, n, visited, out);
        if !is_last {
            out.push(')');
        }
    }
}

// --- INTERFACE GRAPHIQUE (EGUI) ---

struct MyApp {
    input: String,
    isomers: Vec<Molecule>,
    searching: bool,
}

impl Default for MyApp {
    fn default() -> Self {
        Self { input: "C4H10".to_owned(), isomers: vec![], searching: false }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Générateur d'Isomères Rust");
            ui.horizontal(|ui| {
                ui.label("Formule:");
                ui.text_edit_singleline(&mut self.input);
                if ui.button("Générer").clicked() {
                    let (c, h) = parse_formule(&self.input);
                    let mut generator = Generator::new(c as usize, h as usize);
                    generator.run();
                    self.isomers = generator.results.into_iter().collect();
                }
            });

            ui.separator();

            egui::ScrollArea::vertical().show(ui, |ui| {
                for (idx, mol) in self.isomers.iter().enumerate() {
                    ui.group(|ui| {
                        ui.label(format!("Isomère #{}", idx + 1));
                        // Ici, on pourrait dessiner avec egui::Painter
                        ui.label(format!("Matrice d'adjacence: {:?}", mol.adj_matrix));
                    });
                }
            });
        });
    }
}

fn main() -> eframe::Result<()> {
    let native_options = eframe::NativeOptions::default();
    eframe::run_native("Molecule Isomer Visualizer", native_options, Box::new(|_| Ok(Box::new(MyApp::default()))))
}