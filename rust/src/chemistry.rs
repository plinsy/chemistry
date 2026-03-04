use petgraph::graph::{Graph, NodeIndex};
use petgraph::Undirected;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub element: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Bond {
    pub start: usize,
    pub end: usize,
    #[serde(rename = "type")]
    pub bond_type: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Isomer {
    pub smiles: String,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

type MolGraph = Graph<String, u8, Undirected>;

/// Parse a chemical formula like "C4H10" into a hashmap of element counts
pub fn parse_formula(formula: &str) -> Result<HashMap<String, usize>, String> {
    let re = Regex::new(r"([A-Z][a-z]?)(\d*)").unwrap();
    let mut atoms = HashMap::new();

    for cap in re.captures_iter(formula) {
        let symbol = cap[1].to_string();
        let count = if cap[2].is_empty() {
            1
        } else {
            cap[2].parse::<usize>().map_err(|e| e.to_string())?
        };

        *atoms.entry(symbol).or_insert(0) += count;
    }

    Ok(atoms)
}

/// Calculate the degree of unsaturation (DoU)
pub fn calculate_unsaturation(carbon: usize, hydrogen: usize, nitrogen: usize, halogen: usize) -> i32 {
    ((2 * carbon + 2 + nitrogen) as i32 - hydrogen as i32 - halogen as i32) / 2
}

/// Generate all possible isomers for a given molecular formula
pub fn generate_all_skeletons(n_carbons: usize, n_hydrogens: Option<usize>) -> Vec<String> {
    let mut smiles_unique = HashSet::new();

    // Generate acyclic structures (trees)
    generate_acyclic_structures(n_carbons, n_hydrogens, &mut smiles_unique);

    // Generate cyclic structures if there's unsaturation
    if let Some(h) = n_hydrogens {
        if n_carbons >= 3 {
            let dou = calculate_unsaturation(n_carbons, h, 0, 0);
            if dou > 0 {
                generate_cyclic_structures(n_carbons, n_hydrogens, &mut smiles_unique);
            }
        }
    }

    smiles_unique.into_iter().collect()
}

fn generate_acyclic_structures(
    n_carbons: usize,
    n_hydrogens: Option<usize>,
    smiles_set: &mut HashSet<String>,
) {
    if n_carbons == 0 {
        return;
    }

    let mut graph = Graph::<String, u8, Undirected>::new_undirected();
    let nodes: Vec<NodeIndex> = (0..n_carbons)
        .map(|_| graph.add_node("C".to_string()))
        .collect();

    explore_tree(&mut graph, &nodes, n_carbons - 1, n_hydrogens, smiles_set);
}

fn explore_tree(
    graph: &mut MolGraph,
    nodes: &[NodeIndex],
    target_edges: usize,
    n_hydrogens: Option<usize>,
    smiles_set: &mut HashSet<String>,
) {
    if graph.edge_count() == target_edges {
        if is_connected(graph) {
            explore_bond_types(graph, n_hydrogens, smiles_set);
        }
        return;
    }

    let n = nodes.len();
    for i in 0..n {
        for j in (i + 1)..n {
            if !graph.contains_edge(nodes[i], nodes[j]) {
                let valences = count_valences(graph, nodes);
                if valences[i] < 4 && valences[j] < 4 {
                    graph.add_edge(nodes[i], nodes[j], 1);
                    explore_tree(graph, nodes, target_edges, n_hydrogens, smiles_set);
                    if let Some(edge) = graph.find_edge(nodes[i], nodes[j]) {
                        graph.remove_edge(edge);
                    }
                }
            }
        }
    }
}

fn generate_cyclic_structures(
    n_carbons: usize,
    n_hydrogens: Option<usize>,
    smiles_set: &mut HashSet<String>,
) {
    if let Some(h) = n_hydrogens {
        let dou = calculate_unsaturation(n_carbons, h, 0, 0);
        let max_edges = if dou >= 3 {
            (2 * n_carbons).min(n_carbons + dou as usize + 5)
        } else {
            ((n_carbons * 3) / 2 + 2).min(n_carbons + dou as usize + 3)
        };

        for n_edges in n_carbons..=max_edges {
            explore_cyclic_graphs(n_carbons, n_edges, n_hydrogens, smiles_set);
        }
    }
}

fn explore_cyclic_graphs(
    n_carbons: usize,
    n_edges: usize,
    n_hydrogens: Option<usize>,
    smiles_set: &mut HashSet<String>,
) {
    let mut graph = Graph::<String, u8, Undirected>::new_undirected();
    let nodes: Vec<NodeIndex> = (0..n_carbons)
        .map(|_| graph.add_node("C".to_string()))
        .collect();

    let mut node_pairs = Vec::new();
    for i in 0..n_carbons {
        for j in (i + 1)..n_carbons {
            node_pairs.push((nodes[i], nodes[j]));
        }
    }

    build_graph(&mut graph, &nodes, &node_pairs, 0, n_edges, n_hydrogens, smiles_set);
}

fn build_graph(
    graph: &mut MolGraph,
    nodes: &[NodeIndex],
    node_pairs: &[(NodeIndex, NodeIndex)],
    pair_idx: usize,
    target_edges: usize,
    n_hydrogens: Option<usize>,
    smiles_set: &mut HashSet<String>,
) {
    let current_edges = graph.edge_count();

    if current_edges == target_edges {
        if is_connected(graph) {
            explore_bond_types(graph, n_hydrogens, smiles_set);
        }
        return;
    }

    if current_edges > target_edges || pair_idx >= node_pairs.len() {
        return;
    }

    let edges_needed = target_edges - current_edges;
    let pairs_left = node_pairs.len() - pair_idx;

    if pairs_left < edges_needed {
        return;
    }

    let (i, j) = node_pairs[pair_idx];

    // Option 1: Add this edge
    let valences = count_valences(graph, nodes);
    let i_idx = nodes.iter().position(|&n| n == i).unwrap();
    let j_idx = nodes.iter().position(|&n| n == j).unwrap();

    if valences[i_idx] < 4 && valences[j_idx] < 4 {
        graph.add_edge(i, j, 1);
        build_graph(graph, nodes, node_pairs, pair_idx + 1, target_edges, n_hydrogens, smiles_set);
        if let Some(edge) = graph.find_edge(i, j) {
            graph.remove_edge(edge);
        }
    }

    // Option 2: Don't add this edge
    build_graph(graph, nodes, node_pairs, pair_idx + 1, target_edges, n_hydrogens, smiles_set);
}

fn explore_bond_types(
    graph: &mut MolGraph,
    n_hydrogens: Option<usize>,
    smiles_set: &mut HashSet<String>,
) {
    let edges: Vec<_> = graph.edge_indices().collect();
    explore_bond_types_recursive(graph, &edges, 0, n_hydrogens, smiles_set);
}

fn explore_bond_types_recursive(
    graph: &mut MolGraph,
    edges: &[petgraph::graph::EdgeIndex],
    edge_idx: usize,
    n_hydrogens: Option<usize>,
    smiles_set: &mut HashSet<String>,
) {
    if edge_idx == edges.len() {
        // Check if all valences are valid
        let nodes: Vec<_> = graph.node_indices().collect();
        let valences = count_valences(graph, &nodes);
        if valences.iter().all(|&v| v <= 4) {
            if let Some(smiles) = graph_to_smiles(graph, n_hydrogens) {
                smiles_set.insert(smiles);
            }
        }
        return;
    }

    let edge = edges[edge_idx];
    let (u, v) = graph.edge_endpoints(edge).unwrap();
    let nodes: Vec<_> = graph.node_indices().collect();

    // Try bond types 1, 2, 3
    for bond_type in [1u8, 2u8, 3u8] {
        let valences = count_valences(graph, &nodes);
        let u_idx = nodes.iter().position(|&n| n == u).unwrap();
        let v_idx = nodes.iter().position(|&n| n == v).unwrap();

        let current_bond = *graph.edge_weight(edge).unwrap();
        let current_u = valences[u_idx] - current_bond as usize;
        let current_v = valences[v_idx] - current_bond as usize;

        if current_u + bond_type as usize <= 4 && current_v + bond_type as usize <= 4 {
            *graph.edge_weight_mut(edge).unwrap() = bond_type;
            explore_bond_types_recursive(graph, edges, edge_idx + 1, n_hydrogens, smiles_set);
        }
    }

    // Restore default
    *graph.edge_weight_mut(edge).unwrap() = 1;
}

fn count_valences(graph: &MolGraph, nodes: &[NodeIndex]) -> Vec<usize> {
    nodes
        .iter()
        .map(|&node| {
            graph
                .edges(node)
                .map(|e| *e.weight() as usize)
                .sum::<usize>()
        })
        .collect()
}

fn is_connected(graph: &MolGraph) -> bool {
    if graph.node_count() == 0 {
        return true;
    }
    
    use petgraph::visit::Dfs;
    let start = graph.node_indices().next().unwrap();
    let mut dfs = Dfs::new(&*graph, start);
    let mut count = 0;
    
    while dfs.next(&*graph).is_some() {
        count += 1;
    }
    
    count == graph.node_count()
}

fn graph_to_smiles(graph: &MolGraph, n_hydrogens: Option<usize>) -> Option<String> {
    // Simplified SMILES generation
    // In a real implementation, you would use a proper chemistry library
    // Here we create a basic SMILES representation
    
    if graph.node_count() == 0 {
        return None;
    }

    let mut smiles = String::new();
    let mut visited = HashSet::new();
    let start = graph.node_indices().next().unwrap();

    dfs_smiles(graph, start, None, &mut visited, &mut smiles, n_hydrogens);

    Some(smiles)
}

fn dfs_smiles(
    graph: &MolGraph,
    node: NodeIndex,
    parent: Option<NodeIndex>,
    visited: &mut HashSet<NodeIndex>,
    smiles: &mut String,
    _n_hydrogens: Option<usize>,
) {
    visited.insert(node);
    smiles.push('C');

    let mut neighbors: Vec<_> = graph.neighbors(node).collect();
    neighbors.retain(|&n| Some(n) != parent);

    let mut first = true;
    for neighbor in neighbors {
        if !visited.contains(&neighbor) {
            if !first {
                smiles.push('(');
            }

            // Add bond notation if not single
            if let Some(edge) = graph.find_edge(node, neighbor) {
                let bond_type = *graph.edge_weight(edge).unwrap();
                match bond_type {
                    2 => smiles.push('='),
                    3 => smiles.push('#'),
                    _ => {}
                }
            }

            dfs_smiles(graph, neighbor, Some(node), visited, smiles, _n_hydrogens);

            if !first {
                smiles.push(')');
            }
            first = false;
        }
    }
}

/// Generate 3D coordinates for visualization (simplified)
pub fn generate_3d_coords(n_atoms: usize) -> Vec<Atom> {
    use std::f64::consts::PI;
    
    (0..n_atoms)
        .map(|i| {
            let angle = 2.0 * PI * (i as f64) / (n_atoms as f64);
            let radius = 2.0;
            Atom {
                element: "C".to_string(),
                x: radius * angle.cos(),
                y: radius * angle.sin(),
                z: 0.0,
            }
        })
        .collect()
}
