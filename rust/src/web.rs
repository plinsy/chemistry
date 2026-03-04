use actix_files::Files;
use actix_web::{web, HttpResponse, Result};
use serde::{Deserialize};
use tera::{Context, Tera};

use crate::chemistry::{generate_all_skeletons, generate_3d_coords, parse_formula, Bond, Isomer};

#[derive(Deserialize)]
pub struct FormulaRequest {
    formula: String,
}

// #[derive(Serialize)]
// pub struct ErrorResponse {
//     error: String,
// }

lazy_static::lazy_static! {
    pub static ref TEMPLATES: Tera = {
        match Tera::new("templates/**/*.html") {
            Ok(t) => t,
            Err(e) => {
                eprintln!("Parsing error(s): {}", e);
                std::process::exit(1);
            }
        }
    };
}

pub fn config(cfg: &mut web::ServiceConfig) {
    cfg.service(web::resource("/").route(web::get().to(index)))
        .service(web::resource("/api/calculate").route(web::post().to(calculate_isomers)))
        .service(Files::new("/static", "static").show_files_listing());
}

async fn index() -> Result<HttpResponse> {
    let context = Context::new();
    let body = TEMPLATES
        .render("pages/index.html", &context)
        .map_err(|e| actix_web::error::ErrorInternalServerError(e))?;

    Ok(HttpResponse::Ok().content_type("text/html").body(body))
}

async fn calculate_isomers(form: web::Form<FormulaRequest>) -> Result<HttpResponse> {
    let formula = form.formula.trim().to_uppercase();

    if formula.is_empty() {
        return render_error("Please enter a chemical formula");
    }

    // Parse the formula
    let atoms = match parse_formula(&formula) {
        Ok(a) => a,
        Err(e) => return render_error(&format!("Invalid formula: {}", e)),
    };

    // Extract carbon and hydrogen counts
    let c = *atoms.get("C").unwrap_or(&0);
    let h = *atoms.get("H").unwrap_or(&0);

    if c == 0 {
        return render_error("Formula must contain carbon atoms");
    }

    // Generate all isomers
    let smiles_list = generate_all_skeletons(c, Some(h));

    // Prepare isomers data
    let mut isomers = Vec::new();
    for smiles in smiles_list {
        // Generate simplified 3D coordinates
        let num_carbons = smiles.chars().filter(|&c| c == 'C').count();
        let atoms_data = generate_3d_coords(num_carbons);
        
        // Generate bonds from structure
        let bonds_data = generate_bonds_from_smiles(&smiles);

        isomers.push(Isomer {
            smiles: smiles.clone(),
            atoms: atoms_data,
            bonds: bonds_data,
        });
    }

    // Render the molecules partial
    let mut context = Context::new();
    context.insert("isomers", &isomers);
    context.insert("formula", &formula);

    let body = TEMPLATES
        .render("partials/molecules.html", &context)
        .map_err(|e| actix_web::error::ErrorInternalServerError(e))?;

    Ok(HttpResponse::Ok().content_type("text/html").body(body))
}

fn render_error(error: &str) -> Result<HttpResponse> {
    let mut context = Context::new();
    context.insert("error", error);

    let body = TEMPLATES
        .render("partials/error.html", &context)
        .map_err(|e| actix_web::error::ErrorInternalServerError(e))?;

    Ok(HttpResponse::BadRequest()
        .content_type("text/html")
        .body(body))
}

fn generate_bonds_from_smiles(smiles: &str) -> Vec<Bond> {
    let mut bonds = Vec::new();
    let mut current_idx = 0usize;
    let mut parent_stack = Vec::new();
    let mut next_bond_type = 1.0;

    for ch in smiles.chars() {
        match ch {
            'C' => {
                if let Some(&parent_idx) = parent_stack.last() {
                    bonds.push(Bond {
                        start: parent_idx,
                        end: current_idx,
                        bond_type: next_bond_type,
                    });
                    next_bond_type = 1.0;
                }
                parent_stack.push(current_idx);
                current_idx += 1;
            }
            '(' => {
                // Branch start - keep current parent
            }
            ')' => {
                // Branch end - pop to previous parent
                if parent_stack.len() > 1 {
                    parent_stack.pop();
                }
            }
            '=' => {
                next_bond_type = 2.0;
            }
            '#' => {
                next_bond_type = 3.0;
            }
            _ => {}
        }
    }

    bonds
}
