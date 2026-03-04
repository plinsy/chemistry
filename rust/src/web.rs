use actix_files::Files;
use actix_web::{web, HttpResponse, Result};
use serde::{Deserialize, Serialize};
use tera::{Context, Tera};

use crate::chemistry::{molecule_to_smiles, parse_formule, Generator};

#[derive(Deserialize)]
pub struct FormulaRequest {
    formula: String,
}

#[derive(Serialize)]
pub struct AtomView {
    pub index: usize,
    pub element: String,
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(Serialize)]
pub struct BondView {
    pub start: usize,
    pub end: usize,
    pub bond_type: u8,
}

#[derive(Serialize)]
pub struct IsomerView {
    pub smiles: String,
    pub atoms: Vec<AtomView>,
    pub bonds: Vec<BondView>,
}

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

    let (c, h) = parse_formule(&formula);

    if c <= 0 {
        return render_error("Formula must contain carbon atoms");
    }

    // Generate all isomers using the backtracking generator from chemistry.rs
    let mut generator = Generator::new(c as usize, h as usize);
    generator.run();

    let isomers: Vec<IsomerView> = generator
        .results
        .into_iter()
        .map(|mol| {
            let n = mol.n_atoms;
            let smiles = molecule_to_smiles(&mol);

            // Simple circular layout for atom positions
            let atoms = (0..n)
                .map(|i| {
                    let angle = 2.0 * std::f32::consts::PI * i as f32 / n as f32;
                    let radius = 1.5 * (n as f32).sqrt();
                    AtomView {
                        index: i,
                        element: "C".to_string(),
                        x: radius * angle.cos(),
                        y: radius * angle.sin(),
                        z: 0.0,
                    }
                })
                .collect();

            // Bonds from the upper triangle of the adjacency matrix
            let bonds = mol
                .adj_matrix
                .iter()
                .enumerate()
                .flat_map(|(i, row)| {
                    row.iter()
                        .enumerate()
                        .filter(move |&(j, &b)| j > i && b > 0)
                        .map(move |(j, &b)| BondView {
                            start: i,
                            end: j,
                            bond_type: b,
                        })
                })
                .collect();

            IsomerView { smiles, atoms, bonds }
        })
        .collect();

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



