mod chemistry;
mod web;

use actix_web::{middleware, App, HttpServer};
use std::io::Result;

#[actix_web::main]
async fn main() -> Result<()> {
    println!("Starting Chemistry Isomer Calculator...");
    println!("Server running at http://127.0.0.1:5000");

    HttpServer::new(|| {
        App::new()
            .wrap(middleware::Logger::default())
            .configure(web::config)
    })
    .bind(("127.0.0.1", 5000))?
    .run()
    .await
}
