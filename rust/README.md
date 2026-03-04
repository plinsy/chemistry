# Chemistry Isomer Calculator - Rust Edition

High-performance molecular isomer calculator built with Rust and Actix-web.

## Quick Start

```bash
# From the rust/ directory
./run.sh
```

Or manually:
```bash
cargo run --release
```

Visit: http://127.0.0.1:5000

## Documentation

- **README_RUST.md** - Complete documentation
- **SETUP.md** - Installation and troubleshooting  
- **QUICKSTART.md** - Quick reference
- **CONVERSION.md** - Python to Rust conversion details
- **SUMMARY.txt** - Visual overview

## Build

```bash
./build.sh
# or
cargo build --release
```

## Directory Structure

```
rust/
├── src/
│   ├── main.rs         # Entry point
│   ├── chemistry.rs    # Isomer generation
│   └── web.rs          # HTTP routes
├── templates/          # Tera templates (Rust-specific)
├── static/            # CSS/JS assets
├── Cargo.toml         # Dependencies
└── target/            # Build output
```

## Requirements

- Rust 1.70+ (install from https://rustup.rs/)

## Performance

~3-5x faster than Python version for complex molecules.

## Comparison

Run `./compare.sh` to see side-by-side comparison with Python version.
