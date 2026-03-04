# Chemistry Isomer Calculator

Molecular isomer generation and visualization tool. Available in both Python (Flask) and Rust (Actix-web) implementations.

## Project Structure

```
chemistry/
├── 🐍 Python Version (Flask)
│   ├── app.py              # Flask application
│   ├── chemistry.py        # Core chemistry logic
│   ├── templates/          # Jinja2 templates
│   └── static/            # CSS/JS assets
│
└── 🦀 Rust Version (Actix-web)
    ├── rust/
    │   ├── src/           # Rust source code
    │   ├── templates/     # Tera templates
    │   ├── static/        # CSS/JS assets
    │   └── README.md      # Rust documentation
```

## Quick Start

### Python Version

```bash
# Install dependencies
pip install -r requirements.txt  # or use pyproject.toml

# Run Flask app
flask run
# or
python app.py
```

Visit: http://127.0.0.1:5000

### Rust Version

```bash
# Navigate to rust folder
cd rust

# Run (will build automatically)
./run.sh
# or
cargo run --release
```

Visit: http://127.0.0.1:5000

## Features

- ✅ Parse chemical formulas (C4H10, C6H6, etc.)
- ✅ Generate all structural isomers
- ✅ Calculate degree of unsaturation
- ✅ Handle single, double, and triple bonds
- ✅ Interactive 3D molecular visualization (Three.js)
- ✅ Real-time HTMX updates
- ✅ Acyclic and cyclic structures

## Version Comparison

| Feature | Python/Flask | Rust/Actix |
|---------|-------------|------------|
| Performance | Baseline | **3-5x faster** |
| Binary Size | ~500 MB | **~5-10 MB** |
| Type Safety | Runtime | **Compile-time** |
| Dependencies | Many | Few |
| Deployment | Python + venv | **Single binary** |

## Documentation

### Python Version
- This README
- Code comments in `app.py` and `chemistry.py`

### Rust Version
- `rust/README.md` - Overview
- `rust/README_RUST.md` - Complete documentation
- `rust/SETUP.md` - Installation guide
- `rust/QUICKSTART.md` - Quick reference
- `rust/CONVERSION.md` - Migration details
- `rust/SUMMARY.txt` - Visual overview

## Testing

Try these formulas:
- `C4H10` - Butane (2 isomers)
- `C5H12` - Pentane (3 isomers)
- `C6H6` - Benzene (complex structures)
- `C7H16` - Heptane (performance test)

## Requirements

### Python
- Python 3.9+
- Flask
- RDKit
- NetworkX
- matplotlib
- rich

### Rust
- Rust 1.70+
- See `rust/Cargo.toml` for dependencies

## Development

Both versions share the same frontend but have version-specific templates:
- **Python**: Uses Jinja2 templates in `templates/`
- **Rust**: Uses Tera templates in `rust/templates/`

Static files (CSS/JS) are duplicated in both locations.

## License

Educational project - free to use and modify.

## Contributing

Both Python and Rust implementations welcome improvements!

---

**Note**: Both versions run on port 5000. Stop one before starting the other, or modify the port in the respective source files.
