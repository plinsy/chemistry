# Setup Guide - Chemistry Isomer Calculator (Rust)

## Quick Setup Instructions

### Step 1: Install Rust

If you don't have Rust installed, run:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Follow the prompts and restart your terminal, or run:
```bash
source $HOME/.cargo/env
```

Verify installation:
```bash
rustc --version
cargo --version
```

### Step 2: Build the Project

From the project directory:

```bash
./build.sh
```

Or manually:
```bash
cargo build --release
```

### Step 3: Run the Application

```bash
./run.sh
```

Or manually:
```bash
cargo run --release
```

### Step 4: Access the Application

Open your browser and go to:
```
http://127.0.0.1:5000
```

## Development Workflow

### Build for Development (faster compile, slower runtime)
```bash
cargo build
cargo run
```

### Build for Production (slower compile, faster runtime)
```bash
cargo build --release
cargo run --release
```

### Watch Mode (auto-rebuild on changes)
```bash
# Install cargo-watch first
cargo install cargo-watch

# Run with auto-reload
cargo watch -x run
```

## Troubleshooting

### "cargo: command not found"

Make sure Rust is installed and in your PATH:
```bash
source $HOME/.cargo/env
# or
export PATH="$HOME/.cargo/bin:$PATH"
```

### Port 5000 already in use

The Python Flask server might still be running. Kill it:
```bash
# Find the process
lsof -ti:5000

# Kill it
kill -9 $(lsof -ti:5000)
```

Or change the port in `src/main.rs`:
```rust
.bind(("127.0.0.1", 8080))?  // Change 5000 to 8080
```

### Templates not found

Make sure you're running from the project root directory where `templates/` folder exists.

### Build errors

Clean and rebuild:
```bash
cargo clean
cargo build --release
```

## Comparing with Python Version

### Running Both Versions

**Python/Flask version:**
```bash
flask run
# or
python app.py
```

**Rust/Actix version:**
```bash
cargo run --release
```

### Performance Comparison

Try generating isomers for `C7H16`:
- Python: ~5-10 seconds
- Rust: ~1-2 seconds

## What's Different?

### Technology Stack

| Component | Python Version | Rust Version |
|-----------|---------------|--------------|
| Web Framework | Flask | Actix-web |
| Template Engine | Jinja2 | Tera |
| Graph Library | NetworkX | petgraph |
| Chemistry | RDKit | Custom implementation |
| 3D Viz | Matplotlib/Tk | Three.js (web) |

### Features Maintained

✅ Formula parsing  
✅ Isomer generation  
✅ Acyclic & cyclic structures  
✅ Single/double/triple bonds  
✅ 3D visualization  
✅ HTMX dynamic updates  
✅ Responsive UI  

### Features Simplified

⚠️ 3D coordinates: Geometric layout instead of RDKit's conformer generation  
⚠️ SMILES: Custom generation instead of RDKit canonical SMILES  
⚠️ No desktop GUI (web-only)  

## Next Steps

1. Test with various formulas: C4H10, C5H12, C6H6, C3H4
2. Compare performance with Python version
3. Explore the code in `src/` directory
4. Customize styling in `static/css/styles.css`
5. Add new features!

## Resources

- [Rust Book](https://doc.rust-lang.org/book/)
- [Actix-web Documentation](https://actix.rs/)
- [petgraph Documentation](https://docs.rs/petgraph/)
- [Tera Template Engine](https://tera.netlify.app/)

## Support

For issues, check:
1. Rust version: `rustc --version` (should be 1.70+)
2. All files present: `ls -la src/ templates/ static/`
3. Clean build: `cargo clean && cargo build --release`

Happy molecule generation! 🧪⚛️🦀
