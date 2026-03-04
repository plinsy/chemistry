# 🦀 RUST VERSION - QUICK REFERENCE

## 🚀 Quick Start (3 commands)
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh  # Install Rust
source $HOME/.cargo/env                                           # Add to PATH
./run.sh                                                          # Start server
```

## 📋 Common Commands
```bash
# Development (fast compile, slow runtime)
cargo run

# Production (slow compile, fast runtime)  
cargo run --release

# Build only
cargo build --release

# Check for errors (no build)
cargo check

# Auto-reload on changes
cargo watch -x run

# Run tests
cargo test

# Clean build artifacts
cargo clean
```

## 📂 Key Files
```
src/main.rs       - Entry point, starts HTTP server
src/chemistry.rs  - Isomer generation algorithm
src/web.rs        - Routes, handlers, templates
Cargo.toml        - Dependencies & configuration
```

## 🔗 URLs
```
http://127.0.0.1:5000         - Main application
http://127.0.0.1:5000/static/ - Static files
POST /api/calculate           - API endpoint
```

## 🧪 Test Formulas
```
C4H10  - Simple (2 isomers)
C5H12  - Moderate (3 isomers)
C6H6   - Complex (many structures)
C7H16  - Performance test (9 isomers)
```

## 📖 Documentation
```bash
cat SUMMARY.txt      # Quick overview
cat README_RUST.md   # Complete docs
cat SETUP.md         # Installation help
cat CONVERSION.md    # Technical details
./compare.sh         # Python vs Rust
```

## 🛠️ Troubleshooting
```bash
# Cargo not found
source $HOME/.cargo/env

# Port in use
kill -9 $(lsof -ti:5000)

# Build errors
cargo clean && cargo build --release

# Check Rust version (need 1.70+)
rustc --version
```

## ⚡ Performance Tips
```bash
# Always use --release for benchmarks
cargo run --release

# Profile build time
cargo build --release --timings

# Show assembly 
cargo rustc --release -- --emit asm
```

## 🔄 Switching Versions
```bash
# Stop Flask
^C (in flask terminal)

# Start Rust
./run.sh

# Or vice versa
^C (in cargo terminal)
flask run
```

## 📦 Dependencies
```toml
actix-web = "4.5"      # Web framework
actix-files = "0.6"    # Static files
tera = "1.19"          # Templates
petgraph = "0.6"       # Graphs
serde = "1.0"          # Serialization
regex = "1.10"         # Pattern matching
```

## 🎯 Next Steps
1. Read README_RUST.md
2. Try example formulas  
3. Compare with Python version
4. Explore src/ code
5. Run benchmarks

---
💡 Tip: Keep SUMMARY.txt open in another terminal for reference!
