#!/bin/bash

echo "╔════════════════════════════════════════════════════════╗"
echo "║  Chemistry Isomer Calculator - Version Comparison     ║"
echo "╚════════════════════════════════════════════════════════╝"
echo ""

echo "📊 Comparing Python vs Rust implementations..."
echo ""

# Check if Python version exists
if [ -f "../app.py" ]; then
    echo "✅ Python version (app.py) found"
    PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2)
    echo "   Python version: $PYTHON_VERSION"
else
    echo "❌ Python version not found"
fi

echo ""

# Check if Rust version exists
if [ -f "src/main.rs" ]; then
    echo "✅ Rust version (src/main.rs) found"
    if command -v cargo &> /dev/null; then
        RUST_VERSION=$(rustc --version 2>&1 | cut -d' ' -f2)
        echo "   Rust version: $RUST_VERSION"
    else
        echo "   ⚠️  Rust not installed - run: curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
    fi
else
    echo "❌ Rust version not found"
fi

echo ""
echo "─────────────────────────────────────────────────────────"
echo ""

echo "📈 Feature Comparison:"
echo ""
printf "%-30s %-10s %-10s\n" "Feature" "Python" "Rust"
echo "─────────────────────────────────────────────────────────"
printf "%-30s %-10s %-10s\n" "Formula parsing" "✅" "✅"
printf "%-30s %-10s %-10s\n" "Isomer generation" "✅" "✅"
printf "%-30s %-10s %-10s\n" "Web interface" "✅" "✅"
printf "%-30s %-10s %-10s\n" "3D visualization" "✅" "✅"
printf "%-30s %-10s %-10s\n" "HTMX integration" "✅" "✅"
printf "%-30s %-10s %-10s\n" "Desktop GUI" "✅" "❌"
printf "%-30s %-10s %-10s\n" "RDKit integration" "✅" "❌"
printf "%-30s %-10s %-10s\n" "Type safety" "❌" "✅"
printf "%-30s %-10s %-10s\n" "Compiled binary" "❌" "✅"

echo ""
echo "─────────────────────────────────────────────────────────"
echo ""

echo "⚡ Performance Estimates (C7H16):"
echo "   Python:  ~5-10 seconds"
echo "   Rust:    ~1-2 seconds"
echo "   Speedup: ~3-5x faster"
echo ""

echo "💾 Size Comparison:"
if [ -d "../venv" ]; then
    VENV_SIZE=$(du -sh ../venv 2>/dev/null | cut -f1)
    echo "   Python venv: $VENV_SIZE"
else
    echo "   Python venv: ~500 MB (with RDKit)"
fi

if [ -f "target/release/chemistry" ]; then
    RUST_SIZE=$(du -sh target/release/chemistry 2>/dev/null | cut -f1)
    echo "   Rust binary: $RUST_SIZE"
else
    echo "   Rust binary: ~5-10 MB (after build)"
fi

echo ""
echo "─────────────────────────────────────────────────────────"
echo ""

echo "🚀 Quick Start Commands:"
echo ""
echo "Python version:"
echo "   flask run"
echo "   # Access at http://127.0.0.1:5000"
echo ""
echo "Rust version:"
echo "   ./run.sh"
echo "   # Access at http://127.0.0.1:5000"
echo ""

echo "📚 Documentation:"
echo "   Python:  README.md"
echo "   Rust:    README_RUST.md, SETUP.md, CONVERSION.md"
echo ""

echo "╔════════════════════════════════════════════════════════╗"
echo "║  Both versions serve on port 5000                      ║"
echo "║  Stop one before starting the other!                   ║"
echo "╚════════════════════════════════════════════════════════╝"
