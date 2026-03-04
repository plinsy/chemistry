#!/bin/bash

echo "🦀 Building Chemistry Isomer Calculator (Rust Edition)..."

# Navigate to rust directory
cd "$(dirname "$0")"

# Check if Rust is installed
if ! command -v cargo &> /dev/null; then
    echo "❌ Cargo not found. Please install Rust from https://rustup.rs/"
    exit 1
fi

# Build the project
echo "📦 Building in release mode..."
cargo build --release

if [ $? -eq 0 ]; then
    echo "✅ Build successful!"
    echo ""
    echo "To run the application:"
    echo "  cd rust && cargo run --release"
    echo ""
    echo "Or run the binary directly:"
    echo "  ./rust/target/release/chemistry"
    echo ""
    echo "Server will be available at http://127.0.0.1:5000"
else
    echo "❌ Build failed. Please check the error messages above."
    exit 1
fi
