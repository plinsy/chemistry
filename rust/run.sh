#!/bin/bash

clear

echo "🧪 Starting Chemistry Isomer Calculator (Rust)..."
echo ""
echo "Server starting at http://127.0.0.1:5000"
echo "Press Ctrl+C to stop"
echo ""

# Navigate to rust directory
cd "$(dirname "$0")"

cargo run --release
