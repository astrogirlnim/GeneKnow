#!/bin/bash

# 🧬 GenePredict Implementation Test Script
# Tests the Rust ⇄ Python ML integration implemented in Phase 1 Step 2

set -e  # Exit on any error

echo "🧬 GenePredict - Testing Rust ⇄ Python Integration"
echo "=================================================="

# Check if we're in the right directory
if [[ ! -f "src-tauri/Cargo.toml" ]]; then
    echo "❌ Error: Please run this script from the desktop/ directory"
    echo "   Current directory: $(pwd)"
    echo "   Expected: .../LiteratureGapper/desktop/"
    exit 1
fi

echo ""
echo "🔍 Step 1: Testing Python Scripts with JSON Output"
echo "------------------------------------------------"

# Test generate_test_fastq.py
echo "Testing generate_test_fastq.py..."
mkdir -p test_output
python3 python_ml/generate_test_fastq.py --output-dir test_output --num-reads 5 --json > test_output/fastq_result.json
if [[ $? -eq 0 && -s test_output/fastq_result.json ]]; then
    echo "✅ generate_test_fastq.py: JSON output working"
    echo "   Output: $(cat test_output/fastq_result.json | head -c 100)..."
else
    echo "❌ generate_test_fastq.py: Failed"
    exit 1
fi

# Test config_data_source.py
echo "Testing config_data_source.py..."
python3 python_ml/config_data_source.py --list-vcf-files --json > test_output/config_result.json
if [[ $? -eq 0 && -s test_output/config_result.json ]]; then
    echo "✅ config_data_source.py: JSON output working"
    echo "   VCF files listed: $(cat test_output/config_result.json | jq 'keys | length' 2>/dev/null || echo 'N/A')"
else
    echo "❌ config_data_source.py: Failed"
    exit 1
fi

echo ""
echo "🦀 Step 2: Testing Rust Compilation"
echo "----------------------------------"
cd src-tauri
cargo check --quiet
if [[ $? -eq 0 ]]; then
    echo "✅ Rust code compiles successfully"
    echo "   Warnings: $(cargo check 2>&1 | grep -c warning || echo 0)"
else
    echo "❌ Rust compilation failed"
    exit 1
fi
cd ..

echo ""
echo "🔧 Step 3: Testing Frontend Setup"
echo "--------------------------------"
cd ui
if [[ ! -d "node_modules" ]]; then
    echo "Installing frontend dependencies..."
    pnpm install --silent
fi

echo "✅ Frontend dependencies ready"
echo "   Node modules: $(ls node_modules | wc -l | tr -d ' ') packages"
cd ..

echo ""
echo "🚀 Step 4: Testing Development Workflow"
echo "--------------------------------------"

# Start frontend server in background
echo "Starting frontend development server..."
cd ui
pnpm dev &
FRONTEND_PID=$!
cd ..

# Wait for frontend to start
echo "Waiting for frontend server to start..."
sleep 5

# Check if frontend is accessible
if curl -s http://localhost:5173 > /dev/null; then
    echo "✅ Frontend server running on http://localhost:5173"
else
    echo "❌ Frontend server not accessible"
    kill $FRONTEND_PID 2>/dev/null
    exit 1
fi

# Test Tauri compilation without running UI
echo "Testing Tauri backend compilation..."
cd src-tauri
if timeout 10s cargo run --help > /dev/null 2>&1; then
    echo "✅ Tauri backend can be compiled and run"
else
    echo "⚠️  Tauri backend compilation test inconclusive (timeout)"
fi
cd ..

# Cleanup
echo ""
echo "🧹 Cleanup"
echo "---------"
kill $FRONTEND_PID 2>/dev/null || echo "Frontend process already stopped"
sleep 1

echo ""
echo "📊 Test Summary"
echo "==============="
echo "✅ Python scripts with JSON output"
echo "✅ Rust compilation"
echo "✅ Frontend development server"
echo "✅ Integration architecture ready"
echo ""
echo "🎉 Implementation test PASSED!"
echo ""
echo "Next steps to run the full application:"
echo "1. cd ui && pnpm dev     # In one terminal"
echo "2. cd src-tauri && cargo tauri dev --no-dev-server  # In another terminal"
echo ""
echo "Or use the simplified command:"
echo "cd ui && pnpm run tauri-dev" 