#!/bin/bash
# test-bundle.sh - Test the bundled Python runtime and server

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLE_DIR="$SCRIPT_DIR/../bundled_resources"

echo "🧪 Testing Bundled Python Runtime"
echo "================================="

# Check if bundle exists
if [ ! -d "$BUNDLE_DIR" ]; then
    echo "❌ Bundle directory not found at: $BUNDLE_DIR"
    echo "   Please run ./bundle-python.sh first"
    exit 1
fi

# Test Python runtime
echo ""
echo "1️⃣ Testing Python executable..."
if [ -f "$BUNDLE_DIR/python_runtime/bin/python3" ]; then
    PYTHON_EXE="$BUNDLE_DIR/python_runtime/bin/python3"
elif [ -f "$BUNDLE_DIR/python_runtime/python.exe" ]; then
    PYTHON_EXE="$BUNDLE_DIR/python_runtime/python.exe"
else
    echo "❌ Python executable not found"
    exit 1
fi

echo "   Python path: $PYTHON_EXE"
$PYTHON_EXE --version

# Test Python packages
echo ""
echo "2️⃣ Testing installed packages..."
$PYTHON_EXE -c "import flask; print(f'   ✅ Flask {flask.__version__}')"
$PYTHON_EXE -c "import flask_socketio; print('   ✅ Flask-SocketIO installed')"
$PYTHON_EXE -c "import langgraph; print('   ✅ LangGraph installed')"
$PYTHON_EXE -c "import pandas; print(f'   ✅ Pandas {pandas.__version__}')"
$PYTHON_EXE -c "import numpy; print(f'   ✅ NumPy {numpy.__version__}')"

# Test pipeline code
echo ""
echo "3️⃣ Testing pipeline imports..."
cd "$BUNDLE_DIR/geneknow_pipeline"
$PYTHON_EXE -c "import graph; print('   ✅ Pipeline graph module loads')"
$PYTHON_EXE -c "import nodes; print('   ✅ Pipeline nodes module loads')"
$PYTHON_EXE -c "from nodes import cadd_scoring; print('   ✅ CADD scoring module loads')"

# Test database
echo ""
echo "4️⃣ Testing database..."
if [ -f "population_variants.db" ]; then
    echo "   ✅ Database exists"
    DB_SIZE=$(du -h population_variants.db | cut -f1)
    echo "   Size: $DB_SIZE"
elif [ -f ".needs_database_init" ]; then
    echo "   ⏳ Database will be created on first run"
else
    echo "   ⚠️  No database found, creating marker for first-run init"
    touch .needs_database_init
fi

# Test startup script
echo ""
echo "5️⃣ Testing startup script..."
if [ -f "$BUNDLE_DIR/start_api_server.sh" ]; then
    echo "   ✅ Unix startup script found"
    ls -la "$BUNDLE_DIR/start_api_server.sh"
fi
if [ -f "$BUNDLE_DIR/start_api_server.bat" ]; then
    echo "   ✅ Windows startup script found"
fi

# Test API server can import
echo ""
echo "6️⃣ Testing API server imports..."
$PYTHON_EXE -c "import enhanced_api_server; print('   ✅ API server module loads successfully')"

echo ""
echo "✅ All tests passed! Bundle is ready for deployment."
echo ""
echo "To test the server manually:"
echo "  cd $BUNDLE_DIR/geneknow_pipeline"
echo "  $PYTHON_EXE enhanced_api_server.py" 