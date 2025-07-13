#!/bin/bash
# bundle-python-optimized.sh - Optimized Python bundling for GeneKnow
# This script reduces app size from 1.9GB to ~400MB by using lightweight dependencies

set -e  # Exit on error

# Configuration
PYTHON_VERSION="3.11.13"
RELEASE_DATE="20250708"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DESKTOP_DIR="$PROJECT_ROOT/desktop"
BUNDLE_DIR="$DESKTOP_DIR/bundled_resources"

# Platform detection
OS_TYPE=$(uname -s)
ARCH_TYPE=$(uname -m)

echo "🚀 GeneKnow Optimized Python Bundling Script"
echo "=============================================="
echo "Platform: $OS_TYPE ($ARCH_TYPE)"
echo "Python Version: $PYTHON_VERSION"
echo "Bundle Directory: $BUNDLE_DIR"
echo "Target Size: ~300MB (vs 1.9GB original)"
echo ""

# Clean previous bundle
if [ -d "$BUNDLE_DIR" ]; then
    echo "🧹 Cleaning previous bundle..."
    rm -rf "$BUNDLE_DIR"
fi

mkdir -p "$BUNDLE_DIR"

# Determine platform and download URL
case "$OS_TYPE-$ARCH_TYPE" in
    "Darwin-x86_64")
        platform_key="darwin-x86_64"
        python_url="https://github.com/astral-sh/python-build-standalone/releases/download/$RELEASE_DATE/cpython-${PYTHON_VERSION}%2B${RELEASE_DATE}-x86_64-apple-darwin-install_only_stripped.tar.gz"
        ;;
    "Darwin-arm64")
        platform_key="darwin-aarch64"
        python_url="https://github.com/astral-sh/python-build-standalone/releases/download/$RELEASE_DATE/cpython-${PYTHON_VERSION}%2B${RELEASE_DATE}-aarch64-apple-darwin-install_only_stripped.tar.gz"
        ;;
    "Linux-x86_64")
        platform_key="linux-x86_64"
        python_url="https://github.com/astral-sh/python-build-standalone/releases/download/$RELEASE_DATE/cpython-${PYTHON_VERSION}%2B${RELEASE_DATE}-x86_64-unknown-linux-gnu-install_only_stripped.tar.gz"
        ;;
    "MINGW"*|"MSYS"*)
        platform_key="windows-x86_64"
        python_url="https://github.com/astral-sh/python-build-standalone/releases/download/$RELEASE_DATE/cpython-${PYTHON_VERSION}%2B${RELEASE_DATE}-x86_64-pc-windows-msvc-install_only_stripped.tar.gz"
        ;;
    *)
        echo "❌ Unsupported platform: $OS_TYPE-$ARCH_TYPE"
        exit 1
        ;;
esac

# Download Python runtime
echo "📥 Downloading Python runtime..."
echo "   URL: $python_url"
python_archive="python-${PYTHON_VERSION}-${platform_key}.tar.gz"

cd "$BUNDLE_DIR"
curl -L --fail -o "$python_archive" "$python_url" || {
    echo "❌ Failed to download Python runtime"
    exit 1
}

echo "📦 Extracting Python runtime..."
tar -xzf "$python_archive"
rm "$python_archive"

# Find Python runtime directory
if [ -d "python" ]; then
    mv python python_runtime
else
    echo "❌ Python runtime directory not found"
    exit 1
fi

# Get Python executable path
if [ -f "python_runtime/bin/python3" ]; then
    PYTHON_EXE="$BUNDLE_DIR/python_runtime/bin/python3"
elif [ -f "python_runtime/python.exe" ]; then
    PYTHON_EXE="$BUNDLE_DIR/python_runtime/python.exe"
else
    echo "❌ Python executable not found"
    exit 1
fi

echo "🐍 Python executable: $PYTHON_EXE"

# Test Python installation
echo "🧪 Testing Python installation..."
"$PYTHON_EXE" --version || {
    echo "❌ Python installation failed"
    exit 1
}

# Install lightweight dependencies
echo "📦 Installing optimized dependencies..."
if [[ "$platform_key" == "windows-x86_64" ]]; then
    echo "   Using requirements-lite-windows.txt for Windows (pysam excluded)"
    REQUIREMENTS_FILE="$PROJECT_ROOT/geneknow_pipeline/requirements-lite-windows.txt"
else
    echo "   Using requirements-lite.txt for ~85% size reduction (TensorFlow removed)"
    REQUIREMENTS_FILE="$PROJECT_ROOT/geneknow_pipeline/requirements-lite.txt"
fi

"$PYTHON_EXE" -m pip install --no-cache-dir -r "$REQUIREMENTS_FILE" || {
    echo "❌ Failed to install Python dependencies"
    exit 1
}

echo "🧹 Cleaning up Python installation..."
# Remove unnecessary files to reduce size
cd "$BUNDLE_DIR/python_runtime"

# Remove test files and documentation
find . -name "test" -type d -exec rm -rf {} + 2>/dev/null || true
find . -name "tests" -type d -exec rm -rf {} + 2>/dev/null || true
find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
find . -name "*.pyc" -delete 2>/dev/null || true
find . -name "*.pyo" -delete 2>/dev/null || true

# Remove development tools
rm -rf lib/python*/site-packages/pip/_internal/operations/build/ 2>/dev/null || true
rm -rf lib/python*/site-packages/setuptools/ 2>/dev/null || true
rm -rf lib/python*/site-packages/wheel/ 2>/dev/null || true

# Remove documentation and examples
find lib/python*/site-packages -name "docs" -type d -exec rm -rf {} + 2>/dev/null || true
find lib/python*/site-packages -name "examples" -type d -exec rm -rf {} + 2>/dev/null || true
find lib/python*/site-packages -name "*.md" -delete 2>/dev/null || true
find lib/python*/site-packages -name "*.rst" -delete 2>/dev/null || true

# Remove unused language packs (keep only English)
find lib/python*/site-packages -name "locale" -type d -exec sh -c '
    for dir in "$1"/*; do
        if [ -d "$dir" ] && [ "$(basename "$dir")" != "en" ] && [ "$(basename "$dir")" != "en_US" ]; then
            rm -rf "$dir"
        fi
    done
' _ {} \; 2>/dev/null || true

echo "📂 Copying pipeline code..."
cp -r "$PROJECT_ROOT/geneknow_pipeline" "$BUNDLE_DIR/"

# Remove unnecessary files from pipeline
cd "$BUNDLE_DIR/geneknow_pipeline"
rm -rf venv/ 2>/dev/null || true
rm -rf __pycache__/ 2>/dev/null || true
rm -rf .pytest_cache/ 2>/dev/null || true
rm -rf *.log 2>/dev/null || true
find . -name "*.pyc" -delete 2>/dev/null || true

echo "🗄️ Setting up database..."
if [ -f "$PROJECT_ROOT/geneknow_pipeline/tcga_data/population_variants.db" ]; then
    echo "   Found existing database, copying..."
    # Database is already included in the pipeline copy
else
    echo "   Creating minimal database..."
    cd "$BUNDLE_DIR/geneknow_pipeline"
    "$PYTHON_EXE" create_population_database.py
fi

echo "🚀 Creating optimized startup wrapper..."
cat > "$BUNDLE_DIR/start_api_server.sh" << 'EOF'
#!/bin/bash
# Optimized startup script for GeneKnow API server

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_EXE="$SCRIPT_DIR/python_runtime/bin/python3"

# Set environment variables for optimization
export PYTHONPATH="$SCRIPT_DIR:$SCRIPT_DIR/geneknow_pipeline:$PYTHONPATH"
export TF_CPP_MIN_LOG_LEVEL=2  # Reduce TensorFlow logging
export PYTHONOPTIMIZE=1        # Enable Python optimizations
export PYTHONDONTWRITEBYTECODE=1  # Don't create .pyc files

# Change to the bundled resources directory (parent of geneknow_pipeline)
cd "$SCRIPT_DIR"

# Start the API server using module syntax
echo "🚀 Starting GeneKnow API Server (Optimized)..."
exec "$PYTHON_EXE" -m geneknow_pipeline.enhanced_api_server "$@"
EOF

chmod +x "$BUNDLE_DIR/start_api_server.sh"

# Create Windows batch wrapper
echo "🚀 Creating Windows startup wrapper..."
cat > "$BUNDLE_DIR/start_api_server.bat" << 'EOF'
@echo off
setlocal

set SCRIPT_DIR=%~dp0
set PYTHON_EXE=%SCRIPT_DIR%python_runtime\python.exe

REM Set environment variables for optimization
set PYTHONPATH=%SCRIPT_DIR%;%SCRIPT_DIR%geneknow_pipeline;%PYTHONPATH%
set TF_CPP_MIN_LOG_LEVEL=2
set PYTHONOPTIMIZE=1
set PYTHONDONTWRITEBYTECODE=1

REM Change to the bundled resources directory
cd /d "%SCRIPT_DIR%"

REM Start the API server using module syntax
echo Starting GeneKnow API Server (Optimized)...
"%PYTHON_EXE%" -m geneknow_pipeline.enhanced_api_server %*
EOF

echo "📋 Creating bundle manifest..."
# Get final size
BUNDLE_SIZE=$(du -sh "$BUNDLE_DIR" | cut -f1)
DB_SIZE=$(du -sh "$BUNDLE_DIR/geneknow_pipeline/tcga_data/population_variants.db" 2>/dev/null | cut -f1 || echo "N/A")

cat > "$BUNDLE_DIR/manifest.json" << EOF
{
  "bundle_version": "1.0.0-optimized",
  "created_at": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "platform": "$platform_key",
  "python_version": "$PYTHON_VERSION",
  "bundle_size": "$BUNDLE_SIZE",
  "database_size": "$DB_SIZE",
  "optimizations": [
    "TensorFlow completely removed (not used in current codebase)",
    "Removed test files and documentation", 
    "Removed development tools",
    "Removed unused language packs",
    "Bytecode optimization enabled"
  ],
  "components": {
    "python_runtime": "python_runtime/",
    "pipeline_code": "geneknow_pipeline/",
    "database": "geneknow_pipeline/tcga_data/population_variants.db",
    "startup_script": "start_api_server.sh"
  }
}
EOF

echo ""
echo "✅ Optimized bundle created successfully!"
echo "   Platform: $platform_key"
echo "   Size: $BUNDLE_SIZE (vs 1.9GB original)"
echo "   Reduction: ~85% smaller (TensorFlow removed)"
echo "   Location: $BUNDLE_DIR"

echo ""
echo "📦 Bundle contents:"
ls -la "$BUNDLE_DIR"

echo ""
echo "🎉 Ready for Tauri packaging!"
echo "   Size reduction achieved by:"
echo "   • TensorFlow completely removed (1.1GB → 0MB)"
echo "   • Removed test files and documentation"
echo "   • Removed development tools"
echo "   • Optimized Python bytecode"
echo ""
echo "💡 Next steps:"
echo "   1. Test with: ./test-bundle.sh"
echo "   2. Build app: cd ui && pnpm run tauri-build" 