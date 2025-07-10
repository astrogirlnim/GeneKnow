#!/bin/bash
# bundle-python.sh - Bundle Python runtime and dependencies for GenePredict
# This script prepares a complete Python environment for distribution with the Tauri app

set -e  # Exit on error

# Configuration
PYTHON_VERSION="3.11.9"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DESKTOP_DIR="$PROJECT_ROOT/desktop"
BUNDLE_DIR="$DESKTOP_DIR/bundled_resources"

# Platform detection
OS_TYPE=$(uname -s)
ARCH_TYPE=$(uname -m)

echo "🐍 GenePredict Python Bundling Script"
echo "===================================="
echo "Platform: $OS_TYPE ($ARCH_TYPE)"
echo "Python Version: $PYTHON_VERSION"
echo "Bundle Directory: $BUNDLE_DIR"
echo ""

# Clean previous bundle
if [ -d "$BUNDLE_DIR" ]; then
    echo "🧹 Cleaning previous bundle..."
    rm -rf "$BUNDLE_DIR"
fi

mkdir -p "$BUNDLE_DIR"

# Function to download Python build
download_python() {
    local platform=$1
    local python_url=""
    local python_archive="python-$PYTHON_VERSION-$platform.tar.gz"
    
    case $platform in
        "darwin-x86_64")
            python_url="https://github.com/indygreg/python-build-standalone/releases/download/20241016/cpython-${PYTHON_VERSION}+20241016-x86_64-apple-darwin-install_only_stripped.tar.gz"
            ;;
        "darwin-aarch64")
            python_url="https://github.com/indygreg/python-build-standalone/releases/download/20241016/cpython-${PYTHON_VERSION}+20241016-aarch64-apple-darwin-install_only_stripped.tar.gz"
            ;;
        "linux-x86_64")
            python_url="https://github.com/indygreg/python-build-standalone/releases/download/20241016/cpython-${PYTHON_VERSION}+20241016-x86_64-unknown-linux-gnu-install_only_stripped.tar.gz"
            ;;
        "windows-x86_64")
            python_url="https://github.com/indygreg/python-build-standalone/releases/download/20241016/cpython-${PYTHON_VERSION}+20241016-x86_64-pc-windows-msvc-install_only_stripped.tar.gz"
            ;;
        *)
            echo "❌ Unsupported platform: $platform"
            exit 1
            ;;
    esac
    
    echo "📥 Downloading Python for $platform..."
    echo "   URL: $python_url"
    
    cd "$BUNDLE_DIR"
    curl -L -o "$python_archive" "$python_url"
    
    echo "📦 Extracting Python..."
    tar -xzf "$python_archive"
    rm "$python_archive"
    
    # Rename to consistent directory
    mv python "$BUNDLE_DIR/python_runtime"
}

# Determine platform
PLATFORM=""
case "$OS_TYPE" in
    Darwin)
        if [ "$ARCH_TYPE" = "arm64" ]; then
            PLATFORM="darwin-aarch64"
        else
            PLATFORM="darwin-x86_64"
        fi
        ;;
    Linux)
        PLATFORM="linux-x86_64"
        ;;
    MINGW*|MSYS*|CYGWIN*)
        PLATFORM="windows-x86_64"
        ;;
    *)
        echo "❌ Unsupported OS: $OS_TYPE"
        exit 1
        ;;
esac

# Download Python
download_python "$PLATFORM"

# Set up Python executable path
if [ "$PLATFORM" = "windows-x86_64" ]; then
    PYTHON_EXE="$BUNDLE_DIR/python_runtime/python.exe"
    PIP_EXE="$BUNDLE_DIR/python_runtime/Scripts/pip.exe"
else
    PYTHON_EXE="$BUNDLE_DIR/python_runtime/bin/python3"
    PIP_EXE="$BUNDLE_DIR/python_runtime/bin/pip3"
fi

echo "✅ Python runtime downloaded"
echo "   Python: $PYTHON_EXE"
echo "   Pip: $PIP_EXE"

# Copy pipeline code
echo ""
echo "📂 Copying pipeline code..."
mkdir -p "$BUNDLE_DIR/geneknow_pipeline"
cp -r "$PROJECT_ROOT/geneknow_pipeline"/* "$BUNDLE_DIR/geneknow_pipeline/"

# Remove unnecessary files
find "$BUNDLE_DIR/geneknow_pipeline" -name "*.pyc" -delete
find "$BUNDLE_DIR/geneknow_pipeline" -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
find "$BUNDLE_DIR/geneknow_pipeline" -name ".pytest_cache" -type d -exec rm -rf {} + 2>/dev/null || true
find "$BUNDLE_DIR/geneknow_pipeline" -name "venv" -type d -exec rm -rf {} + 2>/dev/null || true
rm -rf "$BUNDLE_DIR/geneknow_pipeline/test_*"

echo "✅ Pipeline code copied"

# Install dependencies
echo ""
echo "📦 Installing Python dependencies..."
cd "$BUNDLE_DIR/geneknow_pipeline"

# Create a temporary requirements file without development dependencies
grep -v -E "(pytest|black|flake8|mypy)" requirements.txt > requirements_production.txt || cp requirements.txt requirements_production.txt

# Install dependencies into the bundled Python
"$PIP_EXE" install --no-cache-dir -r requirements_production.txt

echo "✅ Dependencies installed"

# Create or copy database
echo ""
echo "🗄️ Setting up database..."
if [ -f "$PROJECT_ROOT/geneknow_pipeline/population_variants.db" ]; then
    echo "   Found existing database, copying..."
    cp "$PROJECT_ROOT/geneknow_pipeline/population_variants.db" "$BUNDLE_DIR/geneknow_pipeline/"
else
    echo "   No database found, will create on first run"
    # Create a marker file to trigger database creation on first run
    touch "$BUNDLE_DIR/geneknow_pipeline/.needs_database_init"
fi

# Create startup wrapper script
echo ""
echo "🚀 Creating startup wrapper..."
cat > "$BUNDLE_DIR/start_api_server.sh" << 'EOF'
#!/bin/bash
# Startup wrapper for GeneKnow API server

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_RUNTIME="$SCRIPT_DIR/python_runtime"
PIPELINE_DIR="$SCRIPT_DIR/geneknow_pipeline"

# Set Python path
if [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "win32" ]]; then
    PYTHON_EXE="$PYTHON_RUNTIME/python.exe"
else
    PYTHON_EXE="$PYTHON_RUNTIME/bin/python3"
fi

# Check if database initialization is needed
if [ -f "$PIPELINE_DIR/.needs_database_init" ]; then
    echo "🗄️ Initializing database on first run..."
    cd "$PIPELINE_DIR"
    "$PYTHON_EXE" create_population_database.py --cancer-genes-only
    rm -f "$PIPELINE_DIR/.needs_database_init"
fi

# Start the API server
cd "$PIPELINE_DIR"
exec "$PYTHON_EXE" enhanced_api_server.py
EOF

chmod +x "$BUNDLE_DIR/start_api_server.sh"

# Create Windows batch wrapper
cat > "$BUNDLE_DIR/start_api_server.bat" << 'EOF'
@echo off
setlocal

set SCRIPT_DIR=%~dp0
set PYTHON_RUNTIME=%SCRIPT_DIR%python_runtime
set PIPELINE_DIR=%SCRIPT_DIR%geneknow_pipeline
set PYTHON_EXE=%PYTHON_RUNTIME%\python.exe

REM Check if database initialization is needed
if exist "%PIPELINE_DIR%\.needs_database_init" (
    echo Initializing database on first run...
    cd /d "%PIPELINE_DIR%"
    "%PYTHON_EXE%" create_population_database.py --cancer-genes-only
    del "%PIPELINE_DIR%\.needs_database_init"
)

REM Start the API server
cd /d "%PIPELINE_DIR%"
"%PYTHON_EXE%" enhanced_api_server.py
EOF

# Create bundle manifest
echo ""
echo "📋 Creating bundle manifest..."
cat > "$BUNDLE_DIR/manifest.json" << EOF
{
    "bundle_version": "1.0.0",
    "python_version": "$PYTHON_VERSION",
    "platform": "$PLATFORM",
    "created_at": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
    "components": {
        "python_runtime": true,
        "geneknow_pipeline": true,
        "database": $([ -f "$BUNDLE_DIR/geneknow_pipeline/population_variants.db" ] && echo "true" || echo "false")
    }
}
EOF

# Calculate bundle size
BUNDLE_SIZE=$(du -sh "$BUNDLE_DIR" | cut -f1)

echo ""
echo "✅ Bundle created successfully!"
echo "   Platform: $PLATFORM"
echo "   Size: $BUNDLE_SIZE"
echo "   Location: $BUNDLE_DIR"
echo ""
echo "📦 Bundle contents:"
ls -la "$BUNDLE_DIR/"
echo ""
echo "🎉 Ready for Tauri packaging!" 