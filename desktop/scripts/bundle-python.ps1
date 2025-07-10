# bundle-python.ps1 - Bundle Python runtime and dependencies for GenePredict (Windows)
# This script prepares a complete Python environment for distribution with the Tauri app

$ErrorActionPreference = "Stop"

# Configuration
$PYTHON_VERSION = "3.11.9"
$SCRIPT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$PROJECT_ROOT = Split-Path -Parent (Split-Path -Parent $SCRIPT_DIR)
$DESKTOP_DIR = Join-Path $PROJECT_ROOT "desktop"
$BUNDLE_DIR = Join-Path $DESKTOP_DIR "bundled_resources"

Write-Host "🐍 GenePredict Python Bundling Script (Windows)" -ForegroundColor Cyan
Write-Host "=============================================" -ForegroundColor Cyan
Write-Host "Python Version: $PYTHON_VERSION"
Write-Host "Bundle Directory: $BUNDLE_DIR"
Write-Host ""

# Clean previous bundle
if (Test-Path $BUNDLE_DIR) {
    Write-Host "🧹 Cleaning previous bundle..." -ForegroundColor Yellow
    Remove-Item -Path $BUNDLE_DIR -Recurse -Force
}

New-Item -ItemType Directory -Path $BUNDLE_DIR -Force | Out-Null

# Download Python
Write-Host "📥 Downloading Python for Windows..." -ForegroundColor Green
$pythonUrl = "https://github.com/indygreg/python-build-standalone/releases/download/20241016/cpython-$PYTHON_VERSION+20241016-x86_64-pc-windows-msvc-install_only_stripped.tar.gz"
$pythonArchive = Join-Path $BUNDLE_DIR "python-windows.tar.gz"

Write-Host "   URL: $pythonUrl"
Invoke-WebRequest -Uri $pythonUrl -OutFile $pythonArchive

Write-Host "📦 Extracting Python..." -ForegroundColor Green
# Extract using tar (available on Windows 10+)
Push-Location $BUNDLE_DIR
tar -xzf $pythonArchive
Remove-Item $pythonArchive
Rename-Item "python" "python_runtime"
Pop-Location

$PYTHON_EXE = Join-Path $BUNDLE_DIR "python_runtime\python.exe"
$PIP_EXE = Join-Path $BUNDLE_DIR "python_runtime\Scripts\pip.exe"

Write-Host "✅ Python runtime downloaded" -ForegroundColor Green
Write-Host "   Python: $PYTHON_EXE"
Write-Host "   Pip: $PIP_EXE"

# Copy pipeline code
Write-Host ""
Write-Host "📂 Copying pipeline code..." -ForegroundColor Green
$pipelineSource = Join-Path $PROJECT_ROOT "geneknow_pipeline"
$pipelineTarget = Join-Path $BUNDLE_DIR "geneknow_pipeline"
Copy-Item -Path $pipelineSource -Destination $pipelineTarget -Recurse

# Remove unnecessary files
Get-ChildItem -Path $pipelineTarget -Include "*.pyc" -Recurse | Remove-Item
Get-ChildItem -Path $pipelineTarget -Include "__pycache__" -Directory -Recurse | Remove-Item -Recurse -Force
Get-ChildItem -Path $pipelineTarget -Include ".pytest_cache" -Directory -Recurse | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
Get-ChildItem -Path $pipelineTarget -Include "venv" -Directory -Recurse | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
Get-ChildItem -Path $pipelineTarget -Filter "test_*" | Remove-Item -Recurse

Write-Host "✅ Pipeline code copied" -ForegroundColor Green

# Install dependencies
Write-Host ""
Write-Host "📦 Installing Python dependencies..." -ForegroundColor Green
Push-Location $pipelineTarget

# Create production requirements file
$reqFile = "requirements.txt"
$prodReqFile = "requirements_production.txt"
Get-Content $reqFile | Where-Object { $_ -notmatch "(pytest|black|flake8|mypy)" } | Set-Content $prodReqFile

# Install dependencies
& $PIP_EXE install --no-cache-dir -r $prodReqFile

Pop-Location
Write-Host "✅ Dependencies installed" -ForegroundColor Green

# Create or copy database
Write-Host ""
Write-Host "🗄️ Setting up database..." -ForegroundColor Green
$sourceDb = Join-Path $pipelineSource "population_variants.db"
$targetDb = Join-Path $pipelineTarget "population_variants.db"

if (Test-Path $sourceDb) {
    Write-Host "   Found existing database, copying..."
    Copy-Item -Path $sourceDb -Destination $targetDb
} else {
    Write-Host "   No database found, will create on first run"
    # Create a marker file to trigger database creation on first run
    New-Item -Path (Join-Path $pipelineTarget ".needs_database_init") -ItemType File | Out-Null
}

# Create startup batch file (already created by bash script, but ensure it exists)
$startupBat = Join-Path $BUNDLE_DIR "start_api_server.bat"
if (-not (Test-Path $startupBat)) {
    @'
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
'@ | Set-Content $startupBat
}

# Create bundle manifest
Write-Host ""
Write-Host "📋 Creating bundle manifest..." -ForegroundColor Green
$manifest = @{
    bundle_version = "1.0.0"
    python_version = $PYTHON_VERSION
    platform = "windows-x86_64"
    created_at = (Get-Date -Format "yyyy-MM-ddTHH:mm:ssZ")
    components = @{
        python_runtime = $true
        geneknow_pipeline = $true
        database = (Test-Path $targetDb)
    }
}
$manifest | ConvertTo-Json | Set-Content (Join-Path $BUNDLE_DIR "manifest.json")

# Calculate bundle size
$bundleSize = "{0:N2} MB" -f ((Get-ChildItem $BUNDLE_DIR -Recurse | Measure-Object -Property Length -Sum).Sum / 1MB)

Write-Host ""
Write-Host "✅ Bundle created successfully!" -ForegroundColor Green
Write-Host "   Platform: Windows x64"
Write-Host "   Size: $bundleSize"
Write-Host "   Location: $BUNDLE_DIR"
Write-Host ""
Write-Host "📦 Bundle contents:" -ForegroundColor Green
Get-ChildItem $BUNDLE_DIR | Format-Table Name, Length, LastWriteTime
Write-Host ""
Write-Host "🎉 Ready for Tauri packaging!" -ForegroundColor Green 