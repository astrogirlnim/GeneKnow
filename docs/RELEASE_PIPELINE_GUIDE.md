# 🚀 Release Pipeline Guide for GeneKnow

## Overview

GeneKnow uses an automated GitHub Actions release pipeline that builds the desktop application for all major platforms. The pipeline automatically versions, builds, and releases the app when changes are pushed to the main branch.

## 🎯 Key Features

- **Automatic Version Bumping**: No manual tagging required
- **Cross-Platform Builds**: Windows, macOS (Intel + Apple Silicon), Linux
- **Optimized Bundle Size**: 364MB (vs 1.9GB with TensorFlow)
- **Smart Platform Detection**: Single bundling script handles all platforms
- **Artifact Cleanup**: Aggressive storage management to stay within GitHub limits

## 📦 Build Outputs

Each platform produces multiple distribution formats:

| Platform | Formats | Description |
|----------|---------|-------------|
| **macOS Intel** | `.dmg`, `.app` | Universal installer + direct app |
| **macOS ARM64** | `.dmg`, `.app` | Native Apple Silicon build |
| **Windows** | `.msi`, `.exe` | MSI installer + NSIS installer |
| **Linux** | `.deb`, `.AppImage` | Debian package + portable AppImage |

## 🔧 How It Works

### 1. **Trigger**
- Push to `main` branch
- Manual workflow dispatch
- Changes to desktop/, docs/, or Python files

### 2. **Version Management**
```yaml
# Automatic version calculation
- Reads current version from tauri.conf.json
- Bumps patch/minor/major based on input
- Handles tag conflicts automatically
- Updates version in all config files
```

### 3. **Python Bundling**
The optimized bundling script (`bundle-python-optimized.sh`):
- Detects platform automatically
- Downloads appropriate Python 3.11.13 runtime
- Installs dependencies from `requirements-lite.txt` (no TensorFlow)
- Creates platform-specific startup scripts
- Reduces bundle size by ~85%

### 4. **Platform Detection**
```bash
# Automatic platform detection
case "$OS_TYPE-$ARCH_TYPE" in
    "Darwin-x86_64")     # macOS Intel
    "Darwin-arm64")      # macOS Apple Silicon
    "Linux-x86_64")      # Linux x64
    "MINGW*"|"MSYS*")    # Windows
esac
```

### 5. **Build Process**
- Frontend: React + TypeScript + Tailwind
- Backend: Rust + Tauri
- Python: Bundled runtime with genomic pipeline
- All platforms built in parallel

## 🛠️ Configuration Files

### `tauri.conf.json`
```json
{
  "bundle": {
    "active": true,
    "targets": "all",  // Builds all formats
    "resources": ["../bundled_resources/**/*"]
  }
}
```

### `requirements-lite.txt`
Optimized dependencies without TensorFlow:
- Flask + Flask-SocketIO for API
- BioPython for genomic processing
- NumPy/Pandas for data analysis
- SQLite for local database

## 📊 Bundle Size Optimization

| Component | Original | Optimized | Savings |
|-----------|----------|-----------|---------|
| TensorFlow | 1.1GB | 0MB | 100% |
| Python Runtime | 200MB | 150MB | 25% |
| Dependencies | 600MB | 200MB | 67% |
| **Total** | **1.9GB** | **364MB** | **81%** |

## 🚦 Release Workflow

1. **Cleanup** - Delete old artifacts (keeps only latest)
2. **Version** - Calculate and tag new version
3. **Validate** - Run linting and tests
4. **Build** - Create bundles for all platforms
5. **Release** - Upload to GitHub Releases

## 🐛 Troubleshooting

### Bundle Script Issues
```bash
# Test bundling locally
cd desktop
./scripts/bundle-python-optimized.sh

# Check bundle contents
ls -la bundled_resources/
```

### Platform-Specific Issues
- **Windows**: Ensure Git Bash is installed for CI
- **macOS**: Code signing may be required for distribution
- **Linux**: AppImage needs FUSE to run

### Version Conflicts
The pipeline automatically handles version conflicts by:
1. Checking for existing tags
2. Incrementing version if conflict found
3. Using timestamp suffix as last resort

## 🔒 Security Notes

- All processing happens locally on user's machine
- No external API calls or data transmission
- Python runtime is from official python-build-standalone
- Dependencies are pinned to specific versions

## 📝 Manual Release

If needed, you can trigger a release manually:
```bash
# From GitHub Actions tab
1. Go to "🚀 Release Pipeline"
2. Click "Run workflow"
3. Select version type (patch/minor/major)
4. Click "Run workflow"
```

## 🎯 Best Practices

1. **Always test locally first**: `pnpm run tauri-build`
2. **Check bundle size**: Should be ~350-400MB
3. **Verify all platforms**: Download and test each installer
4. **Monitor storage**: GitHub Actions has storage limits
5. **Update docs**: Keep this guide current with changes

## 📚 Related Documentation

- [Deployment Guide](DEPLOYMENT_GUIDE.md)
- [Dynamic Port Solution](DYNAMIC_PORT_SOLUTION.md)
- [Versioning Guidelines](VERSIONING_GUIDELINES.md)
- [Size Optimization Guide](SIZE_OPTIMIZATION_GUIDE.md) 