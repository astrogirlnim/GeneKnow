# 🧬 GeneKnow - AI-Powered Genomic Risk Assessment

> **Local-First Genomic Analysis Platform** • Built with Tauri, React, and Rust for maximum privacy and security

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Platform](https://img.shields.io/badge/platform-macOS%20%7C%20Windows%20%7C%20Linux-lightgrey.svg)
![Version](https://img.shields.io/badge/version-0.1.0-green.svg)

---

## 🎯 Project Overview

GeneKnow is a **privacy-first genomic risk assessment platform** that processes genetic data entirely on your local machine. No data ever leaves your device, ensuring complete privacy and HIPAA compliance for sensitive genetic information.

### 🔬 Core Capabilities (Planned)
- **FASTQ/BAM File Processing**: Native support for genomic data formats
- **AI Risk Scoring**: TensorFlow-powered breast cancer risk assessment  
- **Variant Analysis**: Deep variant calling and interpretation
- **Multi-language Reports**: AI-generated summaries in English, Hindi, Spanish
- **Privacy-Preserving ML**: OpenMined PySyft integration for federated learning
- **Cross-Platform**: Secure desktop app for macOS, Windows, and Linux

---

## 📁 Current Application State

### ✅ **Phase 1: Foundation - COMPLETED**
- **Tauri Environment**: Cross-platform desktop framework initialized
- **React + TypeScript**: Modern UI with Vite build system  
- **Tailwind CSS**: Utility-first styling with production configuration
- **Rust Backend**: Secure native layer for file processing
- **Development Toolchain**: Hot reload, debugging, and build scripts
- **Logging Infrastructure**: Comprehensive logging with `useLogger` hook

### 🎨 **UI Status**
- **Landing Page**: Beautiful gradient design with GeneKnow branding
- **Interactive Components**: Sample analysis counter with state management
- **Responsive Design**: Mobile-first Tailwind implementation
- **Developer Experience**: Hot reload working for both React and Rust

### 🛠️ **Technical Stack**
```
Frontend:  React 19 + TypeScript + Tailwind CSS 4.1
Backend:   Rust 1.88 + Tauri 2.x
Build:     Vite 7.0 + pnpm 10.12
Platform:  Cross-platform desktop (macOS/Windows/Linux)
```

### 🔖 **Version Management**
GeneKnow uses automated semantic versioning with GitHub Actions:
- **Automated Releases**: Push to `main` triggers automatic version bumping and release builds
- **Cross-Platform Builds**: macOS, Windows, and Linux binaries generated automatically
- **No Manual Tags**: The release pipeline handles all version management to prevent conflicts

**For Developers:**
```bash
# Let the pipeline handle versioning (recommended)
git push origin main

# Manual version bump during development only
./scripts/bump-version.sh patch  # 0.1.0 → 0.1.1
./scripts/bump-version.sh minor  # 0.1.0 → 0.2.0
./scripts/bump-version.sh major  # 0.1.0 → 1.0.0

# Install git hooks to prevent conflicts
./scripts/install-git-hooks.sh
```

📖 **Full Guidelines**: See [`docs/VERSIONING_GUIDELINES.md`](docs/VERSIONING_GUIDELINES.md) for complete versioning workflow and troubleshooting.

---

## 🚀 Quick Start Guide

### Prerequisites
- **Node.js 20+** (LTS recommended)
- **Rust 1.77+** (stable toolchain)
- **pnpm 8+** (package manager)
- **Python 3.8+** (for genomic processing)

### 1. Clone & Setup
```bash
git clone https://github.com/astrogirlnim/GeneKnow.git
cd GeneKnow
```

### 2. Install Dependencies
```bash
# Navigate to UI directory (IMPORTANT: Must be in desktop/ui/)
cd desktop/ui

# Install JavaScript dependencies
pnpm install
```


### 3. Test the Implementation
```bash
# Navigate to desktop directory
cd ..

# Run comprehensive test (recommended first time)
./test_implementation.sh
```

---

## 📦 Installing Pre-built Binaries

### 🍎 macOS Installation (Recommended Method)

**Important for macOS users:** As a small open-source project, our macOS application is not yet code-signed by Apple. You'll need to bypass macOS security (Gatekeeper) to install and run GeneKnow.

**⚠️ Critical Installation Note:** Due to the bundled Python runtime and ML models, GeneKnow requires **programmatic installation** via Terminal to ensure all components are properly copied. **DO NOT** use drag-and-drop installation from the DMG as this will result in a non-functional app.

#### Programmatic Installation Steps:

1. **Download** the `.dmg` file from our [releases page](https://github.com/astrogirlnim/GeneKnow/releases)
2. **Install via Terminal** (copy and paste the commands below):

   ```bash
   # Step 1: Stop any running GeneKnow processes
   pkill -f "GeneKnow.app"
   pkill -f "geneknow_pipeline"
   
   # Step 2: Mount the DMG (replace with your actual DMG filename)
   hdiutil attach ~/Downloads/GeneKnow_*.dmg
   
   # Step 3: Remove old installation if it exists
   sudo rm -rf /Applications/GeneKnow.app
   
   # Step 4: Copy the app to Applications (preserves bundled resources)
   sudo cp -R "/Volumes/GeneKnow*/GeneKnow.app" /Applications/
   
   # Step 5: Remove quarantine attributes (required for unsigned apps)
   sudo xattr -r -d com.apple.quarantine /Applications/GeneKnow.app
   
   # Step 6: Unmount the DMG
   hdiutil detach "/Volumes/GeneKnow*"
   
   # Step 7: Launch the app
   open /Applications/GeneKnow.app
   ```

3. **First Launch**: The app should open automatically. If macOS shows a security warning, click "Open" or go to System Preferences → Security & Privacy → General → "Open Anyway"

#### ✅ Verification Steps:

After installation, verify that the app is working correctly:
```bash
# Check that bundled resources are present
ls -la /Applications/GeneKnow.app/Contents/Resources/_up_/bundled_resources/

# You should see:
# - start_api_server.sh
# - python_runtime/ (directory)
# - geneknow_pipeline/ (directory with .pkl and .db files)
```

If the `bundled_resources` directory is missing or empty, the drag-and-drop installation was used instead of the programmatic method. Please remove the app and reinstall using the Terminal commands above.

#### 🚨 Common Installation Issues:

- **"Command not found" errors**: Make sure you're running the Terminal commands, not drag-and-drop
- **App won't start**: Check that bundled resources are present using the verification steps above
- **"Not Available" in confidence check**: This indicates missing ML models - reinstall using the programmatic method
- **Empty Resources directory**: This confirms drag-and-drop was used - remove the app and reinstall via Terminal

#### Why This Happens
- Apple requires developer certificates ($99/year) for automatic installation
- As a small open-source project, we haven't yet obtained code signing certificates  
- The bundled Python runtime and ML models require proper file permissions and structure preservation
- We're working toward official code signing as our project grows

#### Security Assurance
- ✅ **Open Source**: Our entire codebase is [publicly auditable](https://github.com/astrogirlnim/GeneKnow)
- ✅ **Privacy-First**: Your genetic data never leaves your device
- ✅ **No Network Access**: All processing happens locally
- ✅ **Community Verified**: Code changes are reviewed through public pull requests

### 🪟 Windows Installation

Download the `.msi` or `.exe` installer from our [releases page](https://github.com/astrogirlnim/GeneKnow/releases) and run it normally. No additional security steps are required.

### 🐧 Linux Installation

Download the `.deb` package from our [releases page](https://github.com/astrogirlnim/GeneKnow/releases) and install it using:
```bash
sudo dpkg -i GeneKnow_*.deb
```

No additional security steps are required.

---

### 4. Run Full Application (Complete Integration)

**🚀 Recommended: One-Command Startup**
```bash
# From root directory, navigate to desktop/ui (starts both frontend + backend)
cd desktop/ui
pnpm run tauri-dev
```

This single command:
- ✅ Starts the React frontend (Vite dev server)
- ✅ Starts the Rust backend (Tauri + Python integration)  
- ✅ Opens the desktop application automatically
- ✅ Enables hot reload for both frontend and backend

**🌐 Access Your Application:**
- **Desktop App**: Opens automatically in a native window
- **Web Browser**: http://localhost:5173 (for debugging)

**✅ Verify Integration:**
```bash
# Test the Rust-Python integration (in another terminal)
cd desktop
python3 python_ml/config_data_source.py --list-vcf-files --json
```

---

### Alternative Development Options

**Option B: Manual (Two Terminals)**
```bash
# Terminal 1: Frontend (from desktop/ui/)
cd desktop/ui
pnpm dev

# Terminal 2: Backend (from desktop/src-tauri/)
cd desktop/src-tauri
cargo tauri dev --no-dev-server
```

**Option C: Separate Components**
```bash
# Frontend only (from desktop/ui/)
pnpm dev  # Runs on http://localhost:5173

# Backend only (from desktop/src-tauri/)
cargo run  # Runs Rust backend without UI
```

### 5. Build for Production
```bash
# From desktop/ui directory
pnpm run tauri-build
```

### 🧬 Test Genomic Processing
```bash
# From desktop/ directory
python3 python_ml/generate_test_fastq.py --output-dir test_output --num-reads 10 --json
python3 python_ml/config_data_source.py --list-vcf-files --json
```

---

## 🧬 Genomic Parallel Extraction Pipeline

A production-ready pipeline for extracting specific genomic regions from massive VCF files using parallel processing. This component handles the data processing layer of GeneKnow.

### Quick Start (Team Setup)

**1. First time setup (downloads data once per machine):**
```bash
python3 setup_data.py
```

**2. Extract your regions:**
```bash
python3 extract_by_region_precise.py
```

**That's it!** Large files are shared across team members.

### What This Pipeline Does

- ✅ **Parallel processing** across multiple CPU cores
- ✅ **100% precise extraction** (exact boundaries)
- ✅ **Multi-chromosome support** (handles mixed datasets)
- ✅ **Massive file handling** (tested on 11GB+ files)
- ✅ **Team-friendly** (shared data, no repeated downloads)

### Data Management

**Smart caching system:**
- Downloads large files **once per machine**
- Stores in shared locations (`/shared/genomics-data` → `~/genomics-data` → `./data-cache`)
- Creates symlinks in `test-data/` for project use
- Automatically handles indexing

**Git workflow:**
- ✅ Track: Scripts, regions, documentation
- ❌ Ignore: Large VCF files, outputs, indexes

### Supported File Types

| Dataset | Size | Type | Use Case |
|---------|------|------|----------|
| **Phase 1** | 11GB+ | SNPs + Indels + SVs | Comprehensive variant analysis |
| **Phase 3** | 200MB-1.2GB | High-quality SNPs | Population genetics |

### Usage Examples

**Define your regions in `regions.bed`:**
```
22	29121014	29235591	EWSR1
1	35691274	35801992	TP73
```

**Run extraction:**
```bash
python3 extract_by_region_precise.py
# Output: ✅ Extracted EWSR1 from chr22: 3260 variants (100% precise)
```

**Analyze results:**
```bash
ls -lh output_chunks/
# EWSR1.vcf.gz    598K
# TP73.vcf.gz     4.6M
```

### Team Collaboration

**New team member workflow:**
1. `git clone` the repository
2. `python3 setup_data.py` (uses existing shared data if available)
3. `python3 extract_by_region_precise.py` (ready to go!)

**No need to:**
- Download 11GB+ files individually
- Manage complex data paths
- Worry about incomplete downloads (resume capability)

### Performance

**Tested on:**
- ✅ 11GB Phase 1 files (chr1)
- ✅ 3,260 variants extracted in seconds
- ✅ 100% boundary precision
- ✅ Parallel processing across chromosomes

### Requirements

```bash
# System tools
bcftools
bgzip
wget

# Python (standard library only)
python3
```

### Troubleshooting

**"No such file" errors:** Run `python3 setup_data.py` first

**Slow downloads:** Script uses `wget -c` for resume capability

**Permission errors:** Check write access to data directories

### Architecture

```
Large Files (11GB+)     Shared Storage        Project Directory
├── chr1.vcf.gz    →    ~/genomics-data   →   test-data/ (symlinks)
├── chr22.vcf.gz        (cached once)         ├── Scripts
└── indexes                                   └── Output chunks
```

---

## 🏗️ Project Architecture

```
GeneKnow/
├── desktop/                    # Cross-platform desktop application
│   ├── src-tauri/             # Rust backend (Tauri framework)
│   │   ├── src/
│   │   │   ├── main.rs        # Rust entry point
│   │   │   └── lib.rs         # Core business logic
│   │   ├── tauri.conf.json    # Tauri configuration
│   │   └── Cargo.toml         # Rust dependencies
│   └── ui/                    # React frontend
│       ├── src/
│       │   ├── main.tsx       # React entry point
│       │   ├── App.tsx        # Main application component
│       │   └── hooks/
│       │       └── useLogger.ts # Logging utilities
│       ├── package.json       # Node.js dependencies
│       ├── tailwind.config.ts # Tailwind CSS configuration
│       └── vite.config.ts     # Vite build configuration
├── data_processing/           # Data processing pipelines
│   ├── tcga_download/         # TCGA data download utilities
│   └── testing_data_generator/ # Test data generation
├── extract_by_region_precise.py # Genomic region extraction
├── setup_data.py              # Data setup script
├── docs/                      # Project planning & documentation
├── documentation/             # Technical specifications
└── README.md                  # This file
```

---

## 🔧 Development

### Available Scripts (from `desktop/ui/`)
```bash
pnpm run dev          # Start Vite development server only
pnpm run build        # Build React frontend for production
pnpm run tauri-dev    # Start both React + Tauri (full app)
pnpm run tauri-build  # Build complete desktop application
pnpm run lint         # Run ESLint
```

### Python ML Scripts (from `desktop/`)
```bash
python3 python_ml/generate_test_fastq.py --help
python3 python_ml/config_data_source.py --help
python3 python_ml/extract_by_region.py --help
python3 python_ml/fastq_to_vcf_pipeline.py --help
```

### Environment Variables
- `RUST_LOG=debug` - Enable detailed Rust logging
- `NODE_ENV=development` - React development mode

### Debugging
- **React DevTools**: Available in development mode
- **Rust Logging**: Use `RUST_LOG=debug` for detailed backend logs  
- **Browser Console**: Frontend logs via `useLogger` hook
- **Hot Reload**: Automatic for both React and Rust changes
- **Python Integration**: All scripts support `--json` flag for structured output

### Troubleshooting

**"Command not found: tauri"**
```bash
cd desktop/ui
pnpm install  # Reinstall dependencies including @tauri-apps/cli
```

**"tauri-dev fails to start"**
```bash
# Ensure you're in the correct directory
cd desktop/ui
pwd  # Should show .../GeneKnow/desktop/ui

# Check if dependencies are installed
ls node_modules/.bin/tauri  # Should exist
```

**"Port 5173 is in use"**
```bash
# Kill existing processes
pkill -f "vite|cargo|tauri"
# Or change port in vite.config.ts
```

**"No such file or directory"**
```bash
# Ensure you're in the correct directory:
pwd  # Should show .../GeneKnow/desktop/ui for UI commands
```

**Python script errors**
```bash
# Test Python scripts individually:
cd desktop
python3 python_ml/config_data_source.py --help
```

---

## 🛣️ Next Steps (Roadmap)

### Phase 2: Data Layer (In Progress)
- [ ] Python ML service integration via Tauri plugin
- [ ] FASTQ/BAM file parsing with BioPython
- [ ] TensorFlow pipeline for breast cancer prediction
- [ ] DeepVariant integration for variant calling
- [x] Parallel genomic extraction pipeline

### Phase 3: Interface Layer
- [ ] Drag-and-drop file uploader
- [ ] Variant table with filtering and search
- [ ] Risk visualization and heatmaps
- [ ] Loading states and progress indicators

### Phase 4: Reporting & Compliance
- [ ] Llama 3.1 integration for AI report generation
- [ ] Multi-language support (English, Hindi, Spanish)
- [ ] PDF export functionality
- [ ] Privacy-preserving ML with OpenMined PySyft

### Phase 5: Explorer Mode
- [ ] Interactive variant simulation
- [ ] BRCA1/2 impact visualization
- [ ] Comparative analysis tools

---

## 🔒 Privacy & Security

- **Local Processing**: All data remains on your device
- **No External APIs**: Zero external data transmission
- **Secure File Handling**: Rust-powered native file operations
- **HIPAA Compliant**: Designed for medical data privacy
- **Open Source**: Transparent and auditable codebase

---

## 📚 Documentation

- [`docs/Phase1_Foundation_Tauri_Setup_Plan.md`](docs/Phase1_Foundation_Tauri_Setup_Plan.md) - Detailed setup plan
- [`documentation/GenePredict_BrainLift/`](documentation/GenePredict_BrainLift/) - Technical specifications (legacy)
- [`PRD_V2.md`](PRD_V2.md) - Product Requirements Document

---

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Run tests and linting
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🆘 Support

- **Issues**: [GitHub Issues](../../issues)
- **Discussions**: [GitHub Discussions](../../discussions)
- **Documentation**: See [`docs/`](docs/) directory

---

*"In genomics, privacy is not a feature—it's a fundamental right."*

---

**Last Updated**: January 2025 • **Status**: MVP Foundation Complete ✅
