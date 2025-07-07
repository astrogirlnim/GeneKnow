# 🧬 GenePredict - AI-Powered Genomic Risk Assessment

> **Local-First Genomic Analysis Platform** • Built with Tauri, React, and Rust for maximum privacy and security

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Platform](https://img.shields.io/badge/platform-macOS%20%7C%20Windows%20%7C%20Linux-lightgrey.svg)
![Version](https://img.shields.io/badge/version-0.1.0-green.svg)

---

## 🎯 Project Overview

GenePredict is a **privacy-first genomic risk assessment platform** that processes genetic data entirely on your local machine. No data ever leaves your device, ensuring complete privacy and HIPAA compliance for sensitive genetic information.

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
- **Landing Page**: Beautiful gradient design with GenePredict branding
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

---

## 🚀 Quick Start Guide

### Prerequisites
- **Node.js 20+** (LTS recommended)
- **Rust 1.77+** (stable toolchain)
- **pnpm 8+** (package manager)

### 1. Clone & Setup
```bash
git clone <repository-url>
cd LiteratureGapper
```

### 2. Install Dependencies
```bash
# Navigate to UI directory
cd desktop/ui

# Install JavaScript dependencies
pnpm install
```

### 3. Run Development Server
```bash
# Option 1: Run from UI directory
pnpm run tauri-dev

# Option 2: Run from desktop directory  
cd ../
cargo tauri dev
```

### 4. Build for Production
```bash
# From desktop/ directory
cargo tauri build

```

---

## 🏗️ Project Architecture

```
LiteratureGapper/
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
├── docs/                      # Project planning & documentation
├── documentation/             # Technical specifications
└── README.md                  # This file
```

---

## 🔧 Development

### Available Scripts (from `desktop/ui/`)
```bash
pnpm run dev          # Start Vite development server
pnpm run build        # Build for production
pnpm run tauri-dev    # Start Tauri + React development
pnpm run tauri-build  # Build desktop application
pnpm run lint         # Run ESLint
```

### Environment Variables
- `RUST_LOG=debug` - Enable detailed Rust logging
- `NODE_ENV=development` - React development mode
- `TAURI_DEV_URL` - Custom development server URL

### Debugging
- **React DevTools**: Available in development mode
- **Rust Logging**: Use `RUST_LOG=debug` for detailed backend logs
- **Browser Console**: Frontend logs via `useLogger` hook
- **Hot Reload**: Automatic for both React and Rust changes

---

## 🛣️ Next Steps (Roadmap)

### Phase 2: Data Layer (In Progress)
- [ ] Python ML service integration via Tauri plugin
- [ ] FASTQ/BAM file parsing with BioPython
- [ ] TensorFlow pipeline for breast cancer prediction
- [ ] DeepVariant integration for variant calling

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
- [`documentation/GenePredict_BrainLift/`](documentation/GenePredict_BrainLift/) - Technical specifications
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