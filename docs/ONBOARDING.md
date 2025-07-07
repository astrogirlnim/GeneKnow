# GenePredict Developer Onboarding Guide

Welcome to GenePredict! This guide will help you get up and running with the project quickly.

## 🧬 What is GenePredict?

GenePredict is an AI-powered desktop application for genomic risk assessment that prioritizes privacy through local-only processing. Built with Tauri (Rust), React, and Python ML components, it provides clinicians and researchers with tools to analyze genetic variants and assess disease risks.

## 🎯 Project Overview

- **Type**: Cross-platform desktop application (macOS, Ubuntu, Windows)
- **Architecture**: Tauri + React + Rust + Python ML
- **Processing**: 100% local/offline (no cloud dependencies)
- **Privacy**: GDPR/HIPAA compliant with differential privacy
- **File Support**: VCF, BAM, FASTQ genomic file formats

## 📋 Prerequisites

Before you begin, ensure you have the following installed:

### Required Software

1. **Rust** (1.70+)
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source ~/.cargo/env
   ```

2. **Node.js** (18+) & npm
   ```bash
   # macOS with Homebrew
   brew install node
   
   # Ubuntu/Debian
   curl -fsSL https://deb.nodesource.com/setup_18.x | sudo -E bash -
   sudo apt-get install -y nodejs
   ```

3. **Python** (3.9+)
   ```bash
   # macOS with Homebrew
   brew install python@3.11
   
   # Ubuntu/Debian
   sudo apt update && sudo apt install python3.11 python3.11-venv python3.11-dev
   ```

4. **Docker** & Docker Compose
   - Install from [Docker Desktop](https://www.docker.com/products/docker-desktop)

### Optional Tools

- **Git** (version control)
- **VS Code** (recommended editor)
- **Rust Analyzer** (VS Code extension)
- **Python** extension for VS Code

## 🚀 Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/yourorg/genepredict.git
cd genepredict
```

### 2. Environment Setup

```bash
# Copy environment configuration
cp config/env.example .env

# Edit .env with your local paths
nano .env
```

### 3. Install Dependencies

```bash
# Install all dependencies
npm run install:all

# Or install individually:
cd frontend && npm install
cd ../backend/python && pip install -r requirements.txt
cd ../rust && cargo build
```

### 4. Development with Docker (Recommended)

```bash
# Start all services
docker-compose up -d

# View logs
docker-compose logs -f

# Stop services
docker-compose down
```

### 5. Native Development (Alternative)

```bash
# Terminal 1: Start React frontend
cd frontend
npm run dev

# Terminal 2: Start Python backend
cd backend/python
python -m uvicorn genepredict.main:app --reload

# Terminal 3: Start Tauri app
cd backend/rust
cargo tauri dev
```

## 🏗️ Project Structure

```
genepredict/
├── 📁 frontend/                 # React + Tailwind CSS UI
│   ├── src/
│   │   ├── App.tsx             # Main application component
│   │   ├── components/         # Reusable UI components
│   │   └── utils/              # Frontend utilities
│   ├── package.json
│   └── Dockerfile.dev
│
├── 📁 backend/
│   ├── 📁 rust/                # Tauri application + plugin system
│   │   ├── src/
│   │   │   ├── main.rs         # Application entry point
│   │   │   ├── lib.rs          # Core application logic
│   │   │   ├── plugins.rs      # Plugin system
│   │   │   ├── config.rs       # Configuration management
│   │   │   └── error.rs        # Error handling
│   │   ├── Cargo.toml
│   │   └── tauri.conf.json
│   │
│   └── 📁 python/              # ML processing engine
│       ├── genepredict/
│       │   ├── models/         # AI/ML models
│       │   ├── processors/     # File processors (VCF, BAM, FASTQ)
│       │   └── utils/          # Python utilities
│       ├── requirements.txt
│       ├── pyproject.toml
│       └── Dockerfile.dev
│
├── 📁 docs/                    # Documentation
├── 📁 config/                  # Configuration files
├── 📁 data/                    # Local data storage (gitignored)
├── docker-compose.yml          # Development environment
└── README.md
```

## 🔧 Development Workflow

### Daily Development

1. **Start Environment**
   ```bash
   # With Docker (recommended)
   docker-compose up -d
   
   # Or natively
   npm run dev  # Frontend
   python -m uvicorn genepredict.main:app --reload  # Backend
   cargo tauri dev  # Desktop app
   ```

2. **Make Changes**
   - Frontend: Edit files in `frontend/src/`
   - Rust: Edit files in `backend/rust/src/`
   - Python: Edit files in `backend/python/genepredict/`

3. **Test Changes**
   ```bash
   # Frontend tests
   cd frontend && npm test
   
   # Python tests
   cd backend/python && pytest
   
   # Rust tests
   cd backend/rust && cargo test
   ```

4. **Build for Production**
   ```bash
   # Full production build
   cargo tauri build
   ```

### Code Standards

- **Rust**: Follow `rustfmt` and `clippy` recommendations
- **Python**: Use `black`, `isort`, and `mypy`
- **TypeScript**: Follow ESLint and Prettier configuration
- **Commits**: Use conventional commits format

### Pre-commit Hooks

Install pre-commit hooks to ensure code quality:

```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install

# Run manually
pre-commit run --all-files
```

## 🧪 Testing

### Running Tests

```bash
# All tests
npm run test:all

# Frontend tests
cd frontend && npm test

# Python tests
cd backend/python && pytest -v

# Rust tests
cd backend/rust && cargo test
```

### Test File Examples

Create test files for genomic data processing:

```bash
# Create sample VCF file
mkdir -p data/test
echo "##fileformat=VCFv4.2" > data/test/sample.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" >> data/test/sample.vcf
echo "1	100	.	A	T	30	PASS	." >> data/test/sample.vcf
```

## 🔌 Plugin Development

GenePredict uses a plugin system for file processors and ML models:

### Creating a File Processor Plugin

1. Implement the `GenomicProcessor` trait in Rust
2. Add corresponding Python processor class
3. Register plugin in plugin manager

### Creating an ML Model Plugin

1. Extend `BaseGenomicModel` in Python
2. Implement required abstract methods
3. Add model configuration

See `backend/rust/src/plugins.rs` and `backend/python/genepredict/models/` for examples.

## 🔐 Security & Privacy

GenePredict prioritizes security and privacy:

- **Local Processing**: No data leaves your machine
- **Encryption**: Temporary files are encrypted at rest
- **Differential Privacy**: ML models use privacy-preserving techniques
- **Secure Deletion**: Temporary files are securely wiped
- **GDPR/HIPAA**: Compliant design patterns

### Security Guidelines

1. Never log sensitive genomic data
2. Use secure random number generation
3. Implement proper error handling without data leakage
4. Follow OWASP guidelines for web security

## 📊 Data Management

### Reference Data

Download and configure reference datasets:

```bash
# Create data directories
mkdir -p data/reference/{1000genomes,clinvar,hgmd}

# Download 1000 Genomes data (example)
# wget -O data/reference/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
#   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

### Data Paths Configuration

Update your `.env` file with correct paths:

```bash
DATA_BASE_PATH="/path/to/your/genepredict/data"
REFERENCE_DATA_PATH="${DATA_BASE_PATH}/reference"
MODEL_CACHE_PATH="${DATA_BASE_PATH}/models"
```

## 🐛 Debugging

### Common Issues

1. **Rust compilation errors**
   ```bash
   # Update Rust
   rustup update
   
   # Clean build
   cd backend/rust && cargo clean && cargo build
   ```

2. **Python import errors**
   ```bash
   # Install in development mode
   cd backend/python && pip install -e .
   
   # Check PYTHONPATH
   echo $PYTHONPATH
   ```

3. **Frontend not starting**
   ```bash
   # Clear node_modules
   cd frontend && rm -rf node_modules && npm install
   ```

### Logging

Enable detailed logging:

```bash
# Rust logs
RUST_LOG=debug cargo tauri dev

# Python logs
LOG_LEVEL=DEBUG python -m uvicorn genepredict.main:app --reload

# Frontend logs
Check browser developer console
```

## 📚 Resources

### Documentation

- [Tauri Documentation](https://tauri.app/v1/guides/)
- [React Documentation](https://react.dev/)
- [TensorFlow Documentation](https://www.tensorflow.org/)
- [Genomics File Formats](https://samtools.github.io/hts-specs/)

### Community

- **Issues**: Report bugs and feature requests on GitHub
- **Discussions**: Join project discussions
- **Slack**: Internal team channel for development coordination

## 🎯 Next Steps

1. **Explore the codebase**: Start with `frontend/src/App.tsx` and `backend/rust/src/lib.rs`
2. **Run the test suite**: Ensure everything works on your machine
3. **Pick up a task**: Check the GitHub issues for "good first issue" labels
4. **Ask questions**: Don't hesitate to reach out to the team

## 🤝 Contributing

1. Create a feature branch from `main`
2. Make your changes with tests
3. Run the full test suite
4. Submit a pull request with clear description
5. Address review feedback

Welcome to the team! 🚀 