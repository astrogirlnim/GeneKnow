# GenePredict Tech Context

## Technology Stack

### Frontend Technologies
- **React 18.3.1**: Modern UI framework with hooks and concurrent features
- **TypeScript 5.6.3**: Type safety for complex genomic data structures
- **Tailwind CSS 3.4.15**: Utility-first styling for rapid UI development
- **Vite 6.0.3**: Fast development server and build tool
- **React Testing Library**: Component testing framework

### Native Desktop Layer
- **Tauri 2.1.2**: Secure, cross-platform desktop app framework
- **Rust 1.88.0**: Systems programming for file handling and plugin orchestration
- **tokio**: Async runtime for concurrent genomic processing
- **serde**: JSON serialization between Rust and Python/React layers

### AI/ML Technologies
- **Python 3.13.5**: Core ML and genomic processing language
- **NumPy 2.2.1**: Numerical computing for genomic arrays
- **Pandas 2.2.3**: Data manipulation for genomic datasets
- **scikit-learn 1.6.0**: Machine learning algorithms for risk prediction
- **TensorFlow**: Deep learning models (pending Python 3.13 compatibility)
- **Llama 3.1**: Local LLM via HuggingFace for report generation

### Genomic Processing Libraries
- **BioPython**: FASTA/FASTQ sequence analysis
- **pysam**: BAM/SAM file processing
- **pyVCF**: Variant Call Format parsing
- **aiohttp 3.11.15**: Async HTTP client for TCGA API integration
- **aiofiles 24.1.0**: Async file I/O for large genomic datasets

### Development & Infrastructure
- **Node.js 20.19.2**: JavaScript runtime for frontend tooling
- **Docker & Docker Compose**: Containerized development environment
- **GitHub Actions**: CI/CD pipeline for automated testing
- **pre-commit**: Code quality hooks (Black, ESLint, Prettier)
- **pytest**: Python testing framework

## Development Environment Setup

### Prerequisites (Verified Working)
- **Node.js**: v20.19.2 (with npm 10.9.0)
- **Python**: 3.13.5 with venv support
- **Rust**: 1.88.0 with cargo
- **Git**: Version control for collaboration

### Local Development Commands
```bash
# Install all dependencies
npm run install:all

# Start development environment
npm run dev                    # Frontend only
cargo tauri dev               # Full desktop app

# Python virtual environment
python -m venv venv
source venv/bin/activate      # macOS/Linux
venv\Scripts\activate.bat     # Windows
pip install -r requirements.txt

# Testing
npm test                      # Frontend tests
cargo test                    # Rust tests  
pytest                        # Python tests
```

### Project Structure
```
LiteratureGapper/
├── frontend/                 # React + Tailwind UI
├── backend/
│   ├── rust/                 # Tauri native layer
│   └── python/               # ML and genomic processing
├── data/                     # Local TCGA and reference data
├── scripts/                  # Development and testing utilities
├── docs/                     # Technical documentation
└── memory_bank/              # AI assistant memory system
```

## Environment Configuration

### Required Environment Variables
```bash
# TCGA Integration
TCGA_API_BASE_URL=https://api.gdc.cancer.gov
TCGA_DATA_DIR=./data/tcga
TCGA_CACHE_SIZE_GB=10

# Local Model Paths
HUGGINGFACE_CACHE_DIR=./data/models
LLAMA_MODEL_PATH=./data/models/llama-3.1

# Privacy & Security
GENOMIC_TEMP_DIR=./data/temp
ENABLE_AUDIT_LOGGING=true
MAX_MEMORY_USAGE_GB=8
```

### Data Directory Structure
```
data/
├── tcga/
│   ├── raw/                  # Downloaded TCGA datasets
│   ├── processed/            # Processed reference data
│   ├── cache/                # API response cache
│   └── models/               # Trained ML models
├── reference/
│   ├── 1000genomes/          # Population frequency data
│   └── clinvar/              # Clinical variant database
└── temp/                     # Encrypted temporary files
```

## Platform Requirements

### Minimum System Requirements
- **RAM**: 8GB (16GB recommended for large genomic files)
- **Storage**: 50GB free space for reference datasets
- **CPU**: Multi-core processor (4+ cores recommended)
- **Network**: Internet connection for initial TCGA data download

### Supported Platforms
- **macOS**: 11.0+ (Big Sur and later)
- **Windows**: 10/11 (64-bit)
- **Linux**: Ubuntu 20.04+, Fedora 35+, other systemd distributions

### Browser Requirements (for embedded frontend)
- **WebKit**: Latest (embedded in Tauri)
- **Hardware Acceleration**: GPU support for ML operations

## Technical Constraints & Decisions

### Privacy Constraints
- **No Cloud Processing**: All genomic data must remain local
- **No Telemetry**: No usage data transmitted to external servers
- **Compliance**: GDPR/HIPAA requirements drive architecture decisions

### Performance Constraints
- **Memory Management**: Large genomic files (>1GB) processed in streaming chunks
- **CPU Usage**: ML operations limited to available cores, with user-configurable limits
- **Storage**: Intelligent caching prevents unlimited disk usage

### Compatibility Constraints
- **Python 3.13**: TensorFlow not yet compatible, using scikit-learn temporarily
- **File Formats**: Support for FASTQ, BAM, VCF, 23andMe raw data formats
- **Reference Genomes**: GRCh37/hg19 and GRCh38/hg38 support required

### Security Decisions
- **Rust for File I/O**: Memory safety for genomic data handling
- **Local Encryption**: Sensitive temp files encrypted with OS keychain
- **Sandboxed Processing**: Python ML operations isolated from system

## Development Workflow

### Git Workflow
- **Branch Strategy**: Feature branches from main
- **Commit Standards**: Conventional commits with detailed genomic context
- **Pre-commit Hooks**: Code formatting and security checks
- **No Secrets**: All sensitive config in environment variables

### Testing Strategy
- **Unit Tests**: Each genomic processor plugin independently tested
- **Integration Tests**: Full pipeline testing with synthetic genomic data
- **Manual Testing**: TCGA API connectivity and reference data validation
- **No Test Files in Repo**: Temporary test files automatically cleaned up

### Deployment
- **Local Development**: Direct cargo/npm commands
- **CI/CD**: GitHub Actions for cross-platform builds
- **Distribution**: Platform-specific installers via Tauri bundler 