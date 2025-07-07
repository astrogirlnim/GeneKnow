# GenePredict Development Scripts

This directory contains scripts to help with development, testing, and running the GenePredict application.

## 🚀 Quick Start

For the fastest way to get started:

```bash
# Test if your environment is ready
./scripts/test-setup.sh

# Start the application (recommended)
./scripts/quick-start.sh
```

## 📜 Available Scripts

### Main Development Scripts

#### `dev-setup-and-run.sh` - **Comprehensive Development Environment Setup**
The main script that intelligently sets up your development environment and runs the application.

**Features:**
- ✅ Checks prerequisites (Node.js, Python, Rust)
- ✅ Sets up environment files (.env)
- ✅ Creates necessary data directories
- ✅ Installs all dependencies (Frontend, Python, Rust)
- ✅ Supports multiple run modes (Docker, Native, Tauri)
- ✅ Handles CTRL-C gracefully with proper cleanup
- ✅ Comprehensive logging and status reporting

**Usage:**
```bash
# Default: Setup and run with Tauri (recommended)
./scripts/dev-setup-and-run.sh

# Run with Docker Compose
./scripts/dev-setup-and-run.sh -m docker

# Run natively (Python backend + React frontend)
./scripts/dev-setup-and-run.sh -m native

# Skip setup, just run the app
./scripts/dev-setup-and-run.sh -s

# Force reinstall all dependencies
./scripts/dev-setup-and-run.sh -f

# Enable verbose output
./scripts/dev-setup-and-run.sh -v

# Show help
./scripts/dev-setup-and-run.sh --help
```

**Management Commands:**
```bash
# Check status of running services
./scripts/dev-setup-and-run.sh --status

# View recent logs
./scripts/dev-setup-and-run.sh --logs

# Clean up all processes and temporary files
./scripts/dev-setup-and-run.sh --cleanup
```

#### `quick-start.sh` - **Simple One-Click Start**
A wrapper around the main script with sensible defaults for new users.

**Usage:**
```bash
# Interactive start (asks for confirmation)
./scripts/quick-start.sh

# Automatic start (no confirmation)
./scripts/quick-start.sh --yes
```

#### `test-setup.sh` - **Environment Validation**
Tests your development environment setup without running the full application.

**Usage:**
```bash
./scripts/test-setup.sh
```

**What it tests:**
- Prerequisites (Node.js, Python, Rust)
- Project structure
- Dependency installation status
- Configuration files
- Basic compilation tests

### Legacy Testing Scripts

#### `test-cli.sh` - **Comprehensive CLI Testing**
Advanced testing script with multiple test modes and detailed reporting.

**Usage:**
```bash
# Run all tests
./scripts/test-cli.sh

# Run only frontend tests
./scripts/test-cli.sh -f

# Run with verbose output
./scripts/test-cli.sh -v

# Show help
./scripts/test-cli.sh --help
```

#### `test-tcga.py` - **TCGA Integration Testing**
Tests TCGA (The Cancer Genome Atlas) data integration functionality.

**Usage:**
```bash
python scripts/test-tcga.py
```

#### `test-tcga-simple.py` - **Simple TCGA Testing**
Simplified version of TCGA testing for quick validation.

**Usage:**
```bash
python scripts/test-tcga-simple.py
```

## 🎯 Recommended Workflow

### For New Developers

1. **First Time Setup:**
   ```bash
   # Test your environment
   ./scripts/test-setup.sh
   
   # Start the application
   ./scripts/quick-start.sh
   ```

2. **Daily Development:**
   ```bash
   # Quick start (skips setup if already done)
   ./scripts/dev-setup-and-run.sh -s
   ```

3. **When Dependencies Change:**
   ```bash
   # Force reinstall all dependencies
   ./scripts/dev-setup-and-run.sh -f
   ```

### For Different Development Modes

#### Desktop Application Development (Recommended)
```bash
# Tauri desktop app with hot reload
./scripts/dev-setup-and-run.sh -m tauri
```

#### Web Development
```bash
# Native services (Python backend + React frontend)
./scripts/dev-setup-and-run.sh -m native
```

#### Full Stack with Services
```bash
# Docker Compose with all services
./scripts/dev-setup-and-run.sh -m docker
```

## 🔧 Troubleshooting

### Common Issues

#### Prerequisites Missing
```bash
# Check what's missing
./scripts/test-setup.sh

# Install on macOS
brew install node python rust docker

# Install on Ubuntu
sudo apt update
sudo apt install nodejs npm python3 python3-pip docker.io docker-compose
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

#### Dependencies Not Installing
```bash
# Force reinstall everything
./scripts/dev-setup-and-run.sh -f

# Or install manually:
cd frontend && npm install
cd ../backend/python && python3 -m venv venv && source venv/bin/activate && pip install -r requirements.txt
cd ../rust && cargo build
```

#### Application Won't Start
```bash
# Check what's running
./scripts/dev-setup-and-run.sh --status

# Clean up and try again
./scripts/dev-setup-and-run.sh --cleanup
./scripts/dev-setup-and-run.sh
```

#### Port Conflicts
```bash
# Check what's using ports
lsof -i :3000  # Frontend
lsof -i :8000  # Backend
lsof -i :5432  # PostgreSQL
lsof -i :6379  # Redis

# Kill conflicting processes
sudo kill -9 <PID>
```

### Getting Help

1. **View logs:**
   ```bash
   ./scripts/dev-setup-and-run.sh --logs
   ```

2. **Check service status:**
   ```bash
   ./scripts/dev-setup-and-run.sh --status
   ```

3. **Run with verbose output:**
   ```bash
   ./scripts/dev-setup-and-run.sh -v
   ```

4. **Test environment:**
   ```bash
   ./scripts/test-setup.sh
   ```

## 📁 Generated Files

The scripts create several files during operation:

- `.env` - Environment configuration (auto-generated from `config/env.example`)
- `dev-setup.log` - Detailed log of setup and run operations
- `.dev-pids` - Process IDs of running services (for cleanup)
- `data/` - Data directories for genomic files and models

## 🧹 Cleanup

To clean up all generated files and stop all services:

```bash
./scripts/dev-setup-and-run.sh --cleanup
```

This will:
- Stop all running processes
- Clean up temporary files
- Stop Docker services (if running)
- Remove PID files

## 🔐 Security Notes

- All genomic data processing happens locally
- No user data is sent to external services
- Environment files (`.env`) may contain sensitive paths
- Clean up regularly to avoid accumulating temporary files

---

**Need help?** Run any script with `--help` for detailed usage information. 