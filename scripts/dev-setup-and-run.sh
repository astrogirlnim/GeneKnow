#!/bin/bash

# GenePredict Development Environment Setup and Run Script
# This script intelligently sets up the dev environment and runs the application
# Features:
# - Checks for existing setup and skips unnecessary steps
# - Handles CTRL-C gracefully with proper cleanup
# - Comprehensive logging and status reporting
# - Supports multiple run modes (Docker, Native, Tauri)

set -e  # Exit on any error

# =============================================================================
# Configuration and Constants
# =============================================================================

# Colors for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly PURPLE='\033[0;35m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m' # No Color

# Project paths
readonly PROJECT_ROOT=$(pwd)
readonly FRONTEND_DIR="$PROJECT_ROOT/frontend"
readonly RUST_DIR="$PROJECT_ROOT/backend/rust"
readonly PYTHON_DIR="$PROJECT_ROOT/backend/python"
readonly DATA_DIR="$PROJECT_ROOT/data"
readonly CONFIG_DIR="$PROJECT_ROOT/config"
readonly SCRIPTS_DIR="$PROJECT_ROOT/scripts"

# Configuration
readonly ENV_FILE="$PROJECT_ROOT/.env"
readonly LOG_FILE="$PROJECT_ROOT/dev-setup.log"
readonly PID_FILE="$PROJECT_ROOT/.dev-pids"

# Default settings
RUN_MODE="tauri"  # Options: docker, native, tauri
SKIP_SETUP=false
VERBOSE=false
FORCE_REINSTALL=false
ENABLE_MONITORING=false

# Process tracking
declare -a BACKGROUND_PIDS=()

# =============================================================================
# Utility Functions
# =============================================================================

log() {
    local level="$1"
    shift
    local message="$*"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] [$level] $message" >> "$LOG_FILE"
    
    case "$level" in
        "ERROR")   echo -e "${RED}❌ $message${NC}" ;;
        "SUCCESS") echo -e "${GREEN}✅ $message${NC}" ;;
        "WARNING") echo -e "${YELLOW}⚠️  $message${NC}" ;;
        "INFO")    echo -e "${BLUE}ℹ️  $message${NC}" ;;
        "DEBUG")   [[ "$VERBOSE" == "true" ]] && echo -e "${PURPLE}🔍 $message${NC}" ;;
        "STEP")    echo -e "${CYAN}🔧 $message${NC}" ;;
        *)         echo "$message" ;;
    esac
}

print_header() {
    echo ""
    echo -e "${BLUE}============================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}============================================${NC}"
    echo ""
}

print_section() {
    echo ""
    echo -e "${CYAN}--- $1 ---${NC}"
}

show_help() {
    cat << EOF
GenePredict Development Environment Setup and Run Script

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message
    -v, --verbose           Enable verbose output and debugging
    -m, --mode MODE         Run mode: docker, native, or tauri (default: tauri)
    -s, --skip-setup        Skip environment setup, just run the app
    -f, --force             Force reinstall all dependencies
    -M, --monitoring        Enable monitoring services (Prometheus, Grafana)
    --cleanup               Clean up all processes and temporary files
    --status                Show status of running services
    --logs                  Show recent logs
    --test-ml               Test ML library compatibility

RUN MODES:
    docker      Use Docker Compose for all services
    native      Run services natively (Python backend + React frontend)
    tauri       Run Tauri desktop app (default, recommended)

EXAMPLES:
    $0                      # Setup and run with Tauri (recommended)
    $0 -m docker           # Setup and run with Docker Compose
    $0 -m native -v        # Run natively with verbose output
    $0 -s -m tauri         # Skip setup, just run Tauri app
    $0 --cleanup           # Clean up all running processes
    $0 --status            # Check status of services

ENVIRONMENT:
    The script will create a .env file from config/env.example if it doesn't exist.
    Edit .env to customize paths and settings for your environment.

EOF
}

# =============================================================================
# Signal Handling and Cleanup
# =============================================================================

cleanup() {
    log "INFO" "Received interrupt signal, cleaning up..."
    
    # Kill background processes
    if [[ -f "$PID_FILE" ]]; then
        while IFS= read -r pid; do
            if kill -0 "$pid" 2>/dev/null; then
                log "INFO" "Terminating process $pid"
                kill -TERM "$pid" 2>/dev/null || true
                sleep 2
                kill -KILL "$pid" 2>/dev/null || true
            fi
        done < "$PID_FILE"
        rm -f "$PID_FILE"
    fi
    
    # Kill processes by name if PID file is missing
    pkill -f "vite" 2>/dev/null || true
    pkill -f "cargo tauri" 2>/dev/null || true
    pkill -f "uvicorn" 2>/dev/null || true
    
    # Docker cleanup if running
    if [[ "$RUN_MODE" == "docker" ]]; then
        log "INFO" "Stopping Docker services..."
        docker-compose down 2>/dev/null || true
    fi
    
    # Clean up temporary files
    rm -f "$PROJECT_ROOT"/.dev-* 2>/dev/null || true
    
    log "SUCCESS" "Cleanup completed"
    exit 0
}

# Set up signal handlers
trap cleanup SIGINT SIGTERM EXIT

# =============================================================================
# System Checks
# =============================================================================

check_prerequisites() {
    print_section "Checking Prerequisites"
    
    local missing_deps=()
    
    # Check for required commands
    local required_commands=(
        "node:Node.js (18+)"
        "npm:npm package manager"
        "python3:Python (3.9+)"
        "pip3:Python package manager"
        "cargo:Rust toolchain"
        "rustc:Rust compiler"
    )
    
    for cmd_desc in "${required_commands[@]}"; do
        local cmd="${cmd_desc%%:*}"
        local desc="${cmd_desc##*:}"
        
        if ! command -v "$cmd" &> /dev/null; then
            missing_deps+=("$desc")
            log "ERROR" "$desc not found"
        else
            local version=""
            case "$cmd" in
                "node") version=$(node --version) ;;
                "python3") version=$(python3 --version) ;;
                "cargo") version=$(cargo --version | cut -d' ' -f2) ;;
            esac
            log "SUCCESS" "$desc found: $version"
        fi
    done
    
    # Check Docker if needed
    if [[ "$RUN_MODE" == "docker" ]]; then
        if ! command -v docker &> /dev/null; then
            missing_deps+=("Docker")
            log "ERROR" "Docker not found"
        else
            log "SUCCESS" "Docker found: $(docker --version)"
        fi
        
        if ! command -v docker-compose &> /dev/null; then
            missing_deps+=("Docker Compose")
            log "ERROR" "Docker Compose not found"
        else
            log "SUCCESS" "Docker Compose found: $(docker-compose --version)"
        fi
    fi
    
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        log "ERROR" "Missing required dependencies:"
        for dep in "${missing_deps[@]}"; do
            log "ERROR" "  - $dep"
        done
        
        log "INFO" "Please install missing dependencies:"
        log "INFO" "  macOS: brew install node python rust docker"
        log "INFO" "  Ubuntu: apt install nodejs npm python3 python3-pip docker.io docker-compose"
        log "INFO" "  Rust: curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
        exit 1
    fi
    
    log "SUCCESS" "All prerequisites satisfied"
}

# =============================================================================
# Environment Setup
# =============================================================================

setup_environment() {
    print_section "Setting Up Environment"
    
    # Create .env file if it doesn't exist
    if [[ ! -f "$ENV_FILE" ]]; then
        if [[ -f "$CONFIG_DIR/env.example" ]]; then
            log "STEP" "Creating .env file from template"
            cp "$CONFIG_DIR/env.example" "$ENV_FILE"
            
            # Update paths in .env to match current directory
            if [[ "$OSTYPE" == "darwin"* ]]; then
                sed -i '' "s|/Users/youruser/GenePredict|$PROJECT_ROOT|g" "$ENV_FILE"
            else
                sed -i "s|/Users/youruser/GenePredict|$PROJECT_ROOT|g" "$ENV_FILE"
            fi
            
            log "SUCCESS" "Created .env file with project-specific paths"
            log "INFO" "Edit .env to customize settings for your environment"
        else
            log "WARNING" "No env.example found, creating minimal .env"
            cat > "$ENV_FILE" << EOF
# GenePredict Environment Configuration
APP_NAME="GenePredict"
APP_VERSION="0.1.0"
DEBUG=true
DATA_BASE_PATH="$PROJECT_ROOT/data"
LOG_LEVEL=INFO
EOF
        fi
    else
        log "INFO" ".env file already exists, skipping creation"
    fi
    
    # Create necessary directories
    log "STEP" "Creating data directories"
    mkdir -p "$DATA_DIR"/{tcga/{raw,processed,cache,models},reference/{clinvar,1000genomes,cosmic},temp,user}
    
    # Create log file
    touch "$LOG_FILE"
    log "SUCCESS" "Environment setup completed"
}

# =============================================================================
# Dependency Installation
# =============================================================================

install_frontend_dependencies() {
    print_section "Frontend Dependencies"
    
    cd "$FRONTEND_DIR"
    
    if [[ "$FORCE_REINSTALL" == "true" ]] || [[ ! -d "node_modules" ]]; then
        log "STEP" "Installing frontend dependencies..."
        npm install
        log "SUCCESS" "Frontend dependencies installed"
    else
        log "INFO" "Frontend dependencies already installed, skipping"
    fi
    
    cd "$PROJECT_ROOT"
}

install_python_dependencies() {
    print_section "Python Dependencies"
    
    cd "$PYTHON_DIR"
    
    # Check Python version for compatibility
    local python_version=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    log "INFO" "Detected Python version: $python_version"
    
    # Determine which requirements file to use
    local requirements_file="requirements.txt"
    if [[ $(python3 -c "import sys; print(sys.version_info >= (3, 13))") == "True" ]]; then
        log "INFO" "Using Python 3.13+ compatible requirements (JAX instead of TensorFlow)"
        requirements_file="requirements.txt"  # Already updated for Python 3.13
    else
        log "INFO" "Using legacy requirements (includes TensorFlow for Python 3.9-3.12)"
        if [[ -f "requirements-legacy.txt" ]]; then
            requirements_file="requirements-legacy.txt"
        else
            log "WARNING" "Legacy requirements file not found, using default"
        fi
    fi
    
    # Check if virtual environment exists
    if [[ ! -d "venv" ]] || [[ "$FORCE_REINSTALL" == "true" ]]; then
        log "STEP" "Creating Python virtual environment..."
        python3 -m venv venv
        log "SUCCESS" "Python virtual environment created"
    else
        log "INFO" "Python virtual environment already exists"
    fi
    
    # Activate virtual environment
    source venv/bin/activate
    
    # Install/upgrade pip
    log "STEP" "Upgrading pip..."
    pip install --upgrade pip
    
    # Install dependencies with better error handling
    if [[ "$FORCE_REINSTALL" == "true" ]] || ! pip freeze | grep -q -E "(tensorflow|jax|torch)"; then
        log "STEP" "Installing Python dependencies from $requirements_file..."
        
        # Try to install requirements with better error handling
        if pip install -r "$requirements_file"; then
            log "SUCCESS" "Core Python dependencies installed successfully"
        else
            log "ERROR" "Failed to install some dependencies from $requirements_file"
            log "INFO" "Attempting to install minimal requirements..."
            
            # Install minimal requirements that should work
            local minimal_deps=(
                "numpy>=1.24.0"
                "pandas>=2.0.0"
                "scikit-learn>=1.3.0"
                "fastapi>=0.100.0"
                "uvicorn>=0.23.0"
                "pydantic>=2.0.0"
                "python-dotenv>=1.0.0"
            )
            
            for dep in "${minimal_deps[@]}"; do
                if pip install "$dep"; then
                    log "SUCCESS" "Installed: $dep"
                else
                    log "WARNING" "Failed to install: $dep"
                fi
            done
        fi
        
        # Install additional TCGA dependencies (optional)
        log "STEP" "Installing optional TCGA dependencies..."
        local optional_deps=("requests" "aiohttp" "aiofiles")
        for dep in "${optional_deps[@]}"; do
            if pip install "$dep" 2>/dev/null; then
                log "SUCCESS" "Installed optional dependency: $dep"
            else
                log "WARNING" "Failed to install optional dependency: $dep"
            fi
        done
        
        # Try to install gdctools separately as it may have issues
        if pip install gdctools 2>/dev/null; then
            log "SUCCESS" "Installed gdctools for TCGA integration"
        else
            log "WARNING" "Failed to install gdctools - TCGA integration may be limited"
        fi
        
        log "SUCCESS" "Python dependencies installation completed"
    else
        log "INFO" "Python dependencies already installed, skipping"
    fi
    
    # Verify critical packages are installed
    log "STEP" "Verifying critical packages..."
    local critical_packages=("numpy" "pandas" "fastapi")
    local all_critical_installed=true
    
    for package in "${critical_packages[@]}"; do
        if python -c "import $package" 2>/dev/null; then
            log "SUCCESS" "Verified: $package is installed"
        else
            log "ERROR" "Critical package missing: $package"
            all_critical_installed=false
        fi
    done
    
    if [[ "$all_critical_installed" == "false" ]]; then
        log "ERROR" "Some critical packages are missing. The application may not work properly."
        log "INFO" "Try running with --force to reinstall all dependencies"
    fi
    
    cd "$PROJECT_ROOT"
}

install_rust_dependencies() {
    print_section "Rust Dependencies"
    
    cd "$RUST_DIR"
    
    if [[ "$FORCE_REINSTALL" == "true" ]] || [[ ! -d "target" ]]; then
        log "STEP" "Building Rust dependencies..."
        cargo build
        log "SUCCESS" "Rust dependencies built"
    else
        log "INFO" "Rust dependencies already built, skipping"
    fi
    
    cd "$PROJECT_ROOT"
}

install_all_dependencies() {
    if [[ "$SKIP_SETUP" == "true" ]]; then
        log "INFO" "Skipping dependency installation as requested"
        return
    fi
    
    print_header "Installing Dependencies"
    
    install_frontend_dependencies
    install_python_dependencies
    install_rust_dependencies
    
    log "SUCCESS" "All dependencies installed successfully"
}

# =============================================================================
# Application Runners
# =============================================================================

run_with_docker() {
    print_header "Running with Docker Compose"
    
    # Check if Docker is running
    if ! docker info &> /dev/null; then
        log "ERROR" "Docker daemon is not running. Please start Docker and try again."
        exit 1
    fi
    
    log "STEP" "Starting Docker services..."
    
    # Build and start services
    docker-compose up --build -d
    
    if [[ "$ENABLE_MONITORING" == "true" ]]; then
        log "INFO" "Monitoring enabled - Grafana available at http://localhost:3001"
        log "INFO" "Prometheus available at http://localhost:9090"
    fi
    
    log "SUCCESS" "Docker services started successfully"
    log "INFO" "Frontend available at: http://localhost:3000"
    log "INFO" "Backend API available at: http://localhost:8000"
    
    # Follow logs
    log "INFO" "Following logs (Ctrl+C to stop)..."
    docker-compose logs -f
}

run_native() {
    print_header "Running Native Services"
    
    # Start Python backend
    print_section "Starting Python Backend"
    cd "$PYTHON_DIR"
    source venv/bin/activate
    
    log "STEP" "Starting Python backend server..."
    python -m uvicorn genepredict.main:app --reload --host 127.0.0.1 --port 8000 &
    local python_pid=$!
    BACKGROUND_PIDS+=($python_pid)
    echo "$python_pid" >> "$PID_FILE"
    
    # Wait for backend to start
    sleep 3
    
    # Check if backend is running
    if kill -0 "$python_pid" 2>/dev/null; then
        log "SUCCESS" "Python backend started (PID: $python_pid)"
    else
        log "ERROR" "Failed to start Python backend"
        exit 1
    fi
    
    # Start React frontend
    print_section "Starting React Frontend"
    cd "$FRONTEND_DIR"
    
    log "STEP" "Starting React development server..."
    npm run dev &
    local frontend_pid=$!
    BACKGROUND_PIDS+=($frontend_pid)
    echo "$frontend_pid" >> "$PID_FILE"
    
    # Wait for frontend to start
    sleep 3
    
    if kill -0 "$frontend_pid" 2>/dev/null; then
        log "SUCCESS" "React frontend started (PID: $frontend_pid)"
    else
        log "ERROR" "Failed to start React frontend"
        exit 1
    fi
    
    cd "$PROJECT_ROOT"
    
    log "SUCCESS" "Native services started successfully"
    log "INFO" "Frontend available at: http://localhost:3000"
    log "INFO" "Backend API available at: http://localhost:8000"
    
    # Wait for processes
    wait
}

run_tauri() {
    print_header "Running Tauri Desktop Application"
    
    # Start Python backend first
    print_section "Starting Python Backend"
    cd "$PYTHON_DIR"
    source venv/bin/activate
    
    log "STEP" "Starting Python backend server..."
    python -m uvicorn genepredict.main:app --reload --host 127.0.0.1 --port 8000 &
    local python_pid=$!
    BACKGROUND_PIDS+=($python_pid)
    echo "$python_pid" >> "$PID_FILE"
    
    # Wait for backend to start
    sleep 3
    
    if kill -0 "$python_pid" 2>/dev/null; then
        log "SUCCESS" "Python backend started (PID: $python_pid)"
    else
        log "WARNING" "Python backend may not have started properly, but continuing with Tauri..."
    fi
    
    # Start Tauri application
    print_section "Starting Tauri Desktop App"
    cd "$RUST_DIR"
    
    log "STEP" "Starting Tauri development server..."
    log "INFO" "This will open the desktop application window..."
    
    # Run Tauri in foreground (it will handle the frontend automatically)
    cargo tauri dev
    
    cd "$PROJECT_ROOT"
}

# =============================================================================
# Status and Monitoring
# =============================================================================

show_status() {
    print_header "Service Status"
    
    local services_running=false
    
    # Check Docker services
    if docker-compose ps 2>/dev/null | grep -q "Up"; then
        log "INFO" "Docker services status:"
        docker-compose ps
        services_running=true
    fi
    
    # Check native processes
    if [[ -f "$PID_FILE" ]]; then
        log "INFO" "Native processes:"
        while IFS= read -r pid; do
            if kill -0 "$pid" 2>/dev/null; then
                local cmd=$(ps -p "$pid" -o comm= 2>/dev/null || echo "unknown")
                log "SUCCESS" "Process $pid ($cmd) is running"
                services_running=true
            else
                log "WARNING" "Process $pid is not running"
            fi
        done < "$PID_FILE"
    fi
    
    # Check ports
    local ports=(3000 8000 5432 6379 9000)
    for port in "${ports[@]}"; do
        if lsof -i ":$port" &>/dev/null; then
            local process=$(lsof -ti ":$port" | head -1)
            local cmd=$(ps -p "$process" -o comm= 2>/dev/null || echo "unknown")
            log "INFO" "Port $port is in use by $cmd (PID: $process)"
            services_running=true
        fi
    done
    
    if [[ "$services_running" == "false" ]]; then
        log "INFO" "No GenePredict services are currently running"
    fi
}

show_logs() {
    print_header "Recent Logs"
    
    if [[ -f "$LOG_FILE" ]]; then
        tail -50 "$LOG_FILE"
    else
        log "INFO" "No log file found"
    fi
}

test_ml_compatibility() {
    print_header "ML Compatibility Test"
    
    if [[ ! -f "$PYTHON_DIR/test_ml_compatibility.py" ]]; then
        log "ERROR" "ML compatibility test script not found"
        exit 1
    fi
    
    cd "$PYTHON_DIR"
    
    # Check if virtual environment exists
    if [[ -d "venv" ]]; then
        log "INFO" "Using existing virtual environment"
        source venv/bin/activate
    else
        log "WARNING" "No virtual environment found, using system Python"
    fi
    
    log "INFO" "Running ML compatibility test..."
    python3 test_ml_compatibility.py
    
    cd "$PROJECT_ROOT"
}

# =============================================================================
# Main Functions
# =============================================================================

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            -m|--mode)
                RUN_MODE="$2"
                if [[ ! "$RUN_MODE" =~ ^(docker|native|tauri)$ ]]; then
                    log "ERROR" "Invalid run mode: $RUN_MODE. Must be: docker, native, or tauri"
                    exit 1
                fi
                shift 2
                ;;
            -s|--skip-setup)
                SKIP_SETUP=true
                shift
                ;;
            -f|--force)
                FORCE_REINSTALL=true
                shift
                ;;
            -M|--monitoring)
                ENABLE_MONITORING=true
                shift
                ;;
            --cleanup)
                cleanup
                exit 0
                ;;
            --status)
                show_status
                exit 0
                ;;
            --logs)
                show_logs
                exit 0
                ;;
            --test-ml)
                test_ml_compatibility
                exit 0
                ;;
            *)
                log "ERROR" "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done
}

main() {
    # Initialize log file
    echo "GenePredict Development Setup - $(date)" > "$LOG_FILE"
    
    print_header "GenePredict Development Environment"
    log "INFO" "Starting GenePredict development setup and run script"
    log "INFO" "Run mode: $RUN_MODE"
    log "INFO" "Project root: $PROJECT_ROOT"
    log "INFO" "Log file: $LOG_FILE"
    
    # Run setup steps
    check_prerequisites
    setup_environment
    install_all_dependencies
    
    # Run the application based on mode
    case "$RUN_MODE" in
        "docker")
            run_with_docker
            ;;
        "native")
            run_native
            ;;
        "tauri")
            run_tauri
            ;;
        *)
            log "ERROR" "Invalid run mode: $RUN_MODE"
            exit 1
            ;;
    esac
}

# =============================================================================
# Script Entry Point
# =============================================================================

# Parse command line arguments
parse_arguments "$@"

# Run main function
main

# Note: cleanup() will be called automatically via trap on exit 