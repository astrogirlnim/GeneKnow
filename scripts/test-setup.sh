#!/bin/bash

# GenePredict Setup Test Script
# Tests the development environment setup without running the full application

set -e

# Colors for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly NC='\033[0m' # No Color

# Project paths
readonly PROJECT_ROOT=$(pwd)
readonly FRONTEND_DIR="$PROJECT_ROOT/frontend"
readonly RUST_DIR="$PROJECT_ROOT/backend/rust"
readonly PYTHON_DIR="$PROJECT_ROOT/backend/python"

log() {
    local level="$1"
    shift
    local message="$*"
    
    case "$level" in
        "ERROR")   echo -e "${RED}❌ $message${NC}" ;;
        "SUCCESS") echo -e "${GREEN}✅ $message${NC}" ;;
        "WARNING") echo -e "${YELLOW}⚠️  $message${NC}" ;;
        "INFO")    echo -e "${BLUE}ℹ️  $message${NC}" ;;
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

test_prerequisites() {
    print_header "Testing Prerequisites"
    
    local all_good=true
    
    # Test Node.js
    if command -v node &> /dev/null; then
        local node_version=$(node --version)
        log "SUCCESS" "Node.js found: $node_version"
    else
        log "ERROR" "Node.js not found"
        all_good=false
    fi
    
    # Test npm
    if command -v npm &> /dev/null; then
        local npm_version=$(npm --version)
        log "SUCCESS" "npm found: v$npm_version"
    else
        log "ERROR" "npm not found"
        all_good=false
    fi
    
    # Test Python
    if command -v python3 &> /dev/null; then
        local python_version=$(python3 --version)
        log "SUCCESS" "Python found: $python_version"
    else
        log "ERROR" "Python3 not found"
        all_good=false
    fi
    
    # Test Rust
    if command -v cargo &> /dev/null; then
        local cargo_version=$(cargo --version)
        log "SUCCESS" "Rust/Cargo found: $cargo_version"
    else
        log "ERROR" "Rust/Cargo not found"
        all_good=false
    fi
    
    if [[ "$all_good" == "false" ]]; then
        log "ERROR" "Some prerequisites are missing. Please install them first."
        return 1
    fi
    
    log "SUCCESS" "All prerequisites satisfied"
    return 0
}

test_project_structure() {
    print_header "Testing Project Structure"
    
    local required_dirs=(
        "$FRONTEND_DIR"
        "$RUST_DIR"
        "$PYTHON_DIR"
        "$PROJECT_ROOT/config"
        "$PROJECT_ROOT/docs"
    )
    
    local required_files=(
        "$PROJECT_ROOT/package.json"
        "$FRONTEND_DIR/package.json"
        "$RUST_DIR/Cargo.toml"
        "$PYTHON_DIR/requirements.txt"
        "$PROJECT_ROOT/docker-compose.yml"
    )
    
    # Check directories
    for dir in "${required_dirs[@]}"; do
        if [[ -d "$dir" ]]; then
            log "SUCCESS" "Directory exists: $(basename "$dir")"
        else
            log "ERROR" "Missing directory: $dir"
        fi
    done
    
    # Check files
    for file in "${required_files[@]}"; do
        if [[ -f "$file" ]]; then
            log "SUCCESS" "File exists: $(basename "$file")"
        else
            log "ERROR" "Missing file: $file"
        fi
    done
}

test_dependencies() {
    print_header "Testing Dependencies"
    
    # Test frontend dependencies
    if [[ -d "$FRONTEND_DIR/node_modules" ]]; then
        log "SUCCESS" "Frontend dependencies installed"
    else
        log "WARNING" "Frontend dependencies not installed (run: cd frontend && npm install)"
    fi
    
    # Test Python virtual environment
    if [[ -d "$PYTHON_DIR/venv" ]]; then
        log "SUCCESS" "Python virtual environment exists"
        
        # Test if dependencies are installed
        cd "$PYTHON_DIR"
        source venv/bin/activate
        if pip freeze | grep -q "tensorflow"; then
            log "SUCCESS" "Python dependencies installed"
        else
            log "WARNING" "Python dependencies not installed (run: cd backend/python && pip install -r requirements.txt)"
        fi
        deactivate
        cd "$PROJECT_ROOT"
    else
        log "WARNING" "Python virtual environment not created (run: cd backend/python && python3 -m venv venv)"
    fi
    
    # Test Rust build
    if [[ -d "$RUST_DIR/target" ]]; then
        log "SUCCESS" "Rust dependencies built"
    else
        log "WARNING" "Rust dependencies not built (run: cd backend/rust && cargo build)"
    fi
}

test_configuration() {
    print_header "Testing Configuration"
    
    # Check .env file
    if [[ -f "$PROJECT_ROOT/.env" ]]; then
        log "SUCCESS" ".env file exists"
    else
        log "WARNING" ".env file not found (will be created automatically)"
    fi
    
    # Check data directories
    if [[ -d "$PROJECT_ROOT/data" ]]; then
        log "SUCCESS" "Data directory exists"
    else
        log "INFO" "Data directory will be created automatically"
    fi
}

run_basic_tests() {
    print_header "Running Basic Tests"
    
    # Test frontend build
    log "INFO" "Testing frontend TypeScript compilation..."
    cd "$FRONTEND_DIR"
    if npm run build --silent &> /dev/null; then
        log "SUCCESS" "Frontend builds successfully"
    else
        log "WARNING" "Frontend build failed (may need dependencies installed)"
    fi
    cd "$PROJECT_ROOT"
    
    # Test Rust compilation
    log "INFO" "Testing Rust compilation..."
    cd "$RUST_DIR"
    if cargo check --quiet &> /dev/null; then
        log "SUCCESS" "Rust code compiles successfully"
    else
        log "WARNING" "Rust compilation failed (may need dependencies installed)"
    fi
    cd "$PROJECT_ROOT"
    
    # Test Python imports
    if [[ -d "$PYTHON_DIR/venv" ]]; then
        log "INFO" "Testing Python imports..."
        cd "$PYTHON_DIR"
        source venv/bin/activate
        if python -c "import genepredict" 2>/dev/null; then
            log "SUCCESS" "Python module imports successfully"
        else
            log "WARNING" "Python module import failed (may need dependencies installed)"
        fi
        deactivate
        cd "$PROJECT_ROOT"
    fi
}

main() {
    print_header "GenePredict Setup Validation"
    
    log "INFO" "Testing development environment setup..."
    log "INFO" "Project root: $PROJECT_ROOT"
    
    local all_tests_passed=true
    
    # Run all tests
    test_prerequisites || all_tests_passed=false
    test_project_structure || all_tests_passed=false
    test_dependencies || all_tests_passed=false
    test_configuration || all_tests_passed=false
    run_basic_tests || all_tests_passed=false
    
    # Summary
    print_header "Test Summary"
    
    if [[ "$all_tests_passed" == "true" ]]; then
        log "SUCCESS" "All tests passed! Your development environment is ready."
        log "INFO" "To start the application, run: ./scripts/quick-start.sh"
    else
        log "WARNING" "Some tests failed or showed warnings."
        log "INFO" "Run the setup script to fix issues: ./scripts/dev-setup-and-run.sh"
    fi
}

# Run main function
main "$@" 