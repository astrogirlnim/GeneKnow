#!/bin/bash

# GenePredict CLI Test Script
# This script provides comprehensive testing for all components of GenePredict

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
PROJECT_ROOT=$(pwd)
TEST_DATA_DIR="$PROJECT_ROOT/data/test"
FRONTEND_DIR="$PROJECT_ROOT/frontend"
RUST_DIR="$PROJECT_ROOT/backend/rust"
PYTHON_DIR="$PROJECT_ROOT/backend/python"

# Test flags
RUN_FRONTEND_TESTS=true
RUN_RUST_TESTS=true
RUN_PYTHON_TESTS=true
RUN_INTEGRATION_TESTS=true
RUN_PERFORMANCE_TESTS=false
VERBOSE=false

# Print functions
print_header() {
    echo -e "${BLUE}======================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}======================================${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

# Help function
show_help() {
    echo "GenePredict CLI Test Script"
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help                  Show this help message"
    echo "  -v, --verbose               Enable verbose output"
    echo "  -f, --frontend-only         Run only frontend tests"
    echo "  -r, --rust-only             Run only Rust tests"
    echo "  -p, --python-only           Run only Python tests"
    echo "  -i, --integration-only      Run only integration tests"
    echo "  -P, --performance           Include performance tests"
    echo "  --no-frontend              Skip frontend tests"
    echo "  --no-rust                  Skip Rust tests"
    echo "  --no-python                Skip Python tests"
    echo "  --no-integration           Skip integration tests"
    echo "  --setup-only               Only set up test environment"
    echo "  --cleanup-only             Only clean up test environment"
    echo ""
    echo "Examples:"
    echo "  $0                         # Run all tests"
    echo "  $0 -f                      # Run only frontend tests"
    echo "  $0 -v -P                   # Run all tests with verbose output and performance tests"
    echo "  $0 --no-frontend           # Run all tests except frontend"
    echo ""
}

# Parse command line arguments
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
        -f|--frontend-only)
            RUN_FRONTEND_TESTS=true
            RUN_RUST_TESTS=false
            RUN_PYTHON_TESTS=false
            RUN_INTEGRATION_TESTS=false
            shift
            ;;
        -r|--rust-only)
            RUN_FRONTEND_TESTS=false
            RUN_RUST_TESTS=true
            RUN_PYTHON_TESTS=false
            RUN_INTEGRATION_TESTS=false
            shift
            ;;
        -p|--python-only)
            RUN_FRONTEND_TESTS=false
            RUN_RUST_TESTS=false
            RUN_PYTHON_TESTS=true
            RUN_INTEGRATION_TESTS=false
            shift
            ;;
        -i|--integration-only)
            RUN_FRONTEND_TESTS=false
            RUN_RUST_TESTS=false
            RUN_PYTHON_TESTS=false
            RUN_INTEGRATION_TESTS=true
            shift
            ;;
        -P|--performance)
            RUN_PERFORMANCE_TESTS=true
            shift
            ;;
        --no-frontend)
            RUN_FRONTEND_TESTS=false
            shift
            ;;
        --no-rust)
            RUN_RUST_TESTS=false
            shift
            ;;
        --no-python)
            RUN_PYTHON_TESTS=false
            shift
            ;;
        --no-integration)
            RUN_INTEGRATION_TESTS=false
            shift
            ;;
        --setup-only)
            setup_test_environment
            exit 0
            ;;
        --cleanup-only)
            cleanup_test_environment
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Environment setup
setup_test_environment() {
    print_header "Setting up test environment"
    
    # Create test data directory
    mkdir -p "$TEST_DATA_DIR"
    
    # Create sample VCF file
    cat > "$TEST_DATA_DIR/sample.vcf" << EOF
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=242193529>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	100	rs1	A	T	30	PASS	AF=0.5	GT:DP	0/1:20
1	200	rs2	C	G	40	PASS	AF=0.3	GT:DP	1/1:25
2	300	rs3	G	A	50	PASS	AF=0.7	GT:DP	0/0:15
EOF
    
    # Create sample BAM file header (mock)
    cat > "$TEST_DATA_DIR/sample.bam.header" << EOF
@HD	VN:1.6	SO:coordinate
@SQ	SN:1	LN:249250621
@SQ	SN:2	LN:242193529
@PG	ID:bwa	PN:bwa	VN:0.7.17
EOF
    
    # Create sample FASTQ file
    cat > "$TEST_DATA_DIR/sample.fastq" << EOF
@read1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@read3
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
EOF
    
    # Create invalid test files
    echo "This is not a valid VCF file" > "$TEST_DATA_DIR/invalid.vcf"
    echo "This is not a valid BAM file" > "$TEST_DATA_DIR/invalid.bam"
    echo "This is not a valid FASTQ file" > "$TEST_DATA_DIR/invalid.fastq"
    
    # Create large test file for performance testing
    if [[ "$RUN_PERFORMANCE_TESTS" == "true" ]]; then
        print_info "Creating large test files for performance testing..."
        
        # Create large VCF file
        cat > "$TEST_DATA_DIR/large.vcf" << EOF
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
EOF
        
        # Generate 10000 variants
        for i in {1..10000}; do
            echo -e "1\t$((i * 1000))\trs$i\tA\tT\t30\tPASS\tAF=0.5\tGT\t0/1" >> "$TEST_DATA_DIR/large.vcf"
        done
    fi
    
    print_success "Test environment setup complete"
}

# Environment cleanup
cleanup_test_environment() {
    print_header "Cleaning up test environment"
    
    # Remove test files
    if [[ -d "$TEST_DATA_DIR" ]]; then
        rm -rf "$TEST_DATA_DIR"
        print_success "Test data directory removed"
    fi
    
    # Clean up any temporary files
    find "$PROJECT_ROOT" -name "*.tmp" -delete 2>/dev/null || true
    find "$PROJECT_ROOT" -name ".DS_Store" -delete 2>/dev/null || true
    
    print_success "Test environment cleanup complete"
}

# Check prerequisites
check_prerequisites() {
    print_header "Checking prerequisites"
    
    local missing_tools=()
    
    # Check Node.js
    if ! command -v node &> /dev/null; then
        missing_tools+=("node")
    else
        print_success "Node.js $(node --version) found"
    fi
    
    # Check npm
    if ! command -v npm &> /dev/null; then
        missing_tools+=("npm")
    else
        print_success "npm $(npm --version) found"
    fi
    
    # Check Rust
    if ! command -v cargo &> /dev/null; then
        missing_tools+=("cargo")
    else
        print_success "Rust $(rustc --version) found"
    fi
    
    # Check Python
    if ! command -v python3 &> /dev/null; then
        missing_tools+=("python3")
    else
        print_success "Python $(python3 --version) found"
    fi
    
    # Check pip
    if ! command -v pip &> /dev/null && ! command -v pip3 &> /dev/null; then
        missing_tools+=("pip")
    else
        print_success "pip found"
    fi
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        print_error "Missing required tools: ${missing_tools[*]}"
        print_info "Please install missing tools and try again"
        exit 1
    fi
    
    print_success "All prerequisites satisfied"
}

# Run frontend tests
run_frontend_tests() {
    print_header "Running Frontend Tests"
    
    cd "$FRONTEND_DIR"
    
    # Install dependencies if needed
    if [[ ! -d "node_modules" ]]; then
        print_info "Installing frontend dependencies..."
        npm install
    fi
    
    # Run linting
    print_info "Running ESLint..."
    if npm run lint; then
        print_success "Frontend linting passed"
    else
        print_error "Frontend linting failed"
        return 1
    fi
    
    # Run type checking
    print_info "Running TypeScript type check..."
    if npm run type-check; then
        print_success "TypeScript type check passed"
    else
        print_error "TypeScript type check failed"
        return 1
    fi
    
    # Run unit tests
    print_info "Running frontend unit tests..."
    if npm test -- --watchAll=false; then
        print_success "Frontend unit tests passed"
    else
        print_error "Frontend unit tests failed"
        return 1
    fi
    
    # Run build test
    print_info "Testing frontend build..."
    if npm run build; then
        print_success "Frontend build test passed"
    else
        print_error "Frontend build test failed"
        return 1
    fi
    
    cd "$PROJECT_ROOT"
    print_success "All frontend tests passed"
}

# Run Rust tests
run_rust_tests() {
    print_header "Running Rust Tests"
    
    cd "$RUST_DIR"
    
    # Run clippy
    print_info "Running Clippy..."
    if cargo clippy -- -D warnings; then
        print_success "Clippy check passed"
    else
        print_error "Clippy check failed"
        return 1
    fi
    
    # Run formatting check
    print_info "Running rustfmt check..."
    if cargo fmt -- --check; then
        print_success "Rust formatting check passed"
    else
        print_error "Rust formatting check failed"
        return 1
    fi
    
    # Run unit tests
    print_info "Running Rust unit tests..."
    if cargo test; then
        print_success "Rust unit tests passed"
    else
        print_error "Rust unit tests failed"
        return 1
    fi
    
    # Run security audit
    print_info "Running security audit..."
    if cargo audit; then
        print_success "Security audit passed"
    else
        print_warning "Security audit found issues (check manually)"
    fi
    
    cd "$PROJECT_ROOT"
    print_success "All Rust tests passed"
}

# Run Python tests
run_python_tests() {
    print_header "Running Python Tests"
    
    cd "$PYTHON_DIR"
    
    # Install dependencies if needed
    if [[ ! -d "venv" ]]; then
        print_info "Creating Python virtual environment..."
        python3 -m venv venv
    fi
    
    source venv/bin/activate
    
    print_info "Installing Python dependencies..."
    pip install -r requirements.txt
    pip install -e .
    
    # Run linting
    print_info "Running Python linting..."
    if python -m flake8 genepredict/; then
        print_success "Python linting passed"
    else
        print_error "Python linting failed"
        return 1
    fi
    
    # Run type checking
    print_info "Running mypy type check..."
    if python -m mypy genepredict/; then
        print_success "Python type check passed"
    else
        print_warning "Python type check found issues (check manually)"
    fi
    
    # Run unit tests
    print_info "Running Python unit tests..."
    if python -m pytest tests/ -v; then
        print_success "Python unit tests passed"
    else
        print_error "Python unit tests failed"
        return 1
    fi
    
    # Run coverage
    print_info "Running test coverage..."
    if python -m pytest tests/ --cov=genepredict --cov-report=term-missing; then
        print_success "Test coverage analysis complete"
    else
        print_warning "Test coverage analysis failed"
    fi
    
    deactivate
    cd "$PROJECT_ROOT"
    print_success "All Python tests passed"
}

# Run integration tests
run_integration_tests() {
    print_header "Running Integration Tests"
    
    # Test file processing pipeline
    print_info "Testing file processing pipeline..."
    
    # Test VCF processing
    print_info "Testing VCF file processing..."
    if test_vcf_processing; then
        print_success "VCF processing test passed"
    else
        print_error "VCF processing test failed"
        return 1
    fi
    
    # Test BAM processing
    print_info "Testing BAM file processing..."
    if test_bam_processing; then
        print_success "BAM processing test passed"
    else
        print_error "BAM processing test failed"
        return 1
    fi
    
    # Test FASTQ processing
    print_info "Testing FASTQ file processing..."
    if test_fastq_processing; then
        print_success "FASTQ processing test passed"
    else
        print_error "FASTQ processing test failed"
        return 1
    fi
    
    # Test error handling
    print_info "Testing error handling..."
    if test_error_handling; then
        print_success "Error handling test passed"
    else
        print_error "Error handling test failed"
        return 1
    fi
    
    print_success "All integration tests passed"
}

# Test VCF processing
test_vcf_processing() {
    local test_file="$TEST_DATA_DIR/sample.vcf"
    
    # Test Python VCF processor
    cd "$PYTHON_DIR"
    source venv/bin/activate
    
    python -c "
from genepredict.processors.vcf import VCFProcessor
processor = VCFProcessor()
try:
    result = processor.process('$test_file')
    print('VCF processing successful')
    print(f'Found {len(result.variants)} variants')
    exit(0)
except Exception as e:
    print(f'VCF processing failed: {e}')
    exit(1)
"
    
    local result=$?
    deactivate
    cd "$PROJECT_ROOT"
    return $result
}

# Test BAM processing
test_bam_processing() {
    local test_file="$TEST_DATA_DIR/sample.bam.header"
    
    # Test Python BAM processor
    cd "$PYTHON_DIR"
    source venv/bin/activate
    
    python -c "
from genepredict.processors.bam import BAMProcessor
processor = BAMProcessor()
try:
    # Test with header file (mock BAM)
    result = processor.process('$test_file')
    print('BAM processing successful')
    exit(0)
except Exception as e:
    print(f'BAM processing failed: {e}')
    exit(1)
"
    
    local result=$?
    deactivate
    cd "$PROJECT_ROOT"
    return $result
}

# Test FASTQ processing
test_fastq_processing() {
    local test_file="$TEST_DATA_DIR/sample.fastq"
    
    # Test Python FASTQ processor
    cd "$PYTHON_DIR"
    source venv/bin/activate
    
    python -c "
from genepredict.processors.fastq import FASTQProcessor
processor = FASTQProcessor()
try:
    result = processor.process('$test_file')
    print('FASTQ processing successful')
    print(f'Found {len(result.sequences)} sequences')
    exit(0)
except Exception as e:
    print(f'FASTQ processing failed: {e}')
    exit(1)
"
    
    local result=$?
    deactivate
    cd "$PROJECT_ROOT"
    return $result
}

# Test error handling
test_error_handling() {
    local test_file="$TEST_DATA_DIR/invalid.vcf"
    
    # Test Python error handling
    cd "$PYTHON_DIR"
    source venv/bin/activate
    
    python -c "
from genepredict.processors.vcf import VCFProcessor
processor = VCFProcessor()
try:
    result = processor.process('$test_file')
    print('ERROR: Should have failed on invalid file')
    exit(1)
except Exception as e:
    print(f'Error handling successful: {e}')
    exit(0)
"
    
    local result=$?
    deactivate
    cd "$PROJECT_ROOT"
    return $result
}

# Run performance tests
run_performance_tests() {
    print_header "Running Performance Tests"
    
    if [[ ! -f "$TEST_DATA_DIR/large.vcf" ]]; then
        print_warning "Large test file not found, skipping performance tests"
        return 0
    fi
    
    print_info "Testing large file processing performance..."
    
    cd "$PYTHON_DIR"
    source venv/bin/activate
    
    local start_time=$(date +%s)
    
    python -c "
import time
from genepredict.processors.vcf import VCFProcessor

processor = VCFProcessor()
start_time = time.time()
try:
    result = processor.process('$TEST_DATA_DIR/large.vcf')
    end_time = time.time()
    processing_time = end_time - start_time
    print(f'Large file processing completed in {processing_time:.2f} seconds')
    print(f'Processed {len(result.variants)} variants')
    print(f'Processing rate: {len(result.variants)/processing_time:.2f} variants/second')
    exit(0)
except Exception as e:
    print(f'Performance test failed: {e}')
    exit(1)
"
    
    local result=$?
    deactivate
    cd "$PROJECT_ROOT"
    
    if [[ $result -eq 0 ]]; then
        print_success "Performance tests passed"
    else
        print_error "Performance tests failed"
        return 1
    fi
}

# Main test execution
main() {
    print_header "GenePredict Test Suite"
    
    # Check prerequisites
    check_prerequisites
    
    # Setup test environment
    setup_test_environment
    
    # Track test results
    local test_results=()
    
    # Run tests based on flags
    if [[ "$RUN_FRONTEND_TESTS" == "true" ]]; then
        if run_frontend_tests; then
            test_results+=("Frontend: PASSED")
        else
            test_results+=("Frontend: FAILED")
        fi
    fi
    
    if [[ "$RUN_RUST_TESTS" == "true" ]]; then
        if run_rust_tests; then
            test_results+=("Rust: PASSED")
        else
            test_results+=("Rust: FAILED")
        fi
    fi
    
    if [[ "$RUN_PYTHON_TESTS" == "true" ]]; then
        if run_python_tests; then
            test_results+=("Python: PASSED")
        else
            test_results+=("Python: FAILED")
        fi
    fi
    
    if [[ "$RUN_INTEGRATION_TESTS" == "true" ]]; then
        if run_integration_tests; then
            test_results+=("Integration: PASSED")
        else
            test_results+=("Integration: FAILED")
        fi
    fi
    
    if [[ "$RUN_PERFORMANCE_TESTS" == "true" ]]; then
        if run_performance_tests; then
            test_results+=("Performance: PASSED")
        else
            test_results+=("Performance: FAILED")
        fi
    fi
    
    # Print summary
    print_header "Test Results Summary"
    
    local failed_tests=0
    for result in "${test_results[@]}"; do
        if [[ $result == *"PASSED"* ]]; then
            print_success "$result"
        else
            print_error "$result"
            ((failed_tests++))
        fi
    done
    
    # Clean up test environment
    cleanup_test_environment
    
    if [[ $failed_tests -eq 0 ]]; then
        print_success "All tests passed! 🎉"
        exit 0
    else
        print_error "$failed_tests test suite(s) failed"
        exit 1
    fi
}

# Run main function
main "$@" 