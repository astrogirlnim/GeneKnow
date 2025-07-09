# 🧬 CADD Scoring Testing Guide

> Comprehensive testing framework for the CADD (Combined Annotation Dependent Depletion) scoring node implementation in the GeneKnow genomic pipeline.

## Overview

This testing framework validates all aspects of the CADD scoring implementation as specified in `docs/CADD_Model_Node_Implementation_Plan.md`. It includes unit tests, integration tests, performance benchmarks, memory profiling, and configuration testing.

## Test Suite Structure

### 1. **Unit Tests** (`test_cadd_scoring.py`)
Tests core CADD functionality in isolation:
- Chromosome normalization (chr1 → 1)
- CADD score database lookups
- Risk weight calculations from PHRED scores
- Job tracking and status updates
- Main process function with mock data

### 2. **Minimal Tests** (`test_cadd_minimal.py`)
Quick smoke tests without external dependencies:
- Module imports
- Basic function calls
- Database existence checks
- Minimal state processing

### 3. **Integration Tests** (`test_cadd_integration.py`)
Tests CADD scoring within the full pipeline:
- End-to-end pipeline execution with MAF files
- CADD annotation propagation through nodes
- Statistics collection and reporting
- Feature vector builder integration

### 4. **Comprehensive Tests** (`test_cadd_comprehensive.py`)
Advanced test scenarios:
- **Performance Tests**: Validates <200ms per 10k lookups requirement
- **Memory Tests**: Ensures SQLite cache stays under 200MB
- **Edge Cases**: Malformed data, missing databases, etc.
- **Remote Fallback**: Tests Tabix query structure (when implemented)

### 5. **Configuration Tests** (`test_cadd_config.py`)
Environment variable and configuration handling:
- Default vs custom CADD endpoints
- Enabling/disabling remote lookups
- Custom database paths
- Legacy risk model routing
- PHRED threshold validation

## Running the Tests

### Run All Tests
```bash
cd geneknow_pipeline
python run_cadd_tests.py
```

This generates a comprehensive report with:
- Overall pass/fail statistics
- Performance metrics
- Database validation results
- Detailed JSON report file

### Run Individual Test Suites
```bash
# Unit tests only
python test_cadd_scoring.py

# Quick minimal tests
python test_cadd_minimal.py

# Integration tests
python test_cadd_integration.py

# Comprehensive tests (includes performance)
python test_cadd_comprehensive.py

# Configuration tests
python test_cadd_config.py
```

### Run Specific Test Classes
```bash
# Performance tests only
python -m unittest test_cadd_comprehensive.TestCADDPerformance

# Memory tests only  
python -m unittest test_cadd_comprehensive.TestCADDMemory

# Edge case tests
python -m unittest test_cadd_comprehensive.TestCADDEdgeCases
```

## Test Requirements

### Prerequisites
- Python 3.8+
- Required packages: `psutil` (for memory profiling)
- SQLite3
- Population variants database (or tests will use temporary DBs)

### Environment Variables
Tests respect these environment variables:
- `CADD_DB_PATH`: Custom database location
- `CADD_REMOTE_TABIX`: Remote CADD service endpoint
- `USE_REMOTE_CADD`: Enable/disable remote lookups ("true"/"false")
- `USE_LEGACY_RISK`: Use old risk model path ("true"/"false")

## Performance Benchmarks

### Database Lookup Performance
- **Requirement**: <200ms per 10,000 lookups
- **Test Method**: Creates 100k variant DB, times 10k sequential lookups
- **Expected**: ~0.02ms per lookup with SQLite indexes

### Memory Usage
- **Requirement**: <200MB SQLite page cache
- **Test Method**: 50k variants, 10k lookups, measure RSS increase
- **Expected**: <50MB typical increase

## Test Data

### Mock Variants
Tests use realistic variant data including:
- Cancer genes: BRCA1, BRCA2, TP53, EGFR, etc.
- PHRED scores: Full range from benign (<10) to pathogenic (>25)
- Edge cases: Sex chromosomes (X,Y), mitochondrial (MT)
- Malformed data: Missing fields, invalid types

### Risk Weight Mapping
```
PHRED < 10:    0.1 (benign)
PHRED 10-15:   0.1-0.3 (uncertain)  
PHRED 15-20:   0.3-0.6 (damaging)
PHRED 20-25:   0.6-0.8 (pathogenic)
PHRED > 25:    0.8-1.0 (highly pathogenic)
```

## Debugging Failed Tests

### Common Issues

1. **Database Not Found**
   - Check `population_variants.db` exists
   - Verify `cadd_scores` and `cadd_jobs` tables present
   - Run database creation script if needed

2. **Import Errors**
   - Ensure running from `geneknow_pipeline/` directory
   - Check Python path includes parent directory

3. **Performance Test Failures**
   - May fail on slow systems or under high load
   - Check no other processes accessing database
   - Verify SQLite indexes created properly

4. **Memory Test Failures**
   - Close other memory-intensive applications
   - Check system has sufficient RAM (>2GB recommended)

### Verbose Output
```bash
# Run with verbose logging
python test_cadd_comprehensive.py -v

# Enable debug logging
export LOG_LEVEL=DEBUG
python run_cadd_tests.py
```

## CI/CD Integration

### GitHub Actions Example
```yaml
- name: Run CADD Tests
  run: |
    cd geneknow_pipeline
    pip install psutil
    python run_cadd_tests.py
  env:
    USE_REMOTE_CADD: "false"  # Disable remote calls in CI
```

### Pre-commit Hook
```bash
#!/bin/bash
cd geneknow_pipeline
python test_cadd_minimal.py || exit 1
```

## Test Coverage Goals

- **Unit Test Coverage**: >90% of cadd_scoring.py functions
- **Integration Coverage**: All pipeline paths with CADD enabled
- **Edge Case Coverage**: All error conditions handled gracefully
- **Performance Coverage**: Meets all stated requirements

## Future Test Additions

1. **Remote Tabix Tests** (when implemented)
   - Mock HTTP responses
   - Network timeout handling
   - Cache miss → remote → cache hit flow

2. **Concurrent Access Tests**
   - Multiple processes reading database
   - Write lock handling
   - Connection pool behavior

3. **Large-scale Integration**
   - Full genome VCF processing
   - TCGA cohort analysis
   - Production workload simulation

## Maintenance

- Run full test suite before any CADD-related changes
- Update tests when modifying risk weight calculations
- Add new edge cases as discovered in production
- Keep performance benchmarks aligned with requirements

---

*Last Updated: Based on CADD Model Node Implementation Plan v0.1* 