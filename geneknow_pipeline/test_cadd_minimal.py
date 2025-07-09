#!/usr/bin/env python3
"""
Minimal test for CADD scoring functionality.
Tests core functionality without dependencies.
"""
import os
import sys
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)

print("="*60)
print("CADD Scoring Minimal Test")
print("="*60)

# Test imports
try:
    from nodes import cadd_scoring
    print("✓ Successfully imported cadd_scoring module")
except ImportError as e:
    print(f"✗ Failed to import cadd_scoring: {e}")
    # Try direct import
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nodes.cadd_scoring as cadd_scoring
    print("✓ Successfully imported via direct path")

# Test database existence
db_path = os.path.join(os.path.dirname(__file__), "data", "cadd_scores.db")
if os.path.exists(db_path):
    print(f"✓ CADD database exists at: {db_path}")
    print(f"  Size: {os.path.getsize(db_path) / 1024:.1f} KB")
else:
    print(f"✗ CADD database not found at: {db_path}")

# Test basic functions
print("\nTesting basic functions:")

# Test chromosome normalization
print("- Testing normalize_chromosome:")
test_cases = [('1', 'chr1'), ('chr1', 'chr1'), ('X', 'chrX')]
for input_chr, expected in test_cases:
    result = cadd_scoring.normalize_chromosome(input_chr)
    status = "✓" if result == expected else "✗"
    print(f"  {status} normalize_chromosome('{input_chr}') = '{result}' (expected: '{expected}')")

# Test risk weight calculation
print("\n- Testing calculate_risk_weight:")
test_scores = [(5.0, 0.1), (15.0, 0.3), (25.0, 0.8), (35.0, 1.0)]
for score, expected_range in test_scores:
    result = cadd_scoring.calculate_risk_weight(score)
    status = "✓" if abs(result - expected_range) < 0.2 else "✗"
    print(f"  {status} PHRED {score} -> risk weight {result:.2f}")

# Test process function with minimal state
print("\n- Testing process function:")
test_state = {
    "current_node": None,
    "completed_nodes": [],
    "errors": [],
    "filtered_variants": [
        {
            "variant_id": "test_variant_1",
            "chrom": "chr17",
            "pos": 43044295,
            "ref": "A",
            "alt": "T",
            "gene": "BRCA1"
        }
    ]
}

try:
    result = cadd_scoring.process(test_state)
    print("  ✓ Process function executed successfully")
    
    if "cadd_scoring" in result.get("completed_nodes", []):
        print("  ✓ Node marked as completed")
    
    if "cadd_stats" in result:
        stats = result["cadd_stats"]
        print(f"  ✓ Statistics generated:")
        print(f"    - Total variants: {stats.get('total_variants', 0)}")
        print(f"    - Variants scored: {stats.get('variants_scored', 0)}")
        print(f"    - Lookup failures: {stats.get('lookup_missing', 0)}")
    
    if result.get("filtered_variants") and len(result["filtered_variants"]) > 0:
        var = result["filtered_variants"][0]
        if "cadd_phred" in var:
            print(f"  ✓ CADD score found: PHRED={var['cadd_phred']}")
        else:
            print(f"  ℹ No CADD score found (expected if not in test DB)")
            
except Exception as e:
    print(f"  ✗ Process function failed: {e}")

print("\n" + "="*60)
print("Test complete!")
print("="*60) 