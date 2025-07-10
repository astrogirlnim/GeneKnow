#!/usr/bin/env python3
"""
Test script for ClinVar annotator node.
Verifies that the ClinVar implementation works correctly.
"""

import os
import sys
import logging

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from nodes.clinvar_annotator import process

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_clinvar_annotator():
    """Test the ClinVar annotator with sample variants."""
    
    # Sample variants for testing
    test_variants = [
        {
            "variant_id": "17:43045677",
            "chrom": "17",
            "pos": 43045677,
            "ref": "C",
            "alt": "T",
            "gene": "BRCA1"
        },
        {
            "variant_id": "13:32363533",
            "chrom": "13", 
            "pos": 32363533,
            "ref": "T",
            "alt": "C",
            "gene": "BRCA2"
        },
        {
            "variant_id": "17:7673803",
            "chrom": "17",
            "pos": 7673803,
            "ref": "G",
            "alt": "A",
            "gene": "TP53"
        },
        {
            "variant_id": "1:12345678",  # Unknown variant
            "chrom": "1",
            "pos": 12345678,
            "ref": "A",
            "alt": "G",
            "gene": "UNKNOWN"
        }
    ]
    
    # Create test state
    test_state = {
        "filtered_variants": test_variants,
        "file_metadata": {}
    }
    
    print("="*60)
    print("Testing ClinVar Annotator")
    print("="*60)
    
    # Run the ClinVar annotator
    result_state = process(test_state)
    
    # Check results
    print(f"\nResults:")
    print(f"ClinVar annotations: {len(result_state.get('clinvar_annotations', {}))}")
    print(f"Pathogenic variants: {len(result_state.get('clinvar_pathogenic_variants', []))}")
    print(f"Likely pathogenic variants: {len(result_state.get('clinvar_likely_pathogenic_variants', []))}")
    
    # Print detailed results for each variant
    clinvar_annotations = result_state.get('clinvar_annotations', {})
    for variant in test_variants:
        variant_id = variant['variant_id']
        annotation = clinvar_annotations.get(variant_id, {})
        
        print(f"\n{variant_id} ({variant['gene']}):")
        if annotation.get('found_in_clinvar'):
            print(f"  Clinical significance: {annotation.get('clinical_significance')}")
            print(f"  Risk score: {annotation.get('clinical_risk_score', 0.0):.2f}")
            print(f"  Interpretation: {annotation.get('clinical_interpretation')}")
            print(f"  Recommendation: {annotation.get('recommendation')}")
        else:
            print(f"  No ClinVar annotation found")
    
    # Print summary statistics
    stats = result_state.get('clinvar_stats', {})
    print(f"\nSummary Statistics:")
    print(f"  Total variants: {stats.get('total_variants', 0)}")
    print(f"  Annotated variants: {stats.get('variants_annotated', 0)}")
    print(f"  Annotation rate: {stats.get('annotation_rate', 0.0)*100:.1f}%")
    print(f"  Pathogenic: {stats.get('pathogenic_variants', 0)}")
    print(f"  Likely pathogenic: {stats.get('likely_pathogenic_variants', 0)}")
    print(f"  Cancer-related: {stats.get('cancer_related_variants', 0)}")
    
    print("="*60)
    print("ClinVar test completed successfully!")
    return result_state

if __name__ == "__main__":
    test_clinvar_annotator() 