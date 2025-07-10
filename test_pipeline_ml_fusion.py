#!/usr/bin/env python3
"""Test the complete pipeline to verify ML fusion is used."""

import sys
import os
import json
import logging

# Set up logging to see all INFO messages
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Add geneknow_pipeline to path
sys.path.insert(0, 'geneknow_pipeline')

from graph import run_pipeline

# Create a test VCF file with a known pathogenic variant
test_vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr17,length=81195210>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr17	41245466	.	G	A	100	PASS	.
"""

# Write test file
with open('test_ml_fusion.vcf', 'w') as f:
    f.write(test_vcf_content)

# Run pipeline
print("Running pipeline with test VCF...")
result = run_pipeline(
    'test_ml_fusion.vcf',
    user_preferences={
        'patient_data': {'age': 45, 'sex': 'F'},
        'language': 'en'
    }
)

# Check results
print("\n=== Pipeline Results ===")
print(f"Pipeline status: {result.get('pipeline_status')}")
print(f"Completed nodes: {result.get('completed_nodes', [])}")

# Check if ML fusion ran
if 'ml_fusion' in result.get('completed_nodes', []):
    print("\n✅ ML Fusion node executed!")
    
    # Check ML fusion results
    ml_fusion_results = result.get('ml_fusion_results', {})
    if ml_fusion_results.get('processing_successful'):
        print("✅ ML Fusion processing successful!")
        aggregate = ml_fusion_results.get('aggregate_risk_assessment', {})
        print(f"  - Risk score: {aggregate.get('aggregate_risk_score', 0):.3f}")
        print(f"  - Risk category: {aggregate.get('risk_category', 'unknown')}")
        print(f"  - Confidence: {aggregate.get('confidence', 0):.3f}")
    else:
        print(f"❌ ML Fusion failed: {ml_fusion_results.get('error', 'Unknown error')}")
else:
    print("\n❌ ML Fusion node did NOT execute!")

# Check if risk model used ML fusion
if 'ml_risk_assessment' in result:
    print("\n✅ Risk model used ML fusion results!")
    ml_assessment = result['ml_risk_assessment']
    print(f"  - Method: {ml_assessment.get('method', 'unknown')}")
    print(f"  - Risk category: {ml_assessment.get('risk_category', 'unknown')}")
else:
    print("\n⚠️ Risk model did NOT use ML fusion - used simple calculation")

# Check warnings
warnings = result.get('warnings', [])
for warning in warnings:
    if 'simple risk' in str(warning).lower():
        print(f"\n⚠️ WARNING: {warning}")

# Save full results for inspection
with open('test_ml_fusion_results.json', 'w') as f:
    # Convert any datetime objects to strings
    def json_serializable(obj):
        if hasattr(obj, 'isoformat'):
            return obj.isoformat()
        return str(obj)
    
    json.dump(result, f, indent=2, default=json_serializable)
    print("\nFull results saved to test_ml_fusion_results.json")

# Clean up
os.remove('test_ml_fusion.vcf') 