#!/usr/bin/env python3
"""Test the complete pipeline with a MAF file to see all nodes execute."""

import sys
import os
import json
import logging

# Set up comprehensive logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('test_maf_pipeline.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

# Add geneknow_pipeline to path
sys.path.insert(0, 'geneknow_pipeline')

from graph import run_pipeline

# Use the test MAF file mentioned in the user's output
test_maf = 'test_data/tcga_downloads/3d14b1e2-0555-4d6f-a55b-a56065f915e1.wxs.aliquot_ensemble_masked.maf.gz'

# Check if file exists
if not os.path.exists(test_maf):
    print(f"Test MAF file not found: {test_maf}")
    # Create a simple test MAF
    test_maf = 'test_simple.maf'
    maf_content = """Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode
BRCA1\tchr17\t41244936\t41244936\tMissense_Mutation\tSNP\tG\tG\tA\tTEST-SAMPLE-001
TP53\tchr17\t7577121\t7577121\tMissense_Mutation\tSNP\tG\tG\tA\tTEST-SAMPLE-001
"""
    with open(test_maf, 'w') as f:
        f.write(maf_content)

print(f"Running pipeline with MAF file: {test_maf}")
print("=" * 80)

# Run pipeline
result = run_pipeline(
    test_maf,
    user_preferences={
        'patient_data': {'age': 45, 'sex': 'F'},
        'language': 'en',
        'include_technical': True
    }
)

# Analyze results
print("\n=== Pipeline Results ===")
print(f"Pipeline status: {result.get('pipeline_status')}")
print(f"Processing time: {result.get('processing_time_seconds', 0):.2f} seconds")
print(f"\nCompleted nodes ({len(result.get('completed_nodes', []))}):")
for node in result.get('completed_nodes', []):
    print(f"  - {node}")

# Check key results
print("\n=== Key Results ===")

# Check variant processing
print(f"Total variants: {result.get('variant_count', 0)}")
print(f"Filtered variants: {len(result.get('filtered_variants', []))}")

# Check if ML models ran
if 'ml_fusion_results' in result:
    ml_fusion = result['ml_fusion_results']
    print(f"\n✅ ML Fusion Results Found:")
    print(f"  - Processing successful: {ml_fusion.get('processing_successful', False)}")
    if ml_fusion.get('processing_successful'):
        aggregate = ml_fusion.get('aggregate_risk_assessment', {})
        print(f"  - Risk score: {aggregate.get('aggregate_risk_score', 0):.3f}")
        print(f"  - Risk category: {aggregate.get('risk_category', 'unknown')}")
    else:
        print(f"  - Error: {ml_fusion.get('error', 'Unknown')}")
else:
    print("\n❌ No ML Fusion Results Found")

# Check if static models ran
print("\n=== Static Model Results ===")
print(f"TCGA matches: {'✅' if result.get('tcga_matches') is not None else '❌'}")
print(f"CADD stats: {'✅' if result.get('cadd_stats') is not None else '❌'}")
print(f"ClinVar annotations: {'✅' if result.get('clinvar_annotations') is not None else '❌'}")
print(f"PRS results: {'✅' if result.get('prs_results') is not None else '❌'}")
print(f"Pathway burden: {'✅' if result.get('pathway_burden_results') is not None else '❌'}")

# Check risk scores
risk_scores = result.get('risk_scores', {})
if risk_scores:
    print(f"\n=== Risk Scores ===")
    for cancer, score in risk_scores.items():
        print(f"  - {cancer}: {score}%")

# Check for ML risk assessment
if 'ml_risk_assessment' in result:
    print(f"\n✅ ML Risk Assessment Found:")
    ml_assessment = result['ml_risk_assessment']
    print(f"  - Method: {ml_assessment.get('method', 'unknown')}")
    print(f"  - Risk category: {ml_assessment.get('risk_category', 'unknown')}")
else:
    print(f"\n❌ No ML Risk Assessment - Using simple calculation")

# Check warnings
warnings = result.get('warnings', [])
print(f"\n=== Warnings ({len(warnings)}) ===")
for warning in warnings:
    print(f"  - {warning}")

# Check errors
errors = result.get('errors', [])
print(f"\n=== Errors ({len(errors)}) ===")
for error in errors:
    print(f"  - {error}")

# Save full results
with open('test_maf_pipeline_results.json', 'w') as f:
    def json_serializable(obj):
        if hasattr(obj, 'isoformat'):
            return obj.isoformat()
        return str(obj)
    
    json.dump(result, f, indent=2, default=json_serializable)
    print(f"\nFull results saved to test_maf_pipeline_results.json")

# Clean up test file if we created it
if os.path.exists('test_simple.maf'):
    os.remove('test_simple.maf') 