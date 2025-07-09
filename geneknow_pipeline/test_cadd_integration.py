#!/usr/bin/env python3
"""
Integration test for CADD scoring in the full pipeline.
Tests that CADD scores are properly added to variants during pipeline execution.
"""
import json
import os
import sys
from datetime import datetime

# Add current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from graph import run_pipeline

print("🧬 CADD Integration Test - Full Pipeline")
print("=" * 60)

# Test file path
test_file = "../test_data/test_sample.maf"

if not os.path.exists(test_file):
    print(f"❌ Test file not found: {test_file}")
    print("   Please run: gunzip -c test_data/tcga_downloads/*.maf.gz > test_data/test_sample.maf")
    sys.exit(1)

print(f"✓ Test file found: {test_file}")
print(f"  Size: {os.path.getsize(test_file) / 1024:.1f} KB")

# User preferences
preferences = {
    "language": "en",
    "include_technical": True,
    "patient_data": {
        "age": 50,
        "sex": "F"
    }
}

print("\n🚀 Running pipeline with CADD scoring...")
print("-" * 60)

# Run the pipeline
start_time = datetime.now()
result = run_pipeline(test_file, preferences)
end_time = datetime.now()

print("-" * 60)
print(f"\n⏱️  Pipeline completed in {(end_time - start_time).total_seconds():.2f} seconds")

# Check pipeline status
print(f"\n📊 Pipeline Status: {result.get('pipeline_status', 'unknown')}")
print(f"✓ Completed nodes: {', '.join(result.get('completed_nodes', []))}")

# Check for errors
errors = result.get('errors', [])
if errors:
    print(f"\n❌ Errors ({len(errors)}):")
    for error in errors[:5]:  # Show first 5 errors
        print(f"  - {error.get('node', 'unknown')}: {error.get('error', 'unknown error')}")
else:
    print("\n✅ No errors!")

# Check CADD integration
print("\n🔬 CADD Integration Check:")
print("-" * 40)

# Check if CADD node was executed
if 'cadd_scoring' in result.get('completed_nodes', []):
    print("✓ CADD scoring node executed")
else:
    print("❌ CADD scoring node not executed")

# Check CADD statistics
cadd_stats = result.get('cadd_stats', {})
if cadd_stats and 'error' not in cadd_stats:
    print("\n📊 CADD Statistics:")
    print(f"  - Total variants: {cadd_stats.get('total_variants', 0)}")
    print(f"  - Variants scored: {cadd_stats.get('variants_scored', 0)}")
    print(f"  - Lookup failures: {cadd_stats.get('lookup_missing', 0)}")
    print(f"  - Mean PHRED score: {cadd_stats.get('mean_phred', 0):.1f}")
    print(f"  - Max PHRED score: {cadd_stats.get('max_phred', 0):.1f}")
    print(f"  - Variants with PHRED > 20: {cadd_stats.get('variants_gt20', 0)}")
    
    # Show distribution
    dist = cadd_stats.get('phred_distribution', {})
    if dist:
        print("\n  PHRED Distribution:")
        print(f"    - Benign (<10): {dist.get('benign', 0)}")
        print(f"    - Uncertain (10-15): {dist.get('uncertain', 0)}")
        print(f"    - Damaging (15-20): {dist.get('damaging', 0)}")
        print(f"    - Pathogenic (>20): {dist.get('pathogenic', 0)}")
elif 'error' in cadd_stats:
    print(f"❌ CADD scoring error: {cadd_stats['error']}")
else:
    print("⚠️  No CADD statistics found")

# Check enriched variants
print("\n🧬 Variant Enrichment Check:")
filtered_variants = result.get('filtered_variants', [])
if filtered_variants:
    print(f"  Total filtered variants: {len(filtered_variants)}")
    
    # Check first few variants for CADD annotations
    cadd_annotated = 0
    sample_variants = []
    
    for i, variant in enumerate(filtered_variants[:5]):
        if 'cadd_phred' in variant:
            cadd_annotated += 1
            sample_variants.append({
                'id': variant.get('variant_id', f'Variant_{i}'),
                'gene': variant.get('gene', 'Unknown'),
                'cadd_phred': variant.get('cadd_phred'),
                'cadd_risk_weight': variant.get('cadd_risk_weight'),
                'risk_weight': variant.get('risk_weight')
            })
    
    # Count total annotated
    total_cadd_annotated = sum(1 for v in filtered_variants if 'cadd_phred' in v)
    
    print(f"\n  ✓ Variants with CADD scores: {total_cadd_annotated}/{len(filtered_variants)}")
    
    if sample_variants:
        print("\n  Sample CADD-annotated variants:")
        for v in sample_variants[:3]:
            print(f"    - {v['id']} ({v['gene']})")
            if v['cadd_phred'] is not None:
                print(f"      CADD PHRED: {v['cadd_phred']:.1f}")
                print(f"      CADD risk weight: {v['cadd_risk_weight']:.2f}")
            else:
                print(f"      CADD PHRED: Not found in database")
                print(f"      CADD risk weight: {v['cadd_risk_weight']:.2f} (default)")
            print(f"      Final risk weight: {v['risk_weight']:.2f}")
else:
    print("  ⚠️  No filtered variants found")

# Check feature vector builder
print("\n🔧 Feature Vector Builder Check:")
feature_vector = result.get('feature_vector', {})
if feature_vector:
    print(f"  Status: {feature_vector.get('status', 'unknown')}")
    inputs = feature_vector.get('available_inputs', {})
    if inputs:
        print("  Available inputs:")
        for input_type, available in inputs.items():
            status = "✓" if available else "✗"
            print(f"    {status} {input_type}")

# Save results
output_file = f"test_cadd_integration_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(output_file, 'w') as f:
    json.dump(result, f, indent=2, default=str)

print(f"\n💾 Full results saved to: {output_file}")
print("\n" + "=" * 60)
print("✅ CADD Integration Test Complete!")
print("=" * 60) 