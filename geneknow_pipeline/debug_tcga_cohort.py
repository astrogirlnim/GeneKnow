#!/usr/bin/env python3
"""
Debug script to check if tcga_cohort_sizes is being set properly.
"""
from graph import run_pipeline

# Run the pipeline directly
test_file = "test_data/tcga_downloads/3d14b1e2-0555-4d6f-a55b-a56065f915e1.wxs.aliquot_ensemble_masked.maf.gz"
result = run_pipeline(test_file, {"language": "en"})

print("🔍 Checking for tcga_cohort_sizes in pipeline result:")
print(f"   - tcga_cohort_sizes present: {'tcga_cohort_sizes' in result}")

if 'tcga_cohort_sizes' in result:
    print(f"   - Content: {result['tcga_cohort_sizes']}")
else:
    print("   - Not found in result")
    
# Check if it's in file_metadata
if 'file_metadata' in result and 'tcga_summary' in result['file_metadata']:
    tcga_summary = result['file_metadata']['tcga_summary']
    if 'cohort_sizes' in tcga_summary:
        print("\n✅ Found cohort_sizes in file_metadata.tcga_summary:")
        print(f"   - Content: {tcga_summary['cohort_sizes']}")
        
print("\n📋 All top-level keys in result:")
for key in sorted(result.keys()):
    print(f"   - {key}") 