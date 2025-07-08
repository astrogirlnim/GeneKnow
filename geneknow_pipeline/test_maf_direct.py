#!/usr/bin/env python3
"""
Direct test of MAF processing without API server.
"""
from graph import run_pipeline
import json
import os

# Use the test MAF file we created
maf_path = "test_data/test_sample.maf"

if not os.path.exists(maf_path):
    print(f"❌ MAF file not found: {maf_path}")
    print("   Run test_maf_simple.py first to create it")
    exit(1)

print(f"🧬 Testing MAF processing directly...")
print(f"   File: {maf_path}")

try:
    # Run pipeline directly
    result = run_pipeline(
        file_path=maf_path,
        user_preferences={
            "language": "en",
            "include_technical": True
        }
    )
    
    # Check results
    if result["pipeline_status"] == "completed":
        print("✅ Pipeline completed successfully!")
        
        print(f"\n📊 Results:")
        print(f"   Total variants: {result.get('variant_count', 0)}")
        print(f"   Completed nodes: {result.get('completed_nodes', [])}")
        
        # Check MAF-specific metadata
        maf_info = result.get("file_metadata", {}).get("maf_info", {})
        if maf_info:
            print(f"\n   MAF Analysis:")
            print(f"     - Unique genes: {maf_info.get('unique_genes', 0)}")
            print(f"     - Top mutated genes: {list(maf_info.get('top_mutated_genes', {}).keys())[:5]}")
            print(f"     - Variant classifications: {list(maf_info.get('variant_classifications', {}).keys())[:5]}")
        
        # Check risk scores
        risk_scores = result.get("risk_scores", {})
        if risk_scores:
            print(f"\n   Cancer Risk Scores:")
            for cancer, score in risk_scores.items():
                print(f"     - {cancer}: {score:.1f}%")
        
        # Check TCGA matches
        tcga_matches = result.get("tcga_matches", {})
        if tcga_matches:
            print(f"\n   TCGA Analysis:")
            for cancer_type, matches in tcga_matches.items():
                matched_count = sum(1 for v in matches.values() if v.get('frequency', 0) > 0)
                if matched_count > 0:
                    print(f"     - {cancer_type}: {matched_count} variants matched")
        
        # Check if variants were loaded
        filtered_variants = result.get("filtered_variants", [])
        print(f"\n   Filtered variants: {len(filtered_variants)}")
        if filtered_variants:
            print("   First variant:")
            first_var = filtered_variants[0]
            print(f"     - Gene: {first_var.get('gene')}")
            print(f"     - Position: {first_var.get('chrom')}:{first_var.get('pos')}")
            print(f"     - Change: {first_var.get('ref')} → {first_var.get('alt')}")
        
        # Save detailed results
        with open("test_data/maf_direct_results.json", 'w') as f:
            json.dump(result, f, indent=2, default=str)
        print(f"\n💾 Full results saved to: test_data/maf_direct_results.json")
        
    else:
        print(f"❌ Pipeline failed: {result.get('pipeline_status')}")
        errors = result.get('errors', [])
        for error in errors:
            print(f"   Error in {error.get('node')}: {error.get('error')}")
        
except Exception as e:
    print(f"❌ Pipeline execution error: {e}")
    import traceback
    traceback.print_exc() 