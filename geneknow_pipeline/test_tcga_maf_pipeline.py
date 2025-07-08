#!/usr/bin/env python3
"""
Test script to download a TCGA MAF file and process it through the langgraph pipeline.
"""
import sys
import os
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.append(str(PROJECT_ROOT))

from data_processing.tcga_download.src.simple_blood_cancer_client import (
    search_blood_cancer_files, download_file, TCGAFile
)
from geneknow_pipeline.graph import run_pipeline
import json


def main():
    print("🧬 TCGA MAF File Download and Processing Test")
    print("=" * 60)
    
    # Step 1: Search for blood cancer MAF files
    print("\n📍 Step 1: Searching for blood cancer MAF files...")
    
    # Use the simple client to search for MAF files
    maf_files = []
    all_files = search_blood_cancer_files(max_results=20)
    
    # Filter for MAF files only
    for file in all_files:
        if file.file_type == 'MAF' or file.file_name.endswith('.maf') or file.file_name.endswith('.maf.gz'):
            maf_files.append(file)
    
    if not maf_files:
        print("❌ No MAF files found!")
        return
    
    print(f"\n✅ Found {len(maf_files)} MAF files")
    
    # Select the smallest MAF file for quick testing
    maf_files.sort(key=lambda x: x.file_size)
    selected_file = maf_files[0]
    
    print(f"\n📋 Selected file for testing:")
    print(f"   Name: {selected_file.file_name}")
    print(f"   Size: {selected_file.file_size / (1024**2):.2f} MB")
    print(f"   Strategy: {selected_file.experimental_strategy}")
    
    # Step 2: Download the file
    print(f"\n📍 Step 2: Downloading MAF file...")
    download_dir = "test_data/tcga_downloads"
    Path(download_dir).mkdir(parents=True, exist_ok=True)
    
    success = download_file(selected_file, download_dir)
    
    if not success:
        print("❌ Download failed!")
        return
    
    file_path = Path(download_dir) / selected_file.file_name
    print(f"✅ File downloaded to: {file_path}")
    
    # Step 3: Process through langgraph
    print(f"\n📍 Step 3: Processing through langgraph pipeline...")
    
    try:
        # Run the pipeline
        result = run_pipeline(
            file_path=str(file_path),
            user_preferences={
                "language": "en",
                "include_technical": True
            }
        )
        
        # Check results
        if result["pipeline_status"] == "completed":
            print("\n✅ Pipeline completed successfully!")
            
            print(f"\n📊 Results Summary:")
            print(f"   Total variants: {result.get('variant_count', 0)}")
            print(f"   Completed nodes: {result.get('completed_nodes', [])}")
            
            # Check MAF-specific metadata
            maf_info = result.get("file_metadata", {}).get("maf_info", {})
            if maf_info:
                print(f"\n   MAF Analysis:")
                print(f"     - Unique genes: {maf_info.get('unique_genes', 0)}")
                top_genes = list(maf_info.get('top_mutated_genes', {}).keys())[:5]
                if top_genes:
                    print(f"     - Top mutated genes: {', '.join(top_genes)}")
                variant_types = list(maf_info.get('variant_classifications', {}).keys())[:5]
                if variant_types:
                    print(f"     - Variant types: {', '.join(variant_types)}")
            
            # Check risk scores
            risk_scores = result.get("risk_scores", {})
            if risk_scores:
                print(f"\n   Cancer Risk Scores:")
                for cancer, score in sorted(risk_scores.items(), key=lambda x: x[1], reverse=True):
                    print(f"     - {cancer}: {score:.1f}%")
            
            # Check TCGA matches
            tcga_matches = result.get("tcga_matches", {})
            if tcga_matches:
                print(f"\n   TCGA Database Matches:")
                for cancer_type, matches in tcga_matches.items():
                    matched_count = sum(1 for v in matches.values() if v.get('frequency', 0) > 0)
                    if matched_count > 0:
                        print(f"     - {cancer_type}: {matched_count} variants matched")
            
            # Save results
            results_file = Path("test_data") / f"tcga_maf_pipeline_results_{selected_file.file_id[:8]}.json"
            with open(results_file, 'w') as f:
                json.dump(result, f, indent=2, default=str)
            print(f"\n💾 Full results saved to: {results_file}")
            
        else:
            print(f"\n❌ Pipeline failed: {result.get('pipeline_status')}")
            errors = result.get('errors', [])
            for error in errors:
                print(f"   Error in {error.get('node')}: {error.get('error')}")
        
    except Exception as e:
        print(f"\n❌ Pipeline execution error: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n🎉 Test complete!")


if __name__ == "__main__":
    main() 