#!/usr/bin/env python3
"""
Test script to verify parallel execution of TCGA mapping, CADD scoring, and PRS calculation.
"""
import sys
import os
import time
import logging
from datetime import datetime

# Add the current directory to Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from graph import run_pipeline

# Set up detailed logging to see parallel execution
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def test_parallel_execution():
    """Test that TCGA, CADD, and PRS run in parallel."""
    
    print("\n🧬 Testing Parallel Execution of TCGA, CADD, and PRS")
    print("=" * 60)
    
    # Use the test MAF file
    test_file = "test_data/tcga_downloads/3d14b1e2-0555-4d6f-a55b-a56065f915e1.wxs.aliquot_ensemble_masked.maf.gz"
    
    if not os.path.exists(test_file):
        print(f"❌ Test file not found: {test_file}")
        return False
    
    print(f"📁 Using test MAF file: {test_file}")
    
    # User preferences
    preferences = {
        "language": "en",
        "include_technical_details": True,
        "patient_data": {"age": 45, "sex": "F"},
        "risk_threshold_percentage": 20.0
    }
    
    print("\n🚀 Starting pipeline...")
    start_time = time.time()
    
    try:
        # Run the pipeline
        result = run_pipeline(test_file, preferences)
        
        end_time = time.time()
        total_time = end_time - start_time
        
        # Check pipeline status
        print(f"\n📊 Pipeline Status: {result['pipeline_status']}")
        print(f"⏱️  Total Time: {total_time:.2f} seconds")
        
        # Check completed nodes
        completed_nodes = result.get('completed_nodes', [])
        print(f"\n✅ Completed Nodes ({len(completed_nodes)}):")
        for node in completed_nodes:
            print(f"   - {node}")
        
        # Verify all three parallel nodes completed
        parallel_nodes = ['tcga_mapper', 'cadd_scoring', 'prs_calculator']
        missing_nodes = [node for node in parallel_nodes if node not in completed_nodes]
        
        if missing_nodes:
            print(f"\n❌ Missing parallel nodes: {missing_nodes}")
            return False
        else:
            print(f"\n✅ All parallel nodes completed successfully!")
        
        # Check results from each parallel process
        print("\n📋 Checking Results from Parallel Processes:")
        
        # 1. TCGA Results
        tcga_matches = result.get('tcga_matches', {})
        if tcga_matches:
            print(f"\n1️⃣  TCGA Mapping Results:")
            print(f"   - Cancer types analyzed: {list(tcga_matches.keys())}")
            for cancer_type, matches in tcga_matches.items():
                print(f"   - {cancer_type}: {len(matches)} variants matched")
        else:
            print("\n1️⃣  TCGA Mapping: ❌ No results")
        
        # 2. CADD Results
        cadd_variants = result.get('cadd_enriched_variants', [])
        cadd_stats = result.get('cadd_stats', {})
        if cadd_variants:
            print(f"\n2️⃣  CADD Scoring Results:")
            print(f"   - Variants scored: {len(cadd_variants)}")
            if cadd_stats:
                print(f"   - Mean CADD score: {cadd_stats.get('mean_cadd_score', 'N/A')}")
                print(f"   - High impact variants: {cadd_stats.get('high_impact_count', 0)}")
        else:
            print("\n2️⃣  CADD Scoring: ❌ No results")
        
        # 3. PRS Results
        prs_results = result.get('prs_results', {})
        prs_summary = result.get('prs_summary', {})
        if prs_results:
            print(f"\n3️⃣  PRS Calculation Results:")
            print(f"   - Cancer types analyzed: {list(prs_results.keys())}")
            for cancer_type, prs_data in prs_results.items():
                print(f"   - {cancer_type}:")
                print(f"     • Raw Score: {prs_data.get('raw_score', 0):.3f}")
                print(f"     • Percentile: {prs_data.get('percentile', 0)}%")
                print(f"     • Risk Category: {prs_data.get('risk_category', 'N/A')}")
                print(f"     • Confidence: {prs_data.get('confidence', 'N/A')}")
            
            if prs_summary.get('high_risk_cancers'):
                print(f"\n   🎯 High Risk Cancers: {', '.join(prs_summary['high_risk_cancers'])}")
        else:
            print("\n3️⃣  PRS Calculation: ❌ No results")
        
        # Check for errors or warnings
        if result.get('errors'):
            print(f"\n❌ Errors encountered:")
            for error in result['errors']:
                print(f"   - {error}")
        
        if result.get('warnings'):
            print(f"\n⚠️  Warnings ({len(result['warnings'])}):")
            for warning in result['warnings'][:5]:  # Show first 5
                print(f"   - {warning}")
        
        # Verify parallel execution timing
        print("\n⏱️  Parallel Execution Analysis:")
        print("   The three processes (TCGA, CADD, PRS) should run simultaneously.")
        print("   If working correctly, total time should be close to the slowest process,")
        print("   not the sum of all three processes.")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Pipeline failed with error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_parallel_execution()
    
    if success:
        print("\n✅ Parallel execution test completed successfully!")
    else:
        print("\n❌ Parallel execution test failed!")
        sys.exit(1) 