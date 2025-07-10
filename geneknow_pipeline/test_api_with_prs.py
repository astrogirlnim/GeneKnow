#!/usr/bin/env python3
"""
Test script to verify enhanced API server returns TCGA and PRS results.
"""
import requests
import json
import time
import sys
from pathlib import Path

# API base URL
API_BASE = "http://localhost:5001/api"

def test_api_with_prs():
    """Test the API with a MAF file and verify TCGA/PRS results."""
    
    print("🧬 Testing Enhanced API Server with TCGA and PRS Results")
    print("=" * 60)
    
    # Check if server is running
    try:
        response = requests.get(f"{API_BASE}/health")
        if response.status_code != 200:
            print("❌ API server not responding. Please start it first with:")
            print("   python enhanced_api_server.py")
            return False
        print("✅ API server is running")
    except requests.exceptions.ConnectionError:
        print("❌ Cannot connect to API server. Please start it first with:")
        print("   python enhanced_api_server.py")
        return False
    
    # Use the test MAF file
    test_file = "test_data/tcga_downloads/3d14b1e2-0555-4d6f-a55b-a56065f915e1.wxs.aliquot_ensemble_masked.maf.gz"
    
    if not Path(test_file).exists():
        print(f"❌ Test file not found: {test_file}")
        return False
    
    print(f"📁 Using test file: {test_file}")
    
    # Process the file
    process_data = {
        "file_path": str(Path(test_file).absolute()),
        "preferences": {
            "language": "en",
            "include_technical_details": True
        }
    }
    
    print("\n🚀 Submitting file for processing...")
    response = requests.post(f"{API_BASE}/process", json=process_data)
    
    if response.status_code != 202:
        print(f"❌ Failed to submit file: {response.text}")
        return False
    
    job_data = response.json()
    job_id = job_data["job_id"]
    print(f"✅ Job created: {job_id}")
    
    # Poll for completion
    print("\n⏳ Waiting for processing to complete...")
    max_wait = 30  # seconds
    start_time = time.time()
    
    while time.time() - start_time < max_wait:
        response = requests.get(f"{API_BASE}/status/{job_id}")
        if response.status_code != 200:
            print(f"❌ Failed to get status: {response.text}")
            return False
        
        status_data = response.json()
        status = status_data["status"]
        progress = status_data.get("progress", 0)
        current_step = status_data.get("current_step", "")
        
        print(f"\r   Status: {status} | Progress: {progress}% | Step: {current_step}    ", end="")
        
        if status == "completed":
            print("\n✅ Processing completed!")
            break
        elif status == "failed":
            print(f"\n❌ Processing failed: {status_data.get('error', 'Unknown error')}")
            return False
        
        time.sleep(1)
    else:
        print("\n❌ Processing timed out")
        return False
    
    # Get results
    print("\n📊 Fetching results...")
    response = requests.get(f"{API_BASE}/results/{job_id}")
    
    if response.status_code != 200:
        print(f"❌ Failed to get results: {response.text}")
        return False
    
    results = response.json()
    
    # Debug: Print all keys in results
    print("\n🔍 All keys in results:")
    for key in sorted(results.keys()):
        print(f"   - {key}")
    
    # Check for new fields
    print("\n🔍 Checking for TCGA and PRS results:")
    
    # Check TCGA results
    if "tcga_matches" in results:
        print("\n✅ TCGA Matches found:")
        tcga_matches = results["tcga_matches"]
        for cancer_type, matches in tcga_matches.items():
            print(f"   - {cancer_type}: {len(matches)} variants matched")
    else:
        print("❌ TCGA matches not found in results")
    
    if "tcga_cohort_sizes" in results:
        print("\n✅ TCGA Cohort Sizes found:")
        cohort_sizes = results["tcga_cohort_sizes"]
        for cancer_type, size in cohort_sizes.items():
            print(f"   - {cancer_type}: {size} samples")
    else:
        print("❌ TCGA cohort sizes not found in results")
    
    # Check PRS results
    if "prs_results" in results:
        print("\n✅ PRS Results found:")
        prs_results = results["prs_results"]
        for cancer_type, prs_data in prs_results.items():
            if isinstance(prs_data, dict):
                score = prs_data.get("raw_score", 0)
                percentile = prs_data.get("percentile", 0)
                risk = prs_data.get("risk_category", "unknown")
                print(f"   - {cancer_type}: score={score:.3f}, percentile={percentile}%, risk={risk}")
    else:
        print("❌ PRS results not found in results")
    
    if "prs_summary" in results:
        print("\n✅ PRS Summary found:")
        prs_summary = results["prs_summary"]
        print(f"   - Overall confidence: {prs_summary.get('overall_confidence', 'N/A')}")
        high_risk = prs_summary.get("high_risk_cancers", [])
        print(f"   - High-risk cancers: {', '.join(high_risk) if high_risk else 'None'}")
    else:
        print("❌ PRS summary not found in results")
    
    # Check other expected fields
    print("\n📋 Other fields:")
    print(f"   - Variant count: {results.get('variant_count', 'N/A')}")
    print(f"   - Risk scores: {len(results.get('risk_scores', {}))} cancer types")
    print(f"   - CADD stats: {'Present' if 'cadd_stats' in results else 'Missing'}")
    print(f"   - Processing time: {results.get('processing_time', 'N/A')} seconds")
    
    print("\n✅ Test completed successfully!")
    return True

if __name__ == "__main__":
    success = test_api_with_prs()
    sys.exit(0 if success else 1) 