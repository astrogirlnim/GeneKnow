#!/usr/bin/env python3
"""
Test script to verify API returns survival analysis data.
"""

import requests
import json
import time
import sys
import os

# API configuration
API_URL = "http://localhost:5001"


def test_api_survival_analysis():
    """Test that the API returns survival analysis data."""
    print("🧪 Testing API Survival Analysis Integration")
    print("=" * 60)

    # Check if API is running
    try:
        response = requests.get(f"{API_URL}/api/health", timeout=5)
        if response.status_code != 200:
            print("❌ API server not responding")
            return False
        print("✅ API server is running")
    except requests.exceptions.RequestException as e:
        print(f"❌ Cannot connect to API: {e}")
        return False

    # Use a test VCF file if available
    test_file = None
    possible_files = [
        "test_data/test_variants.vcf",
        "test_data/242b87b1-bf7c-4c1a-bed2-bb077e5ccd00.wxs.aliquot_ensemble_masked.maf.gz",
        "../test_data/test_variants.vcf",
    ]

    for file_path in possible_files:
        if os.path.exists(file_path):
            test_file = os.path.abspath(file_path)
            break

    if not test_file:
        print("❌ No test file found")
        print("Available test files should be in test_data/ directory")
        return False

    print(f"✅ Using test file: {test_file}")

    # Submit job
    print("\n📤 Submitting file for processing...")
    try:
        response = requests.post(
            f"{API_URL}/api/process",
            json={
                "file_path": test_file,
                "preferences": {
                    "patient_data": {"age": 45, "sex": "F"},
                    "language": "en",
                    "include_technical": True,
                },
            },
            timeout=30,
        )

        if response.status_code != 202:
            print(f"❌ Failed to submit job: {response.status_code} - {response.text}")
            return False

        job_data = response.json()
        job_id = job_data["job_id"]
        print(f"✅ Job submitted: {job_id}")

    except requests.exceptions.RequestException as e:
        print(f"❌ Error submitting job: {e}")
        return False

    # Poll for completion
    print("\n⏳ Waiting for processing...")
    max_wait = 120  # 2 minutes
    start_time = time.time()

    while time.time() - start_time < max_wait:
        try:
            response = requests.get(f"{API_URL}/api/status/{job_id}", timeout=10)
            if response.status_code != 200:
                print(f"❌ Error getting status: {response.status_code}")
                return False

            status_data = response.json()
            status = status_data["status"]
            progress = status_data.get("progress", 0)

            print(f"⏳ Progress: {progress}% - Status: {status}")

            if status == "completed":
                print("✅ Processing completed")
                break
            elif status == "failed":
                print(
                    f"❌ Processing failed: {status_data.get('error', 'Unknown error')}"
                )
                return False
            elif status == "cancelled":
                print("❌ Processing was cancelled")
                return False

        except requests.exceptions.RequestException as e:
            print(f"❌ Error polling status: {e}")
            return False

        time.sleep(2)
    else:
        print("❌ Timeout waiting for processing")
        return False

    # Get results
    print("\n📥 Getting results...")
    try:
        response = requests.get(f"{API_URL}/api/results/{job_id}", timeout=30)
        if response.status_code != 200:
            print(f"❌ Error getting results: {response.status_code}")
            return False

        results = response.json()
        print("✅ Results retrieved")

        # Check for survival analysis
        structured_json = results.get("structured_json", {})
        survival_analysis = structured_json.get("survival_analysis")

        if not survival_analysis:
            print("❌ No survival analysis in results")
            print("Available keys in structured_json:")
            for key in structured_json.keys():
                print(f"  - {key}")
            return False

        print("✅ Survival analysis found in results")

        # Check survival curves
        survival_curves = survival_analysis.get("survival_curves", {})
        if not survival_curves:
            print("❌ No survival curves in survival analysis")
            return False

        print(f"✅ Found {len(survival_curves)} survival curves")

        # Validate structure
        for cancer_type, curve_data in survival_curves.items():
            required_fields = [
                "time_points",
                "population_survival",
                "patient_survival",
                "median_survival",
                "hazard_ratio",
                "five_year_survival",
            ]

            for field in required_fields:
                if field not in curve_data:
                    print(f"❌ Missing field {field} in {cancer_type} survival curve")
                    return False

            # Check data types
            if not isinstance(curve_data["time_points"], list):
                print(f"❌ time_points is not a list for {cancer_type}")
                return False

            if not isinstance(curve_data["population_survival"], list):
                print(f"❌ population_survival is not a list for {cancer_type}")
                return False

            if not isinstance(curve_data["patient_survival"], list):
                print(f"❌ patient_survival is not a list for {cancer_type}")
                return False

            print(
                f"✅ {cancer_type}: HR={curve_data['hazard_ratio']}, {len(curve_data['time_points'])} points"
            )

        # Check prognostic factors
        prognostic_factors = survival_analysis.get("prognostic_factors")
        if not prognostic_factors:
            print("❌ No prognostic factors found")
            return False

        print("✅ Prognostic factors found")

        # Check clinical interpretation
        clinical_interpretation = survival_analysis.get("clinical_interpretation")
        if not clinical_interpretation:
            print("❌ No clinical interpretation found")
            return False

        print("✅ Clinical interpretation found")

        # Save results for inspection
        with open("api_survival_results.json", "w") as f:
            json.dump(survival_analysis, f, indent=2)

        print("✅ Results saved to api_survival_results.json")

        return True

    except requests.exceptions.RequestException as e:
        print(f"❌ Error getting results: {e}")
        return False
    except json.JSONDecodeError as e:
        print(f"❌ Error parsing JSON: {e}")
        return False


if __name__ == "__main__":
    success = test_api_survival_analysis()
    print("\n" + "=" * 60)
    if success:
        print("🎉 API Survival Analysis Test: PASSED")
        print("Survival analysis is working correctly through the API!")
    else:
        print("❌ API Survival Analysis Test: FAILED")
        print("Check the API server logs for more details.")

    sys.exit(0 if success else 1)
