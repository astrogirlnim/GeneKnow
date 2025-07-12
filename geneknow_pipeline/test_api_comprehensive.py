#!/usr/bin/env python3
"""
Comprehensive API test suite for GeneKnow Pipeline.
Combines all API testing functionality including basic endpoints, enhanced features,
WebSocket support, and PRS/TCGA results verification.
"""
import requests
import time
import sys
import os
from pathlib import Path
from datetime import datetime

# API configuration
BASE_URL = "http://localhost:5001"
API_BASE = f"{BASE_URL}/api"


class APITestSuite:
    """Comprehensive API testing class."""

    def __init__(self):
        self.results = {"total": 0, "passed": 0, "failed": 0, "tests": []}

    def test(self, name, func):
        """Run a test and record results."""
        self.results["total"] += 1
        print(f"\n🧪 {name}")
        print("-" * 60)

        try:
            result = func()
            if result:
                self.results["passed"] += 1
                print("✅ PASSED")
                self.results["tests"].append({"name": name, "status": "passed"})
            else:
                self.results["failed"] += 1
                print("❌ FAILED")
                self.results["tests"].append({"name": name, "status": "failed"})
        except Exception as e:
            self.results["failed"] += 1
            print(f"❌ ERROR: {str(e)}")
            self.results["tests"].append({"name": name, "status": "error", "error": str(e)})

    def print_summary(self):
        """Print test summary."""
        print("\n" + "=" * 60)
        print("📊 TEST SUMMARY")
        print("=" * 60)
        print(f"Total tests: {self.results['total']}")
        print(f"Passed: {self.results['passed']} ({self.results['passed']/self.results['total']*100:.1f}%)")
        print(f"Failed: {self.results['failed']}")

        if self.results["failed"] > 0:
            print("\nFailed tests:")
            for test in self.results["tests"]:
                if test["status"] != "passed":
                    print(f"  - {test['name']}: {test.get('error', 'Failed')}")


# Basic API Tests
def test_health_check():
    """Test the health endpoint."""
    response = requests.get(f"{API_BASE}/health")
    return response.status_code == 200 and response.json()["status"] == "healthy"


def test_pipeline_info():
    """Test pipeline information endpoint."""
    response = requests.get(f"{API_BASE}/pipeline-info")
    if response.status_code != 200:
        return False

    info = response.json()
    # Check top-level fields
    if "version" not in info or "capabilities" not in info:
        return False

    # Check nested fields
    caps = info["capabilities"]
    required_caps = ["pipeline_nodes", "supported_formats"]
    return all(field in caps for field in required_caps)


def test_supported_formats():
    """Test supported formats endpoint."""
    response = requests.get(f"{API_BASE}/supported-formats")
    if response.status_code != 200:
        return False

    formats = response.json()["formats"]
    # Check that we have the expected format types
    extensions = [fmt["extension"].upper().replace(".", "") for fmt in formats]
    expected = ["FASTQ", "BAM", "VCF", "MAF"]
    return all(fmt in extensions for fmt in expected)


# Enhanced API Tests
def test_job_management():
    """Test job creation and management."""
    # Create a job
    test_file = find_test_file()
    if not test_file:
        assert False, "Test failed"

    process_data = {"file_path": str(test_file), "preferences": {"language": "en"}}

    response = requests.post(f"{API_BASE}/process", json=process_data)
    if response.status_code != 202:
        return False

    job_id = response.json()["job_id"]

    # Check job status
    response = requests.get(f"{API_BASE}/status/{job_id}")
    if response.status_code != 200:
        return False

    status = response.json()["status"]
    return status in ["running", "processing", "completed", "failed"]


def test_job_listing():
    """Test listing all jobs."""
    response = requests.get(f"{API_BASE}/jobs")
    if response.status_code != 200:
        return False

    jobs = response.json()["jobs"]
    return isinstance(jobs, list)


# WebSocket Tests
def test_websocket_wrapper():
    """Test Socket.IO connection (simplified check)."""
    # The API uses Socket.IO which is different from plain WebSockets
    # For now, we'll just check if the Socket.IO endpoint exists
    try:
        # Socket.IO handshake is at /socket.io/
        response = requests.get(f"{BASE_URL}/socket.io/", params={"EIO": "4", "transport": "polling"})
        # Socket.IO returns specific status codes/responses for handshake
        return response.status_code in [200, 400]  # 400 is expected without proper handshake
    except Exception as e:
        print(f"Socket.IO check error: {e}")
        return False


# PRS and TCGA Tests
def test_prs_results():
    """Test that PRS results are included in API response."""
    test_file = find_test_file(prefer_maf=True)
    if not test_file:
        assert False, "Test failed"

    # Process file
    job_id = submit_file_for_processing(test_file)
    if not job_id:
        return False

    # Wait for completion
    if not wait_for_job_completion(job_id):
        return False

    # Get results
    response = requests.get(f"{API_BASE}/results/{job_id}")
    if response.status_code != 200:
        return False

    results = response.json()

    # Check for PRS fields
    has_prs = "prs_results" in results
    has_prs_summary = "prs_summary" in results

    if has_prs:
        print(f"  Found PRS results for {len(results['prs_results'])} cancer types")

    return has_prs and has_prs_summary


def test_tcga_results():
    """Test that TCGA matches are included in API response."""
    test_file = find_test_file(prefer_maf=True)
    if not test_file:
        assert False, "Test failed"

    # Process file
    job_id = submit_file_for_processing(test_file)
    if not job_id:
        return False

    # Wait for completion
    if not wait_for_job_completion(job_id):
        return False

    # Get results
    response = requests.get(f"{API_BASE}/results/{job_id}")
    if response.status_code != 200:
        return False

    results = response.json()

    # Check for TCGA fields
    has_tcga = "tcga_matches" in results
    has_cohort = "tcga_cohort_sizes" in results

    if has_tcga:
        total_matches = sum(len(matches) for matches in results["tcga_matches"].values())
        print(f"  Found {total_matches} total TCGA matches")

    return has_tcga and has_cohort


def test_results_download():
    """Test downloading results as a file."""
    test_file = find_test_file()
    if not test_file:
        assert False, "Test failed"

    # Process file
    job_id = submit_file_for_processing(test_file)
    if not job_id:
        return False

    # Wait for completion
    if not wait_for_job_completion(job_id):
        return False

    # Download results
    response = requests.get(f"{API_BASE}/results/{job_id}/download")
    if response.status_code != 200:
        return False

    # Check that we got JSON data
    try:
        data = response.json()
        return "pipeline_status" in data
    except:
        return False


# Utility Functions
def find_test_file(prefer_maf=False):
    """Find a test file to use."""
    test_files = [
        "test_data/tcga_downloads/3d14b1e2-0555-4d6f-a55b-a56065f915e1.wxs.aliquot_ensemble_masked.maf.gz",
        "test_data/test_sample.ma",
        "test_data/test_variants.vc",
        "../test_R1.fastq.gz",
        "../test-data/sample.vc",
    ]

    if prefer_maf:
        # Move MAF files to front
        test_files = [f for f in test_files if ".ma" in f] + [f for f in test_files if ".ma" not in f]

    for file in test_files:
        if os.path.exists(file):
            return Path(file).absolute()

    return None


def submit_file_for_processing(file_path):
    """Submit a file for processing and return job ID."""
    process_data = {"file_path": str(file_path), "preferences": {"language": "en", "include_technical_details": True}}

    response = requests.post(f"{API_BASE}/process", json=process_data)
    if response.status_code == 202:
        return response.json()["job_id"]

    return None


def wait_for_job_completion(job_id, timeout=30):
    """Wait for a job to complete."""
    start_time = time.time()

    while time.time() - start_time < timeout:
        response = requests.get(f"{API_BASE}/status/{job_id}")
        if response.status_code != 200:
            return False

        status = response.json()["status"]

        if status == "completed":
            return True
        elif status == "failed":
            return False

        time.sleep(1)

    return False


def show_curl_examples():
    """Show example curl commands."""
    print("\n📋 Example curl commands:")
    print("=" * 50)

    print("\n1. Health check:")
    print("curl http://localhost:5001/api/health")

    print("\n2. Get pipeline info:")
    print("curl http://localhost:5001/api/pipeline-info")

    print("\n3. Process a file:")
    print(
        """curl -X POST http://localhost:5001/api/process \\
  -H "Content-Type: application/json" \\
  -d '{
    "file_path": "/path/to/file.vc",
    "preferences": {"language": "en"}
  }'"""
    )

    print("\n4. Check job status:")
    print("curl http://localhost:5001/api/status/JOB_ID")

    print("\n5. Get results:")
    print("curl http://localhost:5001/api/results/JOB_ID")

    print("\n6. Download results:")
    print("curl http://localhost:5001/api/results/JOB_ID/download -o results.json")


def main():
    """Run all API tests."""
    print("🧬 GeneKnow Comprehensive API Test Suite")
    print("=" * 60)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Check if API server is running
    try:
        response = requests.get(f"{API_BASE}/health", timeout=2)
        if response.status_code != 200:
            print("\n❌ API server not responding properly")
            print("Please start it with: python enhanced_api_server.py")
            return 1
    except requests.exceptions.ConnectionError:
        print("\n❌ Cannot connect to API server")
        print("Please start it with: python enhanced_api_server.py")
        return 1

    print("✅ API server is running")

    # Run tests
    suite = APITestSuite()

    # Basic tests
    suite.test("Health Check", test_health_check)
    suite.test("Pipeline Info", test_pipeline_info)
    suite.test("Supported Formats", test_supported_formats)

    # Enhanced API tests
    suite.test("Job Management", test_job_management)
    suite.test("Job Listing", test_job_listing)
    suite.test("Results Download", test_results_download)

    # Socket.IO tests
    suite.test("Socket.IO Connection", test_websocket_wrapper)

    # PRS and TCGA tests
    suite.test("PRS Results", test_prs_results)
    suite.test("TCGA Results", test_tcga_results)

    # Print summary
    suite.print_summary()

    # Show examples
    if suite.results["passed"] > 0:
        show_curl_examples()

    return 0 if suite.results["failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
