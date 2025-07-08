"""
Test script for the enhanced GeneKnow Pipeline API.
Demonstrates various API endpoints and WebSocket functionality.
"""
import requests
import json
import time
import os
import sys
from datetime import datetime
import socketio

# Configuration
API_BASE_URL = "http://localhost:5001"

def test_health_check():
    """Test the health check endpoint."""
    print("\n=== Testing Health Check ===")
    response = requests.get(f"{API_BASE_URL}/api/health")
    if response.status_code == 200:
        data = response.json()
        print(f"✅ API is healthy")
        print(f"   Service: {data['service']}")
        print(f"   Version: {data['version']}")
        print(f"   Active jobs: {data['jobs_active']}")
        return True
    else:
        print(f"❌ Health check failed: {response.status_code}")
        return False

def test_pipeline_info():
    """Test the pipeline info endpoint."""
    print("\n=== Testing Pipeline Info ===")
    response = requests.get(f"{API_BASE_URL}/api/pipeline-info")
    if response.status_code == 200:
        data = response.json()
        print(f"✅ Pipeline: {data['name']} v{data['version']}")
        print(f"   Supported formats: {', '.join(data['capabilities']['supported_formats'])}")
        print(f"   Cancer types: {', '.join(data['capabilities']['cancer_types'])}")
        print(f"   Max file size: {data['capabilities']['max_file_size_gb']}GB")
        return True
    else:
        print(f"❌ Failed to get pipeline info: {response.status_code}")
        return False

def test_supported_formats():
    """Test the supported formats endpoint."""
    print("\n=== Testing Supported Formats ===")
    response = requests.get(f"{API_BASE_URL}/api/supported-formats")
    if response.status_code == 200:
        data = response.json()
        print(f"✅ Found {len(data['formats'])} supported formats:")
        for fmt in data['formats']:
            print(f"   - {fmt['extension']}: {fmt['description']}")
            if fmt['compressed']:
                print(f"     Compressed: {', '.join(fmt['compressed'])}")
        return True
    else:
        print(f"❌ Failed to get supported formats: {response.status_code}")
        return False

def test_process_local_file(file_path):
    """Test processing a local file."""
    print(f"\n=== Testing Local File Processing ===")
    print(f"File: {file_path}")
    
    # Start processing
    response = requests.post(
        f"{API_BASE_URL}/api/process",
        json={
            "file_path": file_path,
            "preferences": {
                "language": "en",
                "include_technical": True,
                "patient_data": {
                    "age": 45,
                    "sex": "F"
                }
            }
        }
    )
    
    if response.status_code != 202:
        print(f"❌ Failed to start processing: {response.status_code}")
        print(f"   Error: {response.json().get('error', 'Unknown error')}")
        return None
    
    job_data = response.json()
    job_id = job_data['job_id']
    print(f"✅ Job created: {job_id}")
    
    # Poll for status
    print("\nTracking progress:")
    last_progress = -1
    while True:
        status_response = requests.get(f"{API_BASE_URL}/api/status/{job_id}")
        if status_response.status_code != 200:
            print(f"❌ Failed to get status: {status_response.status_code}")
            return None
        
        status = status_response.json()
        
        # Print progress if changed
        if status['progress'] != last_progress:
            print(f"   [{status['progress']:3d}%] {status['status']} - {status.get('current_step', 'Initializing')}")
            last_progress = status['progress']
        
        if status['status'] in ['completed', 'failed', 'cancelled']:
            break
        
        time.sleep(0.5)
    
    # Get results if completed
    if status['status'] == 'completed':
        print("\n✅ Processing completed!")
        
        results_response = requests.get(f"{API_BASE_URL}/api/results/{job_id}")
        if results_response.status_code == 200:
            results = results_response.json()
            print(f"\nResults:")
            print(f"   Variant count: {results.get('variant_count', 0)}")
            print(f"   Processing time: {results.get('processing_time_seconds', 0):.2f}s")
            
            risk_scores = results.get('risk_scores', {})
            if risk_scores:
                print(f"\n   Risk Scores:")
                for cancer_type, score in risk_scores.items():
                    print(f"     - {cancer_type}: {score:.2%}")
            
            return job_id
        else:
            print(f"❌ Failed to get results: {results_response.status_code}")
            return None
    else:
        print(f"\n❌ Job {status['status']}: {status.get('error', 'Unknown error')}")
        return None

def test_list_jobs():
    """Test listing jobs."""
    print("\n=== Testing Job List ===")
    response = requests.get(f"{API_BASE_URL}/api/jobs?limit=5")
    if response.status_code == 200:
        data = response.json()
        print(f"✅ Found {data['total']} total jobs")
        print(f"   Showing {len(data['jobs'])} most recent:")
        for job in data['jobs']:
            print(f"   - {job['id']}: {job['status']} ({job['created_at']})")
        return True
    else:
        print(f"❌ Failed to list jobs: {response.status_code}")
        return False

def test_websocket_progress(file_path):
    """Test WebSocket real-time progress updates."""
    print(f"\n=== Testing WebSocket Progress Updates ===")
    print(f"File: {file_path}")
    
    # Initialize SocketIO client
    sio = socketio.Client()
    progress_updates = []
    
    @sio.event
    def connect():
        print("✅ Connected to WebSocket server")
    
    @sio.event
    def disconnect():
        print("🔌 Disconnected from WebSocket server")
    
    @sio.event
    def job_progress(data):
        progress_updates.append(data)
        print(f"   📡 [{data['progress']:3d}%] {data['status']} - {data.get('current_step', 'Processing')}")
    
    # Connect to WebSocket server
    try:
        sio.connect(API_BASE_URL)
    except Exception as e:
        print(f"❌ Failed to connect to WebSocket: {e}")
        return False
    
    # Start processing
    response = requests.post(
        f"{API_BASE_URL}/api/process",
        json={
            "file_path": file_path,
            "preferences": {
                "language": "en",
                "include_technical": False
            }
        }
    )
    
    if response.status_code != 202:
        print(f"❌ Failed to start processing: {response.status_code}")
        sio.disconnect()
        return False
    
    job_data = response.json()
    job_id = job_data['job_id']
    print(f"✅ Job created: {job_id}")
    
    # Subscribe to job progress
    sio.emit('subscribe_job', {'job_id': job_id})
    print("📡 Subscribed to job progress updates")
    
    # Wait for completion
    start_time = time.time()
    timeout = 60  # 60 seconds timeout
    
    while True:
        status_response = requests.get(f"{API_BASE_URL}/api/status/{job_id}")
        if status_response.status_code == 200:
            status = status_response.json()
            if status['status'] in ['completed', 'failed', 'cancelled']:
                break
        
        if time.time() - start_time > timeout:
            print("⏱️  Timeout waiting for job completion")
            break
        
        time.sleep(0.5)
    
    # Unsubscribe and disconnect
    sio.emit('unsubscribe_job', {'job_id': job_id})
    sio.disconnect()
    
    print(f"\n📊 Received {len(progress_updates)} WebSocket updates")
    return True

def main():
    """Run all API tests."""
    print("🧬 GeneKnow Pipeline API Test Suite")
    print("=" * 50)
    
    # Check if API server is running
    if not test_health_check():
        print("\n❌ API server is not running. Please start it with:")
        print("   python enhanced_api_server.py")
        sys.exit(1)
    
    # Run basic endpoint tests
    test_pipeline_info()
    test_supported_formats()
    test_list_jobs()
    
    # Test file processing if test file exists
    test_files = [
        "../test_data/test_annotations.tsv",
        "../test-data/sample.vcf",
        "../test_R1.fastq.gz"
    ]
    
    test_file = None
    for file in test_files:
        if os.path.exists(file):
            test_file = os.path.abspath(file)
            break
    
    if test_file:
        # Test regular processing
        job_id = test_process_local_file(test_file)
        
        if job_id:
            # Test downloading results
            print("\n=== Testing Results Download ===")
            response = requests.get(f"{API_BASE_URL}/api/results/{job_id}/download")
            if response.status_code == 200:
                print(f"✅ Downloaded results ({len(response.content)} bytes)")
            else:
                print(f"❌ Failed to download results: {response.status_code}")
        
        # Test WebSocket functionality
        test_websocket_progress(test_file)
    else:
        print("\n⚠️  No test files found. Skipping file processing tests.")
        print("   Create one of these files to test processing:")
        for file in test_files:
            print(f"   - {file}")
    
    print("\n✅ API test suite completed!")

if __name__ == "__main__":
    main() 