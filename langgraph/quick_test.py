#!/usr/bin/env python3
"""
Quick test script to verify the LangGraph pipeline is working.
"""
import subprocess
import time
import requests
import json
import os
import sys

def test_pipeline():
    """Run a quick test of the pipeline."""
    print("🧬 Quick Pipeline Test")
    print("=" * 50)
    
    # 1. Start the API server
    print("\n1️⃣ Starting API server...")
    server_process = subprocess.Popen(
        [sys.executable, "api_server.py"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    # Wait for server to start
    time.sleep(3)
    
    try:
        # 2. Test if server is running
        print("2️⃣ Checking server status...")
        try:
            response = requests.get("http://localhost:5001/health")
            print(f"   ✅ Server is running: {response.json()}")
        except Exception as e:
            print(f"   ❌ Server not responding: {e}")
            return
        
        # 3. Test with FASTQ file
        print("\n3️⃣ Testing FASTQ processing...")
        test_file = "../test_R1.fastq.gz"
        if not os.path.exists(test_file):
            print(f"   ❌ Test file not found: {test_file}")
            return
            
        response = requests.post(
            "http://localhost:5001/analyze_path",
            json={
                "file_path": test_file,
                "preferences": {
                    "patient_data": {"age": 45, "sex": "F"}
                }
            }
        )
        
        if response.status_code == 200:
            result = response.json()
            print("   ✅ FASTQ processing successful!")
            print(f"   📊 Found {len(result.get('variants', []))} variants")
            print(f"   🎯 Risk scores: {json.dumps(result.get('risk_scores', {}), indent=2)}")
        else:
            print(f"   ❌ Processing failed: {response.status_code}")
            print(f"   Error: {response.text}")
        
        # 4. Check key components
        print("\n4️⃣ Checking key components:")
        components = {
            "TCGA Database": "tcga_data/tcga_variants.db",
            "ML Models": "models/breast_model.pkl",
            "Test Reference": "test_reference/test_genome.fa",
            "Test VCF": "test_data/test_variants.vcf"
        }
        
        for name, path in components.items():
            if os.path.exists(path):
                print(f"   ✅ {name}: {path}")
            else:
                print(f"   ❌ {name}: Missing {path}")
                
    finally:
        # 5. Stop server
        print("\n5️⃣ Stopping server...")
        server_process.terminate()
        server_process.wait()
        print("   ✅ Server stopped")
    
    print("\n" + "=" * 50)
    print("✅ Quick test complete!")
    print("\n📝 Next steps:")
    print("  - Run `python test_pipeline.py` for full tests")
    print("  - Run `python test_parallelization.py` for parallel tests")
    print("  - Start server with `python api_server.py` for manual testing")

if __name__ == "__main__":
    test_pipeline() 