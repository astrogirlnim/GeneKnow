"""
Test script for the GeneKnow API server.
Shows how to call the API endpoints.
"""
import requests
import json

# API base URL
BASE_URL = "http://localhost:5001"

def test_health():
    """Test the health endpoint."""
    print("Testing health endpoint...")
    response = requests.get(f"{BASE_URL}/health")
    print(f"Status: {response.status_code}")
    print(f"Response: {json.dumps(response.json(), indent=2)}")
    print()

def test_analyze_path():
    """Test analysis with file path."""
    print("Testing analyze_path endpoint...")
    
    # Prepare request data
    data = {
        "file_path": "test-data/sample.fastq.gz",  # ⚠️ MOCK path
        "preferences": {
            "language": "en",
            "include_technical_details": True,
            "risk_threshold_percentage": 30.0
        }
    }
    
    # Make request
    response = requests.post(
        f"{BASE_URL}/analyze_path",
        json=data,
        headers={"Content-Type": "application/json"}
    )
    
    print(f"Status: {response.status_code}")
    
    if response.status_code == 200:
        result = response.json()
        print(f"Pipeline Status: {result['status']}")
        print(f"Processing Time: {result['processing_time']:.2f}s")
        
        # Show risk scores
        if 'data' in result:
            risk_scores = result['data']['risk_assessment']['scores']
            print("\nRisk Scores:")
            for cancer, score in risk_scores.items():
                print(f"  - {cancer}: {score}%")
        
        # Show first part of report
        if 'report' in result:
            print("\nReport Preview:")
            print(result['report'][:300] + "...")
    else:
        print(f"Error: {response.json()}")
    print()

def show_curl_examples():
    """Show example curl commands."""
    print("Example curl commands:")
    print("=" * 50)
    
    print("\n1. Health check:")
    print("curl http://localhost:5001/health")
    
    print("\n2. Analyze by path:")
    print("""curl -X POST http://localhost:5001/analyze_path \\
  -H "Content-Type: application/json" \\
  -d '{
    "file_path": "test-data/sample.fastq.gz",
    "preferences": {
      "language": "en",
      "risk_threshold_percentage": 30.0
    }
  }'""")
    
    print("\n3. Upload and analyze file:")
    print("""curl -X POST http://localhost:5001/analyze \\
  -F "file=@/path/to/your/file.fastq.gz" \\
  -F 'preferences={"language":"en","risk_threshold_percentage":30.0}'""")
    print()

def show_typescript_example():
    """Show TypeScript/Tauri integration example."""
    print("\nTypeScript/Tauri Integration Example:")
    print("=" * 50)
    print("""
// In your Tauri frontend (genomicProcessing.ts)

export async function analyzeGenomicFile(filePath: string) {
  try {
    const response = await fetch('http://localhost:5001/analyze_path', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        file_path: filePath,
        preferences: {
          language: 'en',
          include_technical_details: true,
          risk_threshold_percentage: 30.0
        }
      })
    });
    
    const result = await response.json();
    
    if (result.status === 'success') {
      // Update UI with results
      console.log('Risk scores:', result.data.risk_assessment.scores);
      console.log('Report:', result.report);
      return result.data;
    } else {
      throw new Error(result.error || 'Analysis failed');
    }
  } catch (error) {
    console.error('API call failed:', error);
    throw error;
  }
}
""")

if __name__ == "__main__":
    print("🧬 GeneKnow API Test Script")
    print("=" * 50)
    print("\nMake sure the API server is running:")
    print("  python api_server.py")
    print("\nPress Enter to run tests...")
    input()
    
    try:
        test_health()
        test_analyze_path()
    except requests.exceptions.ConnectionError:
        print("❌ Error: Could not connect to API server")
        print("Make sure the server is running: python api_server.py")
    
    show_curl_examples()
    show_typescript_example() 