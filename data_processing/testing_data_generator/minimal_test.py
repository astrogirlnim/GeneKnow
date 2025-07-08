#!/usr/bin/env python3
"""
Minimal test to download files from GDC API
"""

import json
import requests
from pathlib import Path

def test_download():
    """Test downloading files from GDC API"""
    
    print("🧬 Starting GDC API download test...")
    
    base_url = "https://api.gdc.cancer.gov"
    
    # Simple query for small files
    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "files.access",
                    "value": ["open"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_format",
                    "value": ["TSV", "TXT"]
                }
            }
        ]
    }
    
    params = {
        "filters": json.dumps(filters),
        "format": "json",
        "size": "2"  # Just 2 small files
    }
    
    try:
        print("📡 Making API request...")
        response = requests.get(f"{base_url}/files", params=params)
        response.raise_for_status()
        
        data = response.json()
        files = data.get('data', {}).get('hits', [])
        
        print(f"✅ Found {len(files)} files")
        
        if files:
            # Download first file
            file_data = files[0]
            file_id = file_data['id']
            file_name = file_data['file_name']
            file_size = file_data['file_size']
            
            print(f"📥 Downloading: {file_name} ({file_size} bytes)")
            
            # Download file
            download_url = f"{base_url}/data/{file_id}"
            download_response = requests.get(download_url, stream=True)
            download_response.raise_for_status()
            
            # Save file
            download_dir = Path("minimal_test_downloads")
            download_dir.mkdir(exist_ok=True)
            
            file_path = download_dir / file_name
            with open(file_path, 'wb') as f:
                for chunk in download_response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            
            # Check file size
            downloaded_size = file_path.stat().st_size
            print(f"✅ Downloaded {file_name}: {downloaded_size} bytes")
            
            if downloaded_size == file_size:
                print("🎉 Download successful!")
                
                # Show file info
                print(f"📁 File saved to: {file_path}")
                print(f"📊 File info:")
                print(f"   ID: {file_id}")
                print(f"   Format: {file_data.get('data_format', 'Unknown')}")
                print(f"   Strategy: {file_data.get('experimental_strategy', 'Unknown')}")
                print(f"   Size: {file_size} bytes")
                
                return True
            else:
                print(f"❌ Size mismatch: expected {file_size}, got {downloaded_size}")
                return False
                
        else:
            print("⚠️  No files found")
            return False
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_download()
    if success:
        print("\n🎉 Test completed successfully!")
    else:
        print("\n❌ Test failed!") 