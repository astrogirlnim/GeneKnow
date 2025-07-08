#!/usr/bin/env python3
"""
Simple test to download any files from GDC API
"""

import json
import requests
import logging
from pathlib import Path
from gdc_api_client import GDCAPIClient, TCGAFile

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def test_download():
    """Test downloading any files from GDC API"""
    
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
        "size": "3"  # Just 3 files for testing
    }
    
    try:
        logging.info("Getting files for download test...")
        response = requests.get(f"{base_url}/files", params=params)
        response.raise_for_status()
        
        data = response.json()
        files = data.get('data', {}).get('hits', [])
        
        logging.info(f"Found {len(files)} files")
        
        if files:
            # Convert to TCGAFile objects
            tcga_files = []
            for file_data in files:
                tcga_file = TCGAFile(
                    file_id=file_data['id'],
                    file_name=file_data['file_name'],
                    file_size=file_data['file_size'],
                    file_type=file_data['data_format'],
                    submitter_id=file_data.get('submitter_id', ''),
                    experimental_strategy=file_data.get('experimental_strategy', ''),
                    cases=file_data.get('cases', []),
                    md5sum=file_data.get('md5sum', '')
                )
                tcga_files.append(tcga_file)
                logging.info(f"Added file: {tcga_file.file_name} ({tcga_file.file_size} bytes)")
            
            # Download files
            client = GDCAPIClient(download_dir="simple_test_downloads", max_workers=1)
            progress = client.download_files(tcga_files)
            
            # Print summary
            summary = client.get_download_summary()
            logging.info(f"Download Summary:")
            for key, value in summary.items():
                logging.info(f"   {key}: {value}")
                
        else:
            logging.warning("No files found")
        
    except Exception as e:
        logging.error(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_download() 