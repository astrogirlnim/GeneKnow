#!/usr/bin/env python3
"""
Simple debug script to understand GDC API response format
"""

import json
import requests
import logging

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def test_gdc_api():
    """Test the GDC API to understand response format"""
    
    base_url = "https://api.gdc.cancer.gov"
    
    # Simple query for any files
    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "files.access",
                    "value": ["open"]
                }
            }
        ]
    }
    
    params = {
        "filters": json.dumps(filters),
        "expand": "cases,cases.diagnoses,cases.project",
        "format": "json",
        "size": "5"  # Just 5 files for testing
    }
    
    try:
        logging.info("Making API request...")
        response = requests.get(f"{base_url}/files", params=params)
        response.raise_for_status()
        
        data = response.json()
        
        logging.info(f"Response keys: {list(data.keys())}")
        
        # Fix: Get files from data.hits
        files = data.get('data', {}).get('hits', [])
        total_files = data.get('data', {}).get('pagination', {}).get('total', 0)
        
        logging.info(f"Number of files returned: {len(files)}")
        logging.info(f"Total files available: {total_files}")
        
        if files:
            # Handle if files is a dict or list
            if isinstance(files, list) and len(files) > 0:
                logging.info(f"First file keys: {list(files[0].keys())}")
                logging.info(f"First file: {json.dumps(files[0], indent=2)}")
            else:
                logging.info(f"Files structure: {files}")
        else:
            logging.info("No files found with current filters")
        
    except Exception as e:
        logging.error(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_gdc_api() 