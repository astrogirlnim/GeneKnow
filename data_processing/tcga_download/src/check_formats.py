#!/usr/bin/env python3
"""
Check available data formats in GDC API
"""

import json
import requests
import logging
from collections import Counter

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def check_formats():
    """Check what data formats are available"""
    
    base_url = "https://api.gdc.cancer.gov"
    
    # Query for any files
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
        "format": "json",
        "size": "100"  # Get 100 files to sample formats
    }
    
    try:
        logging.info("Getting sample of files to check formats...")
        response = requests.get(f"{base_url}/files", params=params)
        response.raise_for_status()
        
        data = response.json()
        files = data.get('data', {}).get('hits', [])
        
        # Count data formats
        formats = Counter()
        strategies = Counter()
        
        for file_data in files:
            formats[file_data.get('data_format', 'Unknown')] += 1
            strategies[file_data.get('experimental_strategy', 'Unknown')] += 1
        
        logging.info(f"Found {len(files)} files")
        logging.info(f"Data formats: {dict(formats)}")
        logging.info(f"Experimental strategies: {dict(strategies)}")
        
        # Look for sequencing files specifically
        sequencing_files = [f for f in files if f.get('experimental_strategy') in ['RNA-Seq', 'WXS', 'WGS']]
        logging.info(f"Found {len(sequencing_files)} sequencing files")
        
        if sequencing_files:
            seq_formats = Counter()
            for f in sequencing_files:
                seq_formats[f.get('data_format', 'Unknown')] += 1
            logging.info(f"Sequencing file formats: {dict(seq_formats)}")
        
    except Exception as e:
        logging.error(f"Error: {e}")

if __name__ == "__main__":
    check_formats() 