#!/usr/bin/env python3
"""
Download TCGA files with well-known cancer mutations using GDC API Client
"""

import logging
import sys
from pathlib import Path
from gdc_api_client import GDCAPIClient

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)

def download_cancer_mutation_files():
    """Download TCGA files with well-known cancer mutations."""
    
    logging.info("🧬 Starting download of cancer files with well-known mutations")
    
    # Initialize GDC client
    client = GDCAPIClient(
        download_dir="cancer_mutation_files",
        max_workers=2
    )
    
    # Modify the search to look for specific cancer types with known mutations
    # We'll override the blood cancer filter to search for all cancer types
    
    # Search for lung cancer files (KRAS, EGFR mutations common)
    logging.info("\n🫁 Searching for lung cancer files (KRAS, EGFR mutations)...")
    lung_files = search_cancer_files(client, 
                                   primary_sites=["Bronchus and lung"],
                                   cancer_name="Lung")
    
    # Search for colorectal cancer files (KRAS, BRAF mutations common)
    logging.info("\n🔴 Searching for colorectal cancer files (KRAS, BRAF mutations)...")
    colon_files = search_cancer_files(client,
                                    primary_sites=["Colon", "Rectum"],
                                    cancer_name="Colorectal")
    
    # Search for melanoma files (BRAF mutations common)
    logging.info("\n🟤 Searching for melanoma files (BRAF V600E mutations)...")
    melanoma_files = search_cancer_files(client,
                                       primary_sites=["Skin"],
                                       cancer_name="Melanoma")
    
    # Search for breast cancer files (PIK3CA, HER2 mutations common)
    logging.info("\n🎀 Searching for breast cancer files (PIK3CA, HER2 mutations)...")
    breast_files = search_cancer_files(client,
                                     primary_sites=["Breast"],
                                     cancer_name="Breast")
    
    # Combine all files
    all_files = lung_files + colon_files + melanoma_files + breast_files
    
    logging.info(f"\n📊 Total files found across all cancer types: {len(all_files)}")
    
    if all_files:
        # Download up to 10 files total
        files_to_download = all_files[:10]
        logging.info(f"📥 Downloading {len(files_to_download)} files...")
        
        progress = client.download_files(files_to_download)
        
        # Print summary
        summary = client.get_download_summary()
        logging.info("\n📊 Download Summary:")
        for key, value in summary.items():
            logging.info(f"   {key}: {value}")
    else:
        logging.warning("⚠️  No files found!")

def search_cancer_files(client, primary_sites, cancer_name, max_results=5):
    """Search for cancer files with specific primary sites."""
    
    import json
    import requests
    
    # Build custom filter for specific cancer types
    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": ["WXS", "RNA-Seq", "WGS"]
                }
            },
            {
                "op": "in", 
                "content": {
                    "field": "files.data_format",
                    "value": ["MAF", "TSV", "TXT"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.access",
                    "value": ["open"]  # Only open access files
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.primary_site",
                    "value": primary_sites
                }
            }
        ]
    }
    
    params = {
        "filters": json.dumps(filters),
        "expand": "cases,cases.diagnoses,cases.project",
        "format": "json",
        "size": str(max_results)
    }
    
    try:
        response = client.session.get(f"{client.base_url}/files", params=params)
        response.raise_for_status()
        
        data = response.json()
        raw_files = data.get('data', {}).get('hits', [])
        
        logging.info(f"   Found {len(raw_files)} {cancer_name} cancer files")
        
        # Convert to TCGAFile objects
        from gdc_api_client import TCGAFile
        tcga_files = []
        
        for file_data in raw_files:
            try:
                # Prioritize MAF files as they contain mutation data
                if file_data['data_format'] == 'MAF' or 'mutation' in file_data['file_name'].lower():
                    tcga_file = TCGAFile(
                        file_id=file_data['id'],
                        file_name=file_data['file_name'],
                        file_size=file_data['file_size'],
                        file_type=file_data['data_format'],
                        submitter_id=file_data['submitter_id'],
                        experimental_strategy=file_data['experimental_strategy'],
                        cases=file_data.get('cases', []),
                        md5sum=file_data.get('md5sum', '')
                    )
                    tcga_files.append(tcga_file)
                    
                    # Log some details about the file
                    case_info = file_data.get('cases', [{}])[0]
                    project_id = case_info.get('project', {}).get('project_id', 'Unknown')
                    logging.info(f"   ✓ {file_data['file_name']} - {project_id}")
                    
            except Exception as e:
                logging.error(f"   Error processing file: {e}")
                continue
        
        return tcga_files
        
    except Exception as e:
        logging.error(f"   Error searching for {cancer_name} files: {e}")
        return []


if __name__ == "__main__":
    download_cancer_mutation_files() 