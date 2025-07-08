#!/usr/bin/env python3
"""
Simplified GDC API Client for TCGA Blood Cancer Data Download

This version uses metadata directly for cancer classification, no LLM needed.
"""

import json
import requests
from pathlib import Path
from typing import List, Dict
from dataclasses import dataclass

@dataclass
class TCGAFile:
    """Represents a TCGA file with its metadata."""
    file_id: str
    file_name: str
    file_size: int
    file_type: str
    experimental_strategy: str
    cases: List[Dict]
    md5sum: str = ""

def is_blood_cancer(case: Dict) -> bool:
    """
    Check if a case represents blood or bone marrow cancer using metadata.
    
    Args:
        case: Case data from GDC API
        
    Returns:
        True if blood/bone marrow cancer, False otherwise
    """
    # Check primary site
    primary_site = case.get('primary_site', '').lower()
    blood_sites = {
        'blood', 'bone marrow', 'lymph nodes', 'spleen', 
        'hematopoietic and reticuloendothelial systems'
    }
    
    if any(site in primary_site for site in blood_sites):
        return True
    
    # Check disease type
    disease_type = case.get('disease_type', '').lower()
    blood_diseases = {
        'acute myeloid leukemia', 'acute lymphoblastic leukemia',
        'chronic lymphocytic leukemia', 'chronic myelogenous leukemia',
        'hodgkin lymphoma', 'non-hodgkin lymphoma', 'multiple myeloma',
        'myelodysplastic syndromes', 'myelofibrosis', 'leukemia',
        'lymphoma', 'myeloma'
    }
    
    if any(disease in disease_type for disease in blood_diseases):
        return True
    
    # Check diagnoses
    for diagnosis in case.get('diagnoses', []):
        primary_diagnosis = diagnosis.get('primary_diagnosis', '').lower()
        blood_diagnoses = {
            'leukemia', 'lymphoma', 'myeloma', 'myelodysplastic',
            'myelofibrosis', 'hodgkin', 'burkitt'
        }
        
        if any(blood_diag in primary_diagnosis for blood_diag in blood_diagnoses):
            return True
    
    return False

def search_blood_cancer_files(max_results: int = 50) -> List[TCGAFile]:
    """
    Search for blood cancer files using simple metadata filtering.
    
    Args:
        max_results: Maximum number of results to return
        
    Returns:
        List of TCGAFile objects for blood cancer samples
    """
    print(f"🔍 Searching for blood cancer files (max {max_results})...")
    
    base_url = "https://api.gdc.cancer.gov"
    
    # Search for files with blood/bone marrow primary sites
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
                    "value": ["TSV", "TXT", "MAF"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": ["RNA-Seq", "WXS", "WGS"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.primary_site",
                    "value": [
                        "Blood", "Bone Marrow", "Lymph Nodes", "Spleen",
                        "Hematopoietic and reticuloendothelial systems"
                    ]
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
        print("📡 Making API request...")
        response = requests.get(f"{base_url}/files", params=params)
        response.raise_for_status()
        
        data = response.json()
        files = data.get('data', {}).get('hits', [])
        
        print(f"✅ Found {len(files)} files from API")
        
        # Filter for blood cancer files
        blood_cancer_files = []
        for file_data in files:
            cases = file_data.get('cases', [])
            
            # Check if any case is blood cancer
            has_blood_cancer = any(is_blood_cancer(case) for case in cases)
            
            if has_blood_cancer:
                tcga_file = TCGAFile(
                    file_id=file_data['id'],
                    file_name=file_data['file_name'],
                    file_size=file_data['file_size'],
                    file_type=file_data['data_format'],
                    experimental_strategy=file_data.get('experimental_strategy', ''),
                    cases=cases,
                    md5sum=file_data.get('md5sum', '')
                )
                blood_cancer_files.append(tcga_file)
                
                # Show what we found
                case_info = []
                for case in cases:
                    if is_blood_cancer(case):
                        site = case.get('primary_site', 'Unknown')
                        disease = case.get('disease_type', 'Unknown')
                        case_info.append(f"{site} - {disease}")
                
                print(f"🩸 Blood cancer file: {tcga_file.file_name}")
                print(f"   Cases: {'; '.join(case_info)}")
        
        print(f"🎯 Found {len(blood_cancer_files)} blood cancer files")
        return blood_cancer_files
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return []

def download_file(tcga_file: TCGAFile, download_dir: str = "downloads") -> bool:
    """
    Download a single TCGA file.
    
    Args:
        tcga_file: TCGAFile object to download
        download_dir: Directory to download to
        
    Returns:
        True if successful, False otherwise
    """
    base_url = "https://api.gdc.cancer.gov"
    download_url = f"{base_url}/data/{tcga_file.file_id}"
    
    # Create download directory
    Path(download_dir).mkdir(exist_ok=True)
    file_path = Path(download_dir) / tcga_file.file_name
    
    try:
        print(f"📥 Downloading: {tcga_file.file_name} ({tcga_file.file_size} bytes)")
        
        response = requests.get(download_url, stream=True)
        response.raise_for_status()
        
        with open(file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        
        # Verify file size
        downloaded_size = file_path.stat().st_size
        if downloaded_size == tcga_file.file_size:
            print(f"✅ Downloaded: {tcga_file.file_name}")
            return True
        else:
            print(f"❌ Size mismatch: expected {tcga_file.file_size}, got {downloaded_size}")
            return False
            
    except Exception as e:
        print(f"❌ Error downloading {tcga_file.file_name}: {e}")
        return False

def main():
    """Main function to test blood cancer file download."""
    print("🧬 TCGA Blood Cancer Data Download Test")
    print("=" * 50)
    
    # Search for blood cancer files
    blood_files = search_blood_cancer_files(max_results=5)
    
    if not blood_files:
        print("⚠️  No blood cancer files found")
        return
    
    print(f"\n📋 Found {len(blood_files)} blood cancer files:")
    for i, file in enumerate(blood_files, 1):
        print(f"{i}. {file.file_name}")
        print(f"   Size: {file.file_size:,} bytes")
        print(f"   Type: {file.file_type}")
        print(f"   Strategy: {file.experimental_strategy}")
    
    # Download first file as test
    if blood_files:
        print(f"\n📥 Downloading first file for testing...")
        success = download_file(blood_files[0], "blood_cancer_downloads")
        
        if success:
            print("🎉 Download test successful!")
            print(f"📁 File saved to: blood_cancer_downloads/{blood_files[0].file_name}")
        else:
            print("❌ Download test failed!")

if __name__ == "__main__":
    main() 