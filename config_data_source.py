#!/usr/bin/env python3
"""
Configuration for data sources - easily switch between local, external SSD, and remote
"""

import os

# Data source selection: 'local', 'external', or 'remote'
DATA_SOURCE = os.environ.get('GENOMIC_DATA_SOURCE', 'remote')

# External SSD path - update this to match your SSD
EXTERNAL_SSD_PATH = "/Volumes/YourSSDName/genomic-datasets"

def get_vcf_files():
    """Return the appropriate VCF file mapping based on configuration"""
    
    if DATA_SOURCE == 'local':
        # Original local test data
        return {
            "1": "test-data/ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz",
            "22": "test-data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        }
    
    elif DATA_SOURCE == 'external':
        # External SSD with full dataset
        base = f"{EXTERNAL_SSD_PATH}/1000genomes/phase3"
        vcf_files = {}
        
        # All chromosomes 1-22 + X
        for chr in list(range(1, 23)) + ['X']:
            vcf_files[str(chr)] = f"{base}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        
        return vcf_files
    
    else:  # remote
        # Full remote access
        vcf_files = {}
        for chr in list(range(1, 23)) + ['X']:
            vcf_files[str(chr)] = f"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        
        return vcf_files

def get_data_info():
    """Return information about the current data source"""
    if DATA_SOURCE == 'local':
        return "Using local test data (limited chromosomes)"
    elif DATA_SOURCE == 'external':
        return f"Using external SSD data at {EXTERNAL_SSD_PATH}"
    else:
        return "Using remote 1000 Genomes data (no local storage needed)"

# Usage example:
if __name__ == "__main__":
    print(f"🔧 Current data source: {DATA_SOURCE}")
    print(f"📍 {get_data_info()}")
    print(f"\n📊 Available chromosomes: {', '.join(get_vcf_files().keys())}")
    print(f"\n💡 To change source, set environment variable:")
    print(f"   export GENOMIC_DATA_SOURCE=local    # Use test data")
    print(f"   export GENOMIC_DATA_SOURCE=external # Use SSD data") 
    print(f"   export GENOMIC_DATA_SOURCE=remote   # Use remote data") 