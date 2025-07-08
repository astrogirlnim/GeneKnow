#!/usr/bin/env python3
"""
Configuration for data sources - easily switch between local, external SSD, and remote
"""

import os
import json
import argparse
import sys

# Data source selection: 'local', 'external', or 'remote'
DATA_SOURCE = os.environ.get('GENOMIC_DATA_SOURCE', 'remote')

# External SSD path - update this to match your SSD
EXTERNAL_SSD_PATH = "/Volumes/YourSSDName/genomic-datasets"

# Chromosome format for each dataset
DATASET_FORMATS = {
    'local': 'without_chr',      # 1000 Genomes Phase 1
    'remote': 'without_chr',     # 1000 Genomes Phase 3
    'external': 'without_chr',   # 1000 Genomes on SSD
    'gnomad': 'with_chr',        # gnomAD v3 uses chr prefix
    'gnomad_v2': 'without_chr',  # gnomAD v2 doesn't use chr prefix
    'tcga': 'with_chr',          # TCGA uses chr prefix
    'user_upload': 'auto_detect' # Auto-detect for user files
}

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
    
    elif DATA_SOURCE == 'gnomad':
        # gnomAD v3.1.2 - Much larger dataset with more variants
        # WARNING: These files are HUGE (4-30GB each)
        return {
            "1": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz",    # ~30GB
            "2": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr2.vcf.bgz",    # ~33GB
            "3": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr3.vcf.bgz",    # ~26GB
            "7": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr7.vcf.bgz",    # ~22GB
            "11": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr11.vcf.bgz",  # ~19GB
            "13": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz",  # ~14GB
            "17": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz",  # ~12GB
            "19": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr19.vcf.bgz",  # ~9GB
            "22": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr22.vcf.bgz",  # ~6GB
            "X": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz",    # ~17GB
            # Total: ~200GB+ of data across these chromosomes!
        }
    
    elif DATA_SOURCE == 'remote':
        # Full remote access - 1000 Genomes Phase 3
        vcf_files = {}
        for chr in list(range(1, 23)):
            vcf_files[str(chr)] = f"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        
        # Special case for chrX - uses v1c instead of v5b
        vcf_files['X'] = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
        
        return vcf_files
    
    elif DATA_SOURCE == 'gnomad_v2':
        # gnomAD v2.1.1 - Uses GRCh37 coordinates (compatible with regions.bed!)
        return {
            "1": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.1.vcf.bgz",
            "17": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.17.vcf.bgz",
            "22": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.22.vcf.bgz",
            # Add more as needed - these use GRCh37 coordinates!
        }
    
    else:
        raise ValueError(f"Unknown data source: {DATA_SOURCE}. Use 'local', 'remote', 'external', 'gnomad', or 'gnomad_v2'")

def get_data_info():
    """Return information about the current data source"""
    if DATA_SOURCE == 'local':
        return "Using local test data (limited chromosomes)"
    elif DATA_SOURCE == 'external':
        return f"Using external SSD data at {EXTERNAL_SSD_PATH}"
    elif DATA_SOURCE == 'gnomad':
        return "Using gnomAD v3.1.2 data (GRCh38 - HUGE files!)"
    elif DATA_SOURCE == 'gnomad_v2':
        return "Using gnomAD v2.1.1 data (GRCh37 - compatible coordinates)"
    else:
        return "Using remote 1000 Genomes data (no local storage needed)"

# Usage example:
def main():
    parser = argparse.ArgumentParser(description='Get genomic data source configuration')
    parser.add_argument('--list-vcf-files', action='store_true', help='List available VCF files')
    parser.add_argument('--get-data-info', action='store_true', help='Get data source information')
    parser.add_argument('--json', action='store_true', help='Output results in JSON format')
    
    args = parser.parse_args()
    
    try:
        if args.list_vcf_files:
            vcf_files = get_vcf_files()
            if args.json:
                print(json.dumps(vcf_files))
            else:
                print(f"🔧 Current data source: {DATA_SOURCE}")
                print(f"📍 {get_data_info()}")
                print(f"\n📊 Available VCF files:")
                for chrom, path in vcf_files.items():
                    print(f"   chr{chrom}: {path}")
        
        elif args.get_data_info:
            info = {
                "data_source": DATA_SOURCE,
                "description": get_data_info(),
                "available_chromosomes": list(get_vcf_files().keys()),
                "chromosome_format": DATASET_FORMATS.get(DATA_SOURCE, 'auto_detect')
            }
            if args.json:
                print(json.dumps(info))
            else:
                print(f"🔧 Current data source: {info['data_source']}")
                print(f"📍 {info['description']}")
                print(f"📊 Available chromosomes: {', '.join(info['available_chromosomes'])}")
                print(f"🧬 Chromosome format: {info['chromosome_format']}")
        
        else:
            # Default behavior - show overview
            if args.json:
                result = {
                    "data_source": DATA_SOURCE,
                    "description": get_data_info(),
                    "available_chromosomes": list(get_vcf_files().keys()),
                    "vcf_files": get_vcf_files(),
                    "chromosome_format": DATASET_FORMATS.get(DATA_SOURCE, 'auto_detect')
                }
                print(json.dumps(result))
            else:
                print(f"🔧 Current data source: {DATA_SOURCE}")
                print(f"📍 {get_data_info()}")
                print(f"\n📊 Available chromosomes: {', '.join(get_vcf_files().keys())}")
                print(f"\n💡 To change source, set environment variable:")
                print(f"   export GENOMIC_DATA_SOURCE=local    # Use test data")
                print(f"   export GENOMIC_DATA_SOURCE=external # Use SSD data") 
                print(f"   export GENOMIC_DATA_SOURCE=remote   # Use 1000 Genomes data (default)")
                print(f"   export GENOMIC_DATA_SOURCE=gnomad   # Use gnomAD data (HUGE files!)")
    
    except Exception as e:
        error_msg = f"Configuration error: {e}"
        if args.json:
            result = {
                "success": False,
                "error": error_msg,
                "data_source": None,
                "available_chromosomes": [],
                "vcf_files": {}
            }
            print(json.dumps(result))
        else:
            print(f"❌ {error_msg}")
        sys.exit(1)

if __name__ == "__main__":
    main() 