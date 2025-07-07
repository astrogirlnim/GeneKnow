#!/usr/bin/env python3
"""
Check sizes of remote genomic data files and estimate transfer speeds
"""

import subprocess
import requests
import time
from config_data_source import get_vcf_files

def get_remote_file_size(url):
    """Get size of remote file using HTTP HEAD request"""
    try:
        response = requests.head(url, allow_redirects=True, timeout=5)
        size = int(response.headers.get('Content-Length', 0))
        return size
    except:
        return None

def format_size(bytes):
    """Convert bytes to human readable format"""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes < 1024.0:
            return f"{bytes:.1f} {unit}"
        bytes /= 1024.0
    return f"{bytes:.1f} TB"

def check_remote_sizes():
    """Check sizes of all remote VCF files"""
    print("🔍 Checking Remote Genomic Data Sizes...")
    print("="*60)
    
    vcf_files = get_vcf_files()
    
    # Filter only remote files
    remote_files = {k: v for k, v in vcf_files.items() if v.startswith('http')}
    
    if not remote_files:
        print("No remote files configured. Set GENOMIC_DATA_SOURCE=remote")
        return
    
    total_size = 0
    file_sizes = {}
    
    print(f"Checking {len(remote_files)} remote files...\n")
    
    for chrom, url in sorted(remote_files.items(), key=lambda x: (x[0].isdigit(), int(x[0]) if x[0].isdigit() else 999)):
        print(f"Chr {chrom}: ", end='', flush=True)
        size = get_remote_file_size(url)
        
        if size:
            file_sizes[chrom] = size
            total_size += size
            print(f"{format_size(size)}")
        else:
            print("Unable to get size")
    
    print("\n" + "="*60)
    print(f"📊 TOTAL DATASET SIZE: {format_size(total_size)}")
    print("="*60)
    
    # Calculate transfer estimates
    print("\n⏱️  ESTIMATED EXTRACTION TIMES:")
    print("-"*40)
    
    # Based on your test results: ~5MB transferred in ~50 seconds
    observed_speed = 5 * 1024 * 1024 / 50  # bytes per second
    
    print(f"Based on your tests: {format_size(observed_speed)}/s effective speed")
    print(f"(Only transferring needed regions, not full files)\n")
    
    # Example gene regions
    example_genes = [
        ("BRCA1", 81_000),    # 81kb region
        ("TP53", 26_000),     # 26kb region
        ("Full 33 genes", 2_000_000),  # ~2MB total
    ]
    
    for name, region_size in example_genes:
        time_est = region_size / observed_speed
        print(f"{name}: ~{time_est:.1f} seconds")
    
    return file_sizes

def test_transfer_speed(chrom="22", start=29121014, end=29235591):
    """Test actual transfer speed for a region"""
    print(f"\n🚀 Testing actual transfer speed for chr{chrom}:{start}-{end}...")
    
    vcf_files = get_vcf_files()
    if chrom not in vcf_files:
        print(f"Chromosome {chrom} not available")
        return
    
    url = vcf_files[chrom]
    
    # Time the extraction
    start_time = time.time()
    
    cmd = [
        "bcftools", "view",
        "-H",  # Headers only, no data
        "-r", f"{chrom}:{start}-{end}",
        url
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    elapsed = time.time() - start_time
    
    if result.returncode == 0:
        # Count variants
        variants = len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
        data_size = len(result.stdout.encode())
        
        print(f"✅ Transferred {format_size(data_size)} in {elapsed:.2f}s")
        print(f"   Speed: {format_size(data_size/elapsed)}/s")
        print(f"   Variants: {variants}")
    else:
        print(f"❌ Test failed: {result.stderr}")

def add_more_datasets():
    """Show how to add more genomic datasets"""
    print("\n📚 ADDING MORE GENOMIC DATASETS")
    print("="*60)
    
    print("1️⃣ gnomAD (Genome Aggregation Database):")
    print("   Larger dataset with more population diversity")
    print("   Example URLs:")
    print("   https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz")
    
    print("\n2️⃣ ClinVar (Clinical Variants):")
    print("   Disease-associated variants")
    print("   https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz")
    
    print("\n3️⃣ ExAC (Exome Aggregation Consortium):")
    print("   60,706 exome sequences")
    print("   https://storage.googleapis.com/gnomad-public/legacy/exac_browser/ExAC.r1.sites.vep.vcf.gz")
    
    print("\n📝 To add these to your pipeline:")
    print("1. Edit config_data_source.py")
    print("2. Add URLs to the VCF_FILES dictionary")
    print("3. Ensure your BED file regions match the reference genome (GRCh37 vs GRCh38)")
    
    print("\n💡 Example adding gnomAD:")
    print('VCF_FILES["17_gnomad"] = "https://gnomad-public-us-east-1.s3.amazonaws.com/..."')

if __name__ == "__main__":
    # Set to remote mode for this check
    import os
    os.environ['GENOMIC_DATA_SOURCE'] = 'remote'
    
    file_sizes = check_remote_sizes()
    test_transfer_speed()
    add_more_datasets() 