#!/usr/bin/env python3
"""
Test remote genomic data extraction without downloading full files.
This demonstrates how to work with 50GB+ datasets using HTTP range requests.
"""

import os
import time
import subprocess

# Set bcftools cache directory for remote index files
os.environ['HTS_CACHE'] = os.path.join(os.getcwd(), 'remote_index_cache')

# Test with remote 1000 Genomes data - no download needed!
REMOTE_VCF_FILES = {
    # Each file is 200MB - 1.8GB compressed
    "1": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "2": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "3": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "17": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "22": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
}

def test_remote_extraction():
    """Test extracting specific regions from remote VCF files"""
    
    # Create test regions - cancer-related genes across different chromosomes
    test_regions = [
        ("17", 41196312, 41277500, "BRCA1"),    # ~81kb region
        ("17", 7565097, 7590868, "TP53"),       # ~26kb region  
        ("2", 215593487, 215674236, "BARD1"),   # ~81kb region
        ("3", 178866311, 178952497, "PIK3CA"),  # ~86kb region
        ("22", 29121014, 29235591, "EWSR1"),    # ~115kb region
        ("1", 35691274, 35801992, "TP73"),      # ~111kb region
    ]
    
    print("🌐 Testing Remote Genomic Data Extraction")
    print(f"📊 Total remote data available: ~15GB across chromosomes")
    print(f"🎯 Extracting {len(test_regions)} cancer-related gene regions")
    print("="*60)
    
    os.makedirs("remote_test_output", exist_ok=True)
    
    total_start = time.time()
    results = []
    
    for chrom, start, end, gene in test_regions:
        if chrom not in REMOTE_VCF_FILES:
            print(f"⚠️  Skipping {gene} - chr{chrom} not in test set")
            continue
            
        print(f"\n🧬 Extracting {gene} (chr{chrom}:{start:,}-{end:,})")
        
        region_start = time.time()
        remote_url = REMOTE_VCF_FILES[chrom]
        output_file = f"remote_test_output/{gene}.vcf.gz"
        
        # Extract region using bcftools with remote URL
        cmd = [
            "bcftools", "view",
            "-r", f"{chrom}:{start}-{end}",
            "-o", output_file,
            "-Oz",
            remote_url
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Get variant count
            count_cmd = ["bcftools", "view", "-H", output_file]
            count_result = subprocess.run(count_cmd, capture_output=True, text=True, check=True)
            variant_count = len(count_result.stdout.strip().split('\n')) if count_result.stdout.strip() else 0
            
            elapsed = time.time() - region_start
            file_size_kb = os.path.getsize(output_file) / 1024
            
            results.append({
                'gene': gene,
                'chrom': chrom,
                'variants': variant_count,
                'size_kb': file_size_kb,
                'time': elapsed,
                'status': 'success'
            })
            
            print(f"  ✅ Success: {variant_count:,} variants in {elapsed:.2f}s ({file_size_kb:.1f}KB)")
            
        except subprocess.CalledProcessError as e:
            results.append({
                'gene': gene,
                'chrom': chrom,
                'status': 'failed',
                'error': str(e)
            })
            print(f"  ❌ Failed: {e}")
    
    # Summary statistics
    total_time = time.time() - total_start
    successful = [r for r in results if r.get('status') == 'success']
    
    print("\n" + "="*60)
    print("📊 REMOTE EXTRACTION SUMMARY:")
    print("="*60)
    print(f"⏱️  Total time: {total_time:.2f} seconds")
    print(f"✅ Successful extractions: {len(successful)}/{len(test_regions)}")
    
    if successful:
        total_variants = sum(r['variants'] for r in successful)
        total_size = sum(r['size_kb'] for r in successful)
        avg_time = sum(r['time'] for r in successful) / len(successful)
        
        print(f"🧬 Total variants extracted: {total_variants:,}")
        print(f"💾 Total output size: {total_size:.1f} KB ({total_size/1024:.1f} MB)")
        print(f"⚡ Average time per region: {avg_time:.2f} seconds")
        print(f"🌐 Data transferred: ~{total_size*2:.1f} KB (estimated)")
        print(f"\n💡 Note: Only the required genomic regions were transferred,")
        print(f"   not the full {sum([1100,1200,1000,397,196])} MB of chromosome data!")

def verify_bcftools():
    """Check if bcftools is installed and supports remote access"""
    try:
        result = subprocess.run(['bcftools', '--version'], capture_output=True, text=True)
        print(f"✅ bcftools installed: {result.stdout.split()[1]}")
        
        # Check for htslib version (needs to support remote access)
        if 'htslib' in result.stdout:
            print("✅ Remote file access supported")
        return True
    except:
        print("❌ bcftools not found! Install with: brew install bcftools")
        return False

if __name__ == "__main__":
    if verify_bcftools():
        test_remote_extraction()
        
        # Clean up any .tbi files
        from clean_index_files import clean_index_files
        clean_index_files() 