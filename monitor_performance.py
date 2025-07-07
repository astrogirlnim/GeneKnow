#!/usr/bin/env python3
"""
Monitor extraction performance in real-time
"""

import time
import os
import subprocess
from datetime import datetime

def monitor_extraction(regions_file="regions.bed"):
    """Monitor extraction performance with detailed metrics"""
    
    print("📊 GENOMIC EXTRACTION PERFORMANCE MONITOR")
    print("="*60)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Regions file: {regions_file}")
    print(f"Data source: {os.environ.get('GENOMIC_DATA_SOURCE', 'remote')}")
    print("="*60)
    
    # Count regions
    with open(regions_file) as f:
        regions = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
    print(f"📍 Regions to extract: {len(regions)}")
    
    # Check output directory
    output_dir = "output_chunks"
    existing_files = len([f for f in os.listdir(output_dir) if f.endswith('.vcf.gz')]) if os.path.exists(output_dir) else 0
    
    print(f"📁 Existing output files: {existing_files}")
    print("\n🚀 Starting extraction...\n")
    
    # Run extraction with real-time monitoring
    start_time = time.time()
    
    process = subprocess.Popen(
        ["python3", "extract_by_region.py"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )
    
    # Monitor output
    completed = 0
    failed = 0
    
    for line in process.stdout:
        line = line.strip()
        if line:
            print(line)
            
            # Track progress
            if "✅" in line:
                completed += 1
                elapsed = time.time() - start_time
                rate = completed / elapsed
                eta = (len(regions) - completed) / rate if rate > 0 else 0
                
                print(f"   📈 Progress: {completed}/{len(regions)} | "
                      f"Rate: {rate:.2f} genes/sec | "
                      f"ETA: {eta:.0f}s")
            elif "❌" in line:
                failed += 1
    
    process.wait()
    total_time = time.time() - start_time
    
    # Final summary
    print("\n" + "="*60)
    print("📊 PERFORMANCE SUMMARY")
    print("="*60)
    
    # Check actual output
    new_files = len([f for f in os.listdir(output_dir) if f.endswith('.vcf.gz')]) if os.path.exists(output_dir) else 0
    created = new_files - existing_files
    
    # Calculate sizes
    total_size = 0
    if os.path.exists(output_dir):
        for f in os.listdir(output_dir):
            if f.endswith('.vcf.gz'):
                size = os.path.getsize(os.path.join(output_dir, f))
                total_size += size
    
    print(f"⏱️  Total time: {total_time:.1f} seconds")
    print(f"✅ Successful: {completed}")
    print(f"❌ Failed: {failed}")
    print(f"📁 Files created: {created}")
    print(f"💾 Total output size: {total_size/1024/1024:.1f} MB")
    print(f"⚡ Average time per gene: {total_time/len(regions):.1f}s")
    print(f"📊 Extraction rate: {completed/total_time:.2f} genes/sec")
    
    # Data source info
    if os.environ.get('GENOMIC_DATA_SOURCE') == 'remote':
        print(f"\n🌐 Remote extraction performance:")
        print(f"   Effective throughput: {total_size/total_time/1024:.1f} KB/s")
        print(f"   Network efficiency: Only transferred needed regions!")

def compare_sources():
    """Compare performance across different data sources"""
    print("\n🔄 COMPARING DATA SOURCES")
    print("="*60)
    
    # Test with 3 genes
    test_bed = "test_performance.bed"
    with open("regions_large.bed") as f:
        test_regions = f.readlines()[:3]
    
    with open(test_bed, 'w') as f:
        f.writelines(test_regions)
    
    sources = ['local', 'remote']
    results = {}
    
    for source in sources:
        if source == 'local' and not os.path.exists("test-data/ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"):
            print(f"⚠️  Skipping {source} - no local data")
            continue
            
        print(f"\n📊 Testing {source.upper()} source...")
        os.environ['GENOMIC_DATA_SOURCE'] = source
        
        start = time.time()
        subprocess.run(["python3", "extract_by_region.py"], 
                      capture_output=True, text=True,
                      env={**os.environ, 'BED_PATH': test_bed})
        elapsed = time.time() - start
        
        results[source] = elapsed
        print(f"   Time: {elapsed:.1f}s")
    
    # Cleanup
    os.remove(test_bed)
    
    if len(results) >= 2:
        speedup = results['remote'] / results['local']
        print(f"\n📈 Local is {speedup:.1f}x faster than remote")

if __name__ == "__main__":
    # Run performance monitoring
    monitor_extraction()
    
    # Optional: compare sources
    if input("\nCompare data sources? (y/n): ").lower() == 'y':
        compare_sources() 