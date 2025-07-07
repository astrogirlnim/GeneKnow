import os
import subprocess
import time
import resource
from multiprocessing import Pool

# Performance Configuration
CONFIG = {
    "max_processes": None,  # None = use all CPU cores, or set to specific number
    "compression_level": 6,  # 1-9, higher = better compression but slower
    "temp_cleanup": True,    # Clean up temporary files
    "verbose_timing": True,  # Show detailed timing breakdown
    "memory_monitoring": True,  # Track memory usage per extraction
    "continue_on_error": True,  # Continue processing other regions if one fails
    "retry_failed": 1,      # Number of retries for failed extractions
}

# Map chromosomes to VCF files (now using symlinked data)
VCF_FILES = {
    "1": "test-data/ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz",
    "22": "test-data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
}

BED_PATH = "regions.bed"
OUTPUT_DIR = "output_chunks"

def extract_region_precise(line):
    """Extract genomic region with retries and error resilience"""
    chrom, start, end, gene = line.strip().split("\t")
    start, end = int(start), int(end)
    
    for attempt in range(CONFIG["retry_failed"] + 1):
        try:
            return _extract_region_attempt(chrom, start, end, gene, attempt)
        except Exception as e:
            if attempt < CONFIG["retry_failed"]:
                print(f"⚠️  {gene} attempt {attempt + 1} failed, retrying: {e}")
                time.sleep(1)  # Brief pause before retry
                continue
            else:
                return f"❌ {gene} failed after {attempt + 1} attempts: {e}"

def _extract_region_attempt(chrom, start, end, gene, attempt=0):
    """Single extraction attempt with performance monitoring"""
    # Performance monitoring
    start_time = time.time()
    initial_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024  # MB
    
    # Get the right VCF file for this chromosome
    if chrom not in VCF_FILES:
        return f"❌ No VCF file for chromosome {chrom}"
    
    vcf_path = VCF_FILES[chrom]
    temp_file = os.path.join(OUTPUT_DIR, f"{gene}_temp.vcf.gz")
    output_file = os.path.join(OUTPUT_DIR, f"{gene}.vcf.gz")
    region = f"{chrom}:{start}-{end}"
    
    try:
        # Step 1: Extract region (might include boundary variants)
        extract_start = time.time()
        cmd1 = [
            "bcftools", "view",
            "-r", region,
            "-o", temp_file,
            "-Oz",
            vcf_path
        ]
        subprocess.run(cmd1, check=True)
        extract_time = time.time() - extract_start
        
        # Step 2: Filter to EXACT boundaries for 100% precision
        filter_start = time.time()
        
        # Step 3: Write header + precisely filtered variants
        with open(output_file.replace('.gz', ''), 'w') as out:
            # Write header
            header_cmd = ["bcftools", "view", "-h", temp_file]
            subprocess.run(header_cmd, stdout=out, check=True)
            
            # Write only variants within exact boundaries
            filter_cmd = ["bash", "-c", 
                f"bcftools view -H {temp_file} | awk '$2 >= {start} && $2 <= {end}'"]
            subprocess.run(filter_cmd, stdout=out, check=True)
        
        filter_time = time.time() - filter_start
        
        # Step 4: Compress final output
        compress_start = time.time()
        subprocess.run(["bgzip", "-f", "-l", str(CONFIG["compression_level"]), output_file.replace('.gz', '')], check=True)
        compress_time = time.time() - compress_start
        
        # Step 5: Clean up temp file
        if CONFIG["temp_cleanup"]:
            os.remove(temp_file)
        
        # Performance metrics
        total_time = time.time() - start_time
        peak_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024  # MB
        memory_used = peak_memory - initial_memory
        
        # Verify precision and get stats
        count_cmd = ["bash", "-c", f"bcftools view -H {output_file} | wc -l"]
        result = subprocess.run(count_cmd, capture_output=True, text=True, check=True)
        variant_count = int(result.stdout.strip())
        
        # Get file size
        file_size = os.path.getsize(output_file) / 1024  # KB
        
        return (f"✅ {gene} (chr{chrom}): {variant_count} variants, {file_size:.1f}KB | "
                f"Time: {total_time:.2f}s (extract: {extract_time:.2f}s, filter: {filter_time:.2f}s, compress: {compress_time:.2f}s) | "
                f"Memory: +{memory_used:.1f}MB")
        
    except subprocess.CalledProcessError as e:
        return f"❌ Failed for {gene}: {e}"

def main():
    """Main extraction process with performance monitoring"""
    print("🧬 Starting parallel genomic extraction with performance monitoring...")
    
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Input validation
    if not os.path.exists(BED_PATH):
        print(f"❌ ERROR: regions.bed file not found at {BED_PATH}")
        return
        
    with open(BED_PATH) as f:
        regions = f.readlines()
    
    # Validate regions
    if not regions:
        print("❌ ERROR: regions.bed file is empty")
        return
        
    # Parse and validate each region
    valid_regions = []
    for i, line in enumerate(regions, 1):
        line = line.strip()
        if not line or line.startswith('#'):  # Skip empty lines and comments
            continue
            
        parts = line.split('\t')
        if len(parts) != 4:
            print(f"⚠️  WARNING: Line {i} has {len(parts)} columns, expected 4: {line}")
            continue
            
        chrom, start, end, gene = parts
        if chrom not in VCF_FILES:
            print(f"⚠️  WARNING: No VCF file available for chromosome {chrom} (gene: {gene})")
            continue
            
        valid_regions.append(line)
    
    if not valid_regions:
        print("❌ ERROR: No valid regions found to process")
        return

    # Report what will be processed
    chromosomes = set(line.split('\t')[0] for line in valid_regions)
    gene_names = [line.split('\t')[3] for line in valid_regions]
    print(f"📊 Input validation: {len(valid_regions)} valid regions across {len(chromosomes)} chromosomes")
    print(f"🧬 Regions to process: {', '.join(gene_names)}")
    print(f"🖥️  Using {CONFIG['max_processes'] or os.cpu_count()}/{os.cpu_count()} CPU cores (compression level: {CONFIG['compression_level']})")
    
    # Overall performance tracking
    pipeline_start = time.time()
    
    with Pool(processes=CONFIG["max_processes"]) as pool:
        results = pool.map(extract_region_precise, valid_regions)

    pipeline_time = time.time() - pipeline_start
    
    # Print results with performance summary
    print("\n" + "="*80)
    print("EXTRACTION RESULTS:")
    print("="*80)
    
    total_variants = 0
    total_size = 0
    successful_extractions = 0
    
    for r in results:
        print(r)
        if "✅" in r:
            successful_extractions += 1
            # Extract metrics from result string
            if "variants" in r:
                variants = int(r.split("variants")[0].split(": ")[-1])
                total_variants += variants
            if "KB" in r:
                size = float(r.split("KB")[0].split(", ")[-1])
                total_size += size
    
    print("\n" + "="*80)
    print("PIPELINE PERFORMANCE SUMMARY:")
    print("="*80)
    print(f"⏱️  Total time: {pipeline_time:.2f} seconds")
    print(f"📊 Successful extractions: {successful_extractions}/{len(valid_regions)}")
    print(f"🧬 Total variants extracted: {total_variants:,}")
    print(f"💾 Total output size: {total_size:.1f} KB ({total_size/1024:.1f} MB)")
    print(f"⚡ Average time per region: {pipeline_time/len(valid_regions):.2f} seconds")
    print(f"🏃 Variants per second: {total_variants/pipeline_time:.0f}")
    print("="*80)

if __name__ == "__main__":
    main() 