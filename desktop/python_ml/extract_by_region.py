#!/usr/bin/env python3
"""
Extract Genomic Regions
Extracts specific genomic regions from VCF files using BED file input
"""

import os
import subprocess
import time
import resource
import json
import argparse
import sys
from multiprocessing import Pool

# Set bcftools cache directory for remote index files
os.environ['HTS_CACHE'] = os.path.join(os.getcwd(), 'remote_index_cache')

# Import configuration
from config_data_source import get_vcf_files, get_data_info, DATASET_FORMATS, DATA_SOURCE

class GenomicExtractor:
    def __init__(self, bed_file, output_dir, max_processes=None, json_output=False):
        self.bed_file = bed_file
        self.output_dir = output_dir
        self.max_processes = max_processes
        self.json_output = json_output
        
        # Performance Configuration
        self.config = {
            "max_processes": max_processes,  # None = use all CPU cores
            "compression_level": 6,  # 1-9, higher = better compression but slower
            "temp_cleanup": True,    # Clean up temporary files
            "verbose_timing": True,  # Show detailed timing breakdown
            "memory_monitoring": True,  # Track memory usage per extraction
            "continue_on_error": True,  # Continue processing other regions if one fails
            "retry_failed": 1,      # Number of retries for failed extractions
        }
        
        # Get VCF files from configuration
        self.vcf_files = get_vcf_files()
    
    def log(self, message):
        """Log a message (only if not in JSON mode)"""
        if not self.json_output:
            print(message)
    
    def normalize_chromosome(self, chrom, target_format):
        """Convert between chr1 and 1 formats as needed."""
        if target_format == 'with_chr':
            return f"chr{chrom}" if not chrom.startswith('chr') else chrom
        elif target_format == 'without_chr':
            return chrom.replace('chr', '')
        else:  # auto_detect or unknown
            return chrom

    def get_chromosome_format(self):
        """Get the expected chromosome format for current data source."""
        return DATASET_FORMATS.get(DATA_SOURCE, 'auto_detect')

    def extract_region_precise(self, line):
        """Extract genomic region with retries and error resilience"""
        chrom, start, end, gene = line.strip().split("\t")
        start, end = int(start), int(end)
        
        for attempt in range(self.config["retry_failed"] + 1):
            try:
                return self._extract_region_attempt(chrom, start, end, gene, attempt)
            except Exception as e:
                if attempt < self.config["retry_failed"]:
                    self.log(f"⚠️  {gene} attempt {attempt + 1} failed, retrying: {e}")
                    time.sleep(1)  # Brief pause before retry
                    continue
                else:
                    return f"❌ {gene} failed after {attempt + 1} attempts: {e}"

    def _extract_region_attempt(self, chrom, start, end, gene, attempt=0):
        """Single extraction attempt with performance monitoring"""
        # Performance monitoring
        start_time = time.time()
        initial_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024  # MB
        
        # Normalize chromosome for VCF file lookup (remove 'chr' prefix)
        chrom_key = chrom.replace('chr', '')
        
        # Get the right VCF file for this chromosome
        if chrom_key not in self.vcf_files:
            return f"❌ No VCF file for chromosome {chrom}"
        
        vcf_path = self.vcf_files[chrom_key]
        temp_file = os.path.join(self.output_dir, f"{gene}_temp.vcf.gz")
        output_file = os.path.join(self.output_dir, f"{gene}.vcf.gz")
        
        # Normalize chromosome name for the query based on dataset format
        target_format = self.get_chromosome_format()
        chrom_for_query = self.normalize_chromosome(chrom, target_format)
        region = f"{chrom_for_query}:{start}-{end}"
        
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
            subprocess.run(["bgzip", "-f", "-l", str(self.config["compression_level"]), output_file.replace('.gz', '')], check=True)
            compress_time = time.time() - compress_start
            
            # Step 5: Clean up temp file
            if self.config["temp_cleanup"]:
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
            
            return (f"✅ {gene} ({chrom}): {variant_count} variants, {file_size:.1f}KB | "
                    f"Time: {total_time:.2f}s (extract: {extract_time:.2f}s, filter: {filter_time:.2f}s, compress: {compress_time:.2f}s) | "
                    f"Memory: +{memory_used:.1f}MB")
            
        except subprocess.CalledProcessError as e:
            return f"❌ Failed for {gene}: {e}"

    def run_extraction(self):
        """Main extraction process with performance monitoring"""
        self.log("🧬 Starting parallel genomic extraction with performance monitoring...")
        self.log(f"📍 Data source: {get_data_info()}")
        
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # Input validation
        if not os.path.exists(self.bed_file):
            error_msg = f"BED file not found: {self.bed_file}"
            if self.json_output:
                result = {
                    "success": False,
                    "error": error_msg,
                    "total_regions": 0,
                    "successful_extractions": 0,
                    "total_variants": 0,
                    "total_size": 0,
                    "execution_time": 0,
                    "results": []
                }
                print(json.dumps(result))
                return
            else:
                self.log(f"❌ ERROR: {error_msg}")
                return
                
        with open(self.bed_file) as f:
            regions = f.readlines()
        
        # Validate regions
        if not regions:
            error_msg = "BED file is empty"
            if self.json_output:
                result = {
                    "success": False,
                    "error": error_msg,
                    "total_regions": 0,
                    "successful_extractions": 0,
                    "total_variants": 0,
                    "total_size": 0,
                    "execution_time": 0,
                    "results": []
                }
                print(json.dumps(result))
                return
            else:
                self.log(f"❌ ERROR: {error_msg}")
                return
            
        # Parse and validate each region
        valid_regions = []
        for i, line in enumerate(regions, 1):
            line = line.strip()
            if not line or line.startswith('#'):  # Skip empty lines and comments
                continue
                
            parts = line.split('\t')
            if len(parts) != 4:
                self.log(f"⚠️  WARNING: Line {i} has {len(parts)} columns, expected 4: {line}")
                continue
                
            chrom, start, end, gene = parts
            # Normalize chromosome for VCF file lookup
            chrom_key = chrom.replace('chr', '')
            if chrom_key not in self.vcf_files:
                self.log(f"⚠️  WARNING: No VCF file available for chromosome {chrom} (gene: {gene})")
                continue
                
            valid_regions.append(line)
        
        if not valid_regions:
            error_msg = "No valid regions found to process"
            if self.json_output:
                result = {
                    "success": False,
                    "error": error_msg,
                    "total_regions": 0,
                    "successful_extractions": 0,
                    "total_variants": 0,
                    "total_size": 0,
                    "execution_time": 0,
                    "results": []
                }
                print(json.dumps(result))
                return
            else:
                self.log(f"❌ ERROR: {error_msg}")
                return

        # Report what will be processed
        chromosomes = set(line.split('\t')[0] for line in valid_regions)
        gene_names = [line.split('\t')[3] for line in valid_regions]
        self.log(f"📊 Input validation: {len(valid_regions)} valid regions across {len(chromosomes)} chromosomes")
        self.log(f"🧬 Regions to process: {', '.join(gene_names)}")
        self.log(f"🖥️  Using {self.config['max_processes'] or os.cpu_count()}/{os.cpu_count()} CPU cores (compression level: {self.config['compression_level']})")
        
        # Overall performance tracking
        pipeline_start = time.time()
        
        with Pool(processes=self.config["max_processes"]) as pool:
            results = pool.map(self.extract_region_precise, valid_regions)

        pipeline_time = time.time() - pipeline_start
        
        # Process results
        total_variants = 0
        total_size = 0
        successful_extractions = 0
        
        for r in results:
            if not self.json_output:
                self.log(r)
            if "✅" in r:
                successful_extractions += 1
                # Extract metrics from result string
                if "variants" in r:
                    variants = int(r.split("variants")[0].split(": ")[-1])
                    total_variants += variants
                if "KB" in r:
                    size = float(r.split("KB")[0].split(", ")[-1])
                    total_size += size
        
        if self.json_output:
            # Return structured JSON result
            result = {
                "success": True,
                "total_regions": len(valid_regions),
                "successful_extractions": successful_extractions,
                "total_variants": total_variants,
                "total_size": int(total_size),  # KB as integer
                "execution_time": pipeline_time,
                "results": results,
                "error": None
            }
            print(json.dumps(result))
        else:
            # Print performance summary
            self.log("\n" + "="*80)
            self.log("EXTRACTION RESULTS:")
            self.log("="*80)
            
            self.log("\n" + "="*80)
            self.log("PIPELINE PERFORMANCE SUMMARY:")
            self.log("="*80)
            self.log(f"⏱️  Total time: {pipeline_time:.2f} seconds")
            self.log(f"📊 Successful extractions: {successful_extractions}/{len(valid_regions)}")
            self.log(f"🧬 Total variants extracted: {total_variants:,}")
            self.log(f"💾 Total output size: {total_size:.1f} KB ({total_size/1024:.1f} MB)")
            self.log(f"⚡ Average time per region: {pipeline_time/len(valid_regions):.2f} seconds")
            self.log(f"🏃 Variants per second: {total_variants/pipeline_time:.0f}")
            self.log("="*80)

def main():
    parser = argparse.ArgumentParser(description='Extract genomic regions from VCF files')
    parser.add_argument('--bed-file', required=True, help='BED file with regions to extract')
    parser.add_argument('--output-dir', required=True, help='Output directory for extracted regions')
    parser.add_argument('--max-processes', type=int, help='Maximum number of processes to use')
    parser.add_argument('--json', action='store_true', help='Output results in JSON format')
    
    args = parser.parse_args()
    
    try:
        extractor = GenomicExtractor(
            bed_file=args.bed_file,
            output_dir=args.output_dir,
            max_processes=args.max_processes,
            json_output=args.json
        )
        extractor.run_extraction()
    except Exception as e:
        error_msg = f"Extraction failed: {e}"
        if args.json:
            result = {
                "success": False,
                "error": error_msg,
                "total_regions": 0,
                "successful_extractions": 0,
                "total_variants": 0,
                "total_size": 0,
                "execution_time": 0,
                "results": []
            }
            print(json.dumps(result))
        else:
            print(f"❌ {error_msg}")
        sys.exit(1)

if __name__ == "__main__":
    main()
    
    # Clean up any .tbi files that ended up in root
    # TODO: Add cleanup function if needed 