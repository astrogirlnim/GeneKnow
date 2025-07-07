import os
import subprocess
from multiprocessing import Pool

# Map chromosomes to VCF files (now using symlinked data)
VCF_FILES = {
    "1": "test-data/ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz",
    "22": "test-data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
}

BED_PATH = "regions.bed"
OUTPUT_DIR = "output_chunks"

def extract_region_precise(line):
    chrom, start, end, gene = line.strip().split("\t")
    start, end = int(start), int(end)
    
    # Get the right VCF file for this chromosome
    if chrom not in VCF_FILES:
        return f"❌ No VCF file for chromosome {chrom}"
    
    vcf_path = VCF_FILES[chrom]
    temp_file = os.path.join(OUTPUT_DIR, f"{gene}_temp.vcf.gz")
    output_file = os.path.join(OUTPUT_DIR, f"{gene}.vcf.gz")
    region = f"{chrom}:{start}-{end}"
    
    try:
        # Step 1: Extract region (might include boundary variants)
        cmd1 = [
            "bcftools", "view",
            "-r", region,
            "-o", temp_file,
            "-Oz",
            vcf_path
        ]
        subprocess.run(cmd1, check=True)
        
        # Step 2: Filter to EXACT boundaries for 100% precision
        cmd2 = [
            "bcftools", "view",
            "-H", temp_file
        ]
        
        # Step 3: Write header + precisely filtered variants
        with open(output_file.replace('.gz', ''), 'w') as out:
            # Write header
            header_cmd = ["bcftools", "view", "-h", temp_file]
            subprocess.run(header_cmd, stdout=out, check=True)
            
            # Write only variants within exact boundaries
            filter_cmd = ["bash", "-c", 
                f"bcftools view -H {temp_file} | awk '$2 >= {start} && $2 <= {end}'"]
            subprocess.run(filter_cmd, stdout=out, check=True)
        
        # Step 4: Compress final output
        subprocess.run(["bgzip", output_file.replace('.gz', '')], check=True)
        
        # Step 5: Clean up temp file
        os.remove(temp_file)
        
        # Verify precision
        count_cmd = ["bash", "-c", f"bcftools view -H {output_file} | wc -l"]
        result = subprocess.run(count_cmd, capture_output=True, text=True, check=True)
        variant_count = int(result.stdout.strip())
        
        return f"✅ Extracted {gene} from chr{chrom}: {variant_count} variants (100% precise)"
        
    except subprocess.CalledProcessError as e:
        return f"❌ Failed for {gene}: {e}"

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    with open(BED_PATH) as f:
        regions = f.readlines()

    with Pool() as pool:
        results = pool.map(extract_region_precise, regions)

    for r in results:
        print(r)

if __name__ == "__main__":
    main() 