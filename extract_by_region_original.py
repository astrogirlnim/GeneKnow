import os
import subprocess
from multiprocessing import Pool

# Map chromosomes to VCF files
VCF_FILES = {
    "1": "test-data/ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz",
    "22": "test-data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
}

BED_PATH = "regions.bed"
OUTPUT_DIR = "output_chunks"

def extract_region(line):
    chrom, start, end, gene = line.strip().split("\t")
    
    # Get the right VCF file for this chromosome
    if chrom not in VCF_FILES:
        return f"❌ No VCF file for chromosome {chrom}"
    
    vcf_path = VCF_FILES[chrom]
    output_file = os.path.join(OUTPUT_DIR, f"{gene}.vcf.gz")
    region = f"{chrom}:{start}-{end}"
    cmd = [
        "bcftools", "view",
        "-r", region,
        "-o", output_file,
        "-Oz",
        vcf_path
    ]
    try:
        subprocess.run(cmd, check=True)
        return f"✅ Extracted {gene} from chr{chrom}"
    except subprocess.CalledProcessError as e:
        return f"❌ Failed for {gene}: {e}"

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    with open(BED_PATH) as f:
        regions = f.readlines()

    with Pool() as pool:
        results = pool.map(extract_region, regions)

    for r in results:
        print(r)

if __name__ == "__main__":
    main()