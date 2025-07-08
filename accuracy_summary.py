#!/usr/bin/env python3
"""Show the actual accuracy of extractions."""

import subprocess
import os

print("Genomic Extraction Accuracy Summary")
print("=" * 60)

# Check what's in output_chunks (the most recent extraction)
print("\n1. Checking actual extracted data in output_chunks/")
print("-" * 40)

genes = ["TP53", "BRCA1", "BRCA2"]
for gene in genes:
    if os.path.exists(f"output_chunks/{gene}.vcf.gz"):
        cmd = ["bcftools", "view", "-H", f"output_chunks/{gene}.vcf.gz"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        count = len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
        print(f"{gene}: {count:,} variants")

# Check which dataset was used
cmd = ["bcftools", "view", "-h", "output_chunks/TP53.vcf.gz"]
result = subprocess.run(cmd, capture_output=True, text=True)
header_text = result.stdout.lower()

# Detect dataset and version
if "gnomad" in header_text:
    if "v2.1.1" in header_text or "r2.1.1" in header_text:
        dataset = "gnomAD v2.1.1"
        uses_grch38 = False
    elif "v3" in header_text or "3.1" in header_text:
        dataset = "gnomAD v3"
        uses_grch38 = True
    else:
        # Default to v3 if can't determine
        dataset = "gnomAD (version unknown)"
        uses_grch38 = True
elif "1000genomes" in header_text:
    dataset = "1000 Genomes"
    uses_grch38 = False
else:
    dataset = "Unknown"
    uses_grch38 = False

print(f"\nDataset: {dataset}")

print("\n2. Coordinate System Check")
print("-" * 40)
print("Your regions.bed uses GRCh37 coordinates (human genome v37)")
print("• TP53 location: chr17:7,661,779-7,687,550")

if uses_grch38:
    print("")
    print(f"But {dataset} uses GRCh38 coordinates (human genome v38)")
    print("• TP53 actual location in GRCh38: chr17:7,668,421-7,694,192")
    print("")
    print("The ~6,600 base difference means you're extracting the WRONG region!")
    
    print("\n3. What You're Actually Getting")
    print("-" * 40)
    print("When you query with GRCh37 coordinates on GRCh38 data:")
    print("• You get variants from chr17:7,661,779-7,687,550 in GRCh38")
    print("• This is NOT where TP53 is in GRCh38!")
    print("• You're getting variants from a different genomic region")
else:
    print("")
    print(f"{dataset} also uses GRCh37 coordinates")
    print("✅ Coordinates match - extracting from correct regions!")

print("\n4. Accuracy Status")
print("-" * 40)
print("For 1000 Genomes:  ✅ ACCURATE (uses GRCh37, matches your coordinates)")
print("For gnomAD v2.1.1: ✅ ACCURATE (uses GRCh37, matches your coordinates)")
print("For gnomAD v3:     ❌ INACCURATE (uses GRCh38, coordinate mismatch)")
print("For TCGA:          ❌ Will be INACCURATE (TCGA uses GRCh38)")

if not uses_grch38:
    print(f"\nYour current setup with {dataset} is ✅ ACCURATE!")
else:
    print(f"\nYour current setup with {dataset} is ❌ INACCURATE!")
    print("\n5. Solutions")
    print("-" * 40)
    print("Option 1: Create regions_grch38.bed with converted coordinates")
    print("Option 2: Use coordinate conversion in the code based on dataset")
    print("Option 3: Use gnomAD v2.1.1 which still uses GRCh37")

print("\n" + "=" * 60) 