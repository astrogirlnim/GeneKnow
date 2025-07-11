# FASTQ to VCF Integration Guide for GeneKnow

## Overview

This guide shows how to integrate FASTQ processing into your genomic pipeline, allowing you to process raw sequencing data alongside the existing VCF extraction capabilities.

## Pipeline Flow

```
Patient FASTQ → Alignment → Variant Calling → VCF → Risk Analysis → Report
                                                  ↑
                                        Existing Pipeline
```

## Quick Start

### 1. Install Required Tools

```bash
# Required tools
brew install bwa samtools bcftools

# Optional aligners (choose one)
brew install minimap2  # Faster, good for long reads
brew install bowtie2   # Alternative to BWA
```

### 2. Test Data Sources (UPDATED - Working URLs)

#### A. Create Synthetic Test Data (Simplest Option)
```bash
# Create test directory
mkdir -p test-fastq && cd test-fastq

# Generate synthetic FASTQ for quick testing
cat > generate_test_fastq.py << 'EOF'
import random
import gzip

def generate_fastq(filename, num_reads=1000, read_length=150):
    """Generate a simple test FASTQ file"""
    bases = ['A', 'T', 'G', 'C']
    with gzip.open(filename, 'wt') as f:
        for i in range(num_reads):
            # Generate random sequence
            seq = ''.join(random.choice(bases) for _ in range(read_length))
            # Generate quality scores (all high quality for simplicity)
            qual = 'I' * read_length
            # Write FASTQ entry
            f.write(f"@read_{i}\n{seq}\n+\n{qual}\n")

# Generate paired-end reads
generate_fastq('test_R1.fastq.gz')
generate_fastq('test_R2.fastq.gz')
print("Generated test FASTQ files!")
EOF

python3 generate_test_fastq.py

# Download a small reference (E. coli K-12, ~4.6MB)
curl -o ecoli_ref.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip ecoli_ref.fa.gz
mv GCF_000005845.2_ASM584v2_genomic.fna ecoli_ref.fa
```

#### B. Download from SRA using SRA Toolkit (Most Reliable)
```bash
# Install SRA toolkit if not already installed
brew install sra-tools

# Download a small E. coli dataset (SRR390728, ~35MB total)
fastq-dump --split-files --gzip SRR390728 --maxSpotId 100000

# This creates (limited to first 100k reads for speed):
# - SRR390728_1.fastq.gz
# - SRR390728_2.fastq.gz
```

#### C. Direct Download - Verified Working URLs
```bash
# From EBI (European Bioinformatics Institute)
# Small E. coli dataset
wget https://www.ebi.ac.uk/ena/browser/api/fasta/SRR390728?download=true -O SRR390728.fasta
# Note: This gives FASTA not FASTQ, but useful for reference

# Alternative: Use ENA Browser API for FASTQ
curl -X GET "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR390728&result=read_run&fields=fastq_ftp&format=tsv&download=true" -H "accept: */*"
# This returns FTP URLs for the FASTQ files
```

#### D. GitHub Test Data Repositories
```bash
# Option 1: From viral genomes repo (small files)
wget https://raw.githubusercontent.com/zuherJahshan/common_viruses/main/sars-cov-2.fna
# Note: This is FASTA format, not FASTQ

# Option 2: Generate from existing test VCF data
# You already have VCF files - you could simulate reads from them
```

### 3. Run FASTQ to VCF Pipeline

```bash
# Basic usage
python3 fastq_to_vcf_pipeline.py \
  -r reference.fa \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o patient_001 \
  -t 8

# This produces: patient_001.vcf
```

### 4. Quick Test with E. coli

```bash
# Complete test workflow
cd test-fastq

# Run pipeline
python3 ../fastq_to_vcf_pipeline.py \
  -r ecoli_ref.fa \
  -1 ecoli_1.fastq.gz \
  -2 ecoli_2.fastq.gz \
  -o ecoli_test \
  -t 4

# Check output
bcftools stats ecoli_test.vcf
``` 