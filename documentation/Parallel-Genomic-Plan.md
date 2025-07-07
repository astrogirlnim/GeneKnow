# GenePredict: Parallel Genomic Region Extraction Pipeline

"""
This document outlines the architecture and plan for building a parallelized genomic data processor in Cursor. The goal is to split large TCGA VCF files by region (e.g. cancer-relevant genes) and process each section in parallel. We aim to do this entirely locally using Python scripting, and call external tools like `bcftools` via subprocess.
"""

# --------------------------------------
# 🛠 External Setup (Required Once)
# --------------------------------------

"""
Install the following CLI tools locally:

macOS:
    brew install bcftools samtools

Linux:
    sudo apt install bcftools samtools

These will be called via Python subprocess, so ensure they are available in your PATH.
"""

# --------------------------------------
# 📁 Project Structure in Cursor
# --------------------------------------

"""
project/
├── extract_by_region.py         # Main script (entrypoint)
├── regions.bed                  # BED file with gene regions (chr, start, end, gene)
├── utils.py                     # Helper functions (e.g., prioritization, logging)
└── output_chunks/               # Folder for generated region-specific VCF files
"""

# --------------------------------------
# 🧠 System Components
# --------------------------------------

"""
1. BED Region Loader
   - Load and parse `regions.bed`
   - Optionally sort/prioritize regions by gene importance

2. Region Extractor
   - For each region, call:
       bcftools view -r chr:start-end -o output.vcf.gz -Oz input.vcf.gz

3. Parallel Processor
   - Use Python multiprocessing or concurrent.futures to parallelize extraction

4. Optional Post-Processing
   - Add downstream logic to parse each chunk or hand off to ML/LLM components

5. CLI Script Runner
   - Accepts args like `--input`, `--regions`, `--priority`, `--output-dir`
"""

# --------------------------------------
# ✅ Summary for Cursor Agent
# --------------------------------------

"""
We are building a Python-based pipeline to:
1. Read a BED file with gene regions
2. Split a VCF file into smaller region-based chunks using `bcftools`
3. Run all extractions in parallel to improve performance
4. Prioritize certain genes based on cancer associations (e.g. BRCA1, TP53)

Cursor should:
- Scaffold a Python project with the above layout
- Generate functions to extract regions from a VCF
- Add multiprocessing logic for parallel execution
- Handle edge cases: invalid regions, missing inputs, tool not found
- Add optional CLI argument parsing for flexibility
"""

# --------------------------------------
# 🧩 Detailed Implementation Plan
# --------------------------------------

## Files and Responsibilities

### extract_by_region.py
- [ ] Parse `regions.bed`
- [ ] For each line: extract `chr`, `start`, `end`, `gene`
- [ ] Generate output path: `output_chunks/{gene}.vcf.gz`
- [ ] Call `bcftools view -r chr:start-end -o out -Oz input.vcf.gz` via subprocess
- [ ] Use multiprocessing Pool to do all in parallel
- [ ] Accept CLI args: --input, --regions, --output-dir, --priority

### utils.py
- [ ] Load BED file and return list of (chr, start, end, gene)
- [ ] Prioritize list using provided gene priority map
- [ ] Optional: validate that bcftools is installed and callable

## Input/Output Expectations
- Input: VCF file (bgzipped: `.vcf.gz`) and `regions.bed`
- Output: One compressed VCF file per region, e.g.:
  ```
  output_chunks/
  ├── BRCA1.vcf.gz
  ├── BRCA2.vcf.gz
  └── TP53.vcf.gz
  ```

## External Commands Used
- `bcftools view -r {chr}:{start}-{end} -o {output} -Oz {input}`

## Cursor Agent Instructions
- Scaffold the project layout
- Implement functions as described
- Use `multiprocessing.Pool` to parallelize region extraction
- Add error handling and CLI parsing
- Ensure all files are written to the `output_chunks/` directory
- Validate bcftools is installed before execution