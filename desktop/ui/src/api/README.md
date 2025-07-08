# Genomic Processing API

This directory contains the Python scripts for genomic data processing and the TypeScript API wrappers to call them from the Tauri desktop application.

## Available Functions

### 1. FASTQ to VCF Conversion
Converts raw sequencing reads (FASTQ) to variant calls (VCF).

```typescript
import { convertFastqToVcf } from './api/genomicProcessing';

const result = await convertFastqToVcf({
  reference: 'path/to/reference.fa',
  fastq1: 'path/to/reads_R1.fastq.gz',
  fastq2: 'path/to/reads_R2.fastq.gz', // optional for paired-end
  outputPrefix: 'output/sample',
  threads: 4,
  aligner: 'bwa' // or 'minimap2', 'bowtie2'
});
```

### 2. Extract Genomic Regions
Extracts specific genomic regions from VCF files based on a BED file.

```typescript
import { extractGenomicRegions } from './api/genomicProcessing';

const result = await extractGenomicRegions({
  bedFile: 'path/to/regions.bed',
  outputDir: 'output/extracted_regions',
  maxProcesses: 4 // optional, for parallel processing
});
```

### 3. Check Dependencies
Verifies that required bioinformatics tools are installed.

```typescript
import { checkDependencies } from './api/genomicProcessing';

const deps = await checkDependencies();
// Returns: { samtools: true, bcftools: true, bwa: true, ... }
```

### 4. Generate Test FASTQ Files
Creates synthetic FASTQ files for testing.

```typescript
import { generateTestFastqFiles } from './api/genomicProcessing';

const files = await generateTestFastqFiles('/tmp/test', 1000);
// Returns: { success: true, file1: '...', file2: '...' }
```

## Python Scripts

- **fastq_to_vcf_pipeline.py**: Complete pipeline for variant calling
- **extract_by_region.py**: Parallel extraction of genomic regions
- **config_data_source.py**: Configuration for VCF data sources

## Prerequisites

The following tools must be installed on the system:
- Python 3.x
- samtools
- bcftools
- At least one aligner: bwa, minimap2, or bowtie2

Install on macOS:
```bash
brew install samtools bcftools bwa
```

## Architecture

1. **Frontend (TypeScript)**: Provides type-safe API functions
2. **Tauri Backend (Rust)**: Handles IPC and executes Python scripts
3. **Python Scripts**: Perform the actual genomic processing

The communication flow:
```
React Component → TypeScript API → Tauri IPC → Rust Handler → Python Script → Results
``` 