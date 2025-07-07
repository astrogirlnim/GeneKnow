# Genomic Parallel Extraction Pipeline

A production-ready pipeline for extracting specific genomic regions from massive VCF files using parallel processing.

## Quick Start (Team Setup)

**1. First time setup (downloads data once per machine):**
```bash
python3 setup_data.py
```

**2. Extract your regions:**
```bash
python3 extract_by_region_precise.py
```

**That's it!** Large files are shared across team members.

## What This Pipeline Does

- ✅ **Parallel processing** across multiple CPU cores
- ✅ **100% precise extraction** (exact boundaries)
- ✅ **Multi-chromosome support** (handles mixed datasets)
- ✅ **Massive file handling** (tested on 11GB+ files)
- ✅ **Team-friendly** (shared data, no repeated downloads)

## Data Management

**Smart caching system:**
- Downloads large files **once per machine**
- Stores in shared locations (`/shared/genomics-data` → `~/genomics-data` → `./data-cache`)
- Creates symlinks in `test-data/` for project use
- Automatically handles indexing

**Git workflow:**
- ✅ Track: Scripts, regions, documentation
- ❌ Ignore: Large VCF files, outputs, indexes

## Supported File Types

| Dataset | Size | Type | Use Case |
|---------|------|------|----------|
| **Phase 1** | 11GB+ | SNPs + Indels + SVs | Comprehensive variant analysis |
| **Phase 3** | 200MB-1.2GB | High-quality SNPs | Population genetics |

## Usage Examples

**Define your regions in `regions.bed`:**
```
22	29121014	29235591	EWSR1
1	35691274	35801992	TP73
```

**Run extraction:**
```bash
python3 extract_by_region_precise.py
# Output: ✅ Extracted EWSR1 from chr22: 3260 variants (100% precise)
```

**Analyze results:**
```bash
ls -lh output_chunks/
# EWSR1.vcf.gz    598K
# TP73.vcf.gz     4.6M
```

## Team Collaboration

**New team member workflow:**
1. `git clone` the repository
2. `python3 setup_data.py` (uses existing shared data if available)
3. `python3 extract_by_region_precise.py` (ready to go!)

**No need to:**
- Download 11GB+ files individually
- Manage complex data paths
- Worry about incomplete downloads (resume capability)

## Performance

**Tested on:**
- ✅ 11GB Phase 1 files (chr1)
- ✅ 3,260 variants extracted in seconds
- ✅ 100% boundary precision
- ✅ Parallel processing across chromosomes

## Requirements

```bash
# System tools
bcftools
bgzip
wget

# Python (standard library only)
python3
```

## Troubleshooting

**"No such file" errors:** Run `python3 setup_data.py` first

**Slow downloads:** Script uses `wget -c` for resume capability

**Permission errors:** Check write access to data directories

## Architecture

```
Large Files (11GB+)     Shared Storage        Project Directory
├── chr1.vcf.gz    →    ~/genomics-data   →   test-data/ (symlinks)
├── chr22.vcf.gz        (cached once)         ├── Scripts
└── indexes                                   └── Output chunks
``` 