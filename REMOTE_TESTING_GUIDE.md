# Remote Genomic Data Testing Guide

## Overview

You can test your genomic extraction pipeline with datasets of 50GB+ without downloading them! This uses HTTP range requests to fetch only the specific genomic regions you need.

## Available Remote Datasets

### 1. 1000 Genomes Project (EBI)
- **Total Size**: ~15GB compressed per release
- **Access**: Free, no authentication required
- **URL Pattern**: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{N}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`

### 2. AWS S3 1000 Genomes
- **Total Size**: 200TB+ across all versions
- **Access**: Free via HTTP
- **Example**: `https://1000genomes-dragen.s3.amazonaws.com/data/CHR/CHR{N}/`

### 3. gnomAD (Genome Aggregation Database)
- **Total Size**: Multiple TB
- **Access**: Free via Google Cloud
- **Note**: Requires special setup for streaming

## Quick Start

### 1. Test Remote Extraction
```bash
python test_remote_extraction.py
```

This will:
- Extract 6 cancer-related genes from remote VCF files
- Show performance metrics
- Demonstrate HTTP range request efficiency

### 2. Use with Main Pipeline
To use remote data with your main pipeline:

```python
# In extract_by_region.py, uncomment the remote VCF_FILES section
# Or set VCF_FILES directly to remote URLs
```

### 3. Test with Large Gene Set
```bash
# Use the larger BED file with 33 cancer genes
cp regions_large.bed regions.bed
python extract_by_region.py
```

## Performance Characteristics

### Remote vs Local Speed
- **Local SSD**: ~100-500MB/s read speed
- **Remote HTTP**: ~10-50MB/s (depends on connection)
- **But**: You only transfer the regions you need!

### Example Savings
For extracting BRCA1 (81kb region) from chromosome 17:
- **Full chromosome download**: 397MB
- **Region-only transfer**: ~2MB
- **Savings**: 99.5% less data transfer!

## Advanced Usage

### Parallel Remote Extraction
The pipeline automatically parallelizes remote extractions:
```python
CONFIG = {
    "max_processes": 8,  # Use 8 parallel connections
}
```

### Custom Remote Sources
Add your own remote VCF sources:
```python
VCF_FILES = {
    "1": "https://your-server.com/chr1.vcf.gz",
    "2": "s3://your-bucket/chr2.vcf.gz",  # S3 URLs work too!
}
```

### Performance Tuning
1. **Increase parallelism** for better throughput
2. **Use indexed VCFs** (.tbi files) for faster seeks
3. **Cache frequently accessed regions** locally

## Troubleshooting

### "Connection refused" errors
- Check your internet connection
- Some servers limit concurrent connections
- Reduce `max_processes` if needed

### Slow performance
- Remote servers may throttle connections
- Try different mirror sites
- Consider time of day (servers can be busy)

### bcftools errors
Ensure bcftools is compiled with remote file support:
```bash
bcftools --version  # Should show htslib version
```

## Example Output
```
🌐 Testing Remote Genomic Data Extraction
📊 Total remote data available: ~15GB across chromosomes
🎯 Extracting 6 cancer-related gene regions
============================================================

🧬 Extracting BRCA1 (chr17:41,196,312-41,277,500)
  ✅ Success: 3,421 variants in 4.23s (823.1KB)

🧬 Extracting TP53 (chr17:7,565,097-7,590,868)
  ✅ Success: 1,234 variants in 2.11s (234.5KB)

============================================================
📊 REMOTE EXTRACTION SUMMARY:
============================================================
⏱️  Total time: 18.34 seconds
✅ Successful extractions: 6/6
🧬 Total variants extracted: 12,345
💾 Total output size: 2,134.5 KB (2.1 MB)
⚡ Average time per region: 3.06 seconds
🌐 Data transferred: ~4,269.0 KB (estimated)

💡 Note: Only the required genomic regions were transferred,
   not the full 3,893 MB of chromosome data!
``` 