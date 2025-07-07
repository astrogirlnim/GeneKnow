# External SSD Setup for Genomic Demos

## Overview

For demos and production work, using local data on an external SSD provides:
- **10x faster extraction** (0.6s vs 6s per gene)
- **Reliable performance** (no network dependencies)
- **Full dataset access** (all chromosomes)

## Setup Instructions

### 1. Prepare Your External SSD

Edit `setup_external_datasets.sh` and update the path:
```bash
EXTERNAL_SSD="/Volumes/YourSSDName"  # Change this to your SSD name
```

### 2. Download Priority Datasets

Run the setup script to download cancer-relevant chromosomes:
```bash
./setup_external_datasets.sh
```

This downloads ~8GB of priority chromosomes (1,2,3,7,9,10,11,12,13,17,19,20,22,X).

### 3. Configure Your Pipeline

Update `config_data_source.py` with your SSD path:
```python
EXTERNAL_SSD_PATH = "/Volumes/YourSSDName/genomic-datasets"
```

## Usage

### Switch Between Data Sources

```bash
# Use external SSD (fast, for demos)
export GENOMIC_DATA_SOURCE=external
python3 extract_by_region.py

# Use remote data (slow, no storage needed)
export GENOMIC_DATA_SOURCE=remote
python3 extract_by_region.py

# Use test data (limited chromosomes)
export GENOMIC_DATA_SOURCE=local
python3 extract_by_region.py
```

### Check Current Configuration

```bash
python3 config_data_source.py
```

## Storage Requirements

| Dataset | Size | Description |
|---------|------|-------------|
| Priority chromosomes | ~8GB | Cancer-relevant chromosomes |
| Full genome | ~15GB | All chromosomes 1-22 + X |
| With indexes | +5MB | Small .tbi files |

## Performance Comparison

| Operation | External SSD | Remote |
|-----------|--------------|--------|
| Extract BRCA1 | ~0.5s | ~5s |
| 33 cancer genes | ~20s | ~3-5 min |
| Full genome scan | ~2 min | ~20 min |

## Tips for Demos

1. **Pre-download datasets** the night before
2. **Test extraction** to ensure SSD is mounted
3. **Keep remote as backup** in case of SSD issues
4. **Cache results** in `output_chunks/` for instant re-display

## Download Additional Datasets

### Full 1000 Genomes Phase 3
```bash
# Download all chromosomes (~15GB)
for i in {1..22} X; do
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
done
```

### gnomAD (Larger, more variants)
```bash
# Example: gnomAD v3.1.2 chr17 (~5GB)
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz
```

## Troubleshooting

### "File not found" errors
- Check SSD is mounted: `ls /Volumes/`
- Verify path in config: `python3 config_data_source.py`

### Slow performance on SSD
- Check SSD connection (USB 3.0+ recommended)
- Ensure .tbi index files are present
- Consider SSD health/fragmentation

### Switching sources not working
- Check environment variable: `echo $GENOMIC_DATA_SOURCE`
- Restart Python/terminal session
- Verify paths in `config_data_source.py` 