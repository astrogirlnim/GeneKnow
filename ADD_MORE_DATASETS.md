# Adding More Genomic Datasets

## Current Dataset: 1000 Genomes (13.8 GB)
- 2,504 individuals from 26 populations
- ~84 million variants
- Reference: GRCh37/hg19

## 🚀 Quick Start: Add a New Dataset

### 1. Edit `config_data_source.py`:

```python
# In the get_vcf_files() function, add:

# Example: Add gnomAD (more variants, 76k genomes)
elif DATA_SOURCE == 'gnomad':
    return {
        "17": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz",
        # Add more chromosomes...
    }
```

### 2. Set environment variable:
```bash
export GENOMIC_DATA_SOURCE=gnomad
python3 extract_by_region.py
```

## 📊 Available Public Datasets

### 1. gnomAD v3.1 (Huge - 744GB total)
- **Size**: ~30GB per chromosome
- **Samples**: 76,156 genomes
- **Variants**: 600+ million
- **Best for**: Population frequency, rare variants

```python
# Example chr17 (4.7GB compressed):
"17": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz"
```

### 2. ClinVar (Small - 300MB)
- **Size**: Single file, all chromosomes
- **Focus**: Disease-associated variants
- **Best for**: Clinical relevance

```python
"all": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
```

### 3. dbSNP (Medium - 10GB)
- **Size**: Single file
- **Variants**: All known SNPs
- **Best for**: Variant annotation

```python
"all": "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
```

### 4. UK Biobank (Restricted)
- **Size**: 500TB+
- **Samples**: 500,000 individuals
- **Note**: Requires application/approval

## ⚠️ Important Considerations

### Reference Genome Version
Make sure your BED file coordinates match the dataset's reference:
- **GRCh37/hg19**: 1000 Genomes, ExAC
- **GRCh38/hg38**: gnomAD v3, ClinVar, dbSNP

### Convert coordinates if needed:
```bash
# Use UCSC liftOver tool
# https://genome.ucsc.edu/cgi-bin/hgLiftOver
```

### File Formats
- `.vcf.gz` - Standard compressed VCF
- `.vcf.bgz` - Block-gzipped (same as .vcf.gz)
- `.bcf` - Binary VCF (also supported)

## 🔧 Performance Tips

### For Large Datasets (gnomAD):
1. **Test with one chromosome first**
   ```python
   VCF_FILES = {"17": "https://...chr17.vcf.bgz"}
   ```

2. **Increase timeout for slow servers**
   ```python
   # In extract_by_region.py
   subprocess.run(cmd, timeout=300)  # 5 minutes
   ```

3. **Use parallel downloads**
   ```python
   CONFIG["max_processes"] = 4  # Limit parallelism
   ```

## 📝 Example: Add gnomAD + ClinVar

```python
# In config_data_source.py
elif DATA_SOURCE == 'multi':
    return {
        # gnomAD for main analysis
        "1": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz",
        "17": "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz",
        
        # ClinVar for clinical annotations
        "clinvar": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    }
```

## 🌐 Finding More Datasets

1. **NCBI**: https://www.ncbi.nlm.nih.gov/variation/tools/1000genomes/
2. **Ensembl**: https://www.ensembl.org/info/data/ftp/index.html
3. **UCSC**: https://hgdownload.soe.ucsc.edu/downloads.html
4. **AWS Open Data**: https://registry.opendata.aws/

## 💡 Quick Test

After adding a new dataset:
```bash
# Check if it works
python3 check_remote_sizes.py

# Test extraction
head -1 regions.bed > test_region.bed
python3 extract_by_region.py
``` 