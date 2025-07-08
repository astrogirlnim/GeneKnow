# TCGA Testing Data Generator

This directory contains scripts for downloading TCGA (The Cancer Genome Atlas) data from the GDC (Genomic Data Commons) API for testing purposes with the GenePredict project.

## Overview

These scripts allow you to download real genomic data files from TCGA, specifically focusing on blood and bone marrow cancer samples for local testing and development.

## 🚀 Quick Start

1. **Set up Python environment:**
   ```bash
   cd data_processing/testing_data_generator
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install requests
   ```

2. **Download blood cancer data:**
   ```bash
   python simple_blood_cancer_client.py
   ```

## 📁 Scripts Overview

### 1. `simple_blood_cancer_client.py` - **Main Script**
**Purpose:** Download TCGA files specifically for blood and bone marrow cancers

**What it does:**
- Searches GDC API for blood cancer samples using metadata filters
- Identifies files from primary sites: Blood, Bone Marrow, Lymph Nodes, Spleen
- Downloads processed data files (TSV, MAF, TXT formats)
- Uses simple metadata-based classification (no AI/LLM required)

**Usage:**
```bash
python simple_blood_cancer_client.py
```

**Output:**
- Lists found blood cancer files with details
- Downloads first file as test to `blood_cancer_downloads/` directory
- Shows file metadata including cancer types found

**Example output:**
```
🩸 Blood cancer file: sample.rna_seq.gene_counts.tsv
   Cases: Lymph nodes - Mature T- and NK-Cell Lymphomas
```

### 2. `minimal_test.py` - **Basic Download Test**
**Purpose:** Simple test to verify GDC API connectivity and download functionality

**What it does:**
- Downloads any small files from GDC API
- Tests basic download mechanism
- Verifies file integrity

**Usage:**
```bash
python minimal_test.py
```

**Output:**
- Downloads 1 small file to `minimal_test_downloads/`
- Shows download progress and file information

### 3. `check_formats.py` - **Format Discovery**
**Purpose:** Explore what data formats are available in GDC API

**What it does:**
- Queries GDC API for available file formats
- Shows counts of different data formats (TSV, MAF, TXT, etc.)
- Lists experimental strategies (RNA-Seq, WXS, WGS, etc.)

**Usage:**
```bash
python check_formats.py
```

**Output:**
```
Data formats: {'TSV': 23, 'MAF': 8, 'TXT': 12}
Experimental strategies: {'RNA-Seq': 8, 'WXS': 8, 'WGS': 27}
```

### 4. `debug_api.py` - **API Response Inspector**
**Purpose:** Debug and understand GDC API response structure

**What it does:**
- Makes simple API calls to GDC
- Shows raw response structure
- Helps understand data organization

**Usage:**
```bash
python debug_api.py
```

**Output:**
- Shows API response keys and structure
- Displays sample file metadata

## 📊 Data Types Available

The scripts download **processed** genomic data files (not raw FASTQ):

### File Formats:
- **TSV** - Tab-separated values (gene expression counts, copy number, etc.)
- **MAF** - Mutation Annotation Format (variant calls)
- **TXT** - Text files (segmentation data, etc.)

### Experimental Strategies:
- **RNA-Seq** - Gene expression data
- **WXS** - Whole Exome Sequencing (variants)
- **WGS** - Whole Genome Sequencing (variants, copy number)

### Cancer Types Found:
- Acute Myeloid Leukemia
- Acute Lymphoblastic Leukemia
- Chronic Lymphocytic Leukemia
- Hodgkin/Non-Hodgkin Lymphoma
- Multiple Myeloma
- T-Cell/NK-Cell Lymphomas
- Plasma Cell Tumors

## 🔧 Technical Details

### Dependencies:
- **Python 3.7+**
- **requests** library for HTTP calls

### API Endpoint:
- Base URL: `https://api.gdc.cancer.gov`
- Uses public/open access data only
- No authentication required for open data

### File Classification:
Scripts use metadata-based filtering:
- `cases.primary_site` - Anatomical location
- `cases.disease_type` - Disease classification  
- `cases.diagnoses.primary_diagnosis` - Specific diagnosis

### Download Verification:
- File size validation
- Optional MD5 checksum verification (when available)

## 📁 Directory Structure

```
testing_data_generator/
├── README.md                           # This file
├── simple_blood_cancer_client.py       # Main blood cancer downloader
├── minimal_test.py                     # Basic download test
├── check_formats.py                    # Format discovery
├── debug_api.py                        # API debugging
├── venv/                               # Python virtual environment (gitignored)
├── *_downloads/                        # Download directories (gitignored)
└── *.log                               # Log files (gitignored)
```

## 🚫 What's Not Included

- **Raw FASTQ files** - Not available in open access
- **Controlled access data** - Requires dbGaP approval
- **Patient identifiers** - All data is de-identified
- **LLM/AI classification** - Uses simple metadata filtering

## 🎯 Use Cases

### For GenePredict Development:
1. **Test file parsing** - Get real genomic file formats
2. **Validate pipelines** - Use actual TCGA data structure
3. **Demo data** - Blood cancer samples for demonstrations
4. **Integration testing** - Real-world data for testing

### For Research:
1. **Format exploration** - Understand TCGA data organization
2. **Sample discovery** - Find relevant cancer types
3. **Metadata analysis** - Explore available annotations

## 🔒 Privacy & Compliance

- Uses only **open access** data from GDC
- All data is **de-identified** by TCGA
- No patient-specific information downloaded
- Complies with GDC data usage policies

## 🐛 Troubleshooting

### Common Issues:

**"ModuleNotFoundError: No module named 'requests'"**
```bash
pip install requests
```

**"No files found"**
- Check internet connection
- GDC API may be temporarily unavailable
- Try `debug_api.py` to test basic connectivity

**"Download failed"**
- Large files may timeout
- Check available disk space
- Try downloading smaller files first

### Getting Help:
1. Run `debug_api.py` to test API connectivity
2. Check GDC status: https://portal.gdc.cancer.gov/
3. Review GDC API documentation: https://docs.gdc.cancer.gov/API/

## 📝 Notes

- Download speeds depend on file size and network connection
- Large files (>100MB) may take several minutes
- Downloaded files are suitable for local testing and development
- All scripts include extensive logging and progress reporting 