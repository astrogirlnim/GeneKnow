# GeneKnow LangGraph Pipeline

This directory contains the LangGraph implementation for genomic risk assessment pipeline.

## ⚠️ Current Status: MOCK IMPLEMENTATION

All nodes currently use **MOCK DATA** for testing the pipeline flow. Real implementations need to be added for:
- File validation (BioPython)
- Alignment (BWA-MEM2)
- Variant calling (DeepVariant)
- TCGA database queries
- TensorFlow risk models
- LLM report generation (Llama 3.1)

## 📁 Structure

```
langgraph/
├── graph.py              # Main pipeline definition
├── state.py              # State schema
├── nodes/                # Individual processing nodes
│   ├── file_input.py     # ⚠️ MOCK - Validates input files
│   ├── preprocess.py     # ⚠️ MOCK - FASTQ alignment
│   ├── variant_calling.py # ⚠️ MOCK - DeepVariant wrapper
│   ├── qc_filter.py      # ✅ Real logic (needs VCF data)
│   ├── tcga_mapper.py    # ⚠️ MOCK - TCGA matching
│   ├── risk_model.py     # ⚠️ MOCK - ML predictions
│   ├── formatter.py      # ✅ Real formatting logic
│   └── report_writer.py  # ⚠️ MOCK - LLM generation
└── test_pipeline.py      # Test script
```

## 🚀 Quick Start

```bash
# Install dependencies
pip install langgraph langchain tensorflow biopython

# Run test pipeline
cd langgraph
python test_pipeline.py
```

## 🔧 Replacing Mock Implementations

Each node has TODO comments marking what needs to be implemented:

### 1. File Input (`file_input.py`)
```python
# TODO: Real implementation should:
# 1. Check file exists and is readable
# 2. Validate file format (check headers)
# 3. Extract read count for FASTQ or alignment stats for BAM
```

### 2. Preprocessing (`preprocess.py`)
```python
# TODO: Real implementation should:
# 1. Check BWA-MEM2 is installed
# 2. Download/locate reference genome (hg38)
# 3. Run: bwa-mem2 mem -t 8 ref.fa reads.fq | samtools sort -o aligned.bam
```

### 3. Variant Calling (`variant_calling.py`)
```python
# TODO: Real implementation should:
# 1. Check DeepVariant is installed (Docker/Singularity)
# 2. Run DeepVariant with appropriate settings
# 3. Parse resulting VCF file
```

## 📊 Pipeline Flow

```
FASTQ/BAM → Validate → Align (if needed) → Call Variants → QC Filter → 
TCGA Match → Risk Model → Format JSON → Generate Report → PDF
```

## 🧪 Testing

The current mock implementation allows testing:
- Pipeline flow and state management
- Error handling and recovery
- JSON formatting
- Report structure

## 🔄 Integration with Tauri

The pipeline can be called from Tauri backend:

```python
from langgraph.graph import run_pipeline

# In your Tauri command handler
result = run_pipeline(file_path, user_preferences)
return result["structured_json"]  # Return to frontend
```

## 📝 Next Steps

1. **Replace mock file validation** with BioPython
2. **Integrate BWA-MEM2** via subprocess
3. **Set up DeepVariant** Docker container
4. **Create TCGA SQLite database** from MAF files
5. **Train TensorFlow models** on TCGA data
6. **Integrate Llama 3.1** for report generation 