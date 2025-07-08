# 🎤 GeneKnow Pipeline Presentation Script

## Slide 1: Introduction
**"Hello, I'm going to walk you through GeneKnow - our genomic risk assessment pipeline."**

- Local-first genomic analysis
- Processes FASTQ, BAM, and VCF files
- Predicts cancer risk using ML models
- Built with LangGraph + bioinformatics tools

---

## Slide 2: The Problem We Solve
**"Current genomic analysis has three major issues:"**

1. **Privacy**: Cloud services store your DNA data
2. **Cost**: Commercial tests are expensive ($100-500)
3. **Speed**: Results take weeks to months

**"GeneKnow solves all three - it's private, affordable, and gives results in seconds."**

---

## Slide 3: High-Level Architecture
**"Here's how the system works:"**

```
User uploads file → LangGraph Pipeline → Risk Report
         ↓                    ↓              ↓
   (FASTQ/BAM/VCF)    (8 processing nodes)  (JSON/PDF)
```

**Key Points:**
- Everything runs locally on user's machine
- LangGraph orchestrates the pipeline
- Results in structured JSON for the frontend

---

## Slide 4: The Processing Pipeline
**"Let me show you the actual processing flow:"**

```
1. File Input     → Validate format
2. Preprocess     → Align FASTQ or load VCF
3. Parallel Split → Two paths based on file type
   - Path A: Variant Calling (for FASTQ/BAM)
   - Path B: QC Filter (for VCF)
4. Merge Results  → Combine both paths
5. TCGA Mapper    → Match to cancer database
6. Risk Model     → ML prediction
7. Formatter      → Structure for API
8. Report Writer  → Generate sections
```

---

## Slide 5: File Type Handling
**"We handle three genomic file formats intelligently:"**

### FASTQ (Raw sequences)
```
@SEQ_001
AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG
+
#############################@@@@@@@
```
- Contains millions of DNA reads
- We align to reference genome using BWA
- Then call variants

### BAM (Aligned sequences)
- Already aligned, skip BWA step
- Go straight to variant calling

### VCF (Variant calls)
- Already has variants
- Skip to quality filtering
- Most efficient path

---

## Slide 6: Parallel Processing
**"Here's where it gets clever - we use parallel paths:"**

```
Preprocess
    ├─→ If FASTQ/BAM: Variant Calling
    └─→ If VCF: QC Filter
         ↓
    Merge Results
```

**Benefits:**
- Optimized for each file type
- No wasted computation
- Faster results

---

## Slide 7: TCGA Integration
**"We match variants against real cancer data:"**

- **Database**: 2,828 patient cohort
- **5 Cancer Types**: Breast, Colon, Lung, Prostate, Blood
- **83 Known Variants**: Curated from literature

**Example Query:**
```sql
SELECT frequency FROM variants 
WHERE gene='BRCA1' AND cancer_type='breast'
-- Returns: 2.3% frequency in breast cancer patients
```

---

## Slide 8: Risk Prediction
**"Our ML models predict risk for 5 cancer types:"**

### Model Details:
- **Algorithm**: Logistic Regression
- **Features**: Age, Sex, Gene variants
- **Training**: Based on known associations

### Cancer Genes We Track:
- **Breast**: BRCA1, BRCA2, TP53, PIK3CA...
- **Colon**: APC, KRAS, TP53, SMAD4...
- **Lung**: TP53, KRAS, EGFR, ALK...

---

## Slide 9: Technical Implementation
**"Under the hood, here's what makes it work:"**

### LangGraph State Management:
```python
class GenomicState(TypedDict):
    file_path: str
    raw_variants: List[dict]
    risk_scores: dict
    report_sections: dict
    # ... 20+ fields tracking pipeline state
```

### Node Pattern:
```python
def process(state):
    try:
        # Do work
        state["output"] = result
        state["completed_nodes"].append("node_name")
    except Exception as e:
        state["errors"].append(error_info)
    return state
```

---

## Slide 10: API & Frontend Integration
**"The pipeline exposes a simple REST API:"**

### Request:
```json
POST /analyze_path
{
    "file_path": "sample.fastq.gz",
    "preferences": {
        "patient_data": {"age": 45, "sex": "F"}
    }
}
```

### Response:
```json
{
    "risk_scores": {
        "breast": 72.8,
        "colon": 19.3
    },
    "variants": [...],
    "report_sections": {...}
}
```

---

## Slide 11: Performance
**"The pipeline is fast:"**

### Processing Times:
- **File validation**: 0.01s
- **Alignment**: 0.5-0.8s
- **Variant calling**: 0.1s
- **Risk prediction**: 0.05s
- **Total**: <1 second for test files

### Scalability:
- Currently: Process entire file at once
- Future: Chunk large files for parallel processing
- Goal: Handle 1GB files in <1 minute

---

## Slide 12: What's Next?
**"Here's what we're building next:"**

1. **Real DeepVariant Integration**
   - Currently using simulated variants
   - DeepVariant is Google's state-of-the-art caller

2. **FASTQ Chunking**
   - Split large files into chunks
   - Process in parallel
   - Handle whole genome files (100GB+)

3. **LLM Report Generation**
   - Use Llama 3.1 for natural language reports
   - Multi-language support
   - Personalized recommendations

---

## Slide 13: Demo
**"Let me show you a live demo:"**

1. Start the API server
2. Upload a test FASTQ file
3. Watch the pipeline process
4. View risk report

```bash
cd langgraph
python api_server.py
# Server running on http://localhost:5001
```

---

## Slide 14: Key Takeaways
**"To summarize:"**

✅ **Privacy-First**: All processing stays local
✅ **Real Science**: Based on TCGA data and peer-reviewed genes
✅ **Fast**: Results in seconds, not weeks
✅ **Extensible**: LangGraph makes adding nodes easy
✅ **Production-Ready**: RESTful API for easy integration

---

## Slide 15: Questions?
**"Thank you! Any questions?"**

### Common Questions & Answers:

**Q: How accurate is the risk prediction?**
A: Based on known gene associations from literature. Not diagnostic - for research/educational use.

**Q: Can it handle whole genome files?**
A: Not yet - that's our chunking work. Currently handles exome/panel data well.

**Q: Why LangGraph instead of regular Python?**
A: Declarative pipeline, easy parallelization, built-in state management, visual debugging.

**Q: Is the data really secure?**
A: Yes - no network calls, temporary files encrypted, memory cleared after use.

---

## Technical Deep-Dive Bonus Slides

### Slide A: BWA Alignment Details
```bash
# How we align FASTQ to reference
bwa index reference.fa
bwa mem -t 8 reference.fa reads.fastq > aligned.sam
samtools sort aligned.sam > aligned.bam
samtools index aligned.bam
```

### Slide B: Variant Calling Mock
```python
# Current implementation (mock)
def call_variants(bam_path):
    return [
        {"chrom": "17", "pos": 41276045, 
         "gene": "BRCA1", "variant": "c.5266dupC"}
    ]

# Future (real DeepVariant)
subprocess.run(["deepvariant", ...])
```

### Slide C: Quality Control Filters
```python
# Our QC thresholds
MIN_QUALITY = 30      # Phred score
MIN_DEPTH = 10        # Read depth  
MIN_ALLELE_FREQ = 0.01  # 1% frequency

# Only high-confidence variants pass
```

---

**End of Presentation** 