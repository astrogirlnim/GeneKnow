# 🔗 GeneKnow Risk Fusion Architecture

This document outlines how GeneKnow ties together static variant models, machine learning, and LangGraph-based orchestration to produce cancer risk predictions.

---

## 🧬 Overview

GeneKnow’s architecture combines static annotated variant signals (like PRS, ClinVar, and CADD) with a trained TensorFlow model to output a patient-specific risk score. All of this is orchestrated using a modular LangGraph pipeline.

---

## 🧩 Components

### 1. Variant Annotation

* Tools: VEP, bcftools, vcfanno, custom Rust wrappers
* Adds:

  * ClinVar significance tags
  * CADD scores
  * Population frequency (gnomAD)
  * Tumor overlap (TCGA)

### 2. Feature Extraction

* Converts VCF → JSON or CSV row per patient
* Features include:

  * PRS percentile
  * Count of pathogenic ClinVar hits
  * Max CADD score
  * Gene-level burden counts
  * Number of TCGA-enriched variants

### 3. Fusion Layer (ML)

* TensorFlow model predicts cancer risk probability from feature vector
* Can override or adjust score based on:

  * Known pathogenic mutations (ClinVar)
  * Rule-based thresholds (e.g. PRS > 99th percentile)

### 4. Narrative Report

* LLM (e.g. LLaMA 3.1 or similar) generates:

  * Explanation of which features influenced the score
  * Highlight of known mutations
  * Confidence/uncertainty flags

### 5. LangGraph Orchestration

* Each of the above stages is modeled as a LangGraph node:

  * Annotation → Feature Extraction → Fusion → Report Writing
  * Intermediate state (e.g. variant flags) passed via graph memory

---

## ⚖️ Conflict Resolution Strategy

To prevent “model conflict” or double-counting:

* Partition variants:

  * PRS = common SNPs only
  * ClinVar = rare, pathogenic knowns
  * CADD = only non-ClnVar novel variants
* Weighted feature fusion in ML model
* Rule-based overrides for ClinVar flags
* Transparency in narrative output

---

## 🔐 Privacy and Compliance

* Entire pipeline runs locally
* No internet connection required
* All precomputed models (CADD, ClinVar) bundled with build or updated from secure endpoints

---

## 📦 Summary Flow

```
Patient VCF
  ↓
[Annotation Node] → VCF + annotations
  ↓
[Feature Extraction Node] → Feature vector JSON/CSV
  ↓
[Fusion Node] → Cancer risk %
  ↓
[Report Node] → Narrative + chart-ready data
```

This system architecture ensures that GeneKnow provides interpretable, robust, multi-signal cancer risk predictions suitable for both researchers and clinicians.
