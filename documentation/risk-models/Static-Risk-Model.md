# 📘 GeneKnow: Static Annotation Models Overview

This document explains the static models used in the GeneKnow cancer risk analysis system. These models are not trained or learned during runtime — they are external sources of precomputed or curated information that provide powerful, interpretable signals about genomic variants.

---

## 🧬 1. Polygenic Risk Score (PRS)

### 🔍 What it is:

A weighted sum of small-effect SNPs associated with increased disease risk.

### 📦 How it works:

* Each SNP has a known effect size (β) from GWAS
* Patient genotype (0, 1, 2 risk alleles) is multiplied by β
* Sum of all contributions = PRS

### 🧠 Why we use it:

* Captures inherited background risk
* Quantitative and fast to compute
* Based on large human cohorts

### ⚠️ Limitations:

* Doesn’t detect rare, high-penetrance mutations
* Can be population-biased if not properly calibrated

---

## 🟥 2. ClinVar

### 🔍 What it is:

A curated public database of known variants annotated as "Pathogenic", "Likely pathogenic", "Benign", etc. In GeneKnow, this is implemented as a local SQLite database with cancer-focused annotations.

### 📦 How it works:

* Each variant is matched to ClinVar entries via direct database lookup
* Local SQLite database contains clinical significance, review status, and condition information
* Variants are assigned numerical risk scores based on clinical significance:
  * Pathogenic: 1.0
  * Likely pathogenic: 0.8
  * Benign/Likely benign: 0.0
  * Uncertain significance: 0.1
  * Drug response: 0.1, Risk factor: 0.3
* Cancer-related conditions receive 20% scoring bonus
* Review status modifiers (practice guidelines get 10% bonus)

### 🧠 Why we use it:

* Direct clinical interpretation with expert-curated pathogenicity
* Provides numerical risk scores for quantitative analysis
* Cancer-specific focus with enhanced scoring for cancer-related conditions
* Actionable clinical recommendations based on clinical significance
* Confidence scores based on review status and submission quality

### ⚠️ Limitations:

* Limited to variants in our curated database (focused on cancer genes)
* Database requires periodic updates from ClinVar releases
* Some rare variants may not be represented

---

## 🟨 3. CADD (Combined Annotation Dependent Depletion)

### 🔍 What it is:

A machine learning-based score estimating how damaging a variant is to biological function.

### 📦 How it works:

* Each SNV and many indels are pre-scored by CADD
* Scores are PHRED-scaled (e.g. CADD 20 = top 1% most likely deleterious)
* Annotated by lookup (e.g. via vcfanno or VEP)

### 🧠 Why we use it:

* Covers nearly all possible SNVs (even novel ones)
* Scores regulatory, coding, and structural impact

### ⚠️ Limitations:

* Not cancer-specific
* Not interpretable as direct probability of disease
* Doesn’t distinguish gain vs. loss of function

---

## 🧬 4. TCGA Frequency Match (Tumor-Enrichment)

### 🔍 What it is:

Comparison of patient variants to frequencies seen in cancer tumors (e.g. TCGA BRCA dataset)

### 📦 How it works:

* For each patient variant, check if it occurs in tumor cohorts
* Calculate enrichment (e.g. seen in 18% of BRCA tumors vs 0.2% of controls)

### 🧠 Why we use it:

* Adds real-world tumor relevance to rare variants
* Helps confirm that suspicious variants occur in actual cancer patients

### ⚠️ Limitations:

* May reflect passenger mutations
* Requires good filtering to avoid false enrichment

---

## 🟪 5. Gene/Pathway Burden Model

### 🔍 What it is:

Counts damaging variants in key genes or pathways (e.g. DNA repair, TP53, mismatch repair)

### 📦 How it works:

* Define gene sets by cancer relevance
* Count # of rare high-CADD/ClinVar variants in each gene
* Normalize by gene length or patient genome coverage

### 🧠 Why we use it:

* Captures multi-hit polygenic risk not visible through PRS
* Great for exploring biological themes (e.g. DNA repair pathway failure)

### ⚠️ Limitations:

* Requires well-defined gene sets
* May be noisy without good variant filtering

---

These static models form the “feature generators” for our downstream machine learning layer. Each one contributes a different kind of biological signal that, when combined, gives us a more complete picture of patient-specific cancer risk.
