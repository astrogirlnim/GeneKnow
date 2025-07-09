# GeneKnow ML Risk Fusion Model

This document describes how we use machine learning — specifically TensorFlow — to combine static model outputs into a unified cancer risk prediction system.

⸻

## Purpose

To learn how to weight and combine various per-patient risk signals (e.g. PRS, ClinVar hits, CADD burden) into a continuous risk probability (0.0–1.0).

⸻

## Input Format

We convert annotated VCFs into structured, tabular feature vectors. Each row is one patient.

Example Feature Vector:
```
{
  "prs_percentile": 92.4,
  "has_pathogenic_brca1": 1,
  "num_clinvar_pathogenic": 2,
  "max_cadd_score": 33.7,
  "num_high_cadd_variants_in_tp53": 1,
  "num_tcga_matched_variants": 3,
  "cancer_label": 1
}
```

	•	Features are derived from static annotations.
	•	Labels (optional) are used for supervised training.
	•	Can be stored as CSV, JSONL, or TensorFlow Records

⸻

## Model Architecture (MVP)

A simple, transparent feedforward model:

Input (n features) → Dense(32, relu) → Dense(16, relu) → Dense(1, sigmoid)

	•	Loss: Binary Crossentropy
	•	Optimizer: Adam
	•	Metrics: Accuracy, AUC

⸻

## Training Procedure
	•	Normalize input features (min-max or z-score)
	•	Stratify train/test split (preserve class ratio)
	•	Early stopping to prevent overfitting

## Optional Enhancements:
	•	Dropout
	•	Feature importance extraction
	•	Attention over gene/pathway vector inputs
	•	Multi-label classifier (for cancer subtype prediction)

⸻

## Why ML?
	•	Allows us to learn relative importance of PRS vs ClinVar vs TCGA vs CADD
	•	Handles non-linear interactions
	•	Gives smooth probability output (rather than hard thresholds)

⸻

## Considerations
	•	Requires enough labeled samples (real cancer status)
	•	Regularization may be needed to prevent overfitting on small feature sets
	•	Feature engineering and quality is more important than model complexity

⸻

This ML layer is the heart of GeneKnow’s risk fusion engine. It turns curated biological features into actionable risk scores.