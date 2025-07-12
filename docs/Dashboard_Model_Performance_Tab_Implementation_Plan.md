# Dashboard Model Performance Tab Implementation Plan

## Overview
This document details the implementation plan for a new "Model Performance" tab in the Dashboard. The goal is to display meaningful model and analysis metrics for each user job **without using any ground truth labels**. The tab will also provide static reference metrics from model validation (AUC, F1, etc.), but will not attempt to compute sensitivity, specificity, or F1-score for user jobs.

---

## Rationale
- **Why not show sensitivity, specificity, F1-score for user jobs?**
  - These metrics require ground truth (true/false positive/negative) labels, which are not available for user-submitted genomic data.
  - Displaying such metrics for user jobs would be misleading and scientifically invalid.
  - Instead, we focus on metrics that are meaningful for single-sample, label-free analysis: model confidence, risk score distribution, variant metrics, PRS/pathway summary, and risk category breakdown.
- **Model validation metrics (AUC, F1, etc.)**
  - These are computed during model training/validation on public datasets with ground truth.
  - We display them as a static reference, clearly labeled as such, to inform users about the model's general performance.

---

## What Will Be Displayed

### 1. **Your Analysis Metrics** (per user job)
- **Model Confidence:**
  - Mean, min, max confidence (from `metrics.confidence_metrics`)
- **Risk Score Distribution:**
  - Mean, std, max, coefficient of variation (from `metrics.confidence_metrics`)
- **Variant Metrics:**
  - Number of pathogenic, benign, uncertain variants; mean/max CADD score (from `metrics.variant_metrics`)
- **PRS and Pathway Burden:**
  - PRS confidence, high-risk cancers, pathway burden score, high-burden pathways (from `metrics.integration_metrics`)
- **Risk Category Distribution:**
  - Number of predictions in each risk category (low/moderate/high/very high)

### 2. **Model Validation (Reference Only)**
- **AUC-ROC, F1-score, MCC, etc.** (from latest model training/validation)
- **Short explanation:** “These metrics are from model validation on public datasets and do not reflect your individual results.”

---

## Data Sources
- **User Job Metrics:**
  - API: `/api/results/{job_id}`
  - Fields: `metrics.confidence_metrics`, `metrics.variant_metrics`, `metrics.integration_metrics`, `metrics.performance_indicators`
- **Model Validation Metrics:**
  - Static JSON: `geneknow_pipeline/real_data_training_results_FIXED.json` (or hardcoded in frontend)

---

## Relevant Variables
- `metrics.confidence_metrics` (mean_model_confidence, min_model_confidence, max_model_confidence, risk_score_mean, risk_score_std, risk_score_cv, max_risk_score, high_risk_count, ml_fusion_confidence, ml_fusion_risk_category)
- `metrics.variant_metrics` (total_variants, pathogenic_variants, benign_variants, uncertain_variants, pathogenic_ratio, high_cadd_variants, mean_cadd_score, max_cadd_score)
- `metrics.integration_metrics` (prs_confidence, prs_high_risk_cancers, pathway_burden_score, high_burden_pathways)
- `metrics.performance_indicators` (mean_confidence, mean_risk_score, high_confidence_ratio, variant_quality_score, total_predictions)
- **Static validation metrics:** AUC, F1, MCC, etc. (from `real_data_training_results_FIXED.json`)

---

## Relevant Files
- **Frontend:**
  - `desktop/ui/src/pages/DashboardPage.tsx` (add new tab)
  - `desktop/ui/src/components/ModelPerformanceTab.tsx` (new component)
- **Backend:**
  - No changes needed (metrics already available in API)
- **Static Data:**
  - `geneknow_pipeline/real_data_training_results_FIXED.json` (for reference metrics)

---

## Implementation Checklist

### **Backend**
- [x] Ensure all relevant metrics are included in `/api/results/{job_id}` (already implemented)
- [x] Ensure static model validation metrics are available (already in JSON file)

### **Frontend**
- [ ] Add a new tab to `DashboardPage.tsx` labeled “Model Performance”
- [ ] Create `ModelPerformanceTab.tsx` component
    - [ ] Fetch user job metrics from `/api/results/{job_id}`
    - [ ] Display model confidence, risk score stats, variant metrics, PRS/pathway summary, risk category breakdown
    - [ ] Fetch and display static model validation metrics (AUC, F1, etc.)
    - [ ] Add tooltips/info icons for each metric
    - [ ] Add clear disclaimer for static validation metrics
    - [ ] Add comprehensive logging for all data fetches and render steps
- [ ] Update navigation and styling as needed

---

## Example UI Layout

```
[Dashboard Tabs: Report | Model Performance]

---

Model Performance

Your Analysis Metrics
---------------------
- Model Confidence: 0.82 (mean), 0.75 (min), 0.91 (max)
- Risk Score: mean 0.41, std 0.12, max 0.67, CV 0.29
- Pathogenic Variants: 3, Benign: 12, Uncertain: 5
- Mean CADD Score: 23.1, Max CADD: 34.2
- PRS Confidence: high, Pathway Burden Score: 0.81
- Risk Category: 2 high, 4 moderate, 14 low

Model Validation (Reference Only)
---------------------------------
- AUC-ROC: 0.76
- F1-score: 0.63
- MCC: 0.42
- (etc.)
*These metrics are from model validation on public datasets and do not reflect your individual results.*
```

---

## Notes
- **No ground truth metrics (sensitivity, specificity, F1) will be shown for user jobs.**
- **All metrics will be accompanied by tooltips and logs for transparency.**
- **This plan is designed for production and extensibility.** 