# Cancer Risk Prediction Metrics for GenePredict

## 1. Overview
GenePredict aims to predict cancer risk by comparing patient genomic data (FASTQ/BAM files) to TCGA profiles. The "Disease Risk Model" stage uses TensorFlow to generate probability or risk scores for specific cancer types. This report outlines mathematical metrics for prediction, supported by scientific literature, ensuring interpretability, statistical rigor, and suitability for low-resource settings.

## 2. Recommended Mathematical Metrics

### a. Probability-Based Metrics (Classification)
- **Logistic Regression Probability (Softmax for Multi-Class)**  
  - **Description**: Estimates probability of cancer types using logistic regression (binary) or softmax (multi-class).  
  - **Formula**:  
    Binary:  
    \[
    P(y=1|X) = \frac{1}{1 + e^{-(\beta_0 + \beta_1x_1 + \cdots + \beta_nx_n)}}
    \]  
    Multi-class:  
    \[
    P(y=k|X) = \frac{e^{\beta_k^T X}}{\sum_{j=1}^K e^{\beta_j^T X}}
    \]  
  - **Relevance**: Lightweight, interpretable, suitable for low-resource settings.  
  - **Literature**:  
    - Smith & Sheltzer (2022) used Cox models on TCGA, adaptable for classification.  
    - Oncotype DX (Sparano et al., 2018) uses logistic regression for breast cancer risk.

- **Area Under the ROC Curve (AUC-ROC)**  
  - **Description**: Measures ability to distinguish cancer vs. non-cancer.  
  - **Formula**:  
    \[
    \text{AUC} = \int_0^1 \text{TPR}(\text{FPR}) \, d\text{FPR}
    \]  
  - **Relevance**: Standard for imbalanced datasets, interpretable.  
  - **Literature**:  
    - Breast cancer study (2023) achieved high AUC-ROC on TCGA gene expression data.  
    - Elbashir et al. (2022) reported 98.76% precision for CNN-based classification.

### b. Survival-Based Metrics
- **Cox Proportional Hazards Model**  
  - **Description**: Estimates hazard ratio for cancer based on genomic features.  
  - **Formula**:  
    \[
    h(t|X) = h_0(t) \exp(\beta_1x_1 + \cdots + \beta_nx_n)
    \]  
  - **Relevance**: Predicts relative risk, handles TCGA survival data.  
  - **Literature**:  
    - Smith & Sheltzer (2022) identified 100,000+ biomarkers using Cox on TCGA.  
    - Pan-cancer study (2020) validated TCGA clinical endpoints.

- **Kaplan-Meier Survival Curves**  
  - **Description**: Estimates survival probability, stratified by risk groups.  
  - **Formula**:  
    \[
    S(t) = \prod_{t_i \leq t} \left(1 - \frac{d_i}{n_i}\right)
    \]  
  - **Relevance**: Visualizes risk differences for React frontend.  
  - **Literature**:  
    - UCSC Cancer Genomics Browser (2019) used Kaplan-Meier for TCGA survival.  
    - Deep learning study (2021) combined Kaplan-Meier with Cox for 10 cancer types.

### c. Discriminative Metrics
- **Accuracy, Sensitivity, Specificity, F1-Score**  
  - **Description**: Evaluates classification performance.  
  - **Formulas**:  
    \[
    \text{Accuracy} = \frac{\text{TP} + \text{TN}}{\text{TP} + \text{TN} + \text{FP} + \text{FN}}
    \]  
    \[
    \text{Sensitivity} = \frac{\text{TP}}{\text{TP} + \text{FN}}
    \]  
    \[
    \text{Specificity} = \frac{\text{TN}}{\text{TN} + \text{FP}}
    \]  
    \[
    \text{F1-Score} = 2 \cdot \frac{\text{Precision} \cdot \text{Recall}}{\text{Precision} + \text{Recall}}
    \]  
  - **Relevance**: Sensitivity minimizes false negatives, specificity reduces false positives.  
  - **Literature**:  
    - TCGA whole-genome study (2023) reported 98% accuracy.  
    - Elbashir et al. (2022) achieved high F1-scores for breast cancer.

- **Matthews Correlation Coefficient (MCC)**  
  - **Description**: Balanced metric for imbalanced datasets.  
  - **Formula**:  
    \[
    \text{MCC} = \frac{\text{TP} \cdot \text{TN} - \text{FP} \cdot \text{FN}}{\sqrt{(\text{TP} + \text{FP})(\text{TP} + \text{FN})(\text{TN} + \text{FP})(\text{TN} + \text{FN})}}
    \]  
  - **Relevance**: Robust single metric for performance.  
  - **Literature**: Used in TCGA ML studies (2023).

### d. Interpretability Metrics
- **Shapley Values (SHAP)**  
  - **Description**: Assigns contribution scores to genomic features.  
  - **Formula**:  
    \[
    \phi_i = \sum_{S \subseteq N \setminus \{i\}} \frac{|S|!(|N|-|S|-1)!}{|N|!} [f(S \cup \{i\}) - f(S)]
    \]  
  - **Relevance**: Explains variant contributions for interpretable reports.  
  - **Literature**:  
    - Breast cancer study (2023) used SHAP on TCGA data.  
    - ML review (2022) highlights SHAP for genomics.

## 3. Proposed Mathematical Model
- **Feature Selection**: Filter variants, map to TCGA, reduce dimensionality (LASSO/PCA).  
- **Model Training**:  
  - Logistic regression for cancer type probabilities.  
  - Cox regression for hazard ratios.  
  - Ensemble: Weighted average of probabilities and hazard ratios.  
- **Prediction Output**:  
  - Probability scores (0–1) per cancer type.  
  - Hazard ratios or risk percentiles.  
  - SHAP values for variant contributions.  
- **Validation Metrics**: AUC-ROC, sensitivity, specificity, F1-score, MCC, Kaplan-Meier curves, log-rank test.

## 4. Validation Strategy
- **Cross-Validation**: 5/10-fold on TCGA data.  
- **External Validation**: Use 1000 Genomes or ICGC datasets.  
- **Benchmarking**: Compare to TCGA-based models (2021, 2023).  
- **Sensitivity Analysis**: Test robustness to noisy FASTQ/BAM files.  
- **Clinical Relevance**: Correlate with TCGA survival endpoints.

## 5. Literature Support
- Smith & Sheltzer (2022): Cox models for TCGA biomarkers.  
- UCSC Cancer Genomics Browser (2019): Kaplan-Meier for survival.  
- Deep Learning on TCGA (2023): High-accuracy cancer prediction.  
- Breast Cancer Prediction (2023): SHAP for gene ranking.  
- Pan-Cancer Study (2020): Validated TCGA clinical endpoints.  
- ML Review (2022): Summarizes CNNs and SHAP for cancer.

## 6. Implementation Considerations
- **Efficiency**: Use lightweight logistic/Cox models, optimize TensorFlow.  
- **Privacy**: Store TCGA data locally for GDPR/HIPAA compliance.  
- **Scalability**: Batch process FASTQ/BAM, parallelize variant calling.  
- **Frontend**: Visualize probabilities, hazard ratios, SHAP via React/Tailwind.

## 7. Example Output
**JSON Output**:
```json
{
  "patient_id": "P123",
  "cancer_predictions": {
    "breast_cancer": {"probability": 0.85, "hazard_ratio": 2.3},
    "lung_adenocarcinoma": {"probability": 0.12, "hazard_ratio": 1.1}
  },
  "top_variants": [
    {"gene": "TP53", "variant": "p.R175H", "shap_value": 0.45},
    {"gene": "BRCA1", "variant": "c.5266dupC", "shap_value": 0.30}
  ],
  "metrics": {
    "auc_roc": 0.92,
    "sensitivity": 0.89,
    "specificity": 0.87,
    "f1_score": 0.88,
    "mcc": 0.76
  }
}
```
**Visualization**: Heatmaps, Kaplan-Meier curves, probability tables.

## 8. Conclusion
The proposed metrics (logistic regression, Cox regression, AUC-ROC, sensitivity/specificity, F1-score, MCC, Kaplan-Meier, SHAP) ensure accurate, interpretable, and efficient cancer risk prediction. Supported by TCGA studies, the hybrid model balances performance and accessibility for GenePredict’s goals.