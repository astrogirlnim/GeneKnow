# Cancer Risk Prediction Metrics for GenePredict

## 1. Overview
GenePredict aims to predict cancer risk by comparing patient genomic data (FASTQ/BAM files) to TCGA profiles. The "Disease Risk Model" stage uses TensorFlow to generate probability or risk scores for specific cancer types. This report outlines mathematical metrics for prediction, supported by peer-reviewed scientific literature, ensuring interpretability, statistical rigor, and suitability for low-resource settings. A new validation step compares model outputs against TCGA diagnostic data (short-term) and establishes benchmarks for new patient data (long-term), enhancing reliability and generalizability.

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
    - Smith & Sheltzer (2022) used Cox models on TCGA, adaptable for classification. [Link](http://www.tcga-survival.com)
    - Oncotype DX (Sparano et al., 2018) uses logistic regression for breast cancer risk. [Link](https://pubmed.ncbi.nlm.nih.gov/29878226/)

- **Area Under the ROC Curve (AUC-ROC)**  
  - **Description**: Measures ability to distinguish cancer vs. non-cancer.  
  - **Formula**:  
    \[
    \text{AUC} = \int_0^1 \text{TPR}(\text{FPR}) \, d\text{FPR}
    \]  
  - **Relevance**: Standard for imbalanced datasets, interpretable.  
  - **Literature**:  
    - Breast cancer study (2023) achieved high AUC-ROC on TCGA gene expression data. [Link](https://pubmed.ncbi.nlm.nih.gov/36658421/)
    - Elbashir et al. (2022) reported 98.76% precision for CNN-based classification. [Link](https://pubmed.ncbi.nlm.nih.gov/35678312/)

### b. Survival-Based Metrics
- **Cox Proportional Hazards Model**  
  - **Description**: Estimates hazard ratio for cancer based on genomic features.  
  - **Formula**:  
    \[
    h(t|X) = h_0(t) \exp(\beta_1x_1 + \cdots + \beta_nx_n)
    \]  
  - **Relevance**: Predicts relative risk, handles TCGA survival data.  
  - **Literature**:  
    - Smith & Sheltzer (2022) identified 100,000+ biomarkers using Cox on TCGA. [Link](http://www.tcga-survival.com)
    - Pan-cancer study (2020) validated TCGA clinical endpoints. [Link](https://pubmed.ncbi.nlm.nih.gov/32472114/)

- **Kaplan-Meier Survival Curves**  
  - **Description**: Estimates survival probability, stratified by risk groups.  
  - **Formula**:  
    \[
    S(t) = \prod_{t_i \leq t} \left(1 - \frac{d_i}{n_i}\right)
    \]  
  - **Relevance**: Visualizes risk differences for React frontend.  
  - **Literature**:  
    - UCSC Cancer Genomics Browser (2019) used Kaplan-Meier for TCGA survival. [Link](https://pubmed.ncbi.nlm.nih.gov/30936559/)
    - Deep learning study (2021) combined Kaplan-Meier with Cox for 10 cancer types. [Link](https://pubmed.ncbi.nlm.nih.gov/33782624/)

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
    - TCGA whole-genome study (2023) reported 98% accuracy. [Link](https://pubmed.ncbi.nlm.nih.gov/36725512/)
    - Elbashir et al. (2022) achieved high F1-scores for breast cancer. [Link](https://pubmed.ncbi.nlm.nih.gov/35678312/)

- **Matthews Correlation Coefficient (MCC)**  
  - **Description**: Balanced metric for imbalanced datasets.  
  - **Formula**:  
    \[
    \text{MCC} = \frac{\text{TP} \cdot \text{TN} - \text{FP} \cdot \text{FN}}{\sqrt{(\text{TP} + \text{FP})(\text{TP} + \text{FN})(\text{TN} + \text{FP})(\text{TN} + \text{FN})}}
    \]  
  - **Relevance**: Robust single metric for performance.  
  - **Literature**: Used in TCGA ML studies (2023). [Link](https://pubmed.ncbi.nlm.nih.gov/36725512/)

### d. Interpretability Metrics
- **Shapley Values (SHAP)**  
  - **Description**: Assigns contribution scores to genomic features.  
  - **Formula**:  
    \[
    \phi_i = \sum_{S \subseteq N \setminus \{i\}} \frac{|S|!(|N|-|S|-1)!}{|N|!} [f(S \cup \{i\}) - f(S)]
    \]  
  - **Relevance**: Explains variant contributions for interpretable reports.  
  - **Literature**:  
    - Breast cancer study (2023) used SHAP on TCGA data. [Link](https://pubmed.ncbi.nlm.nih.gov/36658421/)
    - ML review (2022) highlights SHAP for genomics. [Link](https://pubmed.ncbi.nlm.nih.gov/35242387/)

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
- Smith & Sheltzer (2022): Cox models for TCGA biomarkers. [Link](http://www.tcga-survival.com)  
- UCSC Cancer Genomics Browser (2019): Kaplan-Meier for survival. [Link](https://pubmed.ncbi.nlm.nih.gov/30936559/)  
- Deep Learning on TCGA (2023): High-accuracy cancer prediction. [Link](https://pubmed.ncbi.nlm.nih.gov/36725512/)  
- Breast Cancer Prediction (2023): SHAP for gene ranking. [Link](https://pubmed.ncbi.nlm.nih.gov/36658421/)  
- Pan-Cancer Study (2020): Validated TCGA clinical endpoints. [Link](https://pubmed.ncbi.nlm.nih.gov/32472114/)  
- ML Review (2022): Summarizes CNNs and SHAP for cancer. [Link](https://pubmed.ncbi.nlm.nih.gov/35242387/)

## 6. Model Output Validation
A new validation step ensures model reliability by comparing predictions with TCGA diagnostic data (short-term) and benchmarking new patient data (long-term). This step is implemented as a new LangGraph node, "Model Output Validation."

### 6.1. Short-term: TCGA Validation
- **Objective**: Compare model predictions (cancer types, hazard ratios) with TCGA clinical diagnoses using preserved patient IDs.
- **Process**:
  - Match patient IDs (e.g., TCGA-XX-XXXX) between model output and TCGA clinical data (from GDC Data Portal: https://portal.gdc.cancer.gov/).
  - Compare predicted cancer types with TCGA `primary_diagnosis`.
  - Compute:
    - **Concordance Rate**: Percentage of matching predictions.
    - **Cohen's Kappa**: Agreement adjusted for chance.
    - **Mean Absolute Error (MAE)**: For hazard ratios vs. survival times.
    - **Log-Rank Test**: For Kaplan-Meier curves of predicted vs. actual risk groups.
- **Tools**:
  - `pandas`: Parse TCGA TSV/JSON clinical data.
  - `scikit-learn`: Compute concordance, Kappa, MAE.
  - `lifelines`: Generate Kaplan-Meier curves and log-rank tests.

### 6.2. Long-term: New Patient Benchmarking
- **Objective**: Establish thresholds for new patient predictions, flagging outliers.
- **Process**:
  - Derive benchmarks from TCGA (e.g., median probabilities, hazard ratios per cancer type).
  - Use external cohorts (ICGC: https://dcc.icgc.org/, 1000 Genomes: https://www.internationalgenome.org/) for population baselines.
  - Flag predictions outside expected ranges (e.g., >2σ from TCGA median probability).
- **Tools**:
  - `numpy`/`scipy`: Statistical thresholding.
  - `pandas`: Manage benchmark data.

### 6.3. Output
- JSON file with validation metrics (concordance, Kappa, MAE, log-rank p-value) and benchmarking flags.
- Visualizations (e.g., concordance table, Kaplan-Meier plots) in the React frontend.

## 7. Updated LangGraph Architecture
The LangGraph flow now includes the "Model Output Validation" node:

```mermaid
flowchart TD
    A["File Input"] --> B["Preprocess Validate + Parse File"]
    B --> C["Variant QC filter low confidence"] & D["Variant Caller DeepVariant Samtools"]
    D --> E["Variant Mapper ← reference genome + TCGA"]
    E --> F["Disease Risk Model"]
    F --> F1["Feature Selection LASSO/PCA"]
    F1 --> F2["Logistic + Cox Models TensorFlow"]
    F2 --> F3["Compute Metrics AUC-ROC, SHAP, Kaplan-Meier"]
    F3 --> G["Interpret Results JSON Output"]
    G --> H["LLM Report Writer markdown pdf json ← Llama 3.1"]
    H --> I["Model Output Validation ← TCGA Clinical Data"]
    I --> J["Visualizations + Frontend charts heatmaps tables ← React + Tailwind"]
    C --> F
```

## 8. Implementation Considerations
- **Efficiency**: Use lightweight logistic/Cox models, optimize TensorFlow with pruning/quantization. Cache TCGA clinical data for fast validation.
- **Privacy**: Store TCGA data and patient outputs locally, encrypted, to comply with GDPR/HIPAA. Avoid external API calls.
- **Scalability**: Batch process FASTQ/BAM files, parallelize variant calling and validation using `joblib`. Handle large TCGA datasets efficiently.
- **Validation Integration**: Ensure JSON output includes validation metrics (concordance, Kappa, MAE) for frontend display. Update React components to show alerts for benchmarked outliers.
- **Frontend**: Visualize probabilities, hazard ratios, SHAP values, and validation metrics (e.g., concordance rate, Kaplan-Meier curves) via React/Tailwind.

## 9. Example Output
**JSON Output** (from Disease Risk Model and Validation):
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
  },
  "validation": {
    "tcga_validation": {
      "concordance_rate": 92.5,
      "cohen_kappa": 0.87,
      "mae_hazard": 0.15,
      "logrank_p": 0.001
    },
    "benchmarking": {
      "breast_cancer": {
        "flags": ["Probability out of range: 0.85"],
        "deviation_score": 0.35
      }
    }
  }
}
```
**Visualization**:
- Heatmaps of SHAP values.
- Kaplan-Meier curves for predicted vs. actual survival.
- Tables of probabilities, hazard ratios, and validation metrics (e.g., 92.5% concordance).

## 10. Conclusion
The proposed metrics (logistic regression, Cox regression, AUC-ROC, sensitivity/specificity, F1-score, MCC, Kaplan-Meier, SHAP) ensure accurate, interpretable, and efficient cancer risk prediction. The new "Model Output Validation" node enhances reliability by validating against TCGA diagnoses and benchmarking new patient data, using tools like `pandas`, `scikit-learn`, and `lifelines`. Supported by TCGA studies, the hybrid model and validation step balance performance and accessibility, aligning with GenePredict's mission for precision medicine in low-resource settings.
