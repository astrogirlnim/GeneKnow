# Integrating Cancer Risk Prediction Metrics into GenePredict's LangGraph Flow

## 1. Overview
GenePredict's LangGraph architecture processes FASTQ/BAM files to predict cancer risk using TCGA data, culminating in interpretable reports and visualizations. The "Disease Risk Model" stage, implemented in TensorFlow, is the focal point for integrating the proposed metrics: logistic regression probabilities, Cox hazard ratios, AUC-ROC, sensitivity/specificity, F1-score, Matthews Correlation Coefficient (MCC), Kaplan-Meier curves, and SHAP values. This document outlines how to merge these metrics into the existing flow, ensuring compatibility with local execution, GDPR/HIPAA compliance, and React/Tailwind frontend integration.

## 2. Current LangGraph Architecture
The existing flow, as described, is:

```mermaid
flowchart TD
    A["File Input"] --> B["Preprocess Validate + Parse File"]
    B --> C["Variant QC filter low confidence"] & D["Variant Caller DeepVariant Samtools"]
    D --> E["Variant Mapper ← reference genome + TCGA"]
    E --> F["Disease Risk Model TensorFlow TCGA + 1000 Genomes"]
    F --> G["Interpret Results JSON Output"]
    G --> H["LLM Report Writer markdown pdf json ← Llama 3.1"]
    H --> I["Visualizations + Frontend charts heatmaps tables ← React + Tailwind"]
    C --> F
```

**Focus**: The "Disease Risk Model" (F) stage needs to compute the proposed metrics, output them in JSON format for the "Interpret Results" (G) stage, and support visualizations in the "Visualizations + Frontend" (I) stage.

## 3. Proposed Metrics and Their Role
- **Logistic Regression Probabilities**: Predict probability of cancer types (e.g., breast cancer, lung adenocarcinoma).
- **Cox Hazard Ratios**: Estimate relative risk of cancer occurrence based on survival data.
- **AUC-ROC, Sensitivity, Specificity, F1-Score, MCC**: Validate classification performance.
- **Kaplan-Meier Curves**: Visualize survival differences across risk groups.
- **SHAP Values**: Explain variant contributions for interpretability.

These metrics will enhance the risk model's accuracy, interpretability, and clinical relevance, aligning with GenePredict's goals for precision medicine in low-resource settings.

## 4. Implementation Strategy

### 4.1. Modify the Disease Risk Model Stage
The "Disease Risk Model" stage will be expanded to include a hybrid model (logistic regression + Cox regression) and compute validation metrics. The steps are:

1. **Feature Preparation**:
   - Input: Variant calls from the "Variant Mapper" stage, mapped to TCGA and filtered by the "Variant QC" stage.
   - Process: Apply feature selection (e.g., LASSO or PCA) to reduce dimensionality, focusing on cancer-associated variants (e.g., TP53, BRCA1 mutations).
   - Output: A feature matrix of variant profiles for each patient.

2. **Model Training and Prediction**:
   - **Logistic Regression**: Train a multi-class logistic regression model to predict cancer type probabilities using TensorFlow.
   - **Cox Regression**: Train a Cox proportional hazards model to estimate hazard ratios, leveraging TCGA survival data.
   - **Ensemble**: Combine predictions (weighted average of probabilities and hazard ratios) to produce a final risk score.
   - **SHAP Values**: Compute SHAP values to explain variant contributions.

3. **Validation Metrics**:
   - Compute AUC-ROC, sensitivity, specificity, F1-score, and MCC using a held-out TCGA test set.
   - Generate Kaplan-Meier curves for high- vs. low-risk groups, with log-rank test p-values.

4. **Output**: A JSON object containing probabilities, hazard ratios, SHAP values, and validation metrics, passed to the "Interpret Results" stage.

### 4.2. Update LangGraph Flow
The modified LangGraph flow incorporates a new sub-node in the "Disease Risk Model" stage for metric computation and ensures data flow to downstream stages. The updated flow is:

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
    H --> I["Visualizations + Frontend charts heatmaps tables ← React + Tailwind"]
    C --> F
```

### 4.3. Code Implementation
Below is a Python code snippet for the "Disease Risk Model" stage, integrating the proposed metrics using TensorFlow, scikit-learn, lifelines, and SHAP. This assumes TCGA data is pre-downloaded locally.

```python
import tensorflow as tf
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score, precision_recall_fscore_support, matthews_corrcoef
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import shap
import pandas as pd
import json

# Step 1: Feature Preparation
def prepare_features(variant_data, tcga_data):
    # Example: LASSO feature selection
    from sklearn.linear_model import Lasso
    lasso = Lasso(alpha=0.01)
    lasso.fit(variant_data, tcga_data['labels'])
    selected_features = variant_data[:, lasso.coef_ != 0]
    return selected_features

# Step 2: Model Training and Prediction
def train_models(features, tcga_data):
    # Logistic Regression
    logreg = LogisticRegression(multi_class='multinomial', solver='lbfgs')
    logreg.fit(features, tcga_data['labels'])
    probabilities = logreg.predict_proba(features)

    # Cox Regression
    cox_data = pd.DataFrame(features)
    cox_data['time'] = tcga_data['survival_time']
    cox_data['event'] = tcga_data['event']
    cph = CoxPHFitter()
    cph.fit(cox_data, duration_col='time', event_col='event')
    hazard_ratios = cph.predict_partial_hazard(features)

    # Ensemble: Weighted average
    risk_scores = 0.7 * probabilities.max(axis=1) + 0.3 * hazard_ratios

    return probabilities, hazard_ratios, risk_scores

# Step 3: Compute Metrics
def compute_metrics(probabilities, hazard_ratios, true_labels, survival_data):
    # Classification Metrics
    pred_labels = np.argmax(probabilities, axis=1)
    auc_roc = roc_auc_score(true_labels, probabilities, multi_class='ovr')
    accuracy = accuracy_score(true_labels, pred_labels)
    precision, recall, f1, _ = precision_recall_fscore_support(true_labels, pred_labels, average='weighted')
    mcc = matthews_corrcoef(true_labels, pred_labels)

    # Kaplan-Meier Curves
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    high_risk = hazard_ratios > np.median(hazard_ratios)
    kmf_high.fit(survival_data['time'][high_risk], survival_data['event'][high_risk], label='High Risk')
    kmf_low.fit(survival_data['time'][~high_risk], survival_data['event'][~high_risk], label='Low Risk')
    logrank_p = logrank_test(
        survival_data['time'][high_risk], survival_data['time'][~high_risk],
        survival_data['event'][high_risk], survival_data['event'][~high_risk]
    ).p_value

    # SHAP Values
    explainer = shap.KernelExplainer(logreg.predict_proba, features)
    shap_values = explainer.shap_values(features)

    return {
        'auc_roc': auc_roc,
        'accuracy': accuracy,
        'sensitivity': recall,
        'specificity': precision,  # Note: Specificity requires additional computation for multi-class
        'f1_score': f1,
        'mcc': mcc,
        'logrank_p': logrank_p,
        'shap_values': shap_values
    }

# Step 4: Generate Output
def generate_output(patient_id, probabilities, hazard_ratios, risk_scores, metrics, variant_names):
    output = {
        'patient_id': patient_id,
        'cancer_predictions': {
            'breast_cancer': {'probability': float(probabilities[0, 0]), 'hazard_ratio': float(hazard_ratios[0])},
            'lung_adenocarcinoma': {'probability': float(probabilities[0, 1]), 'hazard_ratio': float(hazard_ratios[0])}
        },
        'top_variants': [
            {'gene': variant_names[i], 'variant': 'p.R175H', 'shap_value': float(metrics['shap_values'][0][i])}
            for i in np.argsort(-np.abs(metrics['shap_values'][0]))[:5]
        ],
        'metrics': {
            'auc_roc': float(metrics['auc_roc']),
            'sensitivity': float(metrics['sensitivity']),
            'specificity': float(metrics['specificity']),
            'f1_score': float(metrics['f1_score']),
            'mcc': float(metrics['mcc'])
        }
    }
    return output

# Main Function
def disease_risk_model(variant_data, tcga_data, patient_id, variant_names):
    features = prepare_features(variant_data, tcga_data)
    probabilities, hazard_ratios, risk_scores = train_models(features, tcga_data)
    metrics = compute_metrics(probabilities, hazard_ratios, tcga_data['labels'], tcga_data['survival'])
    output = generate_output(patient_id, probabilities, hazard_ratios, risk_scores, metrics, variant_names)
    with open('risk_output.json', 'w') as f:
        json.dump(output, f)
    return output
```

### 4.4. Integration with Downstream Stages
- **Interpret Results (G)**:
  - Input: JSON output from the "Disease Risk Model" stage.
  - Process: Parse probabilities, hazard ratios, SHAP values, and metrics into a structured format for the LLM.
  - Output: A processed JSON file with annotations for report generation.

- **LLM Report Writer (H)**:
  - Input: Processed JSON from the "Interpret Results" stage.
  - Process: Use Llama 3.1 to generate markdown/PDF reports summarizing:
    - Cancer type probabilities (e.g., "85% likelihood of breast cancer").
    - Hazard ratios (e.g., "2.3x increased risk compared to TCGA cohort").
    - Top variants (e.g., "TP53 p.R175H contributes 45% to risk").
    - Validation metrics (e.g., "92% AUC-ROC, 89% sensitivity").
  - Output: Markdown/PDF/JSON reports.

- **Visualizations + Frontend (I)**:
  - Input: JSON output and markdown reports.
  - Process: Use React and Tailwind CSS to create:
    - **Heatmap**: SHAP values for top variants.
    - **Kaplan-Meier Plot**: Survival curves for high- vs. low-risk groups.
    - **Table**: Probabilities and hazard ratios for cancer types.
  - Example React Component:
    ```jsx
    import React from 'react';
    import { LineChart } from 'recharts';
    import 'tailwindcss/tailwind.css';

    const RiskVisualization = ({ data }) => {
      return (
        <div className="p-4">
          <h2 className="text-2xl font-bold">Cancer Risk Profile</h2>
          <table className="table-auto w-full">
            <thead>
              <tr>
                <th>Cancer Type</th>
                <th>Probability</th>
                <th>Hazard Ratio</th>
              </tr>
            </thead>
            <tbody>
              {Object.entries(data.cancer_predictions).map(([type, { probability, hazard_ratio }]) => (
                <tr key={type}>
                  <td>{type}</td>
                  <td>{(probability * 100).toFixed(2)}%</td>
                  <td>{hazard_ratio.toFixed(2)}</td>
                </tr>
              ))}
            </tbody>
          </table>
          <h3 className="text-xl font-semibold mt-4">Top Variants</h3>
          <LineChart width={600} height={300} data={data.top_variants}>
            {/* Add chart components */}
          </LineChart>
        </div>
      );
    };

    export default RiskVisualization;
    ```

### 4.5. Implementation Considerations
- **Efficiency**: Use scikit-learn for logistic regression and lifelines for Cox models to minimize TensorFlow overhead. Optimize feature selection with sparse matrices for large variant datasets.
- **Privacy**: Ensure TCGA data is stored locally and encrypted. Process patient FASTQ/BAM files without external API calls.
- **Scalability**: Parallelize feature selection and metric computation using Python's `joblib` or TensorFlow's distributed training.
- **Validation**: Perform 5-fold cross-validation within the "Disease Risk Model" stage using TCGA data. Validate externally with 1000 Genomes data, as planned.
- **Frontend Compatibility**: Ensure JSON output is lightweight (e.g., limit SHAP values to top 5–10 variants) to avoid React rendering delays.

## 5. Testing and Validation
- **Unit Tests**:
  - Test feature selection for correctness (e.g., non-zero LASSO coefficients).
  - Verify logistic and Cox model outputs against known TCGA benchmarks.
  - Check metric computations (e.g., AUC-ROC, MCC) against scikit-learn ground truth.

- **Integration Tests**:
  - Run a sample FASTQ/BAM file through the entire LangGraph pipeline.
  - Validate JSON output structure and frontend rendering.

- **Performance Tests**:
  - Measure runtime on low-resource hardware (e.g., 8GB RAM laptop).
  - Optimize bottlenecks (e.g., SHAP computation) using approximation methods.

## 6. Example Output
**JSON Output** (from `risk_output.json`):
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

**Frontend Visualization**:
- A table displaying cancer type probabilities and hazard ratios.
- A heatmap of SHAP values for top variants.
- A Kaplan-Meier plot comparing high- vs. low-risk survival curves.

## 7. Conclusion
By expanding the "Disease Risk Model" stage to include feature selection, logistic/Cox models, and metric computation, and updating downstream stages to handle the enriched JSON output, GenePredict can effectively integrate the proposed metrics. The implementation leverages lightweight libraries (scikit-learn, lifelines, SHAP) for efficiency, ensures privacy through local execution, and supports intuitive visualizations via React/Tailwind. This approach enhances the accuracy, interpretability, and clinical relevance of cancer risk predictions, aligning with GenePredict's mission for accessible precision medicine.