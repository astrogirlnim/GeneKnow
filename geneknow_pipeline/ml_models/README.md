# ML Fusion Layer for GeneKnow

This directory contains the **ML Fusion Layer** implementation for GeneKnow's cancer risk assessment pipeline. The fusion layer is the meta-learning component that combines outputs from the 5 static models into a final risk assessment.

## 🎯 Overview

The fusion layer implements the architecture described in `Static-Risk-Model.md`. Instead of predicting variant pathogenicity directly from genomic features, it learns optimal weights for combining pre-computed outputs from:

1. **PRS (Polygenic Risk Score)** - Inherited background risk
2. **ClinVar** - Known pathogenic/benign classifications  
3. **CADD** - Functional impact predictions
4. **TCGA** - Tumor enrichment evidence
5. **Gene Burden** - Pathway-level risk assessment

## 📁 Files

- `fusion_layer.py` - Core fusion layer implementation
- `train_fusion_layer.py` - Training script with synthetic data
- `ml_fusion_node.py` - LangGraph integration node
- `README.md` - This documentation

## 🔧 Architecture

### Input Format
The fusion layer expects 5 numerical/categorical features per variant:

```python
StaticModelInputs(
    prs_score=0.8,                          # 0.0-1.0 (polygenic risk)
    clinvar_classification='pathogenic',     # 'pathogenic', 'benign', 'uncertain', 'not_found'
    cadd_score=25.0,                        # 0.0-50.0 (deleteriousness score)
    tcga_enrichment=3.0,                    # 0.1-20.0 (fold enrichment in tumors)
    gene_burden_score=2.0                   # 0.0-10.0 (damaging variants in key genes)
)
```

### Output Format
The fusion layer produces structured risk assessments:

```python
FusionOutput(
    risk_score=0.75,                        # 0.0-1.0 (final risk score)
    confidence=0.85,                        # 0.0-1.0 (prediction confidence)
    contributing_factors={                  # Feature importance for this prediction
        'clinvar_pathogenic': 0.45,
        'cadd_score': 0.30,
        'prs_score': 0.15,
        'tcga_enrichment': 0.08,
        'gene_burden_score': 0.02
    },
    risk_category='high'                    # 'low', 'moderate', 'high', 'very_high'
)
```

## 🚀 Quick Start

### 1. Train the Fusion Layer

```bash
cd geneknow_pipeline/ml_models
python train_fusion_layer.py
```

This will:
- Generate 5,000 synthetic training samples
- Train gradient boosting, random forest, and logistic regression models
- Save the best model as `best_fusion_model.pkl`
- Create visualization in `training_results.png`

### 2. Test the Fusion Layer

```python
from fusion_layer import FusionLayer, StaticModelInputs

# Load trained model
fusion = FusionLayer()
fusion.load_model('best_fusion_model.pkl')

# Make prediction
test_input = StaticModelInputs(
    prs_score=0.8,
    clinvar_classification='pathogenic',
    cadd_score=25.0,
    tcga_enrichment=3.0,
    gene_burden_score=2.0
)

prediction = fusion.predict(test_input)
print(f"Risk Score: {prediction.risk_score:.3f}")
print(f"Risk Category: {prediction.risk_category}")
```

### 3. Use in Pipeline

The fusion layer integrates into the LangGraph pipeline via `ml_fusion_node.py`:

```python
from nodes.ml_fusion_node import ml_fusion_node

# Add to your pipeline graph
graph = StateGraph(GeneKnowState)
graph.add_node("ml_fusion", ml_fusion_node)
```

## 📊 Training Data

Currently uses synthetic data that simulates realistic static model outputs:

- **PRS scores**: Beta distribution (most people have low inherited risk)
- **ClinVar classifications**: 5% pathogenic, 70% benign, 20% uncertain, 5% not found
- **CADD scores**: Exponential distribution (few high-impact variants)
- **TCGA enrichment**: Log-normal distribution (rare high enrichment)
- **Gene burden**: Poisson distribution (most people have 0-2 damaging variants)

### Real Training Data Integration

To use real training data instead of synthetic:

1. **Collect Pipeline Outputs**: Run your pipeline on known cases and collect the 5 static model outputs
2. **Label with Ground Truth**: Assign risk scores based on known outcomes (cancer diagnosis, family history, etc.)
3. **Replace Synthetic Data**: Modify `train_fusion_layer.py` to load your real data instead of calling `create_synthetic_training_data()`

## 🎯 Model Performance

### Synthetic Data Results
- **Gradient Boosting**: Best performance (lowest MSE)
- **Random Forest**: Good performance, more interpretable
- **Linear Regression**: Fastest inference, adequate performance

### Feature Importance (Typical)
1. **ClinVar pathogenic** (~40-50%) - Known pathogenic variants dominate
2. **CADD score** (~20-30%) - Functional impact is important
3. **PRS score** (~10-20%) - Background genetic risk
4. **TCGA enrichment** (~5-10%) - Tumor evidence supports risk
5. **Gene burden** (~5-10%) - Pathway-level effects

## 🔧 Configuration

### Risk Thresholds
Default risk categorization:
- **Low**: 0.0 - 0.25
- **Moderate**: 0.25 - 0.5  
- **High**: 0.5 - 0.75
- **Very High**: 0.75 - 1.0

### Model Parameters
- **Gradient Boosting**: 100 estimators, depth 3, learning rate 0.1
- **Random Forest**: 100 estimators, depth 5
- **Linear Regression**: Simple linear model, no regularization

## 📈 Pipeline Integration

The fusion layer expects the following data to be available in the pipeline state:

```python
state = {
    'variants': [
        {
            'prs_score': 0.8,                                    # From PRS calculation
            'clinvar': {'clinical_significance': 'Pathogenic'},  # From ClinVar annotation
            'cadd_score': 25.0,                                  # From CADD scoring
            'tcga_enrichment': 3.0,                              # From TCGA frequency matching
            'gene_burden_score': 2.0,                            # From gene burden analysis
            # ... other variant data
        }
    ]
}
```

The fusion node will:
1. Extract static model outputs for each variant
2. Run fusion layer predictions
3. Calculate aggregate risk across all variants
4. Add results to `state['ml_fusion_results']`

## 🧪 Testing

### Unit Tests
```bash
python fusion_layer.py          # Test core functionality
python train_fusion_layer.py    # Test training pipeline
python ml_fusion_node.py        # Test LangGraph integration
```

### Integration Tests
```bash
# Test with mock pipeline state
python -c "
from ml_fusion_node import ml_fusion_node
mock_state = {'variants': [{'prs_score': 0.8, 'clinvar': {'clinical_significance': 'Pathogenic'}, 'cadd_score': 25.0, 'tcga_enrichment': 3.0, 'gene_burden_score': 2.0}]}
result = ml_fusion_node(mock_state)
print(result['ml_fusion_results'])
"
```

## 📋 TODO

### Immediate
- [ ] Replace synthetic training data with real pipeline outputs
- [ ] Add confidence interval estimation
- [ ] Implement ensemble uncertainty quantification

### Future Enhancements
- [ ] Add temporal models (risk changes over time)
- [ ] Include demographic factors (age, sex, ancestry)
- [ ] Add model drift detection
- [ ] Implement active learning for model updates

## 🔍 Troubleshooting

### Common Issues

**Model Not Loading**
```
ML Fusion model not found at geneknow_pipeline/ml_models/best_fusion_model.pkl
```
**Solution**: Run `python train_fusion_layer.py` first to train the model.

**Missing Static Model Outputs**
```
No static model inputs found. Skipping ML fusion.
```
**Solution**: Ensure your pipeline populates `prs_score`, `clinvar`, `cadd_score`, `tcga_enrichment`, and `gene_burden_score` fields.

**Import Errors**
```
ModuleNotFoundError: No module named 'fusion_layer'
```
**Solution**: Make sure you're running from the correct directory and have the proper Python path setup.

## 📚 References

- `Static-Risk-Model.md` - Original architecture specification
- `ML-Risk-Model.md` - ML model requirements
- `Risk-Fusion-Architecture.md` - Fusion layer design

## 🤝 Contributing

When adding new features:
1. Update the `StaticModelInputs` dataclass for new input types
2. Modify the `_extract_static_model_outputs()` method in `ml_fusion_node.py`
3. Retrain the fusion layer with representative data
4. Update tests and documentation

---

**The fusion layer is the brain of GeneKnow's risk assessment** - it intelligently combines multiple sources of evidence to provide the most accurate and interpretable cancer risk predictions possible. 