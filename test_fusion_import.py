#!/usr/bin/env python3
"""Test if fusion layer can be imported and models loaded."""

import sys
import os

# Add paths
sys.path.append('geneknow_pipeline')
sys.path.append('geneknow_pipeline/ml_models')

print("Testing ML Fusion Layer Import...")
print("=" * 60)

# Test 1: Import fusion layer
try:
    from fusion_layer import FusionLayer, StaticModelInputs, FusionOutput
    print("✅ Successfully imported fusion_layer module")
except Exception as e:
    print(f"❌ Failed to import fusion_layer: {e}")
    sys.exit(1)

# Test 2: Create fusion layer instance
try:
    fusion = FusionLayer()
    print("✅ Successfully created FusionLayer instance")
except Exception as e:
    print(f"❌ Failed to create FusionLayer: {e}")
    sys.exit(1)

# Test 3: Check for model files
model_dir = 'geneknow_pipeline/ml_models'
print(f"\nChecking model files in {model_dir}:")
if os.path.exists(model_dir):
    files = os.listdir(model_dir)
    pkl_files = [f for f in files if f.endswith('.pkl')]
    print(f"Found {len(pkl_files)} .pkl files:")
    for f in pkl_files:
        size = os.path.getsize(os.path.join(model_dir, f))
        print(f"  - {f} ({size/1024:.1f} KB)")
else:
    print(f"❌ Model directory not found: {model_dir}")

# Test 4: Try to load model
model_path = os.path.join(model_dir, 'best_fusion_model.pkl')
try:
    fusion.load_model(model_path)
    print(f"\n✅ Successfully loaded model from {model_path}")
except Exception as e:
    print(f"\n❌ Failed to load model: {e}")
    
    # Try alternative models
    alt_models = ['fusion_gradient_boosting.pkl', 'fusion_random_forest.pkl']
    for alt in alt_models:
        alt_path = os.path.join(model_dir, alt)
        if os.path.exists(alt_path):
            try:
                fusion.load_model(alt_path)
                print(f"✅ Successfully loaded alternative model: {alt}")
                break
            except Exception as e2:
                print(f"❌ Failed to load {alt}: {e2}")

# Test 5: Try prediction
print("\nTesting prediction...")
try:
    test_input = StaticModelInputs(
        prs_score=0.8,
        clinvar_classification='pathogenic',
        cadd_score=25.0,
        tcga_enrichment=3.0,
        gene_burden_score=2.0
    )
    
    result = fusion.predict(test_input)
    print("✅ Prediction successful!")
    print(f"  - Risk score: {result.risk_score:.3f}")
    print(f"  - Confidence: {result.confidence:.3f}")
    print(f"  - Risk category: {result.risk_category}")
except Exception as e:
    print(f"❌ Prediction failed: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 60)
print("Test complete!") 