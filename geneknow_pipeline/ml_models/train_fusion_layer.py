#!/usr/bin/env python3
"""
Training Script for the ML Fusion Layer

This script trains the fusion layer using synthetic data that represents
the outputs from the 5 static models (PRS, ClinVar, CADD, TCGA, Gene Burden).

In production, this would be replaced with real training data collected
from your pipeline runs.
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from fusion_layer import FusionLayer, StaticModelInputs, create_synthetic_training_data
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def train_fusion_models(n_samples: int = 5000, model_types: list = None):
    """
    Train fusion layer models with different algorithms.
    
    Args:
        n_samples: Number of synthetic training samples
        model_types: List of model types to train
    """
    if model_types is None:
        model_types = ['gradient_boosting', 'random_forest', 'linear']
    
    logger.info(f"Training fusion layer with {n_samples} samples")
    
    # Generate synthetic training data
    print("📊 Generating synthetic training data...")
    training_data = create_synthetic_training_data(n_samples=n_samples)
    
    # Statistics about the training data
    risk_scores = [item[1] for item in training_data]
    print(f"Training data statistics:")
    print(f"  Mean risk score: {np.mean(risk_scores):.3f}")
    print(f"  Std risk score: {np.std(risk_scores):.3f}")
    print(f"  Min risk score: {np.min(risk_scores):.3f}")
    print(f"  Max risk score: {np.max(risk_scores):.3f}")
    
    # Risk category distribution
    risk_categories = []
    for score in risk_scores:
        if score <= 0.25:
            risk_categories.append('low')
        elif score <= 0.5:
            risk_categories.append('moderate')
        elif score <= 0.75:
            risk_categories.append('high')
        else:
            risk_categories.append('very_high')
    
    from collections import Counter
    category_counts = Counter(risk_categories)
    print(f"Risk category distribution:")
    for category, count in category_counts.items():
        print(f"  {category}: {count} ({count/len(risk_categories)*100:.1f}%)")
    
    # Train different model types
    results = {}
    best_model = None
    best_score = float('inf')
    
    for model_type in model_types:
        print(f"\n🎯 Training {model_type} model...")
        
        fusion = FusionLayer(model_type=model_type)
        training_results = fusion.train(training_data, validation_split=0.2)
        
        results[model_type] = training_results
        
        # Save model
        model_path = f'fusion_{model_type}.pkl'
        fusion.save_model(model_path)
        
        # Track best model
        if training_results['val_mse'] < best_score:
            best_score = training_results['val_mse']
            best_model = model_type
        
        print(f"✅ {model_type} completed:")
        print(f"  Validation MSE: {training_results['val_mse']:.4f}")
        print(f"  CV MSE: {training_results['cv_mse_mean']:.4f} ± {training_results['cv_mse_std']:.4f}")
        
        if training_results['feature_importance']:
            print("📈 Top 5 features:")
            sorted_features = sorted(training_results['feature_importance'].items(), 
                                   key=lambda x: x[1], reverse=True)[:5]
            for feature, importance in sorted_features:
                print(f"  {feature}: {importance:.3f}")
    
    # Save best model as default
    if best_model:
        best_fusion = FusionLayer(model_type=best_model)
        best_fusion.train(training_data)
        best_fusion.save_model('best_fusion_model.pkl')
        print(f"\n🏆 Best model: {best_model} (MSE: {best_score:.4f})")
    
    # Save training results
    training_metadata = {
        'n_samples': n_samples,
        'model_types': model_types,
        'best_model': best_model,
        'best_score': best_score,
        'results': results,
        'risk_distribution': dict(category_counts)
    }
    
    with open('training_metadata.json', 'w') as f:
        json.dump(training_metadata, f, indent=2, default=str)
    
    # Create visualization
    create_training_visualization(results, risk_scores)
    
    return results

def create_training_visualization(results, risk_scores):
    """Create visualization of training results."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Model comparison
    model_names = list(results.keys())
    val_mse = [results[model]['val_mse'] for model in model_names]
    cv_mse = [results[model]['cv_mse_mean'] for model in model_names]
    
    axes[0, 0].bar(model_names, val_mse, alpha=0.7, label='Validation MSE')
    axes[0, 0].bar(model_names, cv_mse, alpha=0.7, label='CV MSE')
    axes[0, 0].set_title('Model Performance Comparison')
    axes[0, 0].set_ylabel('MSE')
    axes[0, 0].legend()
    axes[0, 0].tick_params(axis='x', rotation=45)
    
    # Risk score distribution
    axes[0, 1].hist(risk_scores, bins=30, alpha=0.7, edgecolor='black')
    axes[0, 1].set_title('Risk Score Distribution')
    axes[0, 1].set_xlabel('Risk Score')
    axes[0, 1].set_ylabel('Frequency')
    
    # Feature importance (for best model)
    best_model = min(results.keys(), key=lambda x: results[x]['val_mse'])
    if results[best_model]['feature_importance']:
        features = list(results[best_model]['feature_importance'].keys())
        importances = list(results[best_model]['feature_importance'].values())
        
        # Sort by importance
        sorted_pairs = sorted(zip(features, importances), key=lambda x: x[1], reverse=True)
        features, importances = zip(*sorted_pairs)
        
        axes[1, 0].barh(features, importances, alpha=0.7)
        axes[1, 0].set_title(f'Feature Importance ({best_model})')
        axes[1, 0].set_xlabel('Importance')
    
    # Model complexity comparison
    n_samples = [results[model]['n_samples'] for model in model_names]
    axes[1, 1].scatter(n_samples, val_mse, alpha=0.7)
    for i, model in enumerate(model_names):
        axes[1, 1].annotate(model, (n_samples[i], val_mse[i]))
    axes[1, 1].set_title('Sample Size vs Performance')
    axes[1, 1].set_xlabel('Training Samples')
    axes[1, 1].set_ylabel('Validation MSE')
    
    plt.tight_layout()
    plt.savefig('training_results.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("📊 Training visualization saved to training_results.png")

def test_fusion_predictions():
    """Test the fusion layer with various input scenarios."""
    print("\n🧪 Testing fusion layer predictions...")
    
    # Load best model
    try:
        fusion = FusionLayer()
        fusion.load_model('best_fusion_model.pkl')
    except FileNotFoundError:
        print("❌ Best model not found. Train models first.")
        return
    
    # Test scenarios
    test_cases = [
        {
            'name': 'High Risk Case',
            'inputs': StaticModelInputs(
                prs_score=0.9,
                clinvar_classification='pathogenic',
                cadd_score=30.0,
                tcga_enrichment=8.0,
                gene_burden_score=5.0
            )
        },
        {
            'name': 'Low Risk Case',
            'inputs': StaticModelInputs(
                prs_score=0.1,
                clinvar_classification='benign',
                cadd_score=5.0,
                tcga_enrichment=0.5,
                gene_burden_score=0.0
            )
        },
        {
            'name': 'Moderate Risk Case',
            'inputs': StaticModelInputs(
                prs_score=0.5,
                clinvar_classification='uncertain',
                cadd_score=15.0,
                tcga_enrichment=2.0,
                gene_burden_score=1.0
            )
        },
        {
            'name': 'Unknown Variant Case',
            'inputs': StaticModelInputs(
                prs_score=0.3,
                clinvar_classification='not_found',
                cadd_score=22.0,
                tcga_enrichment=3.0,
                gene_burden_score=2.0
            )
        }
    ]
    
    for test_case in test_cases:
        print(f"\n📋 {test_case['name']}:")
        prediction = fusion.predict(test_case['inputs'])
        
        print(f"  Risk Score: {prediction.risk_score:.3f}")
        print(f"  Risk Category: {prediction.risk_category}")
        print(f"  Confidence: {prediction.confidence:.3f}")
        
        # Show top contributing factors
        if prediction.contributing_factors:
            top_factors = sorted(prediction.contributing_factors.items(), 
                               key=lambda x: x[1], reverse=True)[:3]
            print(f"  Top factors:")
            for factor, contribution in top_factors:
                print(f"    {factor}: {contribution:.3f}")

if __name__ == "__main__":
    print("🔬 GeneKnow Fusion Layer Training")
    print("=" * 50)
    
    # Train models
    results = train_fusion_models(n_samples=5000)
    
    # Test predictions
    test_fusion_predictions()
    
    print("\n✅ Training completed successfully!")
    print("🎯 Models saved in current directory")
    print("📊 Training results saved as training_results.png")
    print("📋 Metadata saved as training_metadata.json") 