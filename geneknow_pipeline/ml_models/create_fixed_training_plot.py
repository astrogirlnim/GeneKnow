#!/usr/bin/env python3
"""
Create training results plot for the fixed (no data leakage) model.
"""
import json
import matplotlib.pyplot as plt
import numpy as np

# Load the results
with open('../real_data_training_results_FIXED.json', 'r') as f:
    results = json.load(f)

# Create the plot
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

# Model comparison
models = list(results['all_models'].keys())
metrics = ['val_mse', 'val_r2', 'val_accuracy', 'val_auc']
metric_labels = ['MSE (↓)', 'R² (↑)', 'Accuracy (↑)', 'AUC (↑)']

for i, (metric, label) in enumerate(zip(metrics, metric_labels)):
    ax = [ax1, ax2, ax3, ax4][i]
    values = [results['all_models'][model][metric] for model in models]
    
    bars = ax.bar(models, values)
    
    # Color best model differently
    best_idx = models.index(results['best_model_type'])
    bars[best_idx].set_color('green')
    bars[best_idx].set_alpha(0.8)
    
    # Add value labels
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{val:.3f}', ha='center', va='bottom')
    
    ax.set_title(f'{label}')
    ax.set_ylabel(label)
    ax.set_xticklabels(models, rotation=45)

# Add overall title
fig.suptitle('ML Fusion Layer Training Results (FIXED - No Data Leakage)\nPredicting ClinVar Pathogenicity from 4 Genomic Features', fontsize=16)

# Add description
plt.figtext(0.5, 0.02, 
           f'Best Model: {results["best_model_type"]} | '
           f'Training Samples: {results["data_stats"]["train_size"]:,} | '
           f'Validation: {results["data_stats"]["val_size"]:,} | '
           f'Test: {results["data_stats"]["test_size"]:,}\n'
           f'Final Test Performance - MSE: {results["final_test_results"]["mse"]:.3f}, '
           f'R²: {results["final_test_results"]["r2"]:.3f}, '
           f'Accuracy: {results["final_test_results"]["accuracy"]:.3f}, '
           f'AUC: {results["final_test_results"]["auc"]:.3f}',
           ha='center', fontsize=12)

plt.tight_layout()
plt.savefig('training_results_FIXED.png', dpi=300, bbox_inches='tight')
print("Created training_results_FIXED.png") 