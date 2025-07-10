#!/usr/bin/env python3
"""
Visualize the Fusion Layer Architecture

Shows how the 5 static models feed into the ML fusion layer
to produce an aggregate risk score.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import numpy as np

def create_fusion_architecture_diagram():
    """Create a visual diagram of the fusion layer architecture."""
    
    fig, ax = plt.subplots(1, 1, figsize=(14, 10))
    
    # Define positions
    static_models = [
        {'name': 'PRS Score', 'desc': 'Polygenic Risk\n(Inherited background)', 'color': '#FF6B6B', 'y': 8},
        {'name': 'ClinVar', 'desc': 'Known Pathogenic\nVariants Database', 'color': '#4ECDC4', 'y': 6.5},
        {'name': 'CADD Score', 'desc': 'Deleteriousness\nPrediction', 'color': '#45B7D1', 'y': 5},
        {'name': 'TCGA Enrichment', 'desc': 'Tumor Frequency\nAnalysis', 'color': '#F7DC6F', 'y': 3.5},
        {'name': 'Gene Burden', 'desc': 'Pathway-level\nDamage Score', 'color': '#BB8FCE', 'y': 2}
    ]
    
    # Draw static model boxes
    for model in static_models:
        # Main box
        box = FancyBboxPatch((0.5, model['y']-0.4), 3, 0.8,
                            boxstyle="round,pad=0.1",
                            facecolor=model['color'],
                            edgecolor='black',
                            linewidth=2,
                            alpha=0.8)
        ax.add_patch(box)
        
        # Model name
        ax.text(2, model['y'], model['name'], 
               fontsize=12, fontweight='bold',
               ha='center', va='center')
        
        # Description
        ax.text(2, model['y']-0.6, model['desc'], 
               fontsize=9, style='italic',
               ha='center', va='top', alpha=0.7)
        
        # Output value example
        if model['name'] == 'PRS Score':
            output = '0.2-0.8'
        elif model['name'] == 'ClinVar':
            output = 'pathogenic/benign'
        elif model['name'] == 'CADD Score':
            output = '0-50'
        elif model['name'] == 'TCGA Enrichment':
            output = '0.1-20x'
        else:
            output = '0-10'
        
        ax.text(4, model['y'], f'→ {output}', 
               fontsize=10, color='darkgreen',
               ha='left', va='center')
    
    # Feature Engineering Box
    feature_box = FancyBboxPatch((5.5, 3), 3, 3,
                                boxstyle="round,pad=0.1",
                                facecolor='lightgray',
                                edgecolor='black',
                                linewidth=2,
                                alpha=0.8)
    ax.add_patch(feature_box)
    
    ax.text(7, 5.5, 'Feature\nEngineering', 
           fontsize=13, fontweight='bold',
           ha='center', va='center')
    
    ax.text(7, 4.5, '• Normalize scores\n• One-hot encode\n• Scale features', 
           fontsize=10,
           ha='center', va='center')
    
    # ML Fusion Layer Box
    fusion_box = FancyBboxPatch((9.5, 2.5), 3.5, 4,
                               boxstyle="round,pad=0.1",
                               facecolor='#2ECC71',
                               edgecolor='black',
                               linewidth=3,
                               alpha=0.9)
    ax.add_patch(fusion_box)
    
    ax.text(11.25, 5.5, 'ML Fusion Layer', 
           fontsize=14, fontweight='bold',
           ha='center', va='center')
    
    ax.text(11.25, 4.5, 'Gradient Boosting\n(Trained on 200K variants)', 
           fontsize=11, style='italic',
           ha='center', va='center')
    
    ax.text(11.25, 3.5, 'Learns optimal\nweights & interactions', 
           fontsize=10,
           ha='center', va='center')
    
    # Draw arrows from static models to feature engineering
    for model in static_models:
        arrow = FancyArrowPatch((3.5, model['y']), (5.5, 4.5),
                               connectionstyle="arc3,rad=0.3",
                               arrowstyle="->",
                               mutation_scale=20,
                               linewidth=2,
                               color='gray',
                               alpha=0.6)
        ax.add_artist(arrow)
    
    # Arrow from feature engineering to fusion
    arrow = FancyArrowPatch((8.5, 4.5), (9.5, 4.5),
                           arrowstyle="->",
                           mutation_scale=25,
                           linewidth=3,
                           color='darkgreen')
    ax.add_artist(arrow)
    
    # Output boxes
    outputs = [
        {'name': 'Risk Score', 'value': '0.0 - 1.0', 'y': 5.5, 'color': '#E74C3C'},
        {'name': 'Confidence', 'value': '0.8 (fixed)', 'y': 4.5, 'color': '#3498DB'},
        {'name': 'Risk Category', 'value': 'low/med/high', 'y': 3.5, 'color': '#9B59B6'}
    ]
    
    for output in outputs:
        # Output arrow
        arrow = FancyArrowPatch((13, 4.5), (14, output['y']),
                               connectionstyle="arc3,rad=0.2",
                               arrowstyle="->",
                               mutation_scale=20,
                               linewidth=2,
                               color=output['color'])
        ax.add_artist(arrow)
        
        # Output box
        out_box = FancyBboxPatch((14, output['y']-0.3), 2.5, 0.6,
                                boxstyle="round,pad=0.05",
                                facecolor=output['color'],
                                alpha=0.3,
                                edgecolor=output['color'],
                                linewidth=2)
        ax.add_patch(out_box)
        
        ax.text(15.25, output['y'], f"{output['name']}\n{output['value']}", 
               fontsize=10,
               ha='center', va='center')
    
    # Add performance metrics
    perf_text = """Model Performance:
• AUC: 0.76
• R²: 0.43
• Comparable to CADD/PolyPhen-2"""
    
    ax.text(11.25, 1.5, perf_text, 
           fontsize=9,
           ha='center', va='center',
           bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.3))
    
    # Title
    ax.text(8.5, 9, 'GeneKnow ML Fusion Layer Architecture', 
           fontsize=18, fontweight='bold',
           ha='center', va='center')
    
    ax.text(8.5, 8.3, 'Aggregating 5 Static Models → Single Risk Score', 
           fontsize=14, style='italic',
           ha='center', va='center', alpha=0.7)
    
    # Input data flow
    ax.text(0.5, 9, 'Patient\nVariant', 
           fontsize=12, fontweight='bold',
           ha='center', va='center',
           bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.5))
    
    # Set axis limits and remove axes
    ax.set_xlim(0, 17)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    plt.tight_layout()
    return fig

def create_example_calculation():
    """Show a concrete example of how fusion works."""
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Example 1: High-risk variant
    ax1.text(0.5, 0.9, 'Example 1: High-Risk Variant', 
            fontsize=16, fontweight='bold',
            transform=ax1.transAxes)
    
    # Input values
    inputs = [
        ('PRS Score:', '0.8', 'High genetic background risk'),
        ('ClinVar:', 'Not Found', 'No prior classification'),
        ('CADD Score:', '30', 'Highly deleterious'),
        ('TCGA Enrichment:', '5x', 'Common in tumors'),
        ('Gene Burden:', '8', 'Many damaged genes')
    ]
    
    y_pos = 0.7
    for name, value, desc in inputs:
        ax1.text(0.1, y_pos, name, fontsize=12, fontweight='bold', transform=ax1.transAxes)
        ax1.text(0.3, y_pos, value, fontsize=12, color='red', transform=ax1.transAxes)
        ax1.text(0.45, y_pos, f'({desc})', fontsize=10, style='italic', alpha=0.7, transform=ax1.transAxes)
        y_pos -= 0.12
    
    # Arrow
    ax1.annotate('', xy=(0.85, 0.4), xytext=(0.75, 0.4),
                arrowprops=dict(arrowstyle='->', lw=3, color='green'),
                transform=ax1.transAxes)
    
    # Output
    ax1.text(0.9, 0.5, 'ML Fusion\nRisk Score:\n0.661', 
            fontsize=14, fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", facecolor='red', alpha=0.3),
            transform=ax1.transAxes)
    
    ax1.text(0.9, 0.2, 'Category:\nHIGH RISK', 
            fontsize=12, ha='center', va='center', color='red',
            transform=ax1.transAxes)
    
    ax1.axis('off')
    
    # Example 2: Low-risk variant
    ax2.text(0.5, 0.9, 'Example 2: Low-Risk Variant', 
            fontsize=16, fontweight='bold',
            transform=ax2.transAxes)
    
    # Input values
    inputs = [
        ('PRS Score:', '0.2', 'Low genetic risk'),
        ('ClinVar:', 'Not Found', 'No prior classification'),
        ('CADD Score:', '5', 'Likely benign'),
        ('TCGA Enrichment:', '0.5x', 'Rare in tumors'),
        ('Gene Burden:', '1', 'Minimal damage')
    ]
    
    y_pos = 0.7
    for name, value, desc in inputs:
        ax2.text(0.1, y_pos, name, fontsize=12, fontweight='bold', transform=ax2.transAxes)
        ax2.text(0.3, y_pos, value, fontsize=12, color='green', transform=ax2.transAxes)
        ax2.text(0.45, y_pos, f'({desc})', fontsize=10, style='italic', alpha=0.7, transform=ax2.transAxes)
        y_pos -= 0.12
    
    # Arrow
    ax2.annotate('', xy=(0.85, 0.4), xytext=(0.75, 0.4),
                arrowprops=dict(arrowstyle='->', lw=3, color='green'),
                transform=ax2.transAxes)
    
    # Output
    ax2.text(0.9, 0.5, 'ML Fusion\nRisk Score:\n0.000', 
            fontsize=14, fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", facecolor='green', alpha=0.3),
            transform=ax2.transAxes)
    
    ax2.text(0.9, 0.2, 'Category:\nLOW RISK', 
            fontsize=12, ha='center', va='center', color='green',
            transform=ax2.transAxes)
    
    ax2.axis('off')
    
    plt.tight_layout()
    return fig

def main():
    """Create and save fusion architecture visualizations."""
    
    print("📊 Creating Fusion Architecture Visualization...")
    
    # Create architecture diagram
    fig1 = create_fusion_architecture_diagram()
    fig1.savefig('fusion_architecture.png', dpi=300, bbox_inches='tight')
    print("✅ Saved fusion_architecture.png")
    
    # Create example calculation
    fig2 = create_example_calculation()
    fig2.savefig('fusion_examples.png', dpi=300, bbox_inches='tight')
    print("✅ Saved fusion_examples.png")
    
    print("\n🎯 Key Points:")
    print("1. The fusion layer takes outputs from 5 independent static models")
    print("2. It learns optimal weights through training on 200K variants")
    print("3. Output is a single risk score (0-1) with confidence")
    print("4. Model achieves AUC=0.76, comparable to published methods")

if __name__ == "__main__":
    main() 