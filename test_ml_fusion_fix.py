#!/usr/bin/env python3
"""Test script to verify ML fusion model is properly loaded and used."""

import os
import sys
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(name)s - %(message)s')

# Add geneknow_pipeline to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'geneknow_pipeline'))

# Direct imports to avoid module issues
import importlib.util

# Load ml_fusion_node directly
ml_fusion_spec = importlib.util.spec_from_file_location(
    "ml_fusion_node", 
    "geneknow_pipeline/nodes/ml_fusion_node.py"
)
ml_fusion_module = importlib.util.module_from_spec(ml_fusion_spec)
ml_fusion_spec.loader.exec_module(ml_fusion_module)
MLFusionNode = ml_fusion_module.MLFusionNode

# Load risk_model directly
risk_model_spec = importlib.util.spec_from_file_location(
    "risk_model",
    "geneknow_pipeline/nodes/risk_model.py"
)
risk_model_module = importlib.util.module_from_spec(risk_model_spec)
risk_model_spec.loader.exec_module(risk_model_module)
risk_model_process = risk_model_module.process

def test_ml_fusion_loading():
    """Test if ML fusion model loads correctly."""
    print("\n=== Testing ML Fusion Model Loading ===")
    
    # Test default model loading
    fusion_node = MLFusionNode()
    print(f"Model loaded: {fusion_node.is_loaded}")
    print(f"Model path: {fusion_node.model_path}")
    
    # Test with mock data
    mock_state = {
        'ml_ready_variants': [
            {
                'gene': 'BRCA1',
                'prs_score': 0.8,
                'clinvar': {'clinical_significance': 'Pathogenic'},
                'cadd_score': 25.0,
                'tcga_enrichment': 3.0,
                'gene_burden_score': 2.0
            },
            {
                'gene': 'TP53',
                'prs_score': 0.2,
                'clinvar': {'clinical_significance': 'Benign'},
                'cadd_score': 5.0,
                'tcga_enrichment': 0.5,
                'gene_burden_score': 0.0
            }
        ]
    }
    
    # Process with fusion node
    result = fusion_node.process(mock_state)
    
    print("\nML Fusion Results:")
    if 'ml_fusion_results' in result:
        fusion_results = result['ml_fusion_results']
        print(f"Processing successful: {fusion_results.get('processing_successful', False)}")
        
        if fusion_results.get('processing_successful'):
            aggregate = fusion_results.get('aggregate_risk_assessment', {})
            print(f"Aggregate risk score: {aggregate.get('aggregate_risk_score', 0):.3f}")
            print(f"Risk category: {aggregate.get('risk_category', 'unknown')}")
            print(f"Confidence: {aggregate.get('confidence', 0):.3f}")
        else:
            print(f"Error: {fusion_results.get('error', 'Unknown')}")
    
    return result

def test_risk_model_integration():
    """Test if risk model properly uses ML fusion results."""
    print("\n\n=== Testing Risk Model Integration ===")
    
    # Create state with ML fusion results
    ml_fusion_results = test_ml_fusion_loading()
    
    state_with_fusion = {
        'filtered_variants': [
            {
                'gene': 'BRCA1',
                'clinical_significance': 'Pathogenic',
                'risk_weight': 0.9
            }
        ],
        'ml_fusion_results': ml_fusion_results.get('ml_fusion_results', {}),
        'ml_ready_variants': [
            {
                'gene': 'BRCA1',
                'prs_score': 0.8,
                'clinvar': {'clinical_significance': 'Pathogenic'},
                'cadd_score': 25.0,
                'tcga_enrichment': 3.0,
                'gene_burden_score': 2.0
            }
        ]
    }
    
    # Process with risk model
    risk_results = risk_model_process(state_with_fusion)
    
    print("\nRisk Model Results:")
    print(f"Risk scores: {risk_results.get('risk_scores', {})}")
    
    # Check if ML fusion was used
    if 'ml_risk_assessment' in risk_results:
        ml_assessment = risk_results['ml_risk_assessment']
        print(f"\n✅ ML Fusion was used!")
        print(f"Method: {ml_assessment.get('method', 'unknown')}")
        print(f"Risk category: {ml_assessment.get('risk_category', 'unknown')}")
        print(f"Confidence: {ml_assessment.get('confidence', 0):.3f}")
    else:
        print("\n⚠️ WARNING: ML Fusion was NOT used - fell back to simple calculation!")

if __name__ == "__main__":
    test_ml_fusion_loading()
    test_risk_model_integration() 