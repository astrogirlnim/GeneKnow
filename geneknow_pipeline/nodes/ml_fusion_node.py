#!/usr/bin/env python3
"""
ML Fusion Node for GeneKnow Pipeline

This node integrates the ML fusion layer into the LangGraph pipeline.
It takes outputs from the 5 static models and produces a final risk assessment.
"""

import os
import sys
import json
import numpy as np
from typing import Dict, List, Any, Optional
from pathlib import Path
import logging

# Add the ml_models directory to the path
sys.path.append(str(Path(__file__).parent.parent / 'ml_models'))

from fusion_layer import FusionLayer, StaticModelInputs, FusionOutput
# Use Dict instead of GeneKnowState for compatibility
from typing import Dict

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MLFusionNode:
    """
    ML Fusion Node that combines static model outputs into final risk assessment.

    This node expects the following to be available in the pipeline state:
    - PRS scores from polygenic risk scoring
    - ClinVar classifications from variant annotation
    - CADD scores from deleteriousness prediction
    - TCGA enrichment from tumor frequency matching
    - Gene burden scores from pathway analysis
    """

    def __init__(self, model_path: str = None):
        """
        Initialize the ML Fusion Node.

        Args:
            model_path: Path to the trained fusion model. If None, uses best model.
        """
        self.fusion_layer = FusionLayer()
        self.model_path = model_path or str(Path(__file__).parent.parent / 'ml_models' / 'best_fusion_model.pkl')
        self.is_loaded = False

        # Load model if it exists
        self._load_model()

    def _load_model(self):
        """Load the trained fusion model."""
        try:
            if os.path.exists(self.model_path):
                self.fusion_layer.load_model(self.model_path)
                self.is_loaded = True
                logger.info(f"ML Fusion model loaded from {self.model_path}")
            else:
                logger.warning(f"ML Fusion model not found at {self.model_path}")
                logger.warning("Use train_fusion_layer.py to train the model first")
        except Exception as e:
            logger.error(f"Error loading fusion model: {e}")
            self.is_loaded = False

    def _extract_static_model_outputs(self, state: Dict[str, Any]) -> List[StaticModelInputs]:
        """
        Extract static model outputs from the pipeline state.

        Args:
            state: Current pipeline state

        Returns:
            List of StaticModelInputs for each variant
        """
        static_inputs = []

        # Get ML-ready variants from feature vector builder
        variants = state.get('ml_ready_variants', [])
        if not variants:
            # Fallback to filtered_variants if ml_ready_variants not available
            variants = state.get('filtered_variants', [])
            if not variants:
                logger.warning("No variants found in pipeline state")
                return static_inputs

        for variant in variants:
            # Extract PRS score (default to population average if not available)
            prs_score = variant.get('prs_score', 0.5)

            # Extract ClinVar classification
            clinvar_classification = 'not_found'
            if 'clinvar' in variant:
                clinvar_data = variant['clinvar']
                if isinstance(clinvar_data, dict):
                    clinical_significance = clinvar_data.get('clinical_significance', '').lower()
                    if 'pathogenic' in clinical_significance:
                        clinvar_classification = 'pathogenic'
                    elif 'benign' in clinical_significance:
                        clinvar_classification = 'benign'
                    elif 'uncertain' in clinical_significance or 'vus' in clinical_significance:
                        clinvar_classification = 'uncertain'
                elif isinstance(clinvar_data, str):
                    clinical_significance = clinvar_data.lower()
                    if 'pathogenic' in clinical_significance:
                        clinvar_classification = 'pathogenic'
                    elif 'benign' in clinical_significance:
                        clinvar_classification = 'benign'
                    elif 'uncertain' in clinical_significance:
                        clinvar_classification = 'uncertain'

            # Extract CADD score
            cadd_score = variant.get('cadd_score', 0.0)
            if isinstance(cadd_score, str):
                try:
                    cadd_score = float(cadd_score)
                except ValueError:
                    cadd_score = 0.0

            # Extract TCGA enrichment
            tcga_enrichment = variant.get('tcga_enrichment', 1.0)
            if isinstance(tcga_enrichment, str):
                try:
                    tcga_enrichment = float(tcga_enrichment)
                except ValueError:
                    tcga_enrichment = 1.0

            # Extract gene burden score
            gene_burden_score = variant.get('gene_burden_score', 0.0)
            if isinstance(gene_burden_score, str):
                try:
                    gene_burden_score = float(gene_burden_score)
                except ValueError:
                    gene_burden_score = 0.0

            # Create StaticModelInputs object
            static_input = StaticModelInputs(
                prs_score=float(prs_score),
                clinvar_classification=clinvar_classification,
                cadd_score=float(cadd_score),
                tcga_enrichment=float(tcga_enrichment),
                gene_burden_score=float(gene_burden_score)
            )

            static_inputs.append(static_input)

        return static_inputs

    def _calculate_aggregate_risk(self, fusion_outputs: List[FusionOutput]) -> Dict[str, Any]:
        """
        Calculate aggregate risk assessment from individual variant predictions.

        Args:
            fusion_outputs: List of fusion outputs for each variant

        Returns:
            Aggregate risk assessment
        """
        if not fusion_outputs:
            return {
                'aggregate_risk_score': 0.0,
                'risk_category': 'low',
                'confidence': 0.0,
                'variant_count': 0,
                'high_risk_variants': 0,
                'contributing_factors': {}
            }

        # Calculate aggregate metrics
        risk_scores = [output.risk_score for output in fusion_outputs]
        confidences = [output.confidence for output in fusion_outputs]

        # Aggregate risk score (weighted average, with higher weights for higher risk)
        weights = np.array(risk_scores)
        weights = weights / np.sum(weights) if np.sum(weights) > 0 else np.ones(len(weights)) / len(weights)
        aggregate_risk_score = np.average(risk_scores, weights=weights)

        # Aggregate confidence (average confidence)
        aggregate_confidence = np.mean(confidences)

        # Count high-risk variants
        high_risk_variants = sum(1 for output in fusion_outputs if output.risk_score > 0.7)

        # Determine aggregate risk category
        if aggregate_risk_score <= 0.25:
            risk_category = 'low'
        elif aggregate_risk_score <= 0.5:
            risk_category = 'moderate'
        elif aggregate_risk_score <= 0.75:
            risk_category = 'high'
        else:
            risk_category = 'very_high'

        # Aggregate contributing factors
        all_factors = {}
        for output in fusion_outputs:
            for factor, contribution in output.contributing_factors.items():
                if factor not in all_factors:
                    all_factors[factor] = []
                all_factors[factor].append(contribution)

        # Average contribution per factor
        aggregated_factors = {
            factor: np.mean(contributions)
            for factor, contributions in all_factors.items()
        }

        return {
            'aggregate_risk_score': float(aggregate_risk_score),
            'risk_category': risk_category,
            'confidence': float(aggregate_confidence),
            'variant_count': len(fusion_outputs),
            'high_risk_variants': high_risk_variants,
            'contributing_factors': aggregated_factors,
            'individual_risk_scores': risk_scores
        }

    def process(self, state: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process the ML fusion step of the pipeline.

        Args:
            state: Current pipeline state

        Returns:
            Updated pipeline state with fusion results
        """
        logger.info("🔬 Starting ML Fusion processing")

        # Check if model is loaded
        if not self.is_loaded:
            logger.error("ML Fusion model not loaded. Cannot process.")
            return {
                'ml_fusion_results': {
                    'error': 'Model not loaded',
                    'aggregate_risk_score': 0.0,
                    'risk_category': 'unknown'
                }
            }

        try:
            # Extract static model outputs from state
            static_inputs = self._extract_static_model_outputs(state)

            if not static_inputs:
                logger.warning("No static model inputs found. Skipping ML fusion.")
                return {
                    'ml_fusion_results': {
                        'error': 'No static model inputs',
                        'aggregate_risk_score': 0.0,
                        'risk_category': 'unknown'
                    }
                }

            # Run fusion layer predictions
            fusion_outputs = []
            for static_input in static_inputs:
                try:
                    fusion_output = self.fusion_layer.predict(static_input)
                    fusion_outputs.append(fusion_output)
                except Exception as e:
                    logger.error(f"Error in fusion prediction: {e}")
                    continue

            # Calculate aggregate risk
            aggregate_results = self._calculate_aggregate_risk(fusion_outputs)

            # Build fusion results
            ml_fusion_results = {
                'aggregate_risk_assessment': aggregate_results,
                'individual_predictions': [output.to_dict() for output in fusion_outputs],
                'static_model_inputs': [inputs.to_dict() for inputs in static_inputs],
                'model_path': self.model_path,
                'processing_successful': True
            }

            # Log results
            logger.info(f"ML Fusion completed successfully:")
            logger.info(f"  Processed {len(fusion_outputs)} variants")
            logger.info(f"  Aggregate risk score: {aggregate_results['aggregate_risk_score']:.3f}")
            logger.info(f"  Risk category: {aggregate_results['risk_category']}")
            logger.info(f"  High-risk variants: {aggregate_results['high_risk_variants']}")

            # Return only the keys this node updates
            return {
                'ml_fusion_results': ml_fusion_results
            }

        except Exception as e:
            logger.error(f"Error in ML Fusion processing: {e}")
            return {
                'ml_fusion_results': {
                    'error': str(e),
                    'processing_successful': False
                }
            }

# LangGraph node function
def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    LangGraph node function for ML fusion processing.

    Args:
        state: Current pipeline state

    Returns:
        Updated pipeline state with ML fusion results
    """
    fusion_node = MLFusionNode()
    return fusion_node.process(state)

# Alternative node function for custom model path
def create_ml_fusion_node(model_path: str = None):
    """
    Create a custom ML fusion node with specific model path.

    Args:
        model_path: Path to the trained fusion model

    Returns:
        Configured ML fusion node function
    """
    def custom_ml_fusion_node(state: Dict[str, Any]) -> Dict[str, Any]:
        fusion_node = MLFusionNode(model_path=model_path)
        return fusion_node.process(state)

    return custom_ml_fusion_node

if __name__ == "__main__":
    # Test the fusion node with mock data
    print("🧪 Testing ML Fusion Node")

    # Create mock state
    mock_state = {
        'variants': [
            {
                'prs_score': 0.8,
                'clinvar': {'clinical_significance': 'Pathogenic'},
                'cadd_score': 25.0,
                'tcga_enrichment': 3.0,
                'gene_burden_score': 2.0
            },
            {
                'prs_score': 0.2,
                'clinvar': {'clinical_significance': 'Benign'},
                'cadd_score': 5.0,
                'tcga_enrichment': 0.5,
                'gene_burden_score': 0.0
            }
        ]
    }

    # Process with fusion node
    result_state = process(mock_state)

    # Print results
    if 'ml_fusion_results' in result_state:
        fusion_results = result_state['ml_fusion_results']
        if fusion_results.get('processing_successful'):
            aggregate = fusion_results['aggregate_risk_assessment']
            print(f"✅ Fusion processing successful")
            print(f"Aggregate risk score: {aggregate['aggregate_risk_score']:.3f}")
            print(f"Risk category: {aggregate['risk_category']}")
            print(f"Variants processed: {aggregate['variant_count']}")
        else:
            print(f"❌ Fusion processing failed: {fusion_results.get('error')}")
    else:
        print("❌ No fusion results found in state")
