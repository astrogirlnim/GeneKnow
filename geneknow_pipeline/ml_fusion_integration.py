#!/usr/bin/env python3
"""
ML Fusion Integration for GeneKnow Pipeline
Integrates trained ML models into the existing LangGraph pipeline.
"""

import logging
import os
import numpy as np
<<<<<<< HEAD
from typing import Dict, List, Any
=======
from typing import Dict, List, Any, Optional
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
from datetime import datetime

from ml_trainer import GenomicMLTrainer
from ml_feature_extractor import GenomicFeatureExtractor

logger = logging.getLogger(__name__)

<<<<<<< HEAD

=======
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
class MLFusionIntegrator:
    """
    Integrates ML fusion models into the GeneKnow pipeline.
    Replaces the stub feature_vector_builder with real ML predictions.
    """
<<<<<<< HEAD

=======
    
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    def __init__(self, model_dir: str = "ml_models"):
        self.model_dir = model_dir
        self.trainer = GenomicMLTrainer()
        self.feature_extractor = GenomicFeatureExtractor()
        self.models_loaded = False
<<<<<<< HEAD

    def load_models(self) -> bool:
        """Load trained ML models."""
        try:
            best_model_path = os.path.join(self.model_dir, "best_model.pkl")
            if os.path.exists(best_model_path):
                self.trainer.load_model(best_model_path)
                self.models_loaded = True
                logger.info(
                    f"ML models loaded successfully: {self.trainer.best_model_name}"
                )
=======
        
    def load_models(self) -> bool:
        """Load trained ML models."""
        try:
            best_model_path = os.path.join(self.model_dir, 'best_model.pkl')
            if os.path.exists(best_model_path):
                self.trainer.load_model(best_model_path)
                self.models_loaded = True
                logger.info(f"ML models loaded successfully: {self.trainer.best_model_name}")
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
                return True
            else:
                logger.warning(f"No trained models found at {best_model_path}")
                return False
        except Exception as e:
            logger.error(f"Failed to load ML models: {e}")
            return False
<<<<<<< HEAD

    def extract_pipeline_features(
        self, variants: List[Dict], state: Dict[str, Any]
    ) -> List[Dict]:
=======
    
    def extract_pipeline_features(self, variants: List[Dict], state: Dict[str, Any]) -> List[Dict]:
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        """
        Extract features from pipeline variants, incorporating all available annotations.
        """
        enriched_variants = []
<<<<<<< HEAD

        # Get additional data from pipeline state
        population_matches = state.get("population_matches", {})
        tcga_matches = state.get("tcga_matches", {})

        for variant in variants:
            # Start with base variant data
            enriched_variant = variant.copy()

=======
        
        # Get additional data from pipeline state
        population_matches = state.get("population_matches", {})
        cadd_stats = state.get("cadd_stats", {})
        tcga_matches = state.get("tcga_matches", {})
        
        for variant in variants:
            # Start with base variant data
            enriched_variant = variant.copy()
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            # Add population frequency data
            variant_id = variant.get("variant_id", "")
            if variant_id in population_matches:
                pop_data = population_matches[variant_id]
<<<<<<< HEAD
                enriched_variant.update(
                    {
                        "population_frequency": pop_data.get(
                            "population_frequency", 0.0
                        ),
                        "clinical_significance": pop_data.get(
                            "clinical_significance", "Unknown"
                        ),
                        "is_pathogenic": pop_data.get("is_pathogenic", 0),
                        "review_status": pop_data.get("review_status", ""),
                    }
                )

            # Add CADD scores if available
            if "cadd_phred" in variant or "cadd_raw" in variant:
                enriched_variant.update(
                    {
                        "cadd_phred": variant.get("cadd_phred", 0.0),
                        "cadd_raw": variant.get("cadd_raw", 0.0),
                    }
                )

=======
                enriched_variant.update({
                    "population_frequency": pop_data.get("population_frequency", 0.0),
                    "clinical_significance": pop_data.get("clinical_significance", "Unknown"),
                    "is_pathogenic": pop_data.get("is_pathogenic", 0),
                    "review_status": pop_data.get("review_status", "")
                })
            
            # Add CADD scores if available
            if "cadd_phred" in variant or "cadd_raw" in variant:
                enriched_variant.update({
                    "cadd_phred": variant.get("cadd_phred", 0.0),
                    "cadd_raw": variant.get("cadd_raw", 0.0)
                })
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            # Add TCGA enrichment data
            tcga_enrichment = 0.0
            for cancer_type, matches in tcga_matches.items():
                if variant_id in matches:
                    match_data = matches[variant_id]
                    if "frequency" in match_data:
                        tcga_enrichment = max(tcga_enrichment, match_data["frequency"])
<<<<<<< HEAD

            enriched_variant["tcga_enrichment_score"] = tcga_enrichment

            # Ensure required fields have defaults
            enriched_variant.setdefault(
                "gnomad_a", enriched_variant.get("population_frequency", 0.0)
            )
            enriched_variant.setdefault(
                "consequence", enriched_variant.get("variant_classification", "unknown")
            )
            enriched_variant.setdefault("gene", "Unknown")

            enriched_variants.append(enriched_variant)

        return enriched_variants

    def predict_pathogenicity(
        self, variants: List[Dict], state: Dict[str, Any]
    ) -> Dict[str, Any]:
=======
            
            enriched_variant["tcga_enrichment_score"] = tcga_enrichment
            
            # Ensure required fields have defaults
            enriched_variant.setdefault("gnomad_af", enriched_variant.get("population_frequency", 0.0))
            enriched_variant.setdefault("consequence", enriched_variant.get("variant_classification", "unknown"))
            enriched_variant.setdefault("gene", "Unknown")
            
            enriched_variants.append(enriched_variant)
        
        return enriched_variants
    
    def predict_pathogenicity(self, variants: List[Dict], state: Dict[str, Any]) -> Dict[str, Any]:
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        """
        Use ML models to predict pathogenicity for variants.
        """
        if not self.models_loaded:
            if not self.load_models():
                logger.warning("ML models not available. Using fallback scoring.")
                return self._fallback_scoring(variants, state)
<<<<<<< HEAD

        try:
            # Extract and enrich features
            enriched_variants = self.extract_pipeline_features(variants, state)

            # Get ML predictions
            ml_results = self.trainer.predict_risk(enriched_variants)

=======
        
        try:
            # Extract and enrich features
            enriched_variants = self.extract_pipeline_features(variants, state)
            
            # Get ML predictions
            ml_results = self.trainer.predict_risk(enriched_variants)
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            # Process results
            ml_predictions = {
                "pathogenicity_scores": ml_results["risk_scores"],
                "pathogenicity_predictions": ml_results["predictions"],
                "model_used": ml_results["model_used"],
<<<<<<< HEAD
                "prediction_confidence": self._calculate_confidence(
                    ml_results["probabilities"]
                ),
            }

=======
                "prediction_confidence": self._calculate_confidence(ml_results["probabilities"])
            }
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            # Update variants with ML scores
            for i, variant in enumerate(variants):
                variant["ml_pathogenicity_score"] = ml_results["risk_scores"][i]
                variant["ml_pathogenicity_prediction"] = ml_results["predictions"][i]
                variant["ml_confidence"] = ml_predictions["prediction_confidence"][i]
<<<<<<< HEAD

            logger.info(f"ML predictions completed for {len(variants)} variants")
            logger.info(f"Model used: {ml_results['model_used']}")

            return ml_predictions

        except Exception as e:
            logger.error(f"ML prediction failed: {e}")
            return self._fallback_scoring(variants, state)

=======
            
            logger.info(f"ML predictions completed for {len(variants)} variants")
            logger.info(f"Model used: {ml_results['model_used']}")
            
            return ml_predictions
            
        except Exception as e:
            logger.error(f"ML prediction failed: {e}")
            return self._fallback_scoring(variants, state)
    
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    def _calculate_confidence(self, probabilities: List[List[float]]) -> List[float]:
        """Calculate prediction confidence based on probability distributions."""
        confidences = []
        for prob_pair in probabilities:
            # Confidence is the maximum probability minus 0.5
            # Higher confidence when probability is closer to 0 or 1
            max_prob = max(prob_pair)
<<<<<<< HEAD
            confidence = abs(max_prob - 0.5) * 2  # Scale to 0 - 1
            confidences.append(confidence)
        return confidences

    def _fallback_scoring(
        self, variants: List[Dict], state: Dict[str, Any]
    ) -> Dict[str, Any]:
=======
            confidence = abs(max_prob - 0.5) * 2  # Scale to 0-1
            confidences.append(confidence)
        return confidences
    
    def _fallback_scoring(self, variants: List[Dict], state: Dict[str, Any]) -> Dict[str, Any]:
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        """
        Fallback scoring when ML models are not available.
        Uses rule-based approach combining available annotations.
        """
        logger.info("Using fallback rule-based scoring")
<<<<<<< HEAD

        pathogenicity_scores = []
        pathogenicity_predictions = []
        confidences = []

        for variant in variants:
            score = 0.0
            confidence = 0.3  # Low confidence for rule-based

=======
        
        pathogenicity_scores = []
        pathogenicity_predictions = []
        confidences = []
        
        for variant in variants:
            score = 0.0
            confidence = 0.3  # Low confidence for rule-based
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            # Clinical significance
            clinical_sig = variant.get("clinical_significance", "").lower()
            if "pathogenic" in clinical_sig and "likely" not in clinical_sig:
                score += 0.8
                confidence = 0.9
            elif "likely_pathogenic" in clinical_sig:
                score += 0.6
                confidence = 0.7
            elif "benign" in clinical_sig:
                score = 0.1
                confidence = 0.8
<<<<<<< HEAD

            # Population frequency
            pop_freq = variant.get("population_frequency", variant.get("gnomad_a", 0.0))
=======
            
            # Population frequency
            pop_freq = variant.get("population_frequency", variant.get("gnomad_af", 0.0))
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            if pop_freq > 0.01:  # Common variant
                score *= 0.3  # Reduce score for common variants
            elif pop_freq < 0.001:  # Rare variant
                score += 0.2
<<<<<<< HEAD

=======
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            # CADD score
            cadd_phred = variant.get("cadd_phred", 0.0)
            if cadd_phred > 20:
                score += 0.3
            elif cadd_phred > 10:
                score += 0.1
<<<<<<< HEAD

            # Consequence severity
            consequence = variant.get("consequence", "").lower()
            if any(
                term in consequence
                for term in ["nonsense", "frameshift", "stop_gained"]
            ):
=======
            
            # Consequence severity
            consequence = variant.get("consequence", "").lower()
            if any(term in consequence for term in ["nonsense", "frameshift", "stop_gained"]):
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
                score += 0.4
            elif "missense" in consequence:
                score += 0.2
            elif "synonymous" in consequence:
                score *= 0.5
<<<<<<< HEAD

=======
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            # TCGA enrichment
            tcga_enrichment = variant.get("tcga_enrichment_score", 0.0)
            if tcga_enrichment > 2.0:
                score += 0.2
<<<<<<< HEAD

            # Normalize score
            score = min(score, 1.0)

            pathogenicity_scores.append(score)
            pathogenicity_predictions.append(1 if score > 0.5 else 0)
            confidences.append(confidence)

=======
            
            # Normalize score
            score = min(score, 1.0)
            
            pathogenicity_scores.append(score)
            pathogenicity_predictions.append(1 if score > 0.5 else 0)
            confidences.append(confidence)
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        return {
            "pathogenicity_scores": pathogenicity_scores,
            "pathogenicity_predictions": pathogenicity_predictions,
            "model_used": "rule_based_fallback",
<<<<<<< HEAD
            "prediction_confidence": confidences,
        }

    def calculate_cancer_risk(
        self,
        variants: List[Dict],
        patient_data: Dict[str, Any],
        ml_predictions: Dict[str, Any],
    ) -> Dict[str, float]:
=======
            "prediction_confidence": confidences
        }
    
    def calculate_cancer_risk(self, variants: List[Dict], patient_data: Dict[str, Any], 
                            ml_predictions: Dict[str, Any]) -> Dict[str, float]:
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        """
        Calculate cancer-specific risk scores using ML pathogenicity predictions.
        """
        # Define cancer-specific gene sets
        cancer_gene_sets = {
<<<<<<< HEAD
            "breast": [
                "BRCA1",
                "BRCA2",
                "PALB2",
                "ATM",
                "CHEK2",
                "TP53",
                "PTEN",
                "CDH1",
            ],
            "colon": [
                "APC",
                "KRAS",
                "TP53",
                "MLH1",
                "MSH2",
                "MSH6",
                "PMS2",
                "SMAD4",
                "BRAF",
            ],
            "lung": ["TP53", "KRAS", "EGFR", "STK11", "KEAP1", "NF1", "RB1", "CDKN2A"],
            "prostate": ["AR", "PTEN", "TP53", "BRCA1", "BRCA2", "ATM", "HOXB13"],
            "blood": ["JAK2", "FLT3", "NPM1", "DNMT3A", "TET2", "IDH1", "IDH2", "TP53"],
        }

=======
            "breast": ["BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "TP53", "PTEN", "CDH1"],
            "colon": ["APC", "KRAS", "TP53", "MLH1", "MSH2", "MSH6", "PMS2", "SMAD4", "BRAF"],
            "lung": ["TP53", "KRAS", "EGFR", "STK11", "KEAP1", "NF1", "RB1", "CDKN2A"],
            "prostate": ["AR", "PTEN", "TP53", "BRCA1", "BRCA2", "ATM", "HOXB13"],
            "blood": ["JAK2", "FLT3", "NPM1", "DNMT3A", "TET2", "IDH1", "IDH2", "TP53"]
        }
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        # Base population risks (lifetime risk percentages)
        base_risks = {
            "breast": 12.9,
            "colon": 4.3,
            "lung": 6.3,
            "prostate": 12.5,
<<<<<<< HEAD
            "blood": 1.8,
        }

        cancer_risks = {}
        pathogenicity_scores = ml_predictions["pathogenicity_scores"]

        for cancer_type, cancer_genes in cancer_gene_sets.items():
            base_risk = base_risks.get(cancer_type, 2.0)
            risk_multiplier = 1.0

            # Count pathogenic variants in cancer-specific genes
            pathogenic_variants = 0
            total_risk_score = 0.0

=======
            "blood": 1.8
        }
        
        cancer_risks = {}
        pathogenicity_scores = ml_predictions["pathogenicity_scores"]
        
        for cancer_type, cancer_genes in cancer_gene_sets.items():
            base_risk = base_risks.get(cancer_type, 2.0)
            risk_multiplier = 1.0
            
            # Count pathogenic variants in cancer-specific genes
            pathogenic_variants = 0
            total_risk_score = 0.0
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            for i, variant in enumerate(variants):
                gene = variant.get("gene", "Unknown")
                if gene in cancer_genes:
                    pathogenicity_score = pathogenicity_scores[i]
                    confidence = ml_predictions["prediction_confidence"][i]
<<<<<<< HEAD

                    # Weight by confidence
                    weighted_score = pathogenicity_score * confidence
                    total_risk_score += weighted_score

                    if pathogenicity_score > 0.5:
                        pathogenic_variants += 1

=======
                    
                    # Weight by confidence
                    weighted_score = pathogenicity_score * confidence
                    total_risk_score += weighted_score
                    
                    if pathogenicity_score > 0.5:
                        pathogenic_variants += 1
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            # Calculate risk multiplier
            if pathogenic_variants > 0:
                # High-confidence pathogenic variants increase risk significantly
                risk_multiplier = 1.0 + (total_risk_score * 2.0)
<<<<<<< HEAD

                # Cap the risk increase
                risk_multiplier = min(risk_multiplier, 10.0)

            # Apply age and sex adjustments
            age = patient_data.get("age", 50)
            sex = patient_data.get("sex", "F")

            # Age adjustment (risk generally increases with age)
            age_factor = 1.0 + (age - 50) * 0.01
            age_factor = max(age_factor, 0.5)

=======
                
                # Cap the risk increase
                risk_multiplier = min(risk_multiplier, 10.0)
            
            # Apply age and sex adjustments
            age = patient_data.get("age", 50)
            sex = patient_data.get("sex", "F")
            
            # Age adjustment (risk generally increases with age)
            age_factor = 1.0 + (age - 50) * 0.01
            age_factor = max(age_factor, 0.5)
            
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            # Sex-specific adjustments
            if cancer_type == "breast" and sex == "M":
                risk_multiplier *= 0.01  # Much lower risk for men
            elif cancer_type == "prostate" and sex == "F":
                risk_multiplier = 0.0  # No risk for women
<<<<<<< HEAD

            # Calculate final risk
            final_risk = base_risk * risk_multiplier * age_factor
            final_risk = min(final_risk, 95.0)  # Cap at 95%

            cancer_risks[cancer_type] = final_risk

=======
            
            # Calculate final risk
            final_risk = base_risk * risk_multiplier * age_factor
            final_risk = min(final_risk, 95.0)  # Cap at 95%
            
            cancer_risks[cancer_type] = final_risk
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        return cancer_risks


def update_feature_vector_builder(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Enhanced feature vector builder that uses ML fusion.
    This replaces the stub in nodes/feature_vector_builder.py
    """
    logger.info("Starting ML-powered feature vector builder")
    state["current_node"] = "feature_vector_builder"
<<<<<<< HEAD

    try:
        # Initialize ML fusion integrator
        integrator = MLFusionIntegrator()

        # Get variants from pipeline
        filtered_variants = state.get("filtered_variants", [])
        patient_data = state.get("patient_data", {})

=======
    
    try:
        # Initialize ML fusion integrator
        integrator = MLFusionIntegrator()
        
        # Get variants from pipeline
        filtered_variants = state.get("filtered_variants", [])
        patient_data = state.get("patient_data", {})
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        if not filtered_variants:
            logger.warning("No variants to process")
            state["feature_vector"] = {"status": "no_variants"}
            return state
<<<<<<< HEAD

        # Get ML pathogenicity predictions
        ml_predictions = integrator.predict_pathogenicity(filtered_variants, state)

        # Calculate cancer-specific risks
        cancer_risks = integrator.calculate_cancer_risk(
            filtered_variants, patient_data, ml_predictions
        )

=======
        
        # Get ML pathogenicity predictions
        ml_predictions = integrator.predict_pathogenicity(filtered_variants, state)
        
        # Calculate cancer-specific risks
        cancer_risks = integrator.calculate_cancer_risk(filtered_variants, patient_data, ml_predictions)
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        # Build comprehensive feature vector
        feature_vector = {
            "status": "ml_fusion_complete",
            "ml_model_used": ml_predictions["model_used"],
            "total_variants": len(filtered_variants),
            "pathogenic_variants": sum(ml_predictions["pathogenicity_predictions"]),
<<<<<<< HEAD
            "average_pathogenicity_score": np.mean(
                ml_predictions["pathogenicity_scores"]
            ),
            "high_confidence_predictions": sum(
                1 for c in ml_predictions["prediction_confidence"] if c > 0.8
            ),
            "cancer_risks": cancer_risks,
            "ml_predictions": ml_predictions,
            "processing_timestamp": datetime.now().isoformat(),
        }

=======
            "average_pathogenicity_score": np.mean(ml_predictions["pathogenicity_scores"]),
            "high_confidence_predictions": sum(1 for c in ml_predictions["prediction_confidence"] if c > 0.8),
            "cancer_risks": cancer_risks,
            "ml_predictions": ml_predictions,
            "processing_timestamp": datetime.now().isoformat()
        }
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        # Update state
        state["feature_vector"] = feature_vector
        state["ml_cancer_risks"] = cancer_risks
        state["ml_pathogenicity_predictions"] = ml_predictions
<<<<<<< HEAD

        # Log summary
        logger.info("ML fusion completed:")
        logger.info(f"  Model used: {ml_predictions['model_used']}")
        logger.info(f"  Total variants: {len(filtered_variants)}")
        logger.info(
            f"  Predicted pathogenic: {sum(ml_predictions['pathogenicity_predictions'])}"
        )
        logger.info(
            f"  Average pathogenicity score: {np.mean(ml_predictions['pathogenicity_scores']):.3f}"
        )
        logger.info(f"  Cancer risks: {cancer_risks}")

        state["completed_nodes"].append("feature_vector_builder")

    except Exception as e:
        logger.error(f"ML fusion feature vector builder failed: {str(e)}")
        state["errors"].append(
            {
                "node": "feature_vector_builder",
                "error": str(e),
                "timestamp": datetime.now(),
            }
        )

=======
        
        # Log summary
        logger.info(f"ML fusion completed:")
        logger.info(f"  Model used: {ml_predictions['model_used']}")
        logger.info(f"  Total variants: {len(filtered_variants)}")
        logger.info(f"  Predicted pathogenic: {sum(ml_predictions['pathogenicity_predictions'])}")
        logger.info(f"  Average pathogenicity score: {np.mean(ml_predictions['pathogenicity_scores']):.3f}")
        logger.info(f"  Cancer risks: {cancer_risks}")
        
        state["completed_nodes"].append("feature_vector_builder")
        
    except Exception as e:
        logger.error(f"ML fusion feature vector builder failed: {str(e)}")
        state["errors"].append({
            "node": "feature_vector_builder",
            "error": str(e),
            "timestamp": datetime.now()
        })
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        # Fallback to basic feature vector
        state["feature_vector"] = {
            "status": "fallback_mode",
            "error": str(e),
            "total_variants": len(state.get("filtered_variants", [])),
<<<<<<< HEAD
            "ml_available": False,
        }

=======
            "ml_available": False
        }
    
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    return state


def main():
    """Test the ML fusion integration."""
    logging.basicConfig(level=logging.INFO)
<<<<<<< HEAD

=======
    
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    # Test with sample data
    test_variants = [
        {
            "variant_id": "17:43044295:A>T",
            "gene": "BRCA1",
            "chrom": "17",
            "pos": 43044295,
<<<<<<< HEAD
            "re": "A",
            "alt": "T",
            "consequence": "missense_variant",
            "quality": 60.0,
            "depth": 50,
=======
            "ref": "A",
            "alt": "T",
            "consequence": "missense_variant",
            "quality": 60.0,
            "depth": 50
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        },
        {
            "variant_id": "17:7579472:G>C",
            "gene": "TP53",
            "chrom": "17",
            "pos": 7579472,
<<<<<<< HEAD
            "re": "G",
            "alt": "C",
            "consequence": "missense_variant",
            "quality": 45.0,
            "depth": 30,
        },
    ]

=======
            "ref": "G",
            "alt": "C",
            "consequence": "missense_variant",
            "quality": 45.0,
            "depth": 30
        }
    ]
    
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    test_state = {
        "filtered_variants": test_variants,
        "patient_data": {"age": 45, "sex": "F"},
        "population_matches": {},
        "tcga_matches": {},
        "cadd_stats": {},
        "errors": [],
<<<<<<< HEAD
        "completed_nodes": [],
    }

    result_state = update_feature_vector_builder(test_state)

    print("ML Fusion Test Results:")
    print(f"Feature vector status: {result_state['feature_vector']['status']}")
    if "ml_cancer_risks" in result_state:
=======
        "completed_nodes": []
    }
    
    # Test ML fusion
    integrator = MLFusionIntegrator()
    result_state = update_feature_vector_builder(test_state)
    
    print("ML Fusion Test Results:")
    print(f"Feature vector status: {result_state['feature_vector']['status']}")
    if 'ml_cancer_risks' in result_state:
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        print(f"Cancer risks: {result_state['ml_cancer_risks']}")


if __name__ == "__main__":
<<<<<<< HEAD
    main()
=======
    main() 
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
