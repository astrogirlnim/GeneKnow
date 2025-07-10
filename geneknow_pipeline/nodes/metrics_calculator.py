"""
Metrics calculator node for GeneKnow pipeline.
Calculates various performance and confidence metrics based on risk model outputs.
Implements metrics from the Cancer Risk Prediction Metrics document.
"""
import logging
import numpy as np
from typing import Dict, Any, List, Optional, Tuple
from datetime import datetime
from collections import defaultdict
import math

logger = logging.getLogger(__name__)


def calculate_confidence_metrics(risk_scores: Dict[str, float], 
                               risk_details: Dict[str, Any],
                               ml_assessment: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate confidence metrics that don't require ground truth.
    These metrics assess the reliability of predictions.
    """
    confidence_metrics = {}
    
    # 1. Model Confidence Score (from risk model)
    model_confidences = []
    for cancer_type, details in risk_details.items():
        if "model_confidence" in details:
            model_confidences.append(details["model_confidence"])
    
    if model_confidences:
        confidence_metrics["mean_model_confidence"] = np.mean(model_confidences)
        confidence_metrics["min_model_confidence"] = np.min(model_confidences)
        confidence_metrics["max_model_confidence"] = np.max(model_confidences)
    
    # 2. ML Fusion Confidence
    if ml_assessment:
        confidence_metrics["ml_fusion_confidence"] = ml_assessment.get("confidence", 0.0)
        confidence_metrics["ml_fusion_risk_category"] = ml_assessment.get("risk_category", "unknown")
    
    # 3. Risk Score Distribution Metrics
    risk_values = list(risk_scores.values())
    if risk_values:
        confidence_metrics["risk_score_mean"] = np.mean(risk_values)
        confidence_metrics["risk_score_std"] = np.std(risk_values)
        confidence_metrics["risk_score_cv"] = np.std(risk_values) / np.mean(risk_values) if np.mean(risk_values) > 0 else 0
        confidence_metrics["max_risk_score"] = np.max(risk_values)
        confidence_metrics["high_risk_count"] = sum(1 for r in risk_values if r > 50)
    
    # 4. Variant Quality Metrics
    # High confidence if we have many high-quality pathogenic variants
    # Low confidence if we have mostly uncertain variants
    
    return confidence_metrics


def calculate_variant_metrics(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate metrics related to variant quality and impact.
    """
    variant_metrics = {}
    
    filtered_variants = state.get("filtered_variants", [])
    ml_ready_variants = state.get("ml_ready_variants", [])
    
    # Basic counts
    variant_metrics["total_variants"] = len(filtered_variants)
    
    # Clinical significance distribution
    clinical_sig_counts = defaultdict(int)
    pathogenic_count = 0
    benign_count = 0
    uncertain_count = 0
    
    for variant in filtered_variants:
        clinical_sig = variant.get("clinical_significance", "Unknown").lower()
        if "pathogenic" in clinical_sig:
            pathogenic_count += 1
            clinical_sig_counts["pathogenic"] += 1
        elif "benign" in clinical_sig:
            benign_count += 1
            clinical_sig_counts["benign"] += 1
        else:
            uncertain_count += 1
            clinical_sig_counts["uncertain"] += 1
    
    variant_metrics["pathogenic_variants"] = pathogenic_count
    variant_metrics["benign_variants"] = benign_count
    variant_metrics["uncertain_variants"] = uncertain_count
    variant_metrics["pathogenic_ratio"] = pathogenic_count / len(filtered_variants) if filtered_variants else 0
    
    # CADD score distribution
    cadd_scores = [v.get("cadd_phred", 0) for v in filtered_variants if "cadd_phred" in v]
    if cadd_scores:
        variant_metrics["mean_cadd_score"] = np.mean(cadd_scores)
        variant_metrics["max_cadd_score"] = np.max(cadd_scores)
        variant_metrics["high_cadd_variants"] = sum(1 for s in cadd_scores if s > 20)
    
    # Gene impact metrics
    genes_affected = set()
    high_impact_genes = set()
    cancer_genes_affected = set()
    
    # Get cancer genes from risk model
    risk_genes = state.get("risk_genes", {})
    all_cancer_genes = set()
    for genes in risk_genes.values():
        all_cancer_genes.update(genes)
    
    for variant in filtered_variants:
        gene = variant.get("gene", "")
        if gene:
            genes_affected.add(gene)
            if gene in all_cancer_genes:
                cancer_genes_affected.add(gene)
            if variant.get("cadd_phred", 0) > 20 or "pathogenic" in variant.get("clinical_significance", "").lower():
                high_impact_genes.add(gene)
    
    variant_metrics["genes_affected"] = len(genes_affected)
    variant_metrics["cancer_genes_affected"] = len(cancer_genes_affected)
    variant_metrics["high_impact_genes"] = len(high_impact_genes)
    
    return variant_metrics


def calculate_prediction_metrics(risk_scores: Dict[str, float],
                               risk_details: Dict[str, Any],
                               ml_assessment: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate metrics about the predictions themselves.
    """
    prediction_metrics = {}
    
    # Risk score statistics by cancer type
    for cancer_type, score in risk_scores.items():
        prediction_metrics[f"{cancer_type}_risk_score"] = score
        
        # Get details for this cancer type
        details = risk_details.get(cancer_type, {})
        if details:
            prediction_metrics[f"{cancer_type}_genes_affected"] = details.get("gene_count", 0)
            prediction_metrics[f"{cancer_type}_pathogenic_genes"] = details.get("pathogenic_count", 0)
    
    # ML fusion aggregate metrics
    if ml_assessment:
        prediction_metrics["aggregate_risk_score"] = ml_assessment.get("aggregate_risk_score", 0.0)
        prediction_metrics["high_risk_variants"] = ml_assessment.get("high_risk_variants", 0)
        
        # Contributing factors
        factors = ml_assessment.get("contributing_factors", {})
        for factor, value in factors.items():
            prediction_metrics[f"factor_{factor}"] = value
    
    return prediction_metrics


def prepare_validation_metrics_structure(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Prepare the structure for validation metrics that would be calculated
    when ground truth data is available (e.g., TCGA validation).
    """
    validation_structure = {
        "validation_ready": False,
        "ground_truth_available": False,
        "metrics_placeholder": {
            # Classification metrics (require ground truth)
            "auc_roc": None,
            "sensitivity": None,
            "specificity": None,
            "f1_score": None,
            "matthews_corrcoef": None,
            "balanced_accuracy": None,
            
            # Regression metrics (for risk scores)
            "mae": None,
            "rmse": None,
            "r2_score": None,
            
            # Agreement metrics
            "concordance_rate": None,
            "cohen_kappa": None,
            
            # Survival metrics (for future implementation)
            "c_index": None,
            "log_rank_p": None
        },
        "validation_note": "These metrics require ground truth labels or TCGA validation data"
    }
    
    # Check if we have any validation data
    if state.get("tcga_validation_data"):
        validation_structure["ground_truth_available"] = True
        validation_structure["validation_ready"] = True
    
    return validation_structure


def calculate_prs_metrics(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate metrics from Polygenic Risk Score results.
    """
    prs_metrics = {}
    
    prs_results = state.get("prs_results", {})
    prs_summary = state.get("prs_summary", {})
    
    if prs_results:
        # PRS score statistics
        prs_scores = []
        prs_percentiles = []
        high_prs_cancers = []
        
        for cancer_type, result in prs_results.items():
            raw_score = result.get("raw_score", 0.0)
            percentile = result.get("percentile", 50)
            confidence = result.get("confidence", "low")
            
            prs_scores.append(raw_score)
            prs_percentiles.append(percentile)
            
            if percentile >= 95:
                high_prs_cancers.append(cancer_type)
            
            prs_metrics[f"prs_{cancer_type}_score"] = raw_score
            prs_metrics[f"prs_{cancer_type}_percentile"] = percentile
            prs_metrics[f"prs_{cancer_type}_confidence"] = confidence
        
        # Aggregate PRS metrics
        if prs_scores:
            prs_metrics["mean_prs_score"] = np.mean(prs_scores)
            prs_metrics["max_prs_percentile"] = np.max(prs_percentiles)
            prs_metrics["high_prs_cancer_count"] = len(high_prs_cancers)
        
        # Overall PRS confidence
        prs_metrics["prs_overall_confidence"] = prs_summary.get("overall_confidence", "low")
    
    return prs_metrics


def calculate_pathway_burden_metrics(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate metrics from pathway burden analysis.
    """
    pathway_metrics = {}
    
    pathway_results = state.get("pathway_burden_results", {})
    pathway_summary = state.get("pathway_burden_summary", {})
    
    if pathway_results:
        burden_scores = []
        high_burden_pathways = []
        
        for pathway_name, result in pathway_results.items():
            burden_score = result.get("burden_score", 0.0)
            burden_scores.append(burden_score)
            
            if burden_score > 0.5:
                high_burden_pathways.append(pathway_name)
            
            pathway_metrics[f"pathway_{pathway_name}_burden"] = burden_score
            pathway_metrics[f"pathway_{pathway_name}_variants"] = result.get("total_variants", 0)
            pathway_metrics[f"pathway_{pathway_name}_damaging"] = result.get("damaging_variants", 0)
        
        # Aggregate pathway metrics
        if burden_scores:
            pathway_metrics["mean_pathway_burden"] = np.mean(burden_scores)
            pathway_metrics["max_pathway_burden"] = np.max(burden_scores)
            pathway_metrics["high_burden_pathway_count"] = len(high_burden_pathways)
        
        # Overall assessment
        pathway_metrics["pathway_risk_level"] = pathway_summary.get("overall_risk_level", "low")
    
    return pathway_metrics


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate comprehensive metrics based on risk model outputs.
    
    This node implements metrics from the Cancer Risk Prediction Metrics document:
    - Confidence metrics (no ground truth needed)
    - Variant quality metrics
    - Prediction statistics
    - Validation structure (for when ground truth is available)
    - Integration of PRS and pathway burden metrics
    """
    logger.info("Starting metrics calculation")
    state["current_node"] = "metrics_calculator"
    
    try:
        # Get risk model outputs
        risk_scores = state.get("risk_scores", {})
        risk_details = state.get("risk_details", {})
        ml_assessment = state.get("ml_risk_assessment", {})
        
        logger.info(f"Calculating metrics for {len(risk_scores)} cancer types")
        
        # Initialize metrics container
        metrics = {
            "timestamp": datetime.now().isoformat(),
            "pipeline_version": state.get("model_version", "1.0.0")
        }
        
        # 1. Calculate confidence metrics
        confidence_metrics = calculate_confidence_metrics(risk_scores, risk_details, ml_assessment)
        metrics["confidence_metrics"] = confidence_metrics
        logger.info(f"Confidence metrics: mean={confidence_metrics.get('mean_model_confidence', 0):.3f}")
        
        # 2. Calculate variant metrics
        variant_metrics = calculate_variant_metrics(state)
        metrics["variant_metrics"] = variant_metrics
        logger.info(f"Variant metrics: {variant_metrics.get('pathogenic_variants', 0)} pathogenic, "
                   f"{variant_metrics.get('uncertain_variants', 0)} uncertain")
        
        # 3. Calculate prediction metrics
        prediction_metrics = calculate_prediction_metrics(risk_scores, risk_details, ml_assessment)
        metrics["prediction_metrics"] = prediction_metrics
        
        # 4. Prepare validation metrics structure
        validation_structure = prepare_validation_metrics_structure(state)
        metrics["validation_metrics"] = validation_structure
        
        # 5. Calculate PRS metrics
        prs_metrics = calculate_prs_metrics(state)
        metrics["prs_metrics"] = prs_metrics
        
        # 6. Calculate pathway burden metrics
        pathway_metrics = calculate_pathway_burden_metrics(state)
        metrics["pathway_metrics"] = pathway_metrics
        
        # 7. Calculate overall risk assessment metrics
        high_risk_cancers = [c for c, score in risk_scores.items() if score > 50]
        metrics["overall_assessment"] = {
            "high_risk_cancer_count": len(high_risk_cancers),
            "high_risk_cancers": high_risk_cancers,
            "max_risk_score": max(risk_scores.values()) if risk_scores else 0,
            "risk_category": ml_assessment.get("risk_category", "unknown"),
            "clinical_action_needed": len(high_risk_cancers) > 0 or 
                                    variant_metrics.get("pathogenic_variants", 0) > 5
        }
        
        # 8. Model performance indicators (without ground truth)
        # These help assess if the model is performing reasonably
        performance_indicators = {
            "variant_coverage": variant_metrics.get("total_variants", 0) > 10,
            "gene_coverage": variant_metrics.get("genes_affected", 0) > 5,
            "model_confidence_adequate": confidence_metrics.get("mean_model_confidence", 0) > 0.7,
            "prs_confidence_adequate": prs_metrics.get("prs_overall_confidence") in ["high", "moderate"],
            "sufficient_evidence": variant_metrics.get("pathogenic_variants", 0) > 0 or 
                                 variant_metrics.get("high_cadd_variants", 0) > 3
        }
        metrics["performance_indicators"] = performance_indicators
        
        # Log summary
        logger.info("Metrics calculation complete:")
        logger.info(f"  Overall risk: {metrics['overall_assessment']['risk_category']}")
        logger.info(f"  High-risk cancers: {', '.join(high_risk_cancers) or 'None'}")
        logger.info(f"  Model confidence: {confidence_metrics.get('mean_model_confidence', 0):.3f}")
        logger.info(f"  Validation ready: {validation_structure['validation_ready']}")
        
        # Update state
        state["metrics"] = metrics
        state["metrics_calculated"] = True
        state["completed_nodes"].append("metrics_calculator")
        
        # Add high-level summary for report
        state["metrics_summary"] = {
            "key_findings": {
                "highest_risk_cancer": max(risk_scores.items(), key=lambda x: x[1])[0] if risk_scores else None,
                "highest_risk_score": max(risk_scores.values()) if risk_scores else 0,
                "pathogenic_variant_count": variant_metrics.get("pathogenic_variants", 0),
                "confidence_level": "high" if performance_indicators.get("model_confidence_adequate") else "moderate"
            },
            "quality_indicators": performance_indicators,
            "validation_status": validation_structure["validation_ready"]
        }
        
    except Exception as e:
        logger.error(f"Metrics calculation failed: {str(e)}")
        state["errors"].append({
            "node": "metrics_calculator",
            "error": str(e),
            "timestamp": datetime.now()
        })
        # Provide minimal metrics on failure
        state["metrics"] = {
            "error": "Metrics calculation failed",
            "timestamp": datetime.now().isoformat()
        }
    
    return state 