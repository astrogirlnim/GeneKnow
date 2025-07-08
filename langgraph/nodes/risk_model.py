"""
Risk model node.
Uses scikit-learn models to predict cancer risk.
"""
import logging
import os
import pickle
import json
import numpy as np
from datetime import datetime
from typing import Dict, Any, List

logger = logging.getLogger(__name__)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate cancer risk scores using ML models.
    
    Updates state with:
    - risk_scores: percentage risk for each cancer type
    - risk_genes: genes contributing to each risk
    """
    logger.info("Starting risk model prediction")
    state["current_node"] = "risk_model"
    
    try:
        filtered_variants = state["filtered_variants"]
        patient_data = state.get("patient_data", {})
        
        # Load model configuration
        models_dir = "models"
        config_path = os.path.join(models_dir, "model_config.json")
        
        if not os.path.exists(config_path):
            # If models don't exist, fallback to simple risk calculation
            logger.warning("ML models not found, using simple risk calculation")
            return _simple_risk_calculation(state, filtered_variants)
        
        with open(config_path, 'r') as f:
            config = json.load(f)
        
        cancer_genes = config["cancer_genes"]
        feature_order = config["feature_order"]
        
        # Extract variant genes
        variant_genes = set()
        gene_to_variant = {}
        for variant in filtered_variants:
            gene = variant.get("gene")
            if gene:
                variant_genes.add(gene)
                gene_to_variant[gene] = variant
        
        logger.info(f"Found variants in genes: {variant_genes}")
        
        # Calculate risk scores for each cancer type
        risk_scores = {}
        risk_genes = {}
        risk_details = {}
        
        for cancer_type, genes in cancer_genes.items():
            logger.info(f"Calculating {cancer_type} cancer risk...")
            
            # Load model and scaler
            model_path = os.path.join(models_dir, f"{cancer_type}_model.pkl")
            scaler_path = os.path.join(models_dir, f"{cancer_type}_scaler.pkl")
            
            if not os.path.exists(model_path) or not os.path.exists(scaler_path):
                logger.warning(f"Model files for {cancer_type} not found")
                continue
            
            with open(model_path, 'rb') as f:
                model = pickle.load(f)
            with open(scaler_path, 'rb') as f:
                scaler = pickle.load(f)
            
            # Prepare feature vector
            features = []
            affected_genes = []
            
            # Gene features (binary: 0 or 1)
            for gene in genes:
                if gene in variant_genes:
                    features.append(1)
                    affected_genes.append(gene)
                else:
                    features.append(0)
            
            # Age feature (normalized 0-1)
            age = patient_data.get("age", 45)  # Default to 45 if not provided
            age_normalized = min(max(age, 20), 80) / 100
            features.append(age_normalized)
            
            # Sex feature (0=M, 1=F)
            sex = patient_data.get("sex", "F")  # Default to F if not provided
            sex_binary = 1 if sex.upper() == "F" else 0
            features.append(sex_binary)
            
            # Make prediction
            features_array = np.array([features])
            features_scaled = scaler.transform(features_array)
            
            # Get probability of high risk
            risk_probability = model.predict_proba(features_scaled)[0, 1]
            risk_percentage = risk_probability * 100
            
            risk_scores[cancer_type] = round(risk_percentage, 1)
            risk_genes[cancer_type] = affected_genes
            
            # Store detailed risk factors
            risk_details[cancer_type] = {
                "affected_genes": affected_genes,
                "gene_count": len(affected_genes),
                "total_genes": len(genes),
                "patient_age": age,
                "patient_sex": sex,
                "model_confidence": round(max(risk_probability, 1-risk_probability), 3)
            }
            
            logger.info(f"{cancer_type} risk: {risk_percentage:.1f}% "
                       f"(genes: {affected_genes}, age: {age}, sex: {sex})")
        
        # Add metadata about the prediction
        state["risk_scores"] = risk_scores
        state["risk_genes"] = risk_genes
        state["risk_details"] = risk_details
        state["model_version"] = config.get("version", "1.0.0")
        
        # Log high-risk findings
        high_risk_cancers = [c for c, r in risk_scores.items() if r > 50]
        if high_risk_cancers:
            logger.warning(f"High risk detected for: {high_risk_cancers}")
            state["high_risk_alert"] = True
            state["high_risk_cancers"] = high_risk_cancers
        
        state["completed_nodes"].append("risk_model")
        logger.info(f"Risk prediction complete: {risk_scores}")
        
    except Exception as e:
        logger.error(f"Risk model failed: {str(e)}")
        state["errors"].append({
            "node": "risk_model",
            "error": str(e),
            "timestamp": datetime.now()
        })
        # Fallback to simple calculation
        return _simple_risk_calculation(state, filtered_variants)
    
    return state


def _simple_risk_calculation(state: Dict[str, Any], 
                           filtered_variants: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Fallback simple risk calculation when ML models aren't available."""
    logger.warning("Using simple risk calculation (ML models not available)")
    
    # High-risk gene weights
    gene_risks = {
        "BRCA1": {"breast": 60, "prostate": 20},
        "BRCA2": {"breast": 45, "prostate": 25},
        "TP53": {"breast": 30, "colon": 25, "lung": 30},
        "APC": {"colon": 50},
        "KRAS": {"colon": 30, "lung": 35},
        "JAK2": {"blood": 40},
        "FLT3": {"blood": 35},
        "EGFR": {"lung": 40},
        "ALK": {"lung": 35}
    }
    
    risk_scores = {"breast": 5.0, "colon": 5.0, "lung": 5.0, "prostate": 5.0, "blood": 5.0}
    risk_genes = {cancer: [] for cancer in risk_scores}
    
    for variant in filtered_variants:
        gene = variant.get("gene")
        if gene in gene_risks:
            for cancer, risk_increase in gene_risks[gene].items():
                risk_scores[cancer] = min(risk_scores[cancer] + risk_increase, 95.0)
                risk_genes[cancer].append(gene)
    
    state["risk_scores"] = risk_scores
    state["risk_genes"] = risk_genes
    state["warnings"].append({
        "node": "risk_model",
        "warning": "Using simple risk calculation - ML models not available",
        "timestamp": datetime.now()
    })
    
    return state 