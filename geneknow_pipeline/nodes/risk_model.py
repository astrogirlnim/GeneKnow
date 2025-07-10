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
from collections import defaultdict

logger = logging.getLogger(__name__)


def extract_enhanced_features(variants: List[Dict[str, Any]], genes: List[str]) -> List[float]:
    """
    Extract comprehensive features for each gene instead of just binary presence.

    For each gene, extracts:
    1. Variant count
    2. Maximum impact score (HIGH=3, MODERATE=2, LOW=1, MODIFIER=0)
    3. Average quality score
    4. Maximum allele frequency
    5. Has frameshift/nonsense mutation
    6. Has homozygous variant
    """
    features = []

    # Create a mapping of genes to their variants
    gene_variants = defaultdict(list)
    for variant in variants:
        gene = variant.get("gene")
        if gene and gene in genes:
            gene_variants[gene].append(variant)

    for gene in genes:
        if gene in gene_variants:
            variants_in_gene = gene_variants[gene]

            # 1. Variant count (capped at 5, normalized)
            variant_count = min(len(variants_in_gene), 5) / 5.0
            features.append(variant_count)

            # 2. Maximum impact score
            impact_scores = {
                "HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0,
                "high": 3, "moderate": 2, "low": 1, "modifier": 0
            }
            max_impact = 0
            for v in variants_in_gene:
                impact = v.get("impact", v.get("consequence", ""))
                score = impact_scores.get(impact, 1)
                max_impact = max(max_impact, score)
            features.append(max_impact / 3.0)  # Normalize to 0-1

            # 3. Average quality score
            qualities = [v.get("quality", 100) for v in variants_in_gene]
            avg_quality = np.mean(qualities) if qualities else 100
            features.append(avg_quality / 100.0)  # Normalize to 0-1

            # 4. Maximum allele frequency
            afs = [v.get("allele_frequency", 0.5) for v in variants_in_gene]
            max_af = max(afs) if afs else 0.5
            features.append(max_af)

            # 5. Has severe consequence (frameshift, nonsense)
            severe_consequences = {"frameshift", "nonsense", "stop_gained", "stop_lost"}
            has_severe = 0
            for v in variants_in_gene:
                if any(sev in str(v.get("consequence", "")).lower() for sev in severe_consequences):
                    has_severe = 1
                    break
            features.append(has_severe)

            # 6. Has homozygous variant
            has_homozygous = 0
            for v in variants_in_gene:
                genotype = v.get("genotype", "")
                if genotype == "1/1" or v.get("allele_frequency", 0) > 0.9:
                    has_homozygous = 1
                    break
            features.append(has_homozygous)
        else:
            # Gene not affected - all features are 0
            features.extend([0, 0, 0, 0, 0, 0])

    return features


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate cancer risk scores using ML fusion results.

    Updates state with:
    - risk_scores: percentage risk for each cancer type
    - risk_genes: genes contributing to each risk
    """
    logger.info("Starting risk model prediction")
    # Note: Don't set current_node to avoid concurrent updates

    # Check if ML fusion results are available
    ml_fusion_results = state.get("ml_fusion_results", {})
    if ml_fusion_results.get("processing_successful"):
        logger.info("Using ML fusion results for risk calculation")
        return _ml_fusion_risk_calculation(state, ml_fusion_results)

    try:
        filtered_variants = state["filtered_variants"]
        patient_data = state.get("patient_data", {})

        # Get PRS results if available
        prs_results = state.get("prs_results", {})
        prs_summary = state.get("prs_summary", {})

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

        # Extract variant genes - now with proper clinical significance
        variant_genes = set()
        pathogenic_genes = set()
        benign_genes = set()
        gene_to_variants = defaultdict(list)

        for variant in filtered_variants:
            gene = variant.get("gene")
            if gene:
                variant_genes.add(gene)
                gene_to_variants[gene].append(variant)

                # Check clinical significance
                clinical_sig = variant.get("clinical_significance", "Unknown").lower()
                risk_weight = variant.get("risk_weight", 0.2)

                if variant.get("is_pathogenic", 0) or "pathogenic" in clinical_sig:
                    pathogenic_genes.add(gene)
                elif "benign" in clinical_sig or risk_weight < 0.15:
                    benign_genes.add(gene)

        logger.info(f"Found variants in {len(variant_genes)} genes: {len(pathogenic_genes)} pathogenic, {len(benign_genes)} benign")

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

            # Check if model expects enhanced features
            expected_features = len(genes) * 6 + 2  # 6 features per gene + age + sex
            legacy_features = len(genes) + 2  # 1 feature per gene + age + sex

            # Extract features based on what the model expects
            # Only count pathogenic/uncertain variants, exclude benign
            risk_contributing_genes = variant_genes - benign_genes

            if hasattr(model, 'n_features_in_') and model.n_features_in_ == expected_features:
                # Use enhanced features (filtered for non-benign variants)
                risk_variants = [v for v in filtered_variants if v.get("gene") not in benign_genes]
                features = extract_enhanced_features(risk_variants, genes)
            else:
                # Use legacy binary features (exclude benign genes)
                features = []
                for gene in genes:
                    if gene in risk_contributing_genes:
                        features.append(1)
                    else:
                        features.append(0)

            affected_genes = [g for g in genes if g in risk_contributing_genes]
            pathogenic_affected = [g for g in affected_genes if g in pathogenic_genes]
            benign_excluded = [g for g in genes if g in benign_genes]

            # Age feature (normalized 0-1)
            age = patient_data.get("age", 45)
            age_normalized = min(max(age, 20), 80) / 100
            features.append(age_normalized)

            # Sex feature (0=M, 1=F)
            sex = patient_data.get("sex", "F")
            sex_binary = 1 if sex.upper() == "F" else 0
            features.append(sex_binary)

            # Make prediction
            features_array = np.array([features])

            # Handle feature mismatch
            if features_array.shape[1] != scaler.n_features_in_:
                logger.warning(f"Feature mismatch for {cancer_type}: got {features_array.shape[1]}, expected {scaler.n_features_in_}")
                # Fallback to simple calculation for this cancer type
                base_risk = 5.0
                gene_risk = len(pathogenic_affected) * 15 + len(affected_genes) * 5
                risk_scores[cancer_type] = min(base_risk + gene_risk, 95.0)
                risk_genes[cancer_type] = affected_genes
                continue

            features_scaled = scaler.transform(features_array)

            # Get probability of high risk
            risk_probability = model.predict_proba(features_scaled)[0, 1]

            # Apply dampening for high variant counts
            if len(filtered_variants) > 100:
                # Reduce confidence when variant count is very high
                dampening_factor = 100 / len(filtered_variants)
                risk_probability = risk_probability * (0.5 + 0.5 * dampening_factor)

            # Boost risk if pathogenic variants found
            if pathogenic_affected:
                risk_probability = min(risk_probability * 1.5, 0.95)

            risk_percentage = risk_probability * 100

            risk_scores[cancer_type] = round(risk_percentage, 1)
            risk_genes[cancer_type] = affected_genes

            # Store detailed risk factors
            risk_details[cancer_type] = {
                "affected_genes": affected_genes,
                "pathogenic_genes": pathogenic_affected,
                "gene_count": len(affected_genes),
                "pathogenic_count": len(pathogenic_affected),
                "total_genes": len(genes),
                "patient_age": age,
                "patient_sex": sex,
                "model_confidence": round(max(risk_probability, 1-risk_probability), 3)
            }

            logger.info(f"{cancer_type} risk: {risk_percentage:.1f}% "
                       f"(genes: {affected_genes}, pathogenic: {pathogenic_affected}, benign excluded: {benign_excluded})")

        # Add metadata about the prediction
        risk_scores_result = risk_scores
        risk_genes_result = risk_genes
        risk_details_result = risk_details
        model_version = config.get("version", "1.0.0")

        # Log high-risk findings
        high_risk_cancers = [c for c, r in risk_scores.items() if r > 50]
        if high_risk_cancers:
            logger.warning(f"High risk detected for: {high_risk_cancers}")

        logger.info(f"Risk prediction complete: {risk_scores}")

        # Return only the keys this node updates
        return {
            "risk_scores": risk_scores_result,
            "risk_genes": risk_genes_result,
            "risk_details": risk_details_result,
            "model_version": model_version,
            "high_risk_alert": True if high_risk_cancers else False,
            "high_risk_cancers": high_risk_cancers
        }

    except Exception as e:
        logger.error(f"Risk model failed: {str(e)}")
        # Fallback to simple calculation
        return _simple_risk_calculation(state, filtered_variants)


def _ml_fusion_risk_calculation(state: Dict[str, Any],
                              ml_fusion_results: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate risk scores from ML fusion results."""
    try:
        # Get aggregate risk assessment
        aggregate = ml_fusion_results.get("aggregate_risk_assessment", {})
        aggregate_risk_score = aggregate.get("aggregate_risk_score", 0.0)
        risk_category = aggregate.get("risk_category", "unknown")
        confidence = aggregate.get("confidence", 0.0)
        high_risk_variants = aggregate.get("high_risk_variants", 0)

        # Get contributing factors
        contributing_factors = aggregate.get("contributing_factors", {})

        # Extract variant genes for risk attribution
        filtered_variants = state.get("filtered_variants", [])
        ml_ready_variants = state.get("ml_ready_variants", [])

        # Group genes by risk level
        high_risk_genes = []
        moderate_risk_genes = []

        # Analyze individual predictions
        individual_predictions = ml_fusion_results.get("individual_predictions", [])
        for i, prediction in enumerate(individual_predictions):
            if i < len(ml_ready_variants):
                variant = ml_ready_variants[i]
                gene = variant.get("gene", "")
                risk_score = prediction.get("risk_score", 0.0)

                if risk_score > 0.7 and gene:
                    high_risk_genes.append(gene)
                elif risk_score > 0.4 and gene:
                    moderate_risk_genes.append(gene)

        # Remove duplicates
        high_risk_genes = list(set(high_risk_genes))
        moderate_risk_genes = list(set(moderate_risk_genes))

        # Calculate cancer-specific risk scores based on ML fusion
        # Base risks adjusted by aggregate risk score
        base_risks = {
            "breast": 15.0,
            "colon": 13.0,
            "lung": 14.0,
            "prostate": 12.0,
            "blood": 8.0,
            "bone": 6.0
        }

        # Scale risks based on ML fusion aggregate score
        risk_multiplier = 1.0 + (aggregate_risk_score * 2.0)  # 1.0 to 3.0 multiplier

        risk_scores = {}
        risk_genes = {}

        for cancer_type, base_risk in base_risks.items():
            # Apply ML fusion multiplier
            adjusted_risk = base_risk * risk_multiplier

            # Further adjust based on specific contributing factors
            if "cadd_score" in contributing_factors and contributing_factors["cadd_score"] > 0.5:
                adjusted_risk *= 1.2
            if "clinvar_classification" in contributing_factors and contributing_factors["clinvar_classification"] > 0.5:
                adjusted_risk *= 1.3
            if "gene_burden_score" in contributing_factors and contributing_factors["gene_burden_score"] > 0.5:
                adjusted_risk *= 1.25

            # Cap at 95%
            risk_scores[cancer_type] = min(round(adjusted_risk, 1), 95.0)

            # Assign genes based on cancer type relevance
            if cancer_type == "breast" and any(g in ["BRCA1", "BRCA2", "PALB2", "ATM"] for g in high_risk_genes):
                risk_genes[cancer_type] = [g for g in high_risk_genes if g in ["BRCA1", "BRCA2", "PALB2", "ATM"]]
            elif cancer_type == "colon" and any(g in ["APC", "KRAS", "MLH1", "MSH2"] for g in high_risk_genes):
                risk_genes[cancer_type] = [g for g in high_risk_genes if g in ["APC", "KRAS", "MLH1", "MSH2"]]
            elif cancer_type == "lung" and any(g in ["EGFR", "KRAS", "ALK", "ROS1"] for g in high_risk_genes):
                risk_genes[cancer_type] = [g for g in high_risk_genes if g in ["EGFR", "KRAS", "ALK", "ROS1"]]
            else:
                # For other cancers, include all high-risk genes
                risk_genes[cancer_type] = high_risk_genes[:3]  # Top 3 genes

        # Add ML fusion metadata
        logger.info(f"ML Fusion Risk Assessment:")
        logger.info(f"  Aggregate risk score: {aggregate_risk_score:.3f}")
        logger.info(f"  Risk category: {risk_category}")
        logger.info(f"  Confidence: {confidence:.3f}")
        logger.info(f"  High-risk variants: {high_risk_variants}")
        logger.info(f"  High-risk genes: {high_risk_genes}")

        for cancer_type, score in risk_scores.items():
            logger.info(f"  {cancer_type}: {score}% (genes: {risk_genes.get(cancer_type, [])})")

        logger.info("Risk model complete (ML fusion)")

        # Return only the keys this node updates
        return {
            "risk_scores": risk_scores,
            "risk_genes": risk_genes,
            "ml_risk_assessment": {
                "method": "ml_fusion",
                "aggregate_risk_score": aggregate_risk_score,
                "risk_category": risk_category,
                "confidence": confidence,
                "high_risk_variants": high_risk_variants,
                "high_risk_genes": high_risk_genes,
                "moderate_risk_genes": moderate_risk_genes,
                "contributing_factors": contributing_factors
            }
        }

    except Exception as e:
        logger.error(f"ML fusion risk calculation failed: {str(e)}")
        # Fallback to simple calculation
        filtered_variants = state.get("filtered_variants", [])
        return _simple_risk_calculation(state, filtered_variants)


def _simple_risk_calculation(state: Dict[str, Any],
                           filtered_variants: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Fallback simple risk calculation when ML models aren't available."""
    logger.warning("Using simple risk calculation (ML models not available)")

    # Get PRS results to incorporate polygenic risk
    prs_results = state.get("prs_results", {})
    prs_summary = state.get("prs_summary", {})

    # High-risk gene weights - now also considering variant impact
    gene_risks = {
        "BRCA1": {"breast": 60, "prostate": 20},
        "BRCA2": {"breast": 45, "prostate": 25},
        "TP53": {"breast": 30, "colon": 25, "lung": 30, "bone": 40, "blood": 35},
        "APC": {"colon": 50},
        "KRAS": {"colon": 30, "lung": 35},
        "JAK2": {"blood": 40},
        "FLT3": {"blood": 35},
        "EGFR": {"lung": 40},
        "ALK": {"lung": 35},
        # Expanded blood cancer genes
        "NPM1": {"blood": 35},
        "DNMT3A": {"blood": 30},
        "IDH1": {"blood": 30, "bone": 35},  # Also in bone (chondrosarcoma)
        "IDH2": {"blood": 30, "bone": 35},  # Also in bone (chondrosarcoma)
        "RUNX1": {"blood": 25},
        "TET2": {"blood": 25},
        "ASXL1": {"blood": 25},
        "CEBPA": {"blood": 20},
        "KIT": {"blood": 20},
        "NRAS": {"blood": 20},
        "CBL": {"blood": 15},
        "EZH2": {"blood": 15},
        "STAT5B": {"blood": 20},  # Important for your test case
        # Bone cancer genes
        "RB1": {"bone": 50},
        "EWSR1": {"bone": 45},
        "FLI1": {"bone": 40},
        "CDKN2A": {"bone": 35},
        "MDM2": {"bone": 30},
        # Additional breast cancer genes
        "PIK3CA": {"breast": 25},
        "PALB2": {"breast": 35},
        "ATM": {"breast": 30},
        "CHEK2": {"breast": 25},
        # Additional colon cancer genes
        "SMAD4": {"colon": 30},
        "BRAF": {"colon": 35},
        "MSH2": {"colon": 40},
        "MLH1": {"colon": 40},
        # Additional lung cancer genes
        "ROS1": {"lung": 35},
        "MET": {"lung": 30},
        # Additional prostate cancer genes
        "AR": {"prostate": 40},
        "PTEN": {"prostate": 35},
        "FOXA1": {"prostate": 25},
        "SPOP": {"prostate": 20}
    }

    # Base risk scores
    risk_scores = {
        "breast": 3.0,
        "colon": 1.7,
        "lung": 2.4,
        "prostate": 2.9,
        "blood": 1.7,
        "bone": 1.9
    }

    # Map PRS cancer types to our risk score keys
    prs_to_risk_mapping = {
        "BRCA": "breast",
        "OVCA": "breast",  # Ovarian cancer shares risk with breast
        "PRAD": "prostate",
        "LUAD": "lung",
        "COAD": "colon",
        "PANCA": "colon"  # Pancreatic shares some risk factors with colon
    }

    # Incorporate PRS scores if available
    if prs_results:
        logger.info("Incorporating PRS scores into risk calculation")
        for prs_cancer, prs_data in prs_results.items():
            if prs_cancer in prs_to_risk_mapping:
                risk_type = prs_to_risk_mapping[prs_cancer]

                # Get PRS percentile and confidence
                percentile = prs_data.get("percentile", 50)
                confidence = prs_data.get("confidence", "low")
                risk_category = prs_data.get("risk_category", "low")

                # Calculate PRS contribution based on percentile
                # High percentile = higher risk
                if percentile >= 95:
                    prs_contribution = 15.0  # High risk
                elif percentile >= 80:
                    prs_contribution = 8.0   # Moderate risk
                elif percentile >= 60:
                    prs_contribution = 3.0   # Slightly elevated
                else:
                    prs_contribution = 0.0   # Average or below

                # Adjust contribution based on confidence
                if confidence == "low":
                    prs_contribution *= 0.3  # Low confidence = reduced impact
                elif confidence == "moderate":
                    prs_contribution *= 0.7
                # High confidence = full contribution

                # Apply PRS contribution
                if prs_contribution > 0:
                    risk_scores[risk_type] = min(risk_scores[risk_type] + prs_contribution, 95.0)
                    logger.info(f"PRS for {prs_cancer} adds {prs_contribution:.1f}% to {risk_type} risk "
                              f"(percentile: {percentile}, confidence: {confidence})")

    risk_genes = {cancer: [] for cancer in risk_scores}
    pathogenic_genes = {cancer: [] for cancer in risk_scores}
    benign_genes = {cancer: [] for cancer in risk_scores}

    # Track how many variants we've processed
    variant_count = len(filtered_variants)
    genes_hit = set()

    # Group variants by gene to avoid double-counting
    gene_variants = defaultdict(list)
    for variant in filtered_variants:
        gene = variant.get("gene")
        if gene:
            gene_variants[gene].append(variant)

    for gene, variants in gene_variants.items():
        if gene not in gene_risks:
            continue

        genes_hit.add(gene)

        # Find the most severe variant in this gene
        max_risk_weight = 0
        most_severe_variant = None

        for variant in variants:
            risk_weight = variant.get("risk_weight", 0.2)  # Default uncertain
            clinical_sig = variant.get("clinical_significance", "Unknown")

            # Only consider the most severe variant per gene
            if risk_weight > max_risk_weight:
                max_risk_weight = risk_weight
                most_severe_variant = variant

        if not most_severe_variant:
            continue

        # Get clinical significance of most severe variant
        is_pathogenic = most_severe_variant.get("is_pathogenic", 0)
        clinical_sig = most_severe_variant.get("clinical_significance", "Unknown")
        risk_weight = most_severe_variant.get("risk_weight", 0.2)

        # Skip benign variants entirely
        if "benign" in clinical_sig.lower() or risk_weight < 0.15:
            for cancer in gene_risks[gene]:
                benign_genes[cancer].append(gene)
            logger.info(f"Skipping benign variant in {gene}: {clinical_sig}")
            continue

        # Calculate risk contribution for each cancer type
        for cancer, base_risk_increase in gene_risks[gene].items():
            # Apply clinical significance weighting
            adjusted_risk = base_risk_increase * risk_weight

            # Additional dampening for high variant counts (tumor samples)
            if variant_count > 100:
                dampening = min(100 / variant_count, 1.0)
                adjusted_risk *= dampening

            # Apply the risk increase
            risk_scores[cancer] = min(risk_scores[cancer] + adjusted_risk, 95.0)

            # Track which genes contributed
            if gene not in risk_genes[cancer]:
                risk_genes[cancer].append(gene)

                if is_pathogenic:
                    pathogenic_genes[cancer].append(gene)
                    logger.warning(f"Pathogenic variant in {gene} increases {cancer} risk by {adjusted_risk:.1f}%")
                else:
                    logger.info(f"Uncertain variant in {gene} increases {cancer} risk by {adjusted_risk:.1f}%")

    # Log summary
    total_benign = sum(len(genes) for genes in benign_genes.values())
    total_pathogenic = sum(len(genes) for genes in pathogenic_genes.values())

    logger.info(f"Simple risk calculation: {len(genes_hit)} cancer genes affected out of {variant_count} variants")
    logger.info(f"  Benign variants: {total_benign} (excluded from risk)")
    logger.info(f"  Pathogenic variants: {total_pathogenic}")

    for cancer, score in risk_scores.items():
        if risk_genes[cancer]:
            logger.info(f"{cancer}: {score:.1f}% (genes: {risk_genes[cancer]}, pathogenic: {pathogenic_genes[cancer]}, benign: {benign_genes[cancer]})")

    # Build result dictionary
    result = {
        "risk_scores": risk_scores,
        "risk_genes": risk_genes,
        "pathogenic_risk_genes": pathogenic_genes,
        "benign_risk_genes": benign_genes,
        "warnings": [{
            "node": "risk_model",
            "warning": "Using simple risk calculation - ML models not available",
            "timestamp": datetime.now()
        }]
    }

    # Add PRS integration metadata if available
    if prs_summary:
        result["risk_integration"] = {
            "prs_high_risk_cancers": prs_summary.get("high_risk_cancers", []),
            "prs_confidence": prs_summary.get("overall_confidence", "not_available"),
            "risk_calculation_method": "integrated_prs_and_variants"
        }
        logger.info(f"Risk calculation integrated PRS data: {prs_summary.get('high_risk_cancers', [])}")

    return result
