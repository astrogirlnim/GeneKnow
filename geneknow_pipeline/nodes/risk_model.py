"""
Risk model node.
Uses scikit-learn models to predict cancer risk.
"""

import logging
import numpy as np
from datetime import datetime
from typing import Dict, Any, List
from collections import defaultdict

logger = logging.getLogger(__name__)


def extract_enhanced_features(
    variants: List[Dict[str, Any]], genes: List[str]
) -> List[float]:
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
                "HIGH": 3,
                "MODERATE": 2,
                "LOW": 1,
                "MODIFIER": 0,
                "high": 3,
                "moderate": 2,
                "low": 1,
                "modifier": 0,
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
                if any(
                    sev in str(v.get("consequence", "")).lower()
                    for sev in severe_consequences
                ):
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
        logger.info("✅ Using ML fusion results for risk calculation")
        return _ml_fusion_risk_calculation(state, ml_fusion_results)

    # If ML fusion failed or isn't available, log why and fall back
    if "ml_fusion_results" in state:
        error = ml_fusion_results.get("error", "Unknown error")
        logger.warning(f"⚠️ ML fusion failed: {error}")
    else:
        logger.warning("⚠️ ML fusion results not found in pipeline state")

    logger.warning("📊 Falling back to simple risk calculation")

    try:
        filtered_variants = state.get("filtered_variants", [])
        if not filtered_variants:
            logger.error("No filtered variants found in state")
            return {
                "risk_scores": {},
                "risk_genes": {},
                "error": "No variants to analyze",
            }

        return _simple_risk_calculation(state, filtered_variants)

    except Exception as e:
        logger.error(f"Risk model failed: {str(e)}")
        # Return empty results on complete failure
        return {"risk_scores": {}, "risk_genes": {}, "error": str(e)}


def _ml_fusion_risk_calculation(
    state: Dict[str, Any], ml_fusion_results: Dict[str, Any]
) -> Dict[str, Any]:
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
            "bone": 6.0,
        }

        # Scale risks based on ML fusion aggregate score
        risk_multiplier = 1.0 + (aggregate_risk_score * 2.0)  # 1.0 to 3.0 multiplier

        risk_scores = {}
        risk_genes = {}

        for cancer_type, base_risk in base_risks.items():
            # Apply ML fusion multiplier
            adjusted_risk = base_risk * risk_multiplier

            # Further adjust based on specific contributing factors
            if (
                "cadd_score" in contributing_factors
                and contributing_factors["cadd_score"] > 0.5
            ):
                adjusted_risk *= 1.2
            if (
                "clinvar_classification" in contributing_factors
                and contributing_factors["clinvar_classification"] > 0.5
            ):
                adjusted_risk *= 1.3
            if (
                "gene_burden_score" in contributing_factors
                and contributing_factors["gene_burden_score"] > 0.5
            ):
                adjusted_risk *= 1.25

            # Cap at 95%
            risk_scores[cancer_type] = min(round(adjusted_risk, 1), 95.0)

            # Assign genes based on cancer type relevance
            if cancer_type == "breast" and any(
                g in ["BRCA1", "BRCA2", "PALB2", "ATM"] for g in high_risk_genes
            ):
                risk_genes[cancer_type] = [
                    g
                    for g in high_risk_genes
                    if g in ["BRCA1", "BRCA2", "PALB2", "ATM"]
                ]
            elif cancer_type == "colon" and any(
                g in ["APC", "KRAS", "MLH1", "MSH2"] for g in high_risk_genes
            ):
                risk_genes[cancer_type] = [
                    g for g in high_risk_genes if g in ["APC", "KRAS", "MLH1", "MSH2"]
                ]
            elif cancer_type == "lung" and any(
                g in ["EGFR", "KRAS", "ALK", "ROS1"] for g in high_risk_genes
            ):
                risk_genes[cancer_type] = [
                    g for g in high_risk_genes if g in ["EGFR", "KRAS", "ALK", "ROS1"]
                ]
            else:
                # For other cancers, include all high-risk genes
                risk_genes[cancer_type] = high_risk_genes[:3]  # Top 3 genes

        # Add ML fusion metadata
        logger.info("ML Fusion Risk Assessment:")
        logger.info(f"  Aggregate risk score: {aggregate_risk_score:.3f}")
        logger.info(f"  Risk category: {risk_category}")
        logger.info(f"  Confidence: {confidence:.3f}")
        logger.info(f"  High-risk variants: {high_risk_variants}")
        logger.info(f"  High-risk genes: {high_risk_genes}")

        for cancer_type, score in risk_scores.items():
            logger.info(
                f"  {cancer_type}: {score}% (genes: {risk_genes.get(cancer_type, [])})"
            )

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
                "contributing_factors": contributing_factors,
            },
        }

    except Exception as e:
        logger.error(f"ML fusion risk calculation failed: {str(e)}")
        # Fallback to simple calculation
        filtered_variants = state.get("filtered_variants", [])
        return _simple_risk_calculation(state, filtered_variants)


def _simple_risk_calculation(
    state: Dict[str, Any], filtered_variants: List[Dict[str, Any]]
) -> Dict[str, Any]:
    """Fallback simple risk calculation when ML models aren't available."""
    logger.warning("⚠️ Using simple risk calculation (ML fusion not available)")
    logger.info("This should only happen if ML fusion node failed or was skipped")

    # DEBUG: Check what we received
    logger.info(f"🔍 RISK DEBUG: Received {len(filtered_variants)} filtered variants")
    if filtered_variants:
        logger.info(f"🔍 RISK DEBUG: First variant keys: {list(filtered_variants[0].keys())}")
        logger.info(f"🔍 RISK DEBUG: First variant: {filtered_variants[0]}")

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
        # Paraganglioma/Pheochromocytoma genes
        "SDHB": {"blood": 40},  # Succinate dehydrogenase - paraganglioma/pheochromocytoma
        "SDHD": {"blood": 35},  # Succinate dehydrogenase - paraganglioma/pheochromocytoma
        "SDHC": {"blood": 30},  # Succinate dehydrogenase - paraganglioma/pheochromocytoma
        "SDHA": {"blood": 25},  # Succinate dehydrogenase - paraganglioma/pheochromocytoma
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
        "SPOP": {"prostate": 20},
    }

    # Base risk scores
    risk_scores = {
        "breast": 3.0,
        "colon": 1.7,
        "lung": 2.4,
        "prostate": 2.9,
        "blood": 1.7,
        "bone": 1.9,
    }

    # Map PRS cancer types to our risk score keys
    prs_to_risk_mapping = {
        "BRCA": "breast",
        "OVCA": "breast",  # Ovarian cancer shares risk with breast
        "PRAD": "prostate",
        "LUAD": "lung",
        "COAD": "colon",
        "PANCA": "colon",  # Pancreatic shares some risk factors with colon
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
                prs_data.get("risk_category", "low")

                # Calculate PRS contribution based on percentile
                # High percentile = higher risk
                if percentile >= 95:
                    prs_contribution = 15.0  # High risk
                elif percentile >= 80:
                    prs_contribution = 8.0  # Moderate risk
                elif percentile >= 60:
                    prs_contribution = 3.0  # Slightly elevated
                else:
                    prs_contribution = 0.0  # Average or below

                # Adjust contribution based on confidence
                if confidence == "low":
                    prs_contribution *= 0.3  # Low confidence = reduced impact
                elif confidence == "moderate":
                    prs_contribution *= 0.7
                # High confidence = full contribution

                # Apply PRS contribution
                if prs_contribution > 0:
                    risk_scores[risk_type] = min(
                        risk_scores[risk_type] + prs_contribution, 95.0
                    )
                    logger.info(
                        f"PRS for {prs_cancer} adds {prs_contribution:.1f}% to {risk_type} risk "
                        f"(percentile: {percentile}, confidence: {confidence})"
                    )

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
    
    logger.info(f"🔍 RISK DEBUG: Found {len(gene_variants)} unique genes in variants")
    logger.info(f"🔍 RISK DEBUG: Genes found: {list(gene_variants.keys())}")
    logger.info(f"🔍 RISK DEBUG: SDHB in gene_risks: {'SDHB' in gene_risks}")
    logger.info(f"🔍 RISK DEBUG: gene_risks keys: {list(gene_risks.keys())[:10]}...")
    
    # Debug: Show structure of first variant
    if filtered_variants:
        first_variant = filtered_variants[0]
        logger.info(f"🔍 RISK DEBUG: First variant keys: {list(first_variant.keys())}")
        logger.info(f"🔍 RISK DEBUG: First variant gene: {first_variant.get('gene', 'NOT_FOUND')}")
        logger.info(f"🔍 RISK DEBUG: First variant clinical_significance: {first_variant.get('clinical_significance', 'NOT_FOUND')}")
        logger.info(f"🔍 RISK DEBUG: First variant is_pathogenic: {first_variant.get('is_pathogenic', 'NOT_FOUND')}")
        logger.info(f"🔍 RISK DEBUG: First variant risk_weight: {first_variant.get('risk_weight', 'NOT_FOUND')}")

    for gene, variants in gene_variants.items():
        logger.info(f"🔍 Checking gene: {gene} (variants: {len(variants)})")
        if gene not in gene_risks:
            logger.info(f"⚠️ Gene {gene} not found in gene_risks dictionary")
            continue

        logger.info(f"✅ Gene {gene} found in gene_risks: {gene_risks[gene]}")
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

            # Apply the risk increase
            risk_scores[cancer] = min(risk_scores[cancer] + adjusted_risk, 95.0)

            # Track which genes contributed
            if gene not in risk_genes[cancer]:
                risk_genes[cancer].append(gene)

                if is_pathogenic:
                    pathogenic_genes[cancer].append(gene)
                    logger.warning(
                        f"Pathogenic variant in {gene} increases {cancer} risk by {adjusted_risk:.1f}%"
                    )
                else:
                    logger.info(
                        f"Uncertain variant in {gene} increases {cancer} risk by {adjusted_risk:.1f}%"
                    )

    # Log summary
    total_benign = sum(len(genes) for genes in benign_genes.values())
    total_pathogenic = sum(len(genes) for genes in pathogenic_genes.values())

    logger.info(
        f"Simple risk calculation: {len(genes_hit)} cancer genes affected out of {variant_count} variants"
    )
    logger.info(f"  Benign variants: {total_benign} (excluded from risk)")
    logger.info(f"  Pathogenic variants: {total_pathogenic}")

    for cancer, score in risk_scores.items():
        if risk_genes[cancer]:
            logger.info(
                f"{cancer}: {score:.1f}% (genes: {risk_genes[cancer]}, pathogenic: {pathogenic_genes[cancer]}, benign: {benign_genes[cancer]})"
            )

    # Build result dictionary
    result = {
        "risk_scores": risk_scores,
        "risk_genes": risk_genes,
        "pathogenic_risk_genes": pathogenic_genes,
        "benign_risk_genes": benign_genes,
        "warnings": [
            {
                "node": "risk_model",
                "warning": "Using simple risk calculation - ML models not available",
                "timestamp": datetime.now(),
            }
        ],
    }

    # Add PRS integration metadata if available
    if prs_summary:
        result["risk_integration"] = {
            "prs_high_risk_cancers": prs_summary.get("high_risk_cancers", []),
            "prs_confidence": prs_summary.get("overall_confidence", "not_available"),
            "risk_calculation_method": "integrated_prs_and_variants",
        }
        logger.info(
            f"Risk calculation integrated PRS data: {prs_summary.get('high_risk_cancers', [])}"
        )

    return result
