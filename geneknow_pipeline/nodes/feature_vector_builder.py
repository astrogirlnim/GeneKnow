"""
Feature vector builder node.
Collects outputs from all static risk model nodes and builds feature vectors for ML fusion.
"""

import logging
from datetime import datetime
from typing import Dict, Any

logger = logging.getLogger(__name__)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Build feature vectors from all static model outputs for ML fusion.

    Prepares variants with scores from:
    - PRS (Polygenic Risk Scores)
    - ClinVar annotations
    - CADD scores
    - TCGA frequency matching
    - Gene/Pathway burden
    """
    logger.info("Building comprehensive feature vectors from static models")
    # Note: Don't set current_node to avoid concurrent updates

    try:
        # Get merged variants with all annotations
        filtered_variants = state.get("filtered_variants", [])

        # Get state-level scores
        prs_results = state.get("prs_results", {})
        state.get("prs_summary", {})
        pathway_burden_results = state.get("pathway_burden_results", {})
        state.get("pathway_burden_summary", {})

        logger.info(
            f"Feature vector builder received {len(filtered_variants)} variants"
        )

        # Prepare enriched variants for ML fusion
        ml_ready_variants = []

        for variant in filtered_variants:
            # Create enhanced variant with all scores needed for ML fusion
            enhanced_variant = variant.copy()

            # 1. PRS Score - use highest cancer-specific PRS as proxy
            # Since PRS is population-level, we'll use the highest risk score
            prs_scores = []
            for cancer_type, prs_data in prs_results.items():
                prs_data.get("raw_score", 0.0)
                percentile = prs_data.get("percentile", 50) / 100.0
                # Use percentile as normalized PRS score
                prs_scores.append(percentile)

            enhanced_variant["prs_score"] = max(prs_scores) if prs_scores else 0.5

            # 2. ClinVar classification
            clinvar_sig = variant.get("clinvar_clinical_significance")
            if clinvar_sig:
                enhanced_variant["clinvar"] = {
                    "clinical_significance": clinvar_sig,
                    "risk_score": variant.get("clinvar_risk_score", 0.0),
                    "interpretation": variant.get("clinvar_interpretation", "unknown"),
                }
            else:
                enhanced_variant["clinvar"] = {
                    "clinical_significance": "not_found",
                    "risk_score": 0.0,
                    "interpretation": "not_found",
                }

            # 3. CADD Score
            cadd_phred = variant.get("cadd_phred", 0.0)
            enhanced_variant["cadd_score"] = float(cadd_phred)

            # 4. TCGA Enrichment
            tcga_relevance = variant.get("tcga_cancer_relevance", 0.0)
            tcga_match = variant.get("tcga_best_match", {})
            if tcga_match:
                # Use enrichment score if available
                enhanced_variant["tcga_enrichment"] = tcga_match.get("enrichment", 1.0)
            else:
                # Otherwise use cancer relevance
                enhanced_variant["tcga_enrichment"] = max(tcga_relevance * 10.0, 1.0)

            # 5. Gene/Pathway Burden Score
            pathway_damage = variant.get("pathway_damage_assessment", {})
            if pathway_damage:
                enhanced_variant["gene_burden_score"] = pathway_damage.get(
                    "damage_score", 0.0
                )
            else:
                # Fallback: calculate from pathway burden results
                gene = variant.get("gene", "")
                burden_score = 0.0
                for pathway_name, pathway_data in pathway_burden_results.items():
                    if gene in pathway_data.get("damaging_genes", []):
                        burden_score = max(
                            burden_score, pathway_data.get("burden_score", 0.0)
                        )
                enhanced_variant["gene_burden_score"] = burden_score

            # Add additional metadata for risk calculation
            enhanced_variant["is_high_impact"] = (
                cadd_phred > 20
                or "pathogenic" in str(clinvar_sig).lower()
                or burden_score > 0.5
            )

            ml_ready_variants.append(enhanced_variant)

        # Log feature extraction summary
        logger.info("=" * 60)
        logger.info("Feature Vector Builder Summary:")
        logger.info(f"  Total variants: {len(ml_ready_variants)}")
        logger.info(
            f"  Variants with PRS data: {sum(1 for v in ml_ready_variants if v.get('prs_score', 0) != 0.5)}"
        )
        logger.info(
            f"  Variants with ClinVar: {sum(1 for v in ml_ready_variants if v['clinvar']['clinical_significance'] != 'not_found')}"
        )
        logger.info(
            f"  Variants with CADD > 20: {sum(1 for v in ml_ready_variants if v.get('cadd_score', 0) > 20)}"
        )
        logger.info(
            f"  Variants with TCGA enrichment > 2: {sum(1 for v in ml_ready_variants if v.get('tcga_enrichment', 1) > 2)}"
        )
        logger.info(
            f"  Variants with pathway burden > 0: {sum(1 for v in ml_ready_variants if v.get('gene_burden_score', 0) > 0)}"
        )
        logger.info(
            f"  High impact variants: {sum(1 for v in ml_ready_variants if v.get('is_high_impact', False))}"
        )

        # Log a sample variant
        if ml_ready_variants:
            sample = ml_ready_variants[0]
            logger.info("Sample ML-ready variant:")
            logger.info(f"  Gene: {sample.get('gene', 'N/A')}")
            logger.info(f"  PRS score: {sample.get('prs_score', 'N/A')}")
            logger.info(f"  ClinVar: {sample['clinvar']['clinical_significance']}")
            logger.info(f"  CADD score: {sample.get('cadd_score', 'N/A')}")
            logger.info(f"  TCGA enrichment: {sample.get('tcga_enrichment', 'N/A')}")
            logger.info(f"  Gene burden: {sample.get('gene_burden_score', 'N/A')}")

        logger.info("=" * 60)

        # Build result with ML-ready variants
        logger.info("Feature vector builder complete - data ready for ML fusion")

        # Return only the keys this node updates
        return {
            "ml_ready_variants": ml_ready_variants,
            "feature_vector": {
                "status": "complete",
                "variant_count": len(ml_ready_variants),
                "features_extracted": [
                    "prs_score",
                    "clinvar",
                    "cadd_score",
                    "tcga_enrichment",
                    "gene_burden_score",
                ],
                "high_impact_count": sum(
                    1 for v in ml_ready_variants if v.get("is_high_impact", False)
                ),
            },
        }

    except Exception as e:
        logger.error(f"Feature vector builder failed: {str(e)}")
        # Return error state updates
        return {
            "ml_ready_variants": [],
            "feature_vector": {"status": "failed", "error": str(e)},
            "errors": [
                {
                    "node": "feature_vector_builder",
                    "error": str(e),
                    "timestamp": datetime.now(),
                }
            ],
        }
