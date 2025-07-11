"""
CADD scoring node for GeneKnow pipeline.
Computes CADD-like deleteriousness scores locally without any external dependencies.
This is the only CADD scoring implementation - designed for offline desktop applications.
No database lookups or remote queries are performed.
"""
import logging
from typing import Dict, Any, List, Optional
from datetime import datetime
import math

logger = logging.getLogger(__name__)

# Variant impact severity scores (CADD-like approach)
IMPACT_SCORES = {
    # Loss of function variants
    "frameshift_variant": 35.0,
    "stop_gained": 34.0,
    "stop_lost": 33.0,
    "start_lost": 32.0,
    "splice_acceptor_variant": 31.0,
    "splice_donor_variant": 30.0,

    # High impact
    "nonsense_mutation": 28.0,
    "nonstop_mutation": 27.0,
    "frameshift_deletion": 26.0,
    "frameshift_insertion": 25.0,

    # Moderate impact
    "missense_variant": 20.0,
    "missense_mutation": 20.0,
    "inframe_deletion": 18.0,
    "inframe_insertion": 17.0,
    "protein_altering_variant": 16.0,

    # Low impact
    "synonymous_variant": 5.0,
    "silent": 5.0,
    "stop_retained_variant": 8.0,

    # Regulatory
    "regulatory_region_variant": 12.0,
    "5_prime_utr_variant": 10.0,
    "3_prime_utr_variant": 9.0,
    "intron_variant": 6.0,
    "intron": 6.0,

    # Splice region
    "splice_region_variant": 15.0,
    "splice_region": 15.0,

    # Other
    "downstream_gene_variant": 3.0,
    "upstream_gene_variant": 3.0,
    "intergenic_variant": 1.0,

    # Default
    "unknown": 10.0
}

# Conservation-based adjustments
GENE_IMPORTANCE_MULTIPLIERS = {
    # Tumor suppressors
    "TP53": 1.5, "BRCA1": 1.5, "BRCA2": 1.5, "APC": 1.4, "PTEN": 1.4,
    "RB1": 1.4, "VHL": 1.3, "MLH1": 1.3, "MSH2": 1.3, "MSH6": 1.3,

    # Oncogenes
    "KRAS": 1.4, "EGFR": 1.4, "BRAF": 1.4, "PIK3CA": 1.3, "MYC": 1.3,
    "ALK": 1.3, "RET": 1.3, "MET": 1.3, "HER2": 1.3, "ERBB2": 1.3,

    # DNA repair genes
    "ATM": 1.3, "ATR": 1.3, "CHEK1": 1.2, "CHEK2": 1.2, "RAD51": 1.2,

    # Cell cycle regulators
    "CDKN2A": 1.3, "CDK4": 1.2, "CCND1": 1.2,

    # Chromatin modifiers
    "ARID1A": 1.2, "SMARCA4": 1.2, "SMARCB1": 1.2, "EZH2": 1.2,

    # Default for unknown genes
    "default": 1.0
}

# Allele frequency penalty (rare variants are more likely deleterious)
def calculate_af_penalty(af: float) -> float:
    """Calculate penalty based on allele frequency (rare = higher score)."""
    if af < 0.0001:  # Ultra-rare
        return 1.5
    elif af < 0.001:  # Very rare
        return 1.3
    elif af < 0.01:   # Rare
        return 1.1
    elif af < 0.05:   # Low frequency
        return 1.0
    else:             # Common
        return 0.8


def calculate_quality_adjustment(quality: float, depth: int) -> float:
    """Adjust score based on variant quality and read depth."""
    quality_factor = min(quality / 100.0, 1.0) if quality > 0 else 0.5
    depth_factor = min(math.log10(depth + 1) / 3.0, 1.0) if depth > 0 else 0.5
    return (quality_factor + depth_factor) / 2.0


def compute_cadd_score(variant: Dict[str, Any]) -> Dict[str, float]:
    """
    Compute CADD-like scores locally based on variant features.

    Returns dict with:
    - raw_score: Raw CADD-like score
    - phred_score: PHRED-scaled score (like CADD PHRED)
    """
    # Get base impact score
    consequence = variant.get("consequence", "").lower()
    impact = variant.get("impact", "").lower()
    variant_class = variant.get("variant_classification", "").lower()

    # Try multiple fields to find the consequence
    base_score = IMPACT_SCORES.get(consequence,
                 IMPACT_SCORES.get(impact,
                 IMPACT_SCORES.get(variant_class,
                 IMPACT_SCORES["unknown"])))

    # Apply gene importance multiplier
    gene = variant.get("gene", "")
    gene_multiplier = GENE_IMPORTANCE_MULTIPLIERS.get(gene,
                      GENE_IMPORTANCE_MULTIPLIERS["default"])

    # Apply allele frequency penalty
    af = variant.get("allele_frequency", 0.5)
    af_penalty = calculate_af_penalty(af)

    # Apply quality adjustment
    quality = variant.get("quality", 100)
    depth = variant.get("depth", 30)
    quality_adj = calculate_quality_adjustment(quality, depth)

    # Calculate final raw score
    raw_score = base_score * gene_multiplier * af_penalty * quality_adj

    # Convert to PHRED scale (similar to CADD)
    # CADD uses -10 * log10(rank/total) but we'll use a simpler mapping
    # Raw scores 0-50 map to PHRED 0-40
    phred_score = min(raw_score * 0.8, 40.0)

    # Add noise for variants with clinical significance
    clinical_sig = variant.get("clinical_significance", "").lower()
    if "pathogenic" in clinical_sig:
        phred_score = max(phred_score, 25.0)
    elif "likely_pathogenic" in clinical_sig:
        phred_score = max(phred_score, 20.0)
    elif "benign" in clinical_sig or "likely_benign" in clinical_sig:
        phred_score = min(phred_score, 10.0)

    return {
        "raw": round(raw_score, 3),
        "phred": round(phred_score, 1)
    }


def calculate_risk_weight(phred_score: float) -> float:
    """
    Calculate normalized risk weight (0-1) from CADD PHRED score.

    Uses thresholds:
    - PHRED < 10: 0.1 (benign)
    - PHRED 10-15: 0.1-0.3 (uncertain)
    - PHRED 15-20: 0.3-0.6 (damaging)
    - PHRED 20-25: 0.6-0.8 (likely pathogenic)
    - PHRED > 25: 0.8-1.0 (pathogenic)
    """
    if phred_score < 10:
        return 0.1
    elif phred_score < 15:
        return 0.1 + 0.2 * (phred_score - 10) / 5
    elif phred_score < 20:
        return 0.3 + 0.3 * (phred_score - 15) / 5
    elif phred_score < 25:
        return 0.6 + 0.2 * (phred_score - 20) / 5
    else:
        return min(0.8 + 0.2 * (phred_score - 25) / 10, 1.0)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compute CADD-like scores locally for all variants.

    This implementation works completely offline without any database lookups
    or remote queries, suitable for desktop applications.

    Updates state with:
    - cadd_enriched_variants: variants with computed CADD scores
    - cadd_stats: summary statistics
    """
    logger.info("Starting offline CADD scoring")
    # Note: Don't set current_node to avoid concurrent updates during parallel execution

    try:
        filtered_variants = state["filtered_variants"]

        # Get cancer genes from risk assessment if available
        risk_genes = state.get("risk_genes", {})
        all_risk_genes = set()
        for cancer_type, genes in risk_genes.items():
            all_risk_genes.update(genes)

        logger.info(f"Computing CADD scores for {len(filtered_variants)} variants")
        logger.info("Using offline scoring algorithm (no database required)")

        # Process variants
        enriched_variants = []
        stats = {
            "job_id": f"cadd_offline_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
            "total_variants": len(filtered_variants),
            "variants_scored": 0,
            "variants_in_cancer_genes": 0,
            "mean_phred": 0.0,
            "max_phred": 0.0,
            "variants_gt20": 0,
            "scoring_method": "offline_algorithm"
        }

        phred_scores = []
        cancer_gene_scores = []

        for i, variant in enumerate(filtered_variants):
            variant_id = variant.get("variant_id", f"{variant['chrom']}:{variant['pos']}")
            gene = variant.get("gene", "Unknown")

            # Track if this is a cancer gene
            is_cancer_gene = gene in all_risk_genes
            if is_cancer_gene:
                stats["variants_in_cancer_genes"] += 1

            # Compute CADD score locally
            logger.debug(f"Computing score for variant {i+1}/{len(filtered_variants)}: {variant_id}")
            cadd_result = compute_cadd_score(variant)

            # Create enriched variant
            enriched_variant = variant.copy()

            phred = cadd_result["phred"]
            raw = cadd_result["raw"]
            risk_weight = calculate_risk_weight(phred)

            enriched_variant["cadd_phred"] = phred
            enriched_variant["cadd_raw"] = raw
            enriched_variant["cadd_risk_weight"] = risk_weight
            enriched_variant["cadd_source"] = "computed_offline"

            # Update existing risk weight if lower
            if "risk_weight" in enriched_variant:
                logger.debug(f"Variant {variant_id}: existing risk_weight={enriched_variant['risk_weight']}, CADD risk_weight={risk_weight}")
                enriched_variant["risk_weight"] = max(enriched_variant["risk_weight"], risk_weight)
            else:
                enriched_variant["risk_weight"] = risk_weight

            # Log scoring info for significant variants
            if phred > 15 or is_cancer_gene:
                log_msg = f"Variant {variant_id}: CADD PHRED={phred:.1f}, raw={raw:.3f}, risk_weight={risk_weight:.2f}"
                if is_cancer_gene:
                    log_msg += f" [CANCER GENE: {gene}]"
                logger.info(log_msg)

            # Update statistics
            stats["variants_scored"] += 1
            phred_scores.append(phred)
            if is_cancer_gene:
                cancer_gene_scores.append(phred)

            if phred > 20:
                stats["variants_gt20"] += 1
                logger.warning(f"High CADD score (>20): {variant_id} in {gene} - PHRED={phred:.1f}")

            enriched_variants.append(enriched_variant)

        # Calculate summary statistics
        if phred_scores:
            stats["mean_phred"] = sum(phred_scores) / len(phred_scores)
            stats["max_phred"] = max(phred_scores)

        # Log cancer gene statistics
        if cancer_gene_scores:
            mean_cancer_phred = sum(cancer_gene_scores) / len(cancer_gene_scores)
            max_cancer_phred = max(cancer_gene_scores)
            logger.info(f"Cancer gene CADD scores: mean={mean_cancer_phred:.1f}, max={max_cancer_phred:.1f}")

        # Log summary
        logger.info("=" * 60)
        logger.info("Offline CADD Scoring Summary:")
        logger.info(f"  Job ID: {stats['job_id']}")
        logger.info(f"  Total variants: {stats['total_variants']}")
        logger.info(f"  Variants scored: {stats['variants_scored']} (100%)")
        logger.info(f"  Cancer gene variants: {stats['variants_in_cancer_genes']}")
        logger.info(f"  Mean PHRED score: {stats['mean_phred']:.1f}")
        logger.info(f"  Max PHRED score: {stats['max_phred']:.1f}")
        logger.info(f"  High-impact variants (>20): {stats['variants_gt20']}")
        logger.info(f"  Scoring method: Offline algorithm (no internet required)")
        logger.info("=" * 60)

        # Note: Don't append to completed_nodes to avoid concurrent updates
        # The merge node will handle tracking completion

        # Return only the keys this node updates
        return {
            "cadd_enriched_variants": enriched_variants,
            "cadd_stats": stats
        }

    except Exception as e:
        logger.error(f"CADD scoring failed: {str(e)}")
        # Don't fail the pipeline, just pass through
        return {
            "cadd_enriched_variants": state["filtered_variants"],
            "cadd_stats": {"error": str(e)},
            "errors": [{
                "node": "cadd_scoring",
                "error": str(e),
                "timestamp": datetime.now()
            }]
        }
