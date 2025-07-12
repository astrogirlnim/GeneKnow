"""
Quality control filtering node.
Filters variants based on quality metrics.
"""

import logging
from datetime import datetime
from typing import Dict, Any, List

logger = logging.getLogger(__name__)

# QC Thresholds (from documentation)
QUAL_THRESHOLD = 30.0
DEPTH_THRESHOLD = 10
ALLELE_FREQ_THRESHOLD = 0.01


def apply_qc_filters(variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Apply QC filters to a list of variants.

    Args:
        variants: List of variant dictionaries

    Returns:
        List of variants passing QC filters
    """
    filtered_variants = []

    for variant in variants:
        # Apply QC filters
        passed_qc = True

        # Quality score filter
        if variant.get("quality", 0) < QUAL_THRESHOLD:
            passed_qc = False

        # Depth filter
        if variant.get("depth", 0) < DEPTH_THRESHOLD:
            passed_qc = False

        # Allele frequency filter
        if variant.get("allele_freq", 0) < ALLELE_FREQ_THRESHOLD:
            passed_qc = False

        if passed_qc:
            filtered_variants.append(variant)

    return filtered_variants


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Filter variants based on quality control metrics.

    Filters:
    - QUAL >= 30
    - Depth >= 10
    - Allele frequency >= 0.01

    Updates state with:
    - filtered_variants: variants passing QC
    """
    logger.info("Starting QC filtering")
    # Note: Don't set current_node to avoid concurrent updates during parallel execution

    try:
        raw_variants = state["raw_variants"]

        # Use the shared apply_qc_filters function
        filtered_variants = apply_qc_filters(raw_variants)

        # Calculate filtering stats
        total_variants = len(raw_variants)
        passed_variants = len(filtered_variants)
        filter_rate = (total_variants - passed_variants) / total_variants * 100 if total_variants > 0 else 0

        logger.info(f"QC filtering complete: {passed_variants}/{total_variants} passed ({filter_rate:.1f}% filtered)")

        # Prepare metadata update
        file_metadata = state.get("file_metadata", {})
        file_metadata["qc_stats"] = {
            "total_variants": total_variants,
            "passed_qc": passed_variants,
            "failed_qc": total_variants - passed_variants,
            "filter_rate_percent": filter_rate,
            "thresholds": {"qual": QUAL_THRESHOLD, "depth": DEPTH_THRESHOLD, "allele_freq": ALLELE_FREQ_THRESHOLD},
        }

        # Note: Don't append to completed_nodes to avoid concurrent updates
        # The merge node will handle tracking completion

        # Return only the keys this node updates
        return {"filtered_variants": filtered_variants, "file_metadata": file_metadata}

    except Exception as e:
        logger.error(f"QC filtering failed: {str(e)}")
        return {
            "errors": [{"node": "qc_filter", "error": str(e), "timestamp": datetime.now()}],
            "pipeline_status": "failed",
        }
