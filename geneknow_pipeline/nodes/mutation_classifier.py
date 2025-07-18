"""
Mutation Classifier Node for GeneKnow pipeline.
Classifies variants into mutation types (SNV, Indel, CNV, Structural).
"""

from typing import Dict
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


def classify_mutation_type(variant: Dict) -> str:
    """Classify variant based on reference and alternate alleles"""
    ref = variant.get("ref", "")
    alt = variant.get("alt", "")

    # Debug logging for classification issues
    variant_id = variant.get("variant_id", "unknown")
    logger.info(f"🔍 CLASSIFYING {variant_id}: ref='{ref}' (len={len(ref)}), alt='{alt}' (len={len(alt)})")

    # Check for indels first (insertions/deletions with "-" symbol)
    if ref == "-" or alt == "-":
        logger.info(f"→ INDEL: {variant_id} (contains '-' symbol)")
        return "indel"
    
    # Single Nucleotide Variant (SNV)
    elif len(ref) == 1 and len(alt) == 1:
        logger.info(f"→ SNV: {variant_id}")
        return "snv"

    # Insertion/Deletion (Indel) - length difference
    elif len(ref) != len(alt):
        logger.info(f"→ INDEL: {variant_id} (length difference)")
        return "indel"

    # Complex variants (could be structural)
    elif len(ref) > 50 or len(alt) > 50:
        return "structural"

    # Default to SNV for other cases
    return "snv"


def process(state: Dict) -> Dict:
    """Process variants and count mutation types"""
    logger.info("Starting mutation classification")
    # Note: Don't set current_node to avoid concurrent updates

    try:
        variants = state.get("filtered_variants", [])

        mutation_counts = {"snv": 0, "indel": 0, "cnv": 0, "structural": 0}

        # Classify each variant
        classified_variants = []
        for variant in variants:
            mutation_type = classify_mutation_type(variant)
            mutation_counts[mutation_type] += 1

            # Add mutation type to variant
            variant = variant.copy()  # Don't modify original
            variant["mutation_type"] = mutation_type
            classified_variants.append(variant)

        logger.info(f"Classified {len(variants)} variants: {mutation_counts}")

        # Return only the fields this node updates
        return {
            "classified_variants": classified_variants,
            "mutation_type_distribution": mutation_counts,
        }

    except Exception as e:
        logger.error(f"Error in mutation classification: {str(e)}")
        # Return error state updates
        return {
            "errors": [{
                "node": "mutation_classifier",
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            }]
        }
