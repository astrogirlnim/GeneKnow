"""
Population frequency mapper node.
Compares variants against gnomAD/ClinVar population frequencies.
"""

import os
import sqlite3
import logging
from datetime import datetime
from typing import Dict, Any

logger = logging.getLogger(__name__)

# Path to population database
POP_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "population_variants.db")


def normalize_chromosome(chrom: str) -> str:
    """Normalize chromosome format (remove 'chr' prefix)."""
    if chrom.startswith("chr"):
        return chrom[3:]
    return chrom


def query_population_database(
    chrom: str, pos: int, ref: str, alt: str
) -> Dict[str, Any]:
    """Query population database for variant frequency."""
    if not os.path.exists(POP_DB_PATH):
        logger.warning(f"Population database not found at {POP_DB_PATH}")
        return {}

    conn = sqlite3.connect(POP_DB_PATH)
    cursor = conn.cursor()

    try:
        # Normalize chromosome format
        norm_chrom = normalize_chromosome(chrom)

        # Try multiple lookup strategies due to data inconsistencies

        # Strategy 1: Exact match with normalized chromosome
        cursor.execute(
            """
        SELECT gene, gnomad_af, clinical_significance, is_pathogenic, consequence
        FROM population_variants
        WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
        """,
            (norm_chrom, pos, ref, alt),
        )

        row = cursor.fetchone()
        if row:
            return {
                "gene": row[0],
                "population_frequency": row[1],
                "clinical_significance": row[2],
                "is_pathogenic": row[3],
                "consequence": row[4],
                "in_population": True,
                "lookup_method": "exact_match",
            }

        # Strategy 2: Position and gene match (since ref/alt might be 'na')
        cursor.execute(
            """
        SELECT gene, gnomad_af, clinical_significance, is_pathogenic, consequence
        FROM population_variants
        WHERE chrom = ? AND pos = ? AND gene IS NOT NULL
        LIMIT 1
        """,
            (norm_chrom, pos),
        )

        row = cursor.fetchone()
        if row:
            logger.info(f"Found variant by position match: {norm_chrom}:{pos}")
            return {
                "gene": row[0],
                "population_frequency": row[1],
                "clinical_significance": row[2],
                "is_pathogenic": row[3],
                "consequence": row[4],
                "in_population": True,
                "lookup_method": "position_match",
            }

        # Strategy 3: Gene-based lookup for known cancer genes
        # Look for any variant in the same gene at nearby positions
        cursor.execute(
            """
        SELECT gene, AVG(gnomad_af) as avg_af, clinical_significance, AVG(is_pathogenic) as avg_pathogenic, consequence
        FROM population_variants
        WHERE gene IN (
            SELECT DISTINCT gene FROM population_variants
            WHERE chrom = ? AND ABS(pos - ?) <= 1000
            AND gene IS NOT NULL AND gene != ''
        )
        AND clinical_significance LIKE '%benign%'
        GROUP BY gene
        LIMIT 1
        """,
            (norm_chrom, pos),
        )

        row = cursor.fetchone()
        if row:
            logger.info(f"Found similar variant in gene {row[0]} near position {pos}")
            return {
                "gene": row[0],
                "population_frequency": row[1] or 0.01,  # Default to 1% if null
                "clinical_significance": row[2],
                "is_pathogenic": int(row[3] > 0.5) if row[3] is not None else 0,
                "consequence": row[4],
                "in_population": True,
                "lookup_method": "gene_vicinity_match",
            }

        # Strategy 4: Not found - return default safe values
        logger.info(
            f"Variant not found in population database: {chrom}:{pos}:{ref}>{alt}"
        )
        return {
            "population_frequency": 0.0,
            "clinical_significance": "Not_in_database",
            "is_pathogenic": 0,
            "in_population": False,
            "lookup_method": "not_found",
        }

    finally:
        conn.close()


def assess_clinical_significance(
    variant: Dict[str, Any], pop_data: Dict[str, Any]
) -> Dict[str, Any]:
    """Assess clinical significance based on multiple factors."""

    # Start with population database info
    clinical_sig = pop_data.get("clinical_significance", "Unknown").lower()
    pop_data.get("is_pathogenic", 0)
    pop_freq = pop_data.get("population_frequency", 0)

    # Override based on variant characteristics
    consequence = variant.get("consequence", "").lower()

    # Benign criteria
    if any(term in clinical_sig for term in ["benign", "likely_benign"]):
        return {
            "final_clinical_significance": "Benign",
            "is_pathogenic": 0,
            "risk_weight": 0.1,  # Very low impact
            "rationale": "Annotated as benign in population database",
        }

    # High frequency = likely benign
    if pop_freq > 0.01:  # >1% population frequency
        return {
            "final_clinical_significance": "Likely_benign_common",
            "is_pathogenic": 0,
            "risk_weight": 0.1,
            "rationale": f"Common variant (AF={pop_freq:.3f})",
        }

    # Synonymous variants are usually benign
    if "synonymous" in consequence or "silent" in consequence:
        return {
            "final_clinical_significance": "Likely_benign_synonymous",
            "is_pathogenic": 0,
            "risk_weight": 0.2,
            "rationale": "Synonymous variant",
        }

    # Pathogenic criteria
    if any(term in clinical_sig for term in ["pathogenic", "likely_pathogenic"]):
        return {
            "final_clinical_significance": "Pathogenic",
            "is_pathogenic": 1,
            "risk_weight": 1.0,
            "rationale": "Annotated as pathogenic",
        }

    # Rare missense variants
    if "missense" in consequence and pop_freq < 0.001:
        return {
            "final_clinical_significance": "Uncertain_significance",
            "is_pathogenic": 0,
            "risk_weight": 0.3,
            "rationale": "Rare missense variant",
        }

    # Default for unknown
    return {
        "final_clinical_significance": "Uncertain_significance",
        "is_pathogenic": 0,
        "risk_weight": 0.2,
        "rationale": "Insufficient evidence",
    }


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compare variants against population frequencies.

    Updates state with:
    - population_matches: frequency of each variant in general population
    - rare_variants: list of variants rare in population (good candidates for pathogenicity)
    """
    logger.info("Starting population frequency mapping")
    # Note: Don't set current_node to avoid concurrent updates

    try:
        filtered_variants = state["filtered_variants"]
        population_matches = {}
        rare_variants = []
        pathogenic_variants = []
        benign_variants = []

        # Process each variant
        for variant in filtered_variants:
            # Query population database
            pop_data = query_population_database(
                variant["chrom"], variant["pos"], variant["ref"], variant["alt"]
            )

            # Assess clinical significance
            clinical_assessment = assess_clinical_significance(variant, pop_data)

            # Store population data
            variant_id = variant["variant_id"]
            population_matches[variant_id] = {**pop_data, **clinical_assessment}

            # Update variant with improved annotations
            variant["population_frequency"] = pop_data.get("population_frequency", 0.0)
            variant["clinical_significance"] = clinical_assessment[
                "final_clinical_significance"
            ]
            variant["is_pathogenic"] = clinical_assessment["is_pathogenic"]
            variant["risk_weight"] = clinical_assessment["risk_weight"]
            variant["assessment_rationale"] = clinical_assessment["rationale"]

            # Categorize variants
            if clinical_assessment["is_pathogenic"]:
                pathogenic_variants.append(variant_id)
                logger.warning(
                    f"⚠️  Pathogenic variant: {variant_id} - {clinical_assessment['rationale']}"
                )
            elif clinical_assessment["final_clinical_significance"] in [
                "Benign",
                "Likely_benign_common",
                "Likely_benign_synonymous",
            ]:
                benign_variants.append(variant_id)
                logger.info(
                    f"✅ Benign variant: {variant_id} - {clinical_assessment['rationale']}"
                )
            elif pop_data.get("population_frequency", 0) < 0.001:
                rare_variants.append(variant_id)
                logger.info(
                    f"🔬 Rare variant: {variant_id} (AF: {pop_data.get('population_frequency', 0):.5f}) - {clinical_assessment['rationale']}"
                )

        # Log summary
        logger.info("Population mapping complete:")
        logger.info(f"  Total variants: {len(filtered_variants)}")
        logger.info(f"  Benign variants: {len(benign_variants)}")
        logger.info(f"  Rare variants: {len(rare_variants)}")
        logger.info(f"  Pathogenic variants: {len(pathogenic_variants)}")

        # Return only the keys this node updates
        return {
            "population_matches": population_matches,
            "rare_variants": rare_variants,
            "pathogenic_variants": pathogenic_variants,
            "benign_variants": benign_variants,
            "filtered_variants": filtered_variants,  # Return updated variants
        }

    except Exception as e:
        logger.error(f"Population mapping failed: {str(e)}")
        return {
            "errors": [
                {
                    "node": "population_mapper",
                    "error": str(e),
                    "timestamp": datetime.now(),
                }
            ]
        }
