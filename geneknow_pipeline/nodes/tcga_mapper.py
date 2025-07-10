"""
TCGA Frequency Mapper Node
Compares patient variants against TCGA tumor frequency database to identify
variants enriched in cancer samples vs normal population.
"""

import os
import sqlite3
import logging
from datetime import datetime
from typing import Dict, Any, List
from collections import defaultdict

logger = logging.getLogger(__name__)

# Path to unified database (contains both population_variants and tcga_variants tables)
TCGA_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "population_variants.db")

def normalize_chromosome(chrom: str) -> str:
    """Normalize chromosome format (remove 'chr' prefix)."""
    if chrom.startswith('chr'):
        return chrom[3:]
    return chrom

def query_tcga_database(chrom: str, pos: int, ref: str, alt: str, cancer_type: str) -> Dict[str, Any]:
    """Query TCGA database for tumor frequency data."""
    if not os.path.exists(TCGA_DB_PATH):
        logger.warning(f"TCGA database not found at {TCGA_DB_PATH}")
        return {}

    conn = sqlite3.connect(TCGA_DB_PATH)
    cursor = conn.cursor()

    try:
        # Normalize chromosome format
        norm_chrom = normalize_chromosome(chrom)

        # Query for exact match
        cursor.execute("""
        SELECT gene, tumor_frequency, normal_frequency, enrichment_score,
               sample_count, total_samples
        FROM tcga_variants
        WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ? AND cancer_type = ?
        """, (norm_chrom, pos, ref, alt, cancer_type))

        row = cursor.fetchone()
        if row:
            return {
                "gene": row[0],
                "tumor_frequency": row[1],
                "normal_frequency": row[2],
                "enrichment_score": row[3],
                "sample_count": row[4],
                "total_samples": row[5],
                "found_in_tcga": True,
                "lookup_method": "exact_match"
            }

        # If no exact match, try gene-based lookup
        cursor.execute("""
        SELECT gene, AVG(tumor_frequency) as avg_tumor_freq,
               AVG(normal_frequency) as avg_normal_freq,
               AVG(enrichment_score) as avg_enrichment,
               SUM(sample_count) as total_samples_with_gene_mutations,
               MAX(total_samples) as cohort_size
        FROM tcga_variants
        WHERE gene IN (
            SELECT DISTINCT gene FROM tcga_variants
            WHERE chrom = ? AND ABS(pos - ?) <= 1000 AND cancer_type = ?
            AND gene IS NOT NULL AND gene != ''
        )
        AND cancer_type = ?
        GROUP BY gene
        LIMIT 1
        """, (norm_chrom, pos, cancer_type, cancer_type))

        row = cursor.fetchone()
        if row and row[0]:  # Make sure gene is not None
            return {
                "gene": row[0],
                "tumor_frequency": row[1] or 0.0,
                "normal_frequency": row[2] or 0.01,
                "enrichment_score": row[3] or 1.0,
                "sample_count": row[4] or 0,
                "total_samples": row[5] or 1000,
                "found_in_tcga": True,
                "lookup_method": "gene_vicinity_match"
            }

        # No match found
        return {
            "found_in_tcga": False,
            "lookup_method": "not_found"
        }

    finally:
        conn.close()

def calculate_enrichment_score(variant: Dict[str, Any], tcga_data: Dict[str, Any]) -> float:
    """Calculate variant enrichment score for cancer risk."""
    if not tcga_data.get("found_in_tcga"):
        return 1.0  # No enrichment if not found in TCGA

    tumor_freq = tcga_data.get("tumor_frequency", 0.0)
    normal_freq = tcga_data.get("normal_frequency", 0.01)

    # Avoid division by zero
    if normal_freq <= 0:
        normal_freq = 0.001

    # Calculate enrichment
    enrichment = tumor_freq / normal_freq

    # Cap enrichment at reasonable levels
    return min(enrichment, 1000.0)

def assess_cancer_relevance(variant: Dict[str, Any], tcga_matches: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Assess overall cancer relevance based on TCGA matches across cancer types."""

    total_enrichment = 0.0
    max_enrichment = 0.0
    best_cancer_match = None
    cancer_types_matched = 0

    for cancer_type, match_data in tcga_matches.items():
        if match_data.get("found_in_tcga"):
            enrichment = match_data.get("enrichment_score", 1.0)
            total_enrichment += enrichment
            cancer_types_matched += 1

            if enrichment > max_enrichment:
                max_enrichment = enrichment
                best_cancer_match = {
                    "cancer_type": cancer_type,
                    "frequency": match_data.get("tumor_frequency", 0.0),
                    "enrichment": enrichment,
                    "sample_count": match_data.get("sample_count", 0),
                    "total_samples": match_data.get("total_samples", 1000)
                }

    # Calculate overall cancer relevance score
    if cancer_types_matched > 0:
        avg_enrichment = total_enrichment / cancer_types_matched
        cancer_relevance = min(avg_enrichment / 10.0, 1.0)  # Normalize to 0-1
    else:
        cancer_relevance = 0.0
        avg_enrichment = 1.0

    return {
        "cancer_relevance_score": cancer_relevance,
        "average_enrichment": avg_enrichment,
        "max_enrichment": max_enrichment,
        "cancer_types_matched": cancer_types_matched,
        "best_match": best_cancer_match
    }

def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Map variants to TCGA tumor frequencies.

    Updates state with:
    - tcga_matches: frequency data for each variant in each cancer type
    - tcga_cohort_sizes: sample sizes for each cancer type
    """
    logger.info("Starting TCGA frequency mapping")
    # Note: Don't set current_node to avoid concurrent updates during parallel execution

    try:
        filtered_variants = state["filtered_variants"]

        # Initialize TCGA matches structure
        cancer_types = ["breast", "colon", "lung", "prostate", "blood"]
        tcga_matches = {cancer_type: {} for cancer_type in cancer_types}

        # Cohort sizes from your documentation
        tcga_cohort_sizes = {
            "breast": 1084,
            "colon": 461,
            "lung": 585,
            "prostate": 498,
            "blood": 200
        }

        # Process each variant
        enriched_variants = []
        total_variants_matched = 0

        for variant in filtered_variants:
            variant_id = variant.get("variant_id", f"{variant.get('chrom')}:{variant.get('pos')}")

            # Query TCGA for each cancer type
            variant_tcga_matches = {}

            for cancer_type in cancer_types:
                tcga_data = query_tcga_database(
                    variant["chrom"],
                    variant["pos"],
                    variant["ref"],
                    variant["alt"],
                    cancer_type
                )

                if tcga_data.get("found_in_tcga"):
                    # Calculate enrichment score
                    enrichment = calculate_enrichment_score(variant, tcga_data)
                    tcga_data["enrichment_score"] = enrichment

                    # Store match data
                    variant_tcga_matches[cancer_type] = tcga_data
                    tcga_matches[cancer_type][variant_id] = tcga_data

                    # Log significant enrichments
                    if enrichment > 10.0:  # More than 10x enriched
                        logger.warning(f"🔥 High enrichment: {variant_id} in {cancer_type} "
                                     f"(tumor: {tcga_data['tumor_frequency']:.1%}, "
                                     f"normal: {tcga_data['normal_frequency']:.1%}, "
                                     f"enrichment: {enrichment:.1f}x)")
                        enriched_variants.append({
                            "variant_id": variant_id,
                            "gene": variant.get("gene"),
                            "cancer_type": cancer_type,
                            "enrichment": enrichment
                        })

            # Assess overall cancer relevance for this variant
            if variant_tcga_matches:
                total_variants_matched += 1
                cancer_assessment = assess_cancer_relevance(variant, variant_tcga_matches)

                # Note: We don't modify the original variant in place to avoid concurrent modification issues
                # The merge node will handle combining TCGA and CADD results
                logger.info(f"🔬 TCGA match: {variant_id} - relevance: {cancer_assessment['cancer_relevance_score']:.2f}, "
                           f"matched in {cancer_assessment['cancer_types_matched']} cancer types")

        # Create summary statistics (will be handled by merge function)
        tcga_summary = {
            "variants_matched": total_variants_matched,
            "total_variants": len(filtered_variants),
            "match_rate": total_variants_matched / len(filtered_variants) if filtered_variants else 0,
            "cancer_types_analyzed": cancer_types,
            "cohort_sizes": tcga_cohort_sizes,
            "highly_enriched_variants": len([v for v in enriched_variants if v["enrichment"] > 50])
        }

        # Log summary
        logger.info(f"TCGA mapping complete:")
        logger.info(f"  Total variants: {len(filtered_variants)}")
        logger.info(f"  Variants matched: {total_variants_matched}")
        logger.info(f"  Match rate: {total_variants_matched/len(filtered_variants)*100:.1f}%" if filtered_variants else "0%")
        logger.info(f"  Highly enriched variants: {len([v for v in enriched_variants if v['enrichment'] > 50])}")

        # Note: Don't append to completed_nodes to avoid concurrent updates
        # The merge node will handle tracking completion

        # Return only the keys this node updates
        return {
            "tcga_matches": tcga_matches,
            "tcga_cohort_sizes": tcga_cohort_sizes,
            "tcga_enriched_variants": enriched_variants,
            "tcga_summary": tcga_summary
        }

    except Exception as e:
        logger.error(f"TCGA mapping failed: {str(e)}")
        return {
            "errors": [{
                "node": "tcga_mapper",
                "error": str(e),
                "timestamp": datetime.now()
            }]
        }
