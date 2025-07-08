"""
TCGA mapper node.
Matches variants against TCGA cancer cohort data.
"""
import os
import sqlite3
import logging
from datetime import datetime
from typing import Dict, Any

logger = logging.getLogger(__name__)

# Path to TCGA database
TCGA_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "tcga_data", "tcga_variants.db")


def query_tcga_database(chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Dict[str, Any]]:
    """
    Query TCGA database for variant frequency across cancer types.
    
    Args:
        chrom: Chromosome
        pos: Position
        ref: Reference allele
        alt: Alternate allele
        
    Returns:
        Dictionary mapping cancer type to frequency data
    """
    if not os.path.exists(TCGA_DB_PATH):
        logger.warning(f"TCGA database not found at {TCGA_DB_PATH}")
        return {}
    
    conn = sqlite3.connect(TCGA_DB_PATH)
    cursor = conn.cursor()
    
    try:
        cursor.execute("""
        SELECT cancer_type, sample_count, total_samples, frequency, gene, clinical_significance
        FROM variants
        WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
        """, (chrom, pos, ref, alt))
        
        results = {}
        for row in cursor.fetchall():
            cancer_type, sample_count, total_samples, frequency, gene, significance = row
            results[cancer_type] = {
                "frequency": frequency,
                "sample_count": sample_count,
                "total_samples": total_samples,
                "gene": gene,
                "clinical_significance": significance
            }
        
        return results
        
    finally:
        conn.close()


def get_cancer_cohort_sizes() -> Dict[str, int]:
    """Get the total sample sizes for each cancer type."""
    if not os.path.exists(TCGA_DB_PATH):
        return {}
    
    conn = sqlite3.connect(TCGA_DB_PATH)
    cursor = conn.cursor()
    
    try:
        cursor.execute("""
        SELECT DISTINCT cancer_type, total_samples
        FROM variants
        """)
        
        cohort_sizes = {}
        for cancer_type, total_samples in cursor.fetchall():
            cohort_sizes[cancer_type] = total_samples
        
        return cohort_sizes
        
    finally:
        conn.close()


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Match variants against TCGA cohorts.
    
    Updates state with:
    - tcga_matches: frequency of each variant in cancer cohorts
    """
    logger.info("Starting TCGA mapping")
    state["current_node"] = "tcga_mapper"
    
    try:
        filtered_variants = state["filtered_variants"]
        
        # Initialize cancer types we're tracking
        # These should match what's in our database
        cancer_types = ["breast", "colon", "lung", "prostate", "blood"]
        tcga_matches = {cancer_type: {} for cancer_type in cancer_types}
        
        # Get cohort sizes
        cohort_sizes = get_cancer_cohort_sizes()
        
        # Process each variant
        matched_count = 0
        for variant in filtered_variants:
            # Query TCGA database
            tcga_results = query_tcga_database(
                variant["chrom"],
                variant["pos"],
                variant["ref"],
                variant["alt"]
            )
            
            if tcga_results:
                matched_count += 1
                logger.debug(f"Found TCGA matches for {variant['variant_id']}: {list(tcga_results.keys())}")
            
            # Add results to appropriate cancer types
            for cancer_type, match_data in tcga_results.items():
                if cancer_type in tcga_matches:
                    tcga_matches[cancer_type][variant["variant_id"]] = {
                        "frequency": match_data["frequency"],
                        "sample_count": match_data["sample_count"],
                        "total_samples": match_data["total_samples"],
                        "clinical_significance": match_data.get("clinical_significance", "Unknown")
                    }
                    
                    # Update variant with TCGA gene if not already set
                    if variant.get("gene") == "Unknown" and match_data.get("gene"):
                        variant["gene"] = match_data["gene"]
        
        # Add unmatched variants with zero frequency
        for variant in filtered_variants:
            variant_id = variant["variant_id"]
            for cancer_type in cancer_types:
                if variant_id not in tcga_matches[cancer_type]:
                    # Variant not found in this cancer type
                    tcga_matches[cancer_type][variant_id] = {
                        "frequency": 0.0,
                        "sample_count": 0,
                        "total_samples": cohort_sizes.get(cancer_type, 0),
                        "clinical_significance": "Not found in TCGA"
                    }
        
        state["tcga_matches"] = tcga_matches
        
        # Log summary statistics
        logger.info(f"TCGA mapping complete: {matched_count}/{len(filtered_variants)} variants found in TCGA")
        
        # Add summary to metadata
        state["file_metadata"]["tcga_summary"] = {
            "variants_matched": matched_count,
            "total_variants": len(filtered_variants),
            "match_rate": matched_count / len(filtered_variants) if filtered_variants else 0,
            "cancer_types_analyzed": list(cancer_types),
            "cohort_sizes": cohort_sizes
        }
        
        state["completed_nodes"].append("tcga_mapper")
        
    except Exception as e:
        logger.error(f"TCGA mapping failed: {str(e)}")
        state["errors"].append({
            "node": "tcga_mapper",
            "error": str(e),
            "timestamp": datetime.now()
        })
        state["pipeline_status"] = "failed"
    
    return state 