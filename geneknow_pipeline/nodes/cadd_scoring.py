"""
CADD scoring node.
Enriches variants with CADD PHRED scores for deleteriousness assessment.
Uses local SQLite cache with remote Tabix fallback.
"""
import os
import sqlite3
import logging
import requests
from typing import Dict, Any, List, Optional, Tuple
from datetime import datetime
import json

logger = logging.getLogger(__name__)

# Configuration
DEFAULT_CADD_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "cadd_scores.db")
CADD_DB_PATH = os.environ.get("CADD_DB_PATH", DEFAULT_CADD_DB_PATH)
CADD_REMOTE_TABIX = os.environ.get("CADD_REMOTE_TABIX", "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/")
USE_REMOTE_CADD = os.environ.get("USE_REMOTE_CADD", "true").lower() == "true"

# Risk weight calculation thresholds
CADD_PHRED_THRESHOLDS = {
    "benign": 10.0,        # < 10: likely benign
    "uncertain": 15.0,     # 10-15: uncertain
    "damaging": 20.0,      # 15-20: possibly damaging  
    "pathogenic": 25.0     # > 25: likely pathogenic
}


def normalize_chromosome(chrom: str) -> str:
    """Normalize chromosome format (chr1 vs 1)."""
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def lookup_cadd_score(chrom: str, pos: int, ref: str, alt: str, 
                     conn: Optional[sqlite3.Connection] = None) -> Optional[Dict[str, float]]:
    """
    Look up CADD score for a variant.
    
    Args:
        chrom: Chromosome (e.g., 'chr1' or '1')
        pos: Position (1-based)
        ref: Reference allele
        alt: Alternative allele
        conn: SQLite connection (optional, will create if not provided)
    
    Returns:
        Dict with 'raw' and 'phred' scores, or None if not found
    """
    logger.debug(f"Looking up CADD score for {chrom}:{pos} {ref}>{alt}")
    
    # Normalize chromosome
    chrom_norm = normalize_chromosome(chrom)
    
    # Try local database first
    close_conn = False
    if conn is None and os.path.exists(CADD_DB_PATH):
        logger.debug(f"Opening CADD database: {CADD_DB_PATH}")
        conn = sqlite3.connect(CADD_DB_PATH, check_same_thread=False)
        close_conn = True
    
    if conn:
        try:
            cursor = conn.cursor()
            cursor.execute(
                "SELECT raw_score, phred_score FROM cadd WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?",
                (chrom_norm, pos, ref, alt)
            )
            result = cursor.fetchone()
            
            if result:
                logger.debug(f"Found CADD score in local DB: raw={result[0]}, phred={result[1]}")
                return {"raw": result[0], "phred": result[1]}
            else:
                logger.debug(f"CADD score not found in local DB for {chrom_norm}:{pos}")
                
        except Exception as e:
            logger.error(f"Database lookup error: {str(e)}")
        finally:
            if close_conn:
                conn.close()
    
    # Try remote Tabix query if enabled and local lookup failed
    if USE_REMOTE_CADD:
        logger.debug(f"Attempting remote CADD lookup for {chrom}:{pos}")
        return query_remote_cadd(chrom_norm, pos, ref, alt)
    
    logger.debug(f"No CADD score found for {chrom}:{pos} {ref}>{alt}")
    return None


def query_remote_cadd(chrom: str, pos: int, ref: str, alt: str) -> Optional[Dict[str, float]]:
    """
    Query remote CADD server via Tabix.
    
    Note: This is a placeholder. In production, implement proper Tabix queries
    or use CADD's REST API.
    """
    logger.warning(f"Remote CADD query not implemented. Would query: {chrom}:{pos} {ref}>{alt}")
    
    # Placeholder for remote query implementation
    # In production:
    # 1. Use pysam.TabixFile with remote URL
    # 2. Or use CADD REST API: https://cadd.gs.washington.edu/api
    # 3. Cache results in local DB for future use
    
    return None


def calculate_risk_weight(phred_score: float) -> float:
    """
    Calculate normalized risk weight (0-1) from CADD PHRED score.
    
    Uses min-max scaling with thresholds:
    - PHRED < 10: 0.1 (benign)
    - PHRED 10-15: 0.1-0.3 (uncertain)
    - PHRED 15-20: 0.3-0.6 (damaging)
    - PHRED 20-25: 0.6-0.8 (likely pathogenic)
    - PHRED > 25: 0.8-1.0 (pathogenic)
    """
    if phred_score < CADD_PHRED_THRESHOLDS["benign"]:
        return 0.1
    elif phred_score < CADD_PHRED_THRESHOLDS["uncertain"]:
        # Scale 10-15 to 0.1-0.3
        return 0.1 + 0.2 * (phred_score - 10) / 5
    elif phred_score < CADD_PHRED_THRESHOLDS["damaging"]:
        # Scale 15-20 to 0.3-0.6
        return 0.3 + 0.3 * (phred_score - 15) / 5
    elif phred_score < CADD_PHRED_THRESHOLDS["pathogenic"]:
        # Scale 20-25 to 0.6-0.8
        return 0.6 + 0.2 * (phred_score - 20) / 5
    else:
        # Scale 25+ to 0.8-1.0 (capped at 1.0)
        return min(0.8 + 0.2 * (phred_score - 25) / 10, 1.0)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Enrich variants with CADD scores.
    
    Updates state with:
    - cadd_enriched_variants: variants with added CADD annotations
    - cadd_stats: summary statistics
    """
    logger.info("Starting CADD scoring enrichment")
    logger.info(f"CADD database path: {CADD_DB_PATH}")
    logger.info(f"Remote CADD enabled: {USE_REMOTE_CADD}")
    state["current_node"] = "cadd_scoring"
    
    try:
        # Get filtered variants from state
        filtered_variants = state.get("filtered_variants", [])
        logger.info(f"Processing {len(filtered_variants)} variants for CADD scoring")
        
        # Initialize statistics
        stats = {
            "total_variants": len(filtered_variants),
            "variants_scored": 0,
            "lookup_missing": 0,
            "mean_phred": 0.0,
            "max_phred": 0.0,
            "variants_gt20": 0,
            "phred_distribution": {
                "benign": 0,        # < 10
                "uncertain": 0,     # 10-15
                "damaging": 0,      # 15-20
                "pathogenic": 0     # > 20
            }
        }
        
        # Open database connection once for all lookups
        conn = None
        if os.path.exists(CADD_DB_PATH):
            conn = sqlite3.connect(CADD_DB_PATH, check_same_thread=False)
            logger.info(f"Connected to CADD database: {CADD_DB_PATH}")
        else:
            logger.warning(f"CADD database not found at {CADD_DB_PATH}")
        
        # Process each variant
        enriched_variants = []
        phred_scores = []
        
        for i, variant in enumerate(filtered_variants):
            variant_id = variant.get("variant_id", f"{variant['chrom']}:{variant['pos']}")
            logger.debug(f"Processing variant {i+1}/{len(filtered_variants)}: {variant_id}")
            
            # Look up CADD score
            cadd_result = lookup_cadd_score(
                variant["chrom"],
                variant["pos"],
                variant["ref"],
                variant["alt"],
                conn
            )
            
            # Create enriched variant
            enriched_variant = variant.copy()
            
            if cadd_result:
                phred = cadd_result["phred"]
                raw = cadd_result["raw"]
                risk_weight = calculate_risk_weight(phred)
                
                enriched_variant["cadd_phred"] = phred
                enriched_variant["cadd_raw"] = raw
                enriched_variant["cadd_risk_weight"] = risk_weight
                
                # Update existing risk weight if lower
                if "risk_weight" in enriched_variant:
                    logger.debug(f"Variant {variant_id}: existing risk_weight={enriched_variant['risk_weight']}, CADD risk_weight={risk_weight}")
                    enriched_variant["risk_weight"] = max(enriched_variant["risk_weight"], risk_weight)
                else:
                    enriched_variant["risk_weight"] = risk_weight
                
                # Log detailed scoring info
                logger.info(f"Variant {variant_id}: CADD PHRED={phred:.1f}, raw={raw:.3f}, risk_weight={risk_weight:.2f}")
                
                # Update statistics
                stats["variants_scored"] += 1
                phred_scores.append(phred)
                
                if phred > 20:
                    stats["variants_gt20"] += 1
                
                # Update distribution
                if phred < CADD_PHRED_THRESHOLDS["benign"]:
                    stats["phred_distribution"]["benign"] += 1
                elif phred < CADD_PHRED_THRESHOLDS["uncertain"]:
                    stats["phred_distribution"]["uncertain"] += 1
                elif phred < CADD_PHRED_THRESHOLDS["damaging"]:
                    stats["phred_distribution"]["damaging"] += 1
                else:
                    stats["phred_distribution"]["pathogenic"] += 1
                    
            else:
                # No CADD score found
                logger.warning(f"No CADD score found for variant {variant_id}")
                enriched_variant["cadd_phred"] = None
                enriched_variant["cadd_raw"] = None
                enriched_variant["cadd_risk_weight"] = 0.2  # Default weight for unknown
                
                # Set risk weight if not already set
                if "risk_weight" not in enriched_variant:
                    enriched_variant["risk_weight"] = 0.2
                
                stats["lookup_missing"] += 1
            
            enriched_variants.append(enriched_variant)
        
        # Calculate summary statistics
        if phred_scores:
            stats["mean_phred"] = sum(phred_scores) / len(phred_scores)
            stats["max_phred"] = max(phred_scores)
        
        # Log summary
        logger.info("=" * 60)
        logger.info("CADD Scoring Summary:")
        logger.info(f"  Total variants: {stats['total_variants']}")
        logger.info(f"  Variants scored: {stats['variants_scored']} ({stats['variants_scored']/stats['total_variants']*100:.1f}%)")
        logger.info(f"  Lookup failures: {stats['lookup_missing']}")
        logger.info(f"  Mean PHRED score: {stats['mean_phred']:.1f}")
        logger.info(f"  Max PHRED score: {stats['max_phred']:.1f}")
        logger.info(f"  Variants with PHRED > 20: {stats['variants_gt20']}")
        logger.info("  PHRED distribution:")
        logger.info(f"    Benign (<10): {stats['phred_distribution']['benign']}")
        logger.info(f"    Uncertain (10-15): {stats['phred_distribution']['uncertain']}")
        logger.info(f"    Damaging (15-20): {stats['phred_distribution']['damaging']}")
        logger.info(f"    Pathogenic (>20): {stats['phred_distribution']['pathogenic']}")
        logger.info("=" * 60)
        
        # Update state
        state["filtered_variants"] = enriched_variants  # Update in place
        state["cadd_enriched_variants"] = enriched_variants
        state["cadd_stats"] = stats
        
        # Close database connection
        if conn:
            conn.close()
            
        state["completed_nodes"].append("cadd_scoring")
        logger.info("CADD scoring enrichment complete")
        
    except Exception as e:
        logger.error(f"CADD scoring failed: {str(e)}")
        state["errors"].append({
            "node": "cadd_scoring",
            "error": str(e),
            "timestamp": datetime.now()
        })
        # Don't fail the pipeline - continue with unenriched variants
        state["cadd_stats"] = {"error": str(e)}
        state["completed_nodes"].append("cadd_scoring")
    
    return state 