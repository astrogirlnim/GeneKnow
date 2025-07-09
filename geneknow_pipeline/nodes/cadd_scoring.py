"""
CADD scoring node.
Enriches variants with CADD PHRED scores for deleteriousness assessment.
Uses cadd_scores table in population_variants.db with job tracking.
"""
import os
import sqlite3
import logging
import requests
import uuid
from typing import Dict, Any, List, Optional, Tuple
from datetime import datetime
import json

logger = logging.getLogger(__name__)

# Configuration - use same database as population_mapper
POP_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "population_variants.db")
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
    """Normalize chromosome format (remove 'chr' prefix to match population_variants)."""
    if chrom.startswith("chr"):
        return chrom[3:]
    return chrom


def lookup_cadd_score(chrom: str, pos: int, ref: str, alt: str, 
                     conn: sqlite3.Connection) -> Optional[Dict[str, Any]]:
    """
    Look up CADD score from local database.
    
    Args:
        chrom: Chromosome
        pos: Position
        ref: Reference allele
        alt: Alternative allele
        conn: Database connection
        
    Returns:
        Dict with raw_score, phred_score, job_id if found, None otherwise
    """
    cursor = conn.cursor()
    
    # Normalize chromosome
    norm_chrom = normalize_chromosome(chrom)
    
    try:
        # Query CADD scores table
        cursor.execute("""
        SELECT raw_score, phred_score, job_id
        FROM cadd_scores
        WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
        """, (norm_chrom, pos, ref, alt))
        
        row = cursor.fetchone()
        if row:
            return {
                "raw": row[0],
                "phred": row[1],
                "job_id": row[2]
            }
            
    except Exception as e:
        logger.error(f"Database error looking up CADD score: {e}")
        
    return None


def query_remote_cadd(chrom: str, pos: int, ref: str, alt: str,
                     job_id: str, conn: sqlite3.Connection) -> Optional[Dict[str, Any]]:
    """
    Query remote CADD service for variant score.
    
    Args:
        chrom: Chromosome
        pos: Position  
        ref: Reference allele
        alt: Alternative allele
        job_id: Job ID for tracking
        conn: Database connection to cache results
        
    Returns:
        Dict with raw and phred scores if successful
    """
    # TODO: Implement actual Tabix remote query
    # For now, just log that we would query
    logger.warning(f"Remote CADD query not implemented. Would query: {chrom}:{pos} {ref}>{alt}")
    
    # In production, this would:
    # 1. Query CADD Tabix service
    # 2. Parse response
    # 3. Cache result in database with job_id
    # 4. Return scores
    
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


def create_job_record(conn: sqlite3.Connection, job_type: str = "pipeline_run") -> str:
    """Create a job tracking record for this CADD scoring run."""
    job_id = f"cadd_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    
    cursor = conn.cursor()
    cursor.execute("""
    INSERT INTO cadd_jobs (job_id, job_type, metadata)
    VALUES (?, ?, ?)
    """, (job_id, job_type, json.dumps({
        "cadd_version": "v1.7",
        "pipeline_node": "cadd_scoring",
        "use_remote": USE_REMOTE_CADD
    })))
    
    conn.commit()
    return job_id


def update_job_status(conn: sqlite3.Connection, job_id: str, 
                     status: str, variant_count: int = 0) -> None:
    """Update job tracking record."""
    cursor = conn.cursor()
    cursor.execute("""
    UPDATE cadd_jobs 
    SET status = ?, completed_at = ?, variant_count = ?
    WHERE job_id = ?
    """, (status, datetime.now() if status in ['completed', 'failed'] else None,
          variant_count, job_id))
    conn.commit()


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Enrich variants with CADD scores for deleteriousness assessment.
    
    Uses metadata from prior steps:
    - Only processes non-benign variants (filtered by population_mapper)
    - Uses gene information already in variants
    - Leverages risk_genes from risk_model if available
    
    Updates state with:
    - cadd_enriched_variants: variants with CADD annotations
    - cadd_stats: summary statistics
    """
    logger.info("Starting CADD scoring enrichment")
    state["current_node"] = "cadd_scoring"
    
    try:
        filtered_variants = state["filtered_variants"]
        
        # Get variants to score based on prior assessments
        # Skip benign variants (already identified by population_mapper)
        benign_variants = set(state.get("benign_variants", []))
        pathogenic_variants = set(state.get("pathogenic_variants", []))
        
        # Get cancer genes from risk assessment if available
        risk_genes = state.get("risk_genes", {})
        all_risk_genes = set()
        for cancer_type, genes in risk_genes.items():
            all_risk_genes.update(genes)
        
        logger.info(f"Processing {len(filtered_variants)} variants")
        logger.info(f"  Skipping {len(benign_variants)} benign variants")
        logger.info(f"  Found {len(all_risk_genes)} cancer risk genes from risk assessment")
        
        # Connect to population database
        if not os.path.exists(POP_DB_PATH):
            logger.warning(f"Population database not found at {POP_DB_PATH}")
            state["cadd_enriched_variants"] = filtered_variants
            state["cadd_stats"] = {"error": "Database not found"}
            state["completed_nodes"].append("cadd_scoring")
            return state
            
        conn = sqlite3.connect(POP_DB_PATH)
        logger.info(f"Connected to population database at {POP_DB_PATH}")
        
        # Create job record
        job_id = create_job_record(conn)
        logger.info(f"Created CADD scoring job: {job_id}")
        
        # Process variants
        enriched_variants = []
        stats = {
            "job_id": job_id,
            "total_variants": len(filtered_variants),
            "variants_scored": 0,
            "variants_skipped_benign": 0,
            "variants_in_cancer_genes": 0,
            "mean_phred": 0.0,
            "max_phred": 0.0,
            "variants_gt20": 0,
            "lookup_missing": 0,
            "remote_lookups": 0
        }
        
        phred_scores = []
        cancer_gene_scores = []
        
        for i, variant in enumerate(filtered_variants):
            variant_id = variant.get("variant_id", f"{variant['chrom']}:{variant['pos']}")
            gene = variant.get("gene", "Unknown")
            
            # Skip if already classified as benign
            if variant_id in benign_variants:
                logger.debug(f"Skipping benign variant: {variant_id}")
                enriched_variants.append(variant)  # Pass through unchanged
                stats["variants_skipped_benign"] += 1
                continue
            
            logger.debug(f"Processing variant {i+1}/{len(filtered_variants)}: {variant_id} in gene {gene}")
            
            # Track if this is a cancer gene
            is_cancer_gene = gene in all_risk_genes
            if is_cancer_gene:
                stats["variants_in_cancer_genes"] += 1
            
            # Look up CADD score
            cadd_result = lookup_cadd_score(
                variant["chrom"],
                variant["pos"],
                variant["ref"],
                variant["alt"],
                conn
            )
            
            # Try remote lookup if not found locally and enabled
            if not cadd_result and USE_REMOTE_CADD:
                cadd_result = query_remote_cadd(
                    variant["chrom"],
                    variant["pos"],
                    variant["ref"],
                    variant["alt"],
                    job_id,
                    conn
                )
                if cadd_result:
                    stats["remote_lookups"] += 1
            
            # Create enriched variant
            enriched_variant = variant.copy()
            
            if cadd_result:
                phred = cadd_result["phred"]
                raw = cadd_result["raw"]
                risk_weight = calculate_risk_weight(phred)
                
                enriched_variant["cadd_phred"] = phred
                enriched_variant["cadd_raw"] = raw
                enriched_variant["cadd_risk_weight"] = risk_weight
                enriched_variant["cadd_job_id"] = cadd_result.get("job_id", job_id)
                
                # Update existing risk weight if lower
                if "risk_weight" in enriched_variant:
                    logger.debug(f"Variant {variant_id}: existing risk_weight={enriched_variant['risk_weight']}, CADD risk_weight={risk_weight}")
                    enriched_variant["risk_weight"] = max(enriched_variant["risk_weight"], risk_weight)
                else:
                    enriched_variant["risk_weight"] = risk_weight
                
                # Log detailed scoring info
                log_msg = f"Variant {variant_id}: CADD PHRED={phred:.1f}, raw={raw:.3f}, risk_weight={risk_weight:.2f}"
                if is_cancer_gene:
                    log_msg += f" [CANCER GENE: {gene}]"
                if variant_id in pathogenic_variants:
                    log_msg += " [PATHOGENIC]"
                    
                logger.info(log_msg)
                
                # Update statistics
                stats["variants_scored"] += 1
                phred_scores.append(phred)
                if is_cancer_gene:
                    cancer_gene_scores.append(phred)
                
                if phred > 20:
                    stats["variants_gt20"] += 1
                    logger.warning(f"High CADD score (>20): {variant_id} in {gene} - PHRED={phred:.1f}")
                    
            else:
                stats["lookup_missing"] += 1
                logger.debug(f"No CADD score found for {variant_id}")
                
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
        
        # Update job status
        update_job_status(conn, job_id, "completed", stats["variants_scored"])
        
        # Close database connection
        conn.close()
        
        # Update state
        state["cadd_enriched_variants"] = enriched_variants
        state["cadd_stats"] = stats
        state["completed_nodes"].append("cadd_scoring")
        
        # Log summary
        logger.info("=" * 60)
        logger.info("CADD Scoring Summary:")
        logger.info(f"  Job ID: {stats['job_id']}")
        logger.info(f"  Total variants: {stats['total_variants']}")
        logger.info(f"  Variants scored: {stats['variants_scored']}")
        logger.info(f"  Benign skipped: {stats['variants_skipped_benign']}")
        logger.info(f"  Cancer gene variants: {stats['variants_in_cancer_genes']}")
        logger.info(f"  Mean PHRED score: {stats['mean_phred']:.1f}")
        logger.info(f"  Max PHRED score: {stats['max_phred']:.1f}")
        logger.info(f"  High-impact variants (>20): {stats['variants_gt20']}")
        logger.info(f"  Missing lookups: {stats['lookup_missing']}")
        logger.info(f"  Remote lookups: {stats['remote_lookups']}")
        logger.info("=" * 60)
        
    except Exception as e:
        logger.error(f"CADD scoring failed: {str(e)}")
        state["errors"].append({
            "node": "cadd_scoring",
            "error": str(e),
            "timestamp": datetime.now()
        })
        # Don't fail the pipeline, just pass through
        state["cadd_enriched_variants"] = state["filtered_variants"]
        state["cadd_stats"] = {"error": str(e)}
        state["completed_nodes"].append("cadd_scoring")
    
    return state 