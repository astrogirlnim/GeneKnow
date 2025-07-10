"""
Polygenic Risk Score (PRS) calculator node.
Calculates weighted sum of small-effect SNPs associated with cancer risk.
"""
import os
import sqlite3
import logging
import math
from datetime import datetime
from typing import Dict, Any, List, Tuple, Optional
from collections import defaultdict

logger = logging.getLogger(__name__)

# Path to PRS database
PRS_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "prs_snps.db")

# Cancer types we support
SUPPORTED_CANCERS = ["BRCA", "OVCA", "PRAD", "LUAD", "COAD", "PANCA"]

# Population groups
POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR", "GLOBAL"]

# Confidence thresholds
MIN_COVERAGE_HIGH_CONF = 0.8
MIN_COVERAGE_MODERATE_CONF = 0.5


def create_prs_database():
    """Create PRS database with example SNPs if it doesn't exist."""
    if os.path.exists(PRS_DB_PATH):
        return

    logger.info("Creating PRS database with example data...")
    conn = sqlite3.connect(PRS_DB_PATH)
    cursor = conn.cursor()

    # Create table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS prs_snps (
        snp_id TEXT PRIMARY KEY,
        chrom TEXT,
        pos INTEGER,
        ref TEXT,
        alt TEXT,
        risk_allele TEXT,
        effect_size REAL,
        cancer_type TEXT,
        p_value REAL,
        source_study TEXT,
        population TEXT,
        maf REAL
    )
    """)

    # Add example PRS SNPs for different cancers
    example_snps = [
        # Breast cancer SNPs
        ("17:43044295:G>A", "17", 43044295, "G", "A", "A", 2.3, "BRCA", 1e-20, "PMID:28108588", "EUR", 0.0001),
        ("13:32911888:A>G", "13", 32911888, "A", "G", "G", 1.8, "BRCA", 1e-15, "PMID:28108588", "EUR", 0.0002),
        ("10:123352317:C>T", "10", 123352317, "C", "T", "T", 0.15, "BRCA", 5e-8, "PMID:28108588", "EUR", 0.23),
        ("5:56031884:T>C", "5", 56031884, "T", "C", "C", 0.12, "BRCA", 3e-9, "PMID:28108588", "EUR", 0.31),

        # Ovarian cancer SNPs (some overlap with BRCA)
        ("17:43044295:G>A", "17", 43044295, "G", "A", "A", 3.1, "OVCA", 1e-25, "PMID:28346442", "EUR", 0.0001),
        ("9:22125503:G>C", "9", 22125503, "G", "C", "C", 0.25, "OVCA", 1e-10, "PMID:28346442", "EUR", 0.18),

        # Prostate cancer SNPs
        ("8:128081119:T>C", "8", 128081119, "T", "C", "C", 0.18, "PRAD", 2e-11, "PMID:29892016", "EUR", 0.42),
        ("17:47440843:A>G", "17", 47440843, "A", "G", "G", 0.22, "PRAD", 8e-13, "PMID:29892016", "EUR", 0.35),

        # Population-specific effect sizes for same SNP
        ("10:123352317:C>T", "10", 123352317, "C", "T", "T", 0.08, "BRCA", 5e-6, "PMID:31427789", "AFR", 0.41),
        ("10:123352317:C>T", "10", 123352317, "C", "T", "T", 0.13, "BRCA", 2e-7, "PMID:31427789", "EAS", 0.19),
    ]

    cursor.executemany("""
    INSERT OR REPLACE INTO prs_snps VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, example_snps)

    conn.commit()
    conn.close()
    logger.info(f"Created PRS database with {len(example_snps)} example SNPs")


def normalize_chromosome(chrom: str) -> str:
    """Normalize chromosome format."""
    if chrom.startswith('chr'):
        return chrom[3:]
    return chrom


def infer_genotype(variant: Dict[str, Any], maf: Optional[float] = None) -> Tuple[int, str]:
    """
    Infer genotype (number of risk alleles) when not directly observed.

    Returns: (num_risk_alleles, inference_method)
    """
    # Check if genotype is provided
    if "genotype" in variant:
        gt = variant["genotype"]
        if gt == "0/0":
            return (0, "observed")
        elif gt in ["0/1", "1/0"]:
            return (1, "observed")
        elif gt == "1/1":
            return (2, "observed")

    # Use population frequency for inference
    pop_freq = variant.get("population_frequency", 0)

    # For rare variants, assume heterozygous
    if pop_freq < 0.01:
        return (1, "rare_variant_assumption")

    # For common variants, use Hardy-Weinberg equilibrium
    if maf and maf > 0:
        # Calculate genotype probabilities
        p = maf  # Risk allele frequency
        q = 1 - p  # Reference allele frequency

        # If we see the variant, it's either het or hom alt
        # P(het | has variant) = 2pq / (2pq + p²)
        prob_het = (2 * p * q) / (2 * p * q + p * p)

        # Most likely heterozygous for common variants
        if prob_het > 0.8:
            return (1, "hardy_weinberg_likely_het")
        else:
            return (2, "hardy_weinberg_possible_hom")

    # Default conservative assumption
    return (1, "default_heterozygous")


def query_prs_database(chrom: str, pos: int, ref: str, alt: str,
                      cancer_type: Optional[str] = None, population: str = "EUR") -> List[Dict[str, Any]]:
    """Query PRS database for effect sizes."""
    if not os.path.exists(PRS_DB_PATH):
        create_prs_database()

    conn = sqlite3.connect(PRS_DB_PATH)
    cursor = conn.cursor()

    results = []
    try:
        norm_chrom = normalize_chromosome(chrom)

        # Build query based on parameters
        query = """
        SELECT chrom, pos, ref, alt, risk_allele, effect_size,
               cancer_type, p_value, source_study, population, maf
        FROM prs_snps
        WHERE chrom = ? AND pos = ?
        """
        params = [norm_chrom, pos]

        if cancer_type:
            query += " AND cancer_type = ?"
            params.append(cancer_type)

        cursor.execute(query, params)

        for row in cursor.fetchall():
            # Check if ref/alt match (considering strand flips)
            db_ref, db_alt = row[2], row[3]
            if (ref == db_ref and alt == db_alt) or \
               (ref == db_alt and alt == db_ref):  # Strand flip

                results.append({
                    "chrom": row[0],
                    "pos": row[1],
                    "ref": row[2],
                    "alt": row[3],
                    "risk_allele": row[4],
                    "effect_size": row[5],
                    "cancer_type": row[6],
                    "p_value": row[7],
                    "source_study": row[8],
                    "population": row[9],
                    "maf": row[10]
                })

    finally:
        conn.close()

    return results


def calculate_cancer_specific_prs(variants: List[Dict[str, Any]],
                                 cancer_type: str,
                                 population: str = "EUR") -> Dict[str, Any]:
    """Calculate PRS for a specific cancer type."""

    # Get all PRS SNPs for this cancer type
    conn = sqlite3.connect(PRS_DB_PATH)
    cursor = conn.cursor()

    cursor.execute("""
    SELECT COUNT(*) FROM prs_snps
    WHERE cancer_type = ? AND population = ?
    """, (cancer_type, population))

    total_snps = cursor.fetchone()[0]

    # If no population-specific data, fall back to EUR or GLOBAL
    if total_snps == 0 and population != "EUR":
        logger.warning(f"No PRS data for {cancer_type} in {population}, using EUR")
        population = "EUR"
        cursor.execute("""
        SELECT COUNT(*) FROM prs_snps
        WHERE cancer_type = ? AND population = ?
        """, (cancer_type, population))
        total_snps = cursor.fetchone()[0]

    conn.close()

    # Calculate PRS
    prs_sum = 0.0
    matched_snps = 0
    contributing_snps = []

    for variant in variants:
        # Query for PRS effect
        prs_hits = query_prs_database(
            variant["chrom"],
            variant["pos"],
            variant["ref"],
            variant["alt"],
            cancer_type,
            population
        )

        for hit in prs_hits:
            # Determine number of risk alleles
            num_risk_alleles, inference_method = infer_genotype(variant, hit["maf"])

            # Calculate contribution
            contribution = hit["effect_size"] * num_risk_alleles
            prs_sum += contribution
            matched_snps += 1

            contributing_snps.append({
                "variant_id": variant["variant_id"],
                "effect_size": hit["effect_size"],
                "risk_alleles": num_risk_alleles,
                "contribution": contribution,
                "inference_method": inference_method,
                "p_value": hit["p_value"],
                "source": hit["source_study"]
            })

    # Calculate coverage
    coverage = matched_snps / max(total_snps, 1)

    # Adjust score based on coverage
    if coverage > 0:
        adjusted_prs = prs_sum / coverage
    else:
        adjusted_prs = 0.0

    # Determine confidence level
    if coverage >= MIN_COVERAGE_HIGH_CONF:
        confidence = "high"
    elif coverage >= MIN_COVERAGE_MODERATE_CONF:
        confidence = "moderate"
    else:
        confidence = "low"

    # Calculate percentile (simplified - would use population distribution in practice)
    # Using standard normal approximation
    z_score = adjusted_prs / 0.5  # Assuming SD of 0.5
    percentile = int(100 * (1 + math.erf(z_score / math.sqrt(2))) / 2)
    percentile = max(1, min(99, percentile))  # Bound between 1-99

    # Risk category
    if percentile >= 95:
        risk_category = "high"
    elif percentile >= 80:
        risk_category = "moderate"
    else:
        risk_category = "low"

    return {
        "cancer_type": cancer_type,
        "raw_score": prs_sum,
        "adjusted_score": adjusted_prs,
        "coverage": coverage,
        "matched_snps": matched_snps,
        "total_snps": total_snps,
        "confidence": confidence,
        "percentile": percentile,
        "risk_category": risk_category,
        "population": population,
        "contributing_snps": contributing_snps
    }


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate Polygenic Risk Scores for multiple cancer types.

    Updates state with:
    - prs_results: comprehensive PRS analysis for each cancer type
    - prs_summary: overall PRS assessment
    """
    logger.info("Starting PRS calculation")
    # Note: Don't set current_node to avoid concurrent updates during parallel execution

    try:
        # Ensure PRS database exists
        create_prs_database()

        filtered_variants = state["filtered_variants"]

        # Determine patient population if available
        patient_population = state.get("patient_population", "EUR")

        # Calculate PRS for each cancer type
        prs_results = {}
        high_risk_cancers = []

        for cancer_type in SUPPORTED_CANCERS:
            logger.info(f"Calculating PRS for {cancer_type}")

            result = calculate_cancer_specific_prs(
                filtered_variants,
                cancer_type,
                patient_population
            )

            prs_results[cancer_type] = result

            # Track high-risk cancers
            if result["risk_category"] == "high":
                high_risk_cancers.append(cancer_type)

            # Log summary
            logger.info(
                f"  {cancer_type} PRS: {result['raw_score']:.3f} "
                f"(coverage: {result['coverage']:.1%}, "
                f"percentile: {result['percentile']}%, "
                f"confidence: {result['confidence']})"
            )

        # Create overall summary
        prs_summary = {
            "high_risk_cancers": high_risk_cancers,
            "primary_concern": high_risk_cancers[0] if high_risk_cancers else None,
            "overall_confidence": "high" if all(
                r["confidence"] == "high" for r in prs_results.values()
            ) else "moderate" if any(
                r["confidence"] == "high" for r in prs_results.values()
            ) else "low",
            "limitations": []
        }

        # Add limitations
        for cancer_type, result in prs_results.items():
            if result["confidence"] == "low":
                prs_summary["limitations"].append(
                    f"Low SNP coverage for {cancer_type} PRS ({result['coverage']:.1%})"
                )

        if patient_population != "EUR":
            prs_summary["limitations"].append(
                f"PRS may be less accurate for {patient_population} population"
            )

        # Log overall summary
        logger.info(f"PRS calculation complete:")
        logger.info(f"  High-risk cancers: {', '.join(high_risk_cancers) or 'None'}")
        logger.info(f"  Overall confidence: {prs_summary['overall_confidence']}")

        if prs_summary["limitations"]:
            logger.warning("  Limitations:")
            for limitation in prs_summary["limitations"]:
                logger.warning(f"    - {limitation}")

        # Note: completed_nodes is updated by the merge function for parallel nodes

        # Return only the keys this node updates
        return {
            "prs_results": prs_results,
            "prs_summary": prs_summary
        }

    except Exception as e:
        logger.error(f"PRS calculation failed: {str(e)}")

        # Return error state updates
        return {
            "prs_results": {},
            "prs_summary": {
                "error": str(e),
                "overall_confidence": "failed"
            },
            "errors": [{
                "node": "prs_calculator",
                "error": str(e),
                "timestamp": datetime.now()
            }]
        }
