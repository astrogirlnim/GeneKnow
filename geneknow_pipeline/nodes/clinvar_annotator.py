"""
ClinVar annotator node for GeneKnow pipeline.
Annotates variants with clinical significance from ClinVar database.
Provides direct clinical interpretation with expert knowledge.
"""

import os
import sqlite3
import logging
from datetime import datetime
from typing import Dict, Any, List, Optional, Tuple
from collections import defaultdict

logger = logging.getLogger(__name__)

# Path to ClinVar database
CLINVAR_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "clinvar_annotations.db")

# Clinical significance mapping to risk scores
CLINICAL_SIGNIFICANCE_SCORES = {
    "pathogenic": 1.0,
    "likely_pathogenic": 0.8,
    "pathogenic/likely_pathogenic": 0.9,
    "drug_response": 0.1,
    "risk_factor": 0.3,
    "affects": 0.2,
    "association": 0.1,
    "protective": -0.1,
    "benign": 0.0,
    "likely_benign": 0.0,
    "benign/likely_benign": 0.0,
    "uncertain_significance": 0.1,
    "conflicting_interpretations_of_pathogenicity": 0.2,
    "other": 0.0,
    "not_provided": 0.0
}

# Cancer-related conditions for enhanced scoring
CANCER_CONDITIONS = [
    "breast cancer", "ovarian cancer", "colorectal cancer", "lung cancer", 
    "prostate cancer", "pancreatic cancer", "melanoma", "lymphoma", "leukemia",
    "hereditary breast and ovarian cancer", "lynch syndrome", "familial adenomatous polyposis",
    "hereditary nonpolyposis colorectal cancer", "li-fraumeni syndrome", "cowden syndrome",
    "peutz-jeghers syndrome", "hereditary diffuse gastric cancer", "multiple endocrine neoplasia",
    "von hippel-lindau disease", "neurofibromatosis", "tuberous sclerosis", "gorlin syndrome"
]

def normalize_chromosome(chrom: str) -> str:
    """Normalize chromosome format (remove 'chr' prefix)."""
    if chrom.startswith('chr'):
        return chrom[3:]
    return chrom

def create_clinvar_database():
    """Create ClinVar database with example data if it doesn't exist."""
    if os.path.exists(CLINVAR_DB_PATH):
        return
    
    logger.info("Creating ClinVar database with example data...")
    conn = sqlite3.connect(CLINVAR_DB_PATH)
    cursor = conn.cursor()
    
    # Create table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS clinvar_variants (
        variation_id INTEGER PRIMARY KEY,
        chrom TEXT,
        pos INTEGER,
        ref TEXT,
        alt TEXT,
        gene TEXT,
        clinical_significance TEXT,
        review_status TEXT,
        condition TEXT,
        molecular_consequence TEXT,
        allele_id INTEGER,
        rs_id TEXT,
        rcv_accession TEXT,
        last_evaluated TEXT,
        guidelines TEXT,
        cancer_related INTEGER DEFAULT 0,
        pathogenicity_score REAL DEFAULT 0.0
    )
    """)
    
    # Create indexes for fast lookups
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_position ON clinvar_variants(chrom, pos, ref, alt)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene ON clinvar_variants(gene)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_clinical_sig ON clinvar_variants(clinical_significance)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_cancer_related ON clinvar_variants(cancer_related)")
    
    # Insert example ClinVar data (based on known pathogenic variants)
    example_variants = [
        # BRCA1 pathogenic variants
        (1, "17", 43045677, "C", "T", "BRCA1", "pathogenic", "criteria_provided_single_submitter", 
         "Hereditary breast and ovarian cancer syndrome", "missense_variant", 12345, "rs80356868", 
         "RCV000077282", "2020-01-15", "ACMG", 1, 1.0),
        
        (2, "17", 43057135, "G", "A", "BRCA1", "pathogenic", "criteria_provided_multiple_submitters", 
         "Hereditary breast and ovarian cancer syndrome", "nonsense", 12346, "rs80357906", 
         "RCV000077283", "2020-02-20", "ACMG", 1, 1.0),
        
        # BRCA2 pathogenic variants
        (3, "13", 32363533, "T", "C", "BRCA2", "pathogenic", "practice_guideline", 
         "Hereditary breast and ovarian cancer syndrome", "missense_variant", 12347, "rs80359550", 
         "RCV000077284", "2020-03-10", "ACMG", 1, 1.0),
        
        (4, "13", 32379617, "G", "T", "BRCA2", "likely_pathogenic", "criteria_provided_multiple_submitters", 
         "Hereditary breast and ovarian cancer syndrome", "splice_acceptor_variant", 12348, "rs80359876", 
         "RCV000077285", "2020-04-05", "ACMG", 1, 0.8),
        
        # TP53 pathogenic variants
        (5, "17", 7673803, "G", "A", "TP53", "pathogenic", "criteria_provided_multiple_submitters", 
         "Li-Fraumeni syndrome", "missense_variant", 12349, "rs28934578", 
         "RCV000077286", "2020-05-12", "ACMG", 1, 1.0),
        
        (6, "17", 7674894, "C", "T", "TP53", "pathogenic", "practice_guideline", 
         "Li-Fraumeni syndrome", "nonsense", 12350, "rs121912651", 
         "RCV000077287", "2020-06-18", "ACMG", 1, 1.0),
        
        # MLH1 pathogenic variants (Lynch syndrome)
        (7, "3", 37092337, "A", "G", "MLH1", "pathogenic", "criteria_provided_multiple_submitters", 
         "Lynch syndrome", "missense_variant", 12351, "rs63750447", 
         "RCV000077288", "2020-07-22", "ACMG", 1, 1.0),
        
        # MSH2 pathogenic variants (Lynch syndrome)
        (8, "2", 47806747, "G", "A", "MSH2", "likely_pathogenic", "criteria_provided_single_submitter", 
         "Lynch syndrome", "splice_donor_variant", 12352, "rs63750323", 
         "RCV000077289", "2020-08-14", "ACMG", 1, 0.8),
        
        # APC pathogenic variants (FAP)
        (9, "5", 112839514, "C", "T", "APC", "pathogenic", "criteria_provided_multiple_submitters", 
         "Familial adenomatous polyposis", "nonsense", 12353, "rs137854420", 
         "RCV000077290", "2020-09-30", "ACMG", 1, 1.0),
        
        # PALB2 pathogenic variants
        (10, "16", 23652178, "G", "A", "PALB2", "pathogenic", "criteria_provided_multiple_submitters", 
         "Hereditary breast and ovarian cancer syndrome", "nonsense", 12354, "rs180177143", 
         "RCV000077291", "2020-10-25", "ACMG", 1, 1.0),
        
        # Benign variants for contrast
        (11, "17", 43045678, "A", "G", "BRCA1", "benign", "criteria_provided_multiple_submitters", 
         "not_provided", "synonymous_variant", 12355, "rs1799950", 
         "RCV000077292", "2020-11-12", "ACMG", 0, 0.0),
        
        (12, "13", 32363534, "C", "T", "BRCA2", "likely_benign", "criteria_provided_single_submitter", 
         "not_provided", "synonymous_variant", 12356, "rs766173", 
         "RCV000077293", "2020-12-05", "ACMG", 0, 0.0),
        
        # Uncertain significance variants
        (13, "17", 43045679, "T", "A", "BRCA1", "uncertain_significance", "criteria_provided_single_submitter", 
         "Hereditary breast and ovarian cancer syndrome", "missense_variant", 12357, "rs1799966", 
         "RCV000077294", "2021-01-20", "ACMG", 1, 0.1),
        
        # Drug response variants
        (14, "19", 15990431, "C", "T", "CYP2D6", "drug_response", "criteria_provided_multiple_submitters", 
         "drug_response", "missense_variant", 12358, "rs16947", 
         "RCV000077295", "2021-02-15", "PharmGKB", 0, 0.1),
        
        # Risk factor variants
        (15, "9", 22125504, "C", "G", "CDKN2A", "risk_factor", "criteria_provided_single_submitter", 
         "Melanoma", "missense_variant", 12359, "rs11515", 
         "RCV000077296", "2021-03-10", "ACMG", 1, 0.3)
    ]
    
    # Insert the example data
    cursor.executemany("""
    INSERT OR IGNORE INTO clinvar_variants 
    (variation_id, chrom, pos, ref, alt, gene, clinical_significance, review_status, 
     condition, molecular_consequence, allele_id, rs_id, rcv_accession, 
     last_evaluated, guidelines, cancer_related, pathogenicity_score)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, example_variants)
    
    conn.commit()
    conn.close()
    logger.info(f"Created ClinVar database with {len(example_variants)} example variants")

def query_clinvar_database(chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
    """Query ClinVar database for clinical significance."""
    if not os.path.exists(CLINVAR_DB_PATH):
        create_clinvar_database()
    
    conn = sqlite3.connect(CLINVAR_DB_PATH)
    cursor = conn.cursor()
    
    try:
        # Normalize chromosome format
        norm_chrom = normalize_chromosome(chrom)
        
        # Query for exact match
        cursor.execute("""
        SELECT variation_id, gene, clinical_significance, review_status, condition, 
               molecular_consequence, allele_id, rs_id, rcv_accession, 
               last_evaluated, guidelines, cancer_related, pathogenicity_score
        FROM clinvar_variants
        WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
        ORDER BY pathogenicity_score DESC, variation_id ASC
        LIMIT 1
        """, (norm_chrom, pos, ref, alt))
        
        row = cursor.fetchone()
        if row:
            return {
                "variation_id": row[0],
                "gene": row[1],
                "clinical_significance": row[2],
                "review_status": row[3],
                "condition": row[4],
                "molecular_consequence": row[5],
                "allele_id": row[6],
                "rs_id": row[7],
                "rcv_accession": row[8],
                "last_evaluated": row[9],
                "guidelines": row[10],
                "cancer_related": bool(row[11]),
                "pathogenicity_score": row[12],
                "found_in_clinvar": True,
                "lookup_method": "exact_match"
            }
        
        # If no exact match, try gene-based lookup for nearby variants
        cursor.execute("""
        SELECT variation_id, gene, clinical_significance, review_status, condition, 
               molecular_consequence, allele_id, rs_id, rcv_accession, 
               last_evaluated, guidelines, cancer_related, pathogenicity_score
        FROM clinvar_variants
        WHERE chrom = ? AND ABS(pos - ?) <= 100 AND gene IS NOT NULL AND gene != ''
        ORDER BY pathogenicity_score DESC, ABS(pos - ?) ASC
        LIMIT 1
        """, (norm_chrom, pos, pos))
        
        row = cursor.fetchone()
        if row:
            return {
                "variation_id": row[0],
                "gene": row[1],
                "clinical_significance": row[2],
                "review_status": row[3],
                "condition": row[4],
                "molecular_consequence": row[5],
                "allele_id": row[6],
                "rs_id": row[7],
                "rcv_accession": row[8],
                "last_evaluated": row[9],
                "guidelines": row[10],
                "cancer_related": bool(row[11]),
                "pathogenicity_score": row[12],
                "found_in_clinvar": True,
                "lookup_method": "gene_vicinity_match"
            }
        
        # No match found
        return {
            "found_in_clinvar": False,
            "lookup_method": "not_found"
        }
        
    finally:
        conn.close()

def calculate_clinical_risk_score(variant: Dict[str, Any], clinvar_data: Dict[str, Any]) -> float:
    """Calculate clinical risk score based on ClinVar annotation."""
    if not clinvar_data.get("found_in_clinvar"):
        return 0.0  # No clinical evidence
    
    clinical_sig = clinvar_data.get("clinical_significance", "").lower()
    base_score = CLINICAL_SIGNIFICANCE_SCORES.get(clinical_sig, 0.0)
    
    # Apply cancer-related bonus
    if clinvar_data.get("cancer_related"):
        base_score *= 1.2  # 20% bonus for cancer-related conditions
    
    # Apply review status modifier
    review_status = clinvar_data.get("review_status", "").lower()
    if "practice_guideline" in review_status:
        base_score *= 1.1  # 10% bonus for practice guidelines
    elif "multiple_submitters" in review_status:
        base_score *= 1.05  # 5% bonus for multiple submitters
    elif "single_submitter" in review_status:
        base_score *= 0.95  # 5% penalty for single submitter
    
    # Ensure score is within bounds
    return max(0.0, min(1.0, base_score))

def assess_clinical_significance(variant: Dict[str, Any], clinvar_data: Dict[str, Any]) -> Dict[str, Any]:
    """Assess clinical significance and provide interpretation."""
    if not clinvar_data.get("found_in_clinvar"):
        return {
            "clinical_interpretation": "no_clinical_evidence",
            "confidence": "none",
            "actionability": "none",
            "recommendation": "No clinical evidence available"
        }
    
    clinical_sig = clinvar_data.get("clinical_significance", "").lower()
    
    # Determine interpretation
    if "pathogenic" in clinical_sig and "likely" not in clinical_sig:
        interpretation = "pathogenic"
        confidence = "high"
        actionability = "high"
        recommendation = "Consider genetic counseling and enhanced screening"
    elif "likely_pathogenic" in clinical_sig:
        interpretation = "likely_pathogenic"
        confidence = "moderate"
        actionability = "moderate"
        recommendation = "Consider genetic counseling and family history evaluation"
    elif "benign" in clinical_sig:
        interpretation = "benign"
        confidence = "high"
        actionability = "low"
        recommendation = "No special action needed"
    elif "uncertain_significance" in clinical_sig:
        interpretation = "uncertain"
        confidence = "low"
        actionability = "low"
        recommendation = "Monitor for updates in clinical interpretation"
    elif "drug_response" in clinical_sig:
        interpretation = "pharmacogenomic"
        confidence = "moderate"
        actionability = "moderate"
        recommendation = "Consider pharmacogenomic testing for drug selection"
    elif "risk_factor" in clinical_sig:
        interpretation = "risk_factor"
        confidence = "moderate"
        actionability = "low"
        recommendation = "Consider as part of overall risk assessment"
    else:
        interpretation = "other"
        confidence = "low"
        actionability = "low"
        recommendation = "Consult with genetics professional"
    
    return {
        "clinical_interpretation": interpretation,
        "confidence": confidence,
        "actionability": actionability,
        "recommendation": recommendation
    }

def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Annotate variants with ClinVar clinical significance.
    
    Updates state with:
    - clinvar_annotations: clinical significance data for each variant
    - clinvar_stats: summary statistics
    - clinvar_pathogenic_variants: list of pathogenic/likely pathogenic variants
    """
    logger.info("Starting ClinVar annotation")
    # Note: Don't set current_node to avoid concurrent updates during parallel execution
    
    try:
        # Ensure ClinVar database exists
        create_clinvar_database()
        
        filtered_variants = state["filtered_variants"]
        
        # Initialize results
        clinvar_annotations = {}
        pathogenic_variants = []
        likely_pathogenic_variants = []
        benign_variants = []
        vus_variants = []
        drug_response_variants = []
        risk_factor_variants = []
        
        # Process each variant
        annotated_variants = []
        total_variants_annotated = 0
        
        for variant in filtered_variants:
            variant_id = variant.get("variant_id", f"{variant.get('chrom')}:{variant.get('pos')}")
            
            # Query ClinVar database
            clinvar_data = query_clinvar_database(
                variant["chrom"],
                variant["pos"],
                variant["ref"],
                variant["alt"]
            )
            
            # Calculate clinical risk score
            clinical_risk_score = calculate_clinical_risk_score(variant, clinvar_data)
            
            # Assess clinical significance
            clinical_assessment = assess_clinical_significance(variant, clinvar_data)
            
            # Combine data
            annotation = {
                **clinvar_data,
                "clinical_risk_score": clinical_risk_score,
                **clinical_assessment
            }
            
            clinvar_annotations[variant_id] = annotation
            
            # Categorize variants
            if clinvar_data.get("found_in_clinvar"):
                total_variants_annotated += 1
                clinical_sig = clinvar_data.get("clinical_significance", "").lower()
                
                if "pathogenic" in clinical_sig and "likely" not in clinical_sig:
                    pathogenic_variants.append({
                        "variant_id": variant_id,
                        "gene": variant.get("gene"),
                        "clinical_significance": clinical_sig,
                        "condition": clinvar_data.get("condition"),
                        "risk_score": clinical_risk_score
                    })
                elif "likely_pathogenic" in clinical_sig:
                    likely_pathogenic_variants.append({
                        "variant_id": variant_id,
                        "gene": variant.get("gene"),
                        "clinical_significance": clinical_sig,
                        "condition": clinvar_data.get("condition"),
                        "risk_score": clinical_risk_score
                    })
                elif "benign" in clinical_sig:
                    benign_variants.append({
                        "variant_id": variant_id,
                        "gene": variant.get("gene"),
                        "clinical_significance": clinical_sig
                    })
                elif "uncertain_significance" in clinical_sig:
                    vus_variants.append({
                        "variant_id": variant_id,
                        "gene": variant.get("gene"),
                        "clinical_significance": clinical_sig
                    })
                elif "drug_response" in clinical_sig:
                    drug_response_variants.append({
                        "variant_id": variant_id,
                        "gene": variant.get("gene"),
                        "clinical_significance": clinical_sig
                    })
                elif "risk_factor" in clinical_sig:
                    risk_factor_variants.append({
                        "variant_id": variant_id,
                        "gene": variant.get("gene"),
                        "clinical_significance": clinical_sig
                    })
                
                # Log significant findings
                if clinical_risk_score > 0.7:
                    logger.warning(f"🔴 High clinical risk: {variant_id} in {variant.get('gene', 'unknown')} "
                                 f"({clinical_sig}, risk score: {clinical_risk_score:.2f})")
                elif clinical_risk_score > 0.3:
                    logger.info(f"🟡 Moderate clinical risk: {variant_id} in {variant.get('gene', 'unknown')} "
                               f"({clinical_sig}, risk score: {clinical_risk_score:.2f})")
        
        # Calculate summary statistics
        clinvar_stats = {
            "total_variants": len(filtered_variants),
            "variants_annotated": total_variants_annotated,
            "annotation_rate": total_variants_annotated / len(filtered_variants) if filtered_variants else 0,
            "pathogenic_variants": len(pathogenic_variants),
            "likely_pathogenic_variants": len(likely_pathogenic_variants),
            "benign_variants": len(benign_variants),
            "vus_variants": len(vus_variants),
            "drug_response_variants": len(drug_response_variants),
            "risk_factor_variants": len(risk_factor_variants),
            "high_risk_variants": len([v for v in pathogenic_variants + likely_pathogenic_variants 
                                     if v["risk_score"] > 0.7]),
            "cancer_related_variants": len([a for a in clinvar_annotations.values() 
                                          if a.get("cancer_related")])
        }
        
        # Update file metadata
        file_metadata = state.get("file_metadata", {})
        file_metadata["clinvar_summary"] = clinvar_stats
        
        # Log summary
        logger.info(f"ClinVar annotation complete:")
        logger.info(f"  Total variants: {len(filtered_variants)}")
        logger.info(f"  Variants annotated: {total_variants_annotated}")
        logger.info(f"  Annotation rate: {total_variants_annotated/len(filtered_variants)*100:.1f}%" if filtered_variants else "0%")
        logger.info(f"  Pathogenic: {len(pathogenic_variants)}")
        logger.info(f"  Likely pathogenic: {len(likely_pathogenic_variants)}")
        logger.info(f"  Benign: {len(benign_variants)}")
        logger.info(f"  VUS: {len(vus_variants)}")
        logger.info(f"  Cancer-related: {clinvar_stats['cancer_related_variants']}")
        
        # Return only the keys this node updates
        return {
            "clinvar_annotations": clinvar_annotations,
            "clinvar_stats": clinvar_stats,
            "clinvar_pathogenic_variants": pathogenic_variants,
            "clinvar_likely_pathogenic_variants": likely_pathogenic_variants,
            "clinvar_benign_variants": benign_variants,
            "clinvar_vus_variants": vus_variants,
            "clinvar_drug_response_variants": drug_response_variants,
            "clinvar_risk_factor_variants": risk_factor_variants,
            "file_metadata": file_metadata
        }
        
    except Exception as e:
        logger.error(f"ClinVar annotation failed: {str(e)}")
        return {
            "errors": [{
                "node": "clinvar_annotator",
                "error": str(e),
                "timestamp": datetime.now()
            }]
        } 