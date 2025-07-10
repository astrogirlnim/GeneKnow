"""
Gene/Pathway Burden Model for GeneKnow pipeline.
Counts damaging variants in key cancer pathways and calculates pathway-specific burden scores.
Integrates with CADD, ClinVar, and population frequency data to identify high-impact variants.
"""

import logging
from typing import Dict, Any, List, Set
from datetime import datetime
from collections import defaultdict

logger = logging.getLogger(__name__)

# Cancer-relevant pathways and their associated genes
CANCER_PATHWAYS = {
    "dna_repair": {
        "genes": ["BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "RAD51", "PALB2", "RAD50", "RAD52", "XRCC1", "XRCC2", "XRCC3"],
        "description": "DNA damage response and homologous recombination",
        "cancer_relevance": 0.9,
        "weight": 1.2
    },
    "mismatch_repair": {
        "genes": ["MLH1", "MSH2", "MSH6", "PMS2", "MSH3", "EPCAM"],
        "description": "DNA mismatch repair pathway",
        "cancer_relevance": 0.95,
        "weight": 1.3
    },
    "tumor_suppressors": {
        "genes": ["TP53", "RB1", "APC", "PTEN", "VHL", "NF1", "NF2", "CDKN2A", "CDKN2B", "STK11"],
        "description": "Classical tumor suppressor genes",
        "cancer_relevance": 0.95,
        "weight": 1.3
    },
    "oncogenes": {
        "genes": ["KRAS", "EGFR", "BRAF", "PIK3CA", "MYC", "ALK", "RET", "MET", "HER2", "ERBB2", "HRAS", "NRAS"],
        "description": "Activated oncogenes",
        "cancer_relevance": 0.85,
        "weight": 1.1
    },
    "cell_cycle": {
        "genes": ["CDKN2A", "CDK4", "CCND1", "RB1", "CDKN1A", "CDKN1B", "E2F1", "MDM2", "MDM4"],
        "description": "Cell cycle regulation",
        "cancer_relevance": 0.8,
        "weight": 1.1
    },
    "chromatin_remodeling": {
        "genes": ["ARID1A", "SMARCA4", "SMARCB1", "EZH2", "SETD2", "KMT2D", "CREBBP", "EP300", "ASXL1"],
        "description": "Chromatin modification and remodeling",
        "cancer_relevance": 0.75,
        "weight": 1.0
    },
    "apoptosis": {
        "genes": ["TP53", "BCL2", "BAX", "APAF1", "CASP8", "CASP9", "PUMA", "BAK1", "BID", "NOXA"],
        "description": "Programmed cell death pathway",
        "cancer_relevance": 0.8,
        "weight": 1.0
    },
    "wnt_signaling": {
        "genes": ["APC", "CTNNB1", "AXIN1", "AXIN2", "GSK3B", "WNT1", "FZD1", "LRP5", "LRP6"],
        "description": "Wnt/beta-catenin signaling",
        "cancer_relevance": 0.7,
        "weight": 0.9
    },
    "pi3k_akt": {
        "genes": ["PIK3CA", "PIK3CB", "PIK3R1", "AKT1", "AKT2", "PTEN", "TSC1", "TSC2", "MTOR"],
        "description": "PI3K/AKT/mTOR signaling pathway",
        "cancer_relevance": 0.8,
        "weight": 1.0
    },
    "ras_mapk": {
        "genes": ["KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3"],
        "description": "RAS/MAPK signaling pathway",
        "cancer_relevance": 0.85,
        "weight": 1.1
    }
}

# Variant impact classifications for burden calculation
DAMAGING_CONSEQUENCES = {
    "high_impact": [
        "frameshift_variant", "stop_gained", "stop_lost", "start_lost", 
        "splice_acceptor_variant", "splice_donor_variant", "nonsense_mutation",
        "frameshift_deletion", "frameshift_insertion"
    ],
    "moderate_impact": [
        "missense_variant", "missense_mutation", "inframe_deletion", "inframe_insertion",
        "protein_altering_variant", "splice_region_variant"
    ],
    "low_impact": [
        "synonymous_variant", "stop_retained_variant", "5_prime_utr_variant",
        "3_prime_utr_variant"
    ]
}

def is_damaging_variant(variant: Dict[str, Any]) -> Dict[str, Any]:
    """
    Determine if a variant is damaging based on multiple criteria.
    
    Returns dict with:
    - is_damaging: boolean
    - damage_score: float (0-1)
    - damage_reasons: list of reasons
    """
    damage_reasons = []
    damage_score = 0.0
    
    # Check CADD score
    cadd_phred = variant.get("cadd_phred", 0)
    if cadd_phred > 25:
        damage_score += 0.4
        damage_reasons.append(f"High CADD score ({cadd_phred:.1f})")
    elif cadd_phred > 15:
        damage_score += 0.2
        damage_reasons.append(f"Moderate CADD score ({cadd_phred:.1f})")
    
    # Check ClinVar annotation
    clinical_sig = variant.get("clinical_significance", "").lower()
    if "pathogenic" in clinical_sig and "likely" not in clinical_sig:
        damage_score += 0.5
        damage_reasons.append("ClinVar pathogenic")
    elif "likely_pathogenic" in clinical_sig:
        damage_score += 0.3
        damage_reasons.append("ClinVar likely pathogenic")
    
    # Check consequence type
    consequence = variant.get("consequence", "").lower()
    impact = variant.get("impact", "").lower()
    variant_class = variant.get("variant_classification", "").lower()
    
    for cons in DAMAGING_CONSEQUENCES["high_impact"]:
        if cons in consequence or cons in impact or cons in variant_class:
            damage_score += 0.4
            damage_reasons.append(f"High impact consequence ({cons})")
            break
    else:
        for cons in DAMAGING_CONSEQUENCES["moderate_impact"]:
            if cons in consequence or cons in impact or cons in variant_class:
                damage_score += 0.2
                damage_reasons.append(f"Moderate impact consequence ({cons})")
                break
    
    # Check allele frequency (rare variants more likely damaging)
    af = variant.get("allele_frequency", 0.5)
    pop_freq = variant.get("population_frequency", af)
    
    if pop_freq < 0.0001:  # Ultra-rare
        damage_score += 0.2
        damage_reasons.append(f"Ultra-rare variant (AF={pop_freq:.6f})")
    elif pop_freq < 0.001:  # Very rare
        damage_score += 0.1
        damage_reasons.append(f"Very rare variant (AF={pop_freq:.5f})")
    elif pop_freq > 0.01:  # Common variants less likely damaging
        damage_score *= 0.5
        damage_reasons.append(f"Common variant reduces damage score (AF={pop_freq:.3f})")
    
    # Check if variant is in a cancer gene (from CADD scoring)
    gene = variant.get("gene", "")
    if any(gene in pathway_data["genes"] for pathway_data in CANCER_PATHWAYS.values()):
        damage_score += 0.1
        damage_reasons.append(f"Variant in cancer gene ({gene})")
    
    # Normalize damage score to 0-1 range
    damage_score = min(damage_score, 1.0)
    
    # Consider damaging if score > 0.3 or has strong evidence
    is_damaging = damage_score > 0.3 or any(reason.startswith("ClinVar pathogenic") for reason in damage_reasons)
    
    return {
        "is_damaging": is_damaging,
        "damage_score": damage_score,
        "damage_reasons": damage_reasons
    }

def calculate_pathway_burden(pathway_name: str, pathway_data: Dict[str, Any], 
                           variants: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Calculate burden score for a specific pathway."""
    
    pathway_genes = set(pathway_data["genes"])
    pathway_variants = []
    damaging_variants = []
    gene_variant_counts = defaultdict(int)
    gene_damaging_counts = defaultdict(int)
    
    # Collect variants in this pathway
    for variant in variants:
        gene = variant.get("gene", "")
        if gene in pathway_genes:
            pathway_variants.append(variant)
            gene_variant_counts[gene] += 1
            
            # Check if variant is damaging
            damage_assessment = is_damaging_variant(variant)
            variant["pathway_damage_assessment"] = damage_assessment
            
            if damage_assessment["is_damaging"]:
                damaging_variants.append(variant)
                gene_damaging_counts[gene] += 1
    
    # Calculate burden metrics
    total_variants = len(pathway_variants)
    total_damaging = len(damaging_variants)
    genes_with_variants = len(gene_variant_counts)
    genes_with_damaging = len(gene_damaging_counts)
    pathway_size = len(pathway_genes)
    
    # Calculate raw burden score
    if total_variants > 0:
        raw_burden = total_damaging / total_variants
    else:
        raw_burden = 0.0
    
    # Apply pathway-specific weighting
    pathway_weight = pathway_data.get("weight", 1.0)
    cancer_relevance = pathway_data.get("cancer_relevance", 0.8)
    
    # Normalize by pathway size (larger pathways get adjusted scores)
    size_adjustment = min(pathway_size / 10.0, 1.0)  # Normalize to max 10 genes
    
    # Calculate final burden score
    burden_score = raw_burden * pathway_weight * cancer_relevance * size_adjustment
    
    # Determine risk level
    if burden_score > 0.6:
        risk_level = "high"
    elif burden_score > 0.3:
        risk_level = "moderate"
    elif burden_score > 0.1:
        risk_level = "low"
    else:
        risk_level = "minimal"
    
    # Find top damaging variant
    top_variant = None
    max_damage_score = 0
    for variant in damaging_variants:
        damage_score = variant.get("pathway_damage_assessment", {}).get("damage_score", 0)
        if damage_score > max_damage_score:
            max_damage_score = damage_score
            top_variant = {
                "variant_id": variant.get("variant_id"),
                "gene": variant.get("gene"),
                "consequence": variant.get("consequence"),
                "damage_score": damage_score,
                "cadd_phred": variant.get("cadd_phred"),
                "clinical_significance": variant.get("clinical_significance")
            }
    
    # Multi-hit genes (genes with multiple damaging variants)
    multi_hit_genes = [gene for gene, count in gene_damaging_counts.items() if count > 1]
    
    return {
        "pathway_name": pathway_name,
        "description": pathway_data["description"],
        "total_variants": total_variants,
        "damaging_variants": total_damaging,
        "burden_score": round(burden_score, 3),
        "raw_burden": round(raw_burden, 3),
        "risk_level": risk_level,
        "genes_in_pathway": pathway_size,
        "genes_with_variants": genes_with_variants,
        "genes_with_damaging": genes_with_damaging,
        "contributing_genes": list(gene_variant_counts.keys()),
        "damaging_genes": list(gene_damaging_counts.keys()),
        "multi_hit_genes": multi_hit_genes,
        "top_variant": top_variant,
        "gene_variant_counts": dict(gene_variant_counts),
        "gene_damaging_counts": dict(gene_damaging_counts),
        "pathway_weight": pathway_weight,
        "cancer_relevance": cancer_relevance
    }

def assess_overall_burden(pathway_results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Assess overall pathway burden across all pathways."""
    
    high_burden_pathways = []
    moderate_burden_pathways = []
    total_damaging_variants = 0
    total_variants = 0
    
    # Aggregate statistics
    for pathway_name, results in pathway_results.items():
        risk_level = results["risk_level"]
        if risk_level == "high":
            high_burden_pathways.append(pathway_name)
        elif risk_level == "moderate":
            moderate_burden_pathways.append(pathway_name)
        
        total_damaging_variants += results["damaging_variants"]
        total_variants += results["total_variants"]
    
    # Calculate overall burden score (weighted average)
    if total_variants > 0:
        overall_burden_score = sum(
            results["burden_score"] * results["total_variants"] 
            for results in pathway_results.values()
        ) / total_variants
    else:
        overall_burden_score = 0.0
    
    # Determine primary concern
    primary_concern = None
    max_burden_score = 0
    for pathway_name, results in pathway_results.items():
        if results["burden_score"] > max_burden_score:
            max_burden_score = results["burden_score"]
            primary_concern = pathway_name
    
    # Identify pathway interactions (genes appearing in multiple pathways)
    gene_pathway_map = defaultdict(set)
    for pathway_name, results in pathway_results.items():
        for gene in results["contributing_genes"]:
            gene_pathway_map[gene].add(pathway_name)
    
    multi_pathway_genes = {gene: list(pathways) for gene, pathways in gene_pathway_map.items() if len(pathways) > 1}
    
    return {
        "overall_burden_score": round(overall_burden_score, 3),
        "high_burden_pathways": high_burden_pathways,
        "moderate_burden_pathways": moderate_burden_pathways,
        "primary_concern": primary_concern,
        "total_damaging_variants": total_damaging_variants,
        "total_variants": total_variants,
        "pathways_analyzed": len(pathway_results),
        "multi_pathway_genes": multi_pathway_genes,
        "pathway_crosstalk": len(multi_pathway_genes) > 0
    }

def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate pathway burden scores for all cancer-relevant pathways.
    
    Updates state with:
    - pathway_burden_results: detailed analysis for each pathway
    - pathway_burden_summary: overall burden assessment
    """
    logger.info("Starting pathway burden analysis")
    # Note: Don't set current_node to avoid concurrent updates during parallel execution
    
    try:
        filtered_variants = state["filtered_variants"]
        
        logger.info(f"Analyzing pathway burden for {len(filtered_variants)} variants across {len(CANCER_PATHWAYS)} pathways")
        
        # Calculate burden for each pathway
        pathway_results = {}
        
        for pathway_name, pathway_data in CANCER_PATHWAYS.items():
            logger.info(f"Analyzing {pathway_name} pathway ({len(pathway_data['genes'])} genes)")
            
            burden_result = calculate_pathway_burden(pathway_name, pathway_data, filtered_variants)
            pathway_results[pathway_name] = burden_result
            
            # Log significant findings
            if burden_result["risk_level"] in ["high", "moderate"]:
                logger.warning(f"🔴 {burden_result['risk_level'].upper()} burden in {pathway_name}: "
                             f"{burden_result['damaging_variants']}/{burden_result['total_variants']} variants "
                             f"(score: {burden_result['burden_score']:.3f})")
                
                if burden_result["multi_hit_genes"]:
                    logger.warning(f"  Multi-hit genes: {', '.join(burden_result['multi_hit_genes'])}")
                
                if burden_result["top_variant"]:
                    top = burden_result["top_variant"]
                    logger.warning(f"  Top variant: {top['variant_id']} in {top['gene']} "
                                 f"(damage score: {top['damage_score']:.2f})")
            
            elif burden_result["total_variants"] > 0:
                logger.info(f"✅ {pathway_name}: {burden_result['damaging_variants']}/{burden_result['total_variants']} variants "
                           f"(score: {burden_result['burden_score']:.3f})")
        
        # Calculate overall burden assessment
        burden_summary = assess_overall_burden(pathway_results)
        
        # Log comprehensive summary
        logger.info("=" * 60)
        logger.info("Pathway Burden Analysis Summary:")
        logger.info(f"  Overall burden score: {burden_summary['overall_burden_score']:.3f}")
        logger.info(f"  High burden pathways: {', '.join(burden_summary['high_burden_pathways']) or 'None'}")
        logger.info(f"  Moderate burden pathways: {', '.join(burden_summary['moderate_burden_pathways']) or 'None'}")
        logger.info(f"  Primary concern: {burden_summary['primary_concern'] or 'None'}")
        logger.info(f"  Total damaging variants: {burden_summary['total_damaging_variants']}")
        logger.info(f"  Pathway crosstalk detected: {burden_summary['pathway_crosstalk']}")
        
        if burden_summary["multi_pathway_genes"]:
            logger.info(f"  Multi-pathway genes: {len(burden_summary['multi_pathway_genes'])}")
            for gene, pathways in burden_summary["multi_pathway_genes"].items():
                logger.info(f"    {gene}: {', '.join(pathways)}")
        
        logger.info("=" * 60)
        
        # Note: Don't append to completed_nodes to avoid concurrent updates
        # The merge node will handle tracking completion
        
        # IMPORTANT: Return the modified filtered_variants with pathway_damage_assessment added
        # Since we modified variants in-place during calculate_pathway_burden, we need to
        # return them to ensure the modifications are preserved after parallel execution
        # 
        # Also create a separate key for pathway-enriched variants to avoid conflicts
        # with other parallel nodes that might also be modifying filtered_variants
        return {
            "pathway_burden_results": pathway_results,
            "pathway_burden_summary": burden_summary,
            "pathway_enriched_variants": filtered_variants  # Use a unique key to avoid conflicts
        }
        
    except Exception as e:
        logger.error(f"Pathway burden analysis failed: {str(e)}")
        # Still return empty results so the merge node can detect completion
        return {
            "pathway_burden_results": {},
            "pathway_burden_summary": {
                "error": str(e),
                "overall_burden_score": 0.0,
                "high_burden_pathways": [],
                "moderate_burden_pathways": [],
                "primary_concern": None,
                "total_damaging_variants": 0,
                "total_variants": 0,
                "pathways_analyzed": 0,
                "multi_pathway_genes": {},
                "pathway_crosstalk": False
            },
            "pathway_enriched_variants": state.get("filtered_variants", []),  # Return unchanged variants
            "errors": state.get("errors", []) + [{
                "node": "pathway_burden",
                "error": str(e),
                "timestamp": datetime.now()
            }]
        } 