"""
Feature vector builder node (STUB).
Collects outputs from all static risk model nodes and builds a unified feature vector.
Currently a passthrough - will be implemented when all 5 models are ready.
"""
import logging
from datetime import datetime
from typing import Dict, Any, List

logger = logging.getLogger(__name__)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Feature vector builder node (STUB).
    
    Future: Will build feature vectors from all 5 static models:
    - PRS (Polygenic Risk Scores)
    - ClinVar annotations  
    - CADD scores
    - TCGA frequency matching
    - Gene/Pathway burden
    
    For now: Analyzes available data without modifying filtered_variants
    (since merge node already combined parallel results).
    """
    logger.info("Starting feature vector builder (STUB)")
    state["current_node"] = "feature_vector_builder"
    
    try:
        # Use the merged variants from the parallel processing
        # Don't modify filtered_variants since merge node already handled this
        filtered_variants = state.get("filtered_variants", [])
        
        # Debug logging to understand the actual variant content
        logger.info(f"Feature vector builder received {len(filtered_variants)} variants")
        if filtered_variants and len(filtered_variants) > 0:
            # Log details of first few variants
            for i, variant in enumerate(filtered_variants[:3]):
                logger.info(f"Variant {i+1} sample:")
                logger.info(f"  Keys: {sorted(variant.keys())}")
                logger.info(f"  Gene: {variant.get('gene', 'N/A')}")
                logger.info(f"  CADD PHRED: {variant.get('cadd_phred', 'MISSING')}")
                logger.info(f"  TCGA relevance: {variant.get('tcga_cancer_relevance', 'MISSING')}")
                logger.info(f"  ClinVar significance: {variant.get('clinvar_clinical_significance', 'MISSING')}")
                logger.info(f"  Has pathway_damage_assessment: {'pathway_damage_assessment' in variant}")
        
        # Log what we have so far
        logger.info("=" * 60)
        logger.info("Feature Vector Builder - Current Inputs:")
        
        # CADD scores
        if "cadd_stats" in state:
            logger.info(f"✓ CADD scores available: {state['cadd_stats'].get('variants_scored', 0)} variants")
            logger.info(f"  Mean PHRED: {state['cadd_stats'].get('mean_phred', 0):.1f}")
            logger.info(f"  Max PHRED: {state['cadd_stats'].get('max_phred', 0):.1f}")
            logger.info(f"  Cancer gene variants: {state['cadd_stats'].get('variants_in_cancer_genes', 0)}")
        else:
            logger.info("✗ CADD scores not available")
        
        # PRS scores  
        if "prs_results" in state:
            logger.info(f"✓ PRS scores available: {len(state.get('prs_results', {}))} cancer types")
            prs_summary = state.get('prs_summary', {})
            if prs_summary:
                logger.info(f"  Overall PRS confidence: {prs_summary.get('overall_confidence', 'unknown')}")
        else:
            logger.info("✗ PRS scores not available")
        
        # ClinVar annotations 
        if "clinvar_stats" in state:
            logger.info(f"✓ ClinVar annotations available: {state['clinvar_stats'].get('variants_annotated', 0)} variants")
            logger.info(f"  Pathogenic variants: {state['clinvar_stats'].get('pathogenic_variants', 0)}")
            logger.info(f"  Likely pathogenic variants: {state['clinvar_stats'].get('likely_pathogenic_variants', 0)}")
            logger.info(f"  Cancer-related variants: {state['clinvar_stats'].get('cancer_related_variants', 0)}")
        else:
            logger.info("✗ ClinVar annotations not available")
        
        # TCGA frequency matches (existing, need to integrate)
        if "tcga_matches" in state:
            logger.info(f"✓ TCGA matches available: {len(state.get('tcga_matches', {}))} cancer types")
        else:
            logger.info("✗ TCGA matches not available")
        
        # Gene/Pathway burden
        if "pathway_burden_results" in state:
            pathway_results = state.get("pathway_burden_results", {})
            pathway_summary = state.get("pathway_burden_summary", {})
            logger.info(f"✓ Pathway burden analysis available: {len(pathway_results)} pathways")
            if pathway_summary:
                logger.info(f"  Overall burden score: {pathway_summary.get('overall_burden_score', 0):.3f}")
                high_burden = pathway_summary.get('high_burden_pathways', [])
                if high_burden:
                    logger.info(f"  High burden pathways: {', '.join(high_burden)}")
        else:
            logger.info("✗ Pathway burden analysis not available")
        logger.info("=" * 60)
        
        # Count variants with each type of annotation
        variants_with_cadd = sum(1 for v in filtered_variants if "cadd_phred" in v)
        variants_with_population = sum(1 for v in filtered_variants if "population_frequency" in v)
        variants_with_pathogenic = sum(1 for v in filtered_variants if v.get("is_pathogenic"))
        variants_with_tcga = sum(1 for v in filtered_variants if v.get("tcga_cancer_relevance"))
        variants_with_clinvar = sum(1 for v in filtered_variants if v.get("clinvar_clinical_significance"))
        variants_with_pathway_burden = sum(1 for v in filtered_variants if "pathway_damage_assessment" in v)
        
        logger.info(f"Variant annotation summary:")
        logger.info(f"  Total variants: {len(filtered_variants)}")
        logger.info(f"  With CADD scores: {variants_with_cadd}")
        logger.info(f"  With population data: {variants_with_population}")
        logger.info(f"  With TCGA relevance: {variants_with_tcga}")
        logger.info(f"  With ClinVar annotations: {variants_with_clinvar}")
        logger.info(f"  With pathway burden assessment: {variants_with_pathway_burden}")
        logger.info(f"  Pathogenic: {variants_with_pathogenic}")
        
        # For now, just pass through
        # In future, build actual feature vector here
        state["feature_vector"] = {
            "status": "stub",
            "message": "Feature vector builder not yet implemented",
            "available_inputs": {
                "cadd": "cadd_stats" in state,
                "prs": "prs_results" in state,
                "clinvar": "clinvar_stats" in state,
                "tcga": "tcga_matches" in state,
                "pathway_burden": "pathway_burden_results" in state
            },
            "annotation_summary": {
                "total_variants": len(filtered_variants),
                "with_cadd": variants_with_cadd,
                "with_population": variants_with_population,
                "with_tcga": variants_with_tcga,
                "with_clinvar": variants_with_clinvar,
                "with_pathway_burden": variants_with_pathway_burden,
                "pathogenic": variants_with_pathogenic
            }
        }
        
        state["completed_nodes"].append("feature_vector_builder")
        logger.info("Feature vector builder complete (passthrough)")
        
    except Exception as e:
        logger.error(f"Feature vector builder failed: {str(e)}")
        state["errors"].append({
            "node": "feature_vector_builder",
            "error": str(e),
            "timestamp": datetime.now()
        })
    
    return state 