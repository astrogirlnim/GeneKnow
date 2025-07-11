"""
Mutation Classifier Node for GeneKnow pipeline.
Classifies variants into mutation types (SNV, Indel, CNV, Structural).
"""
from typing import Dict, List
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


def classify_mutation_type(variant: Dict) -> str:
    """Classify variant based on reference and alternate alleles"""
    ref = variant.get('ref', '')
    alt = variant.get('alt', '')
    
    # Single Nucleotide Variant (SNV)
    if len(ref) == 1 and len(alt) == 1:
        return 'snv'
    
    # Insertion/Deletion (Indel)
    elif len(ref) != len(alt):
        return 'indel'
    
    # Complex variants (could be structural)
    elif len(ref) > 50 or len(alt) > 50:
        return 'structural'
    
    # Default to SNV for other cases
    return 'snv'


def process(state: Dict) -> Dict:
    """Process variants and count mutation types"""
    logger.info("Starting mutation classification")
    state["current_node"] = "mutation_classifier"
    
    try:
        variants = state.get('filtered_variants', [])
        
        mutation_counts = {
            'snv': 0,
            'indel': 0,
            'cnv': 0,
            'structural': 0
        }
        
        # Classify each variant
        classified_variants = []
        for variant in variants:
            mutation_type = classify_mutation_type(variant)
            mutation_counts[mutation_type] += 1
            
            # Add mutation type to variant
            variant = variant.copy()  # Don't modify original
            variant['mutation_type'] = mutation_type
            classified_variants.append(variant)
        
        # Update state
        state['classified_variants'] = classified_variants
        state['mutation_types'] = mutation_counts
        
        # Add to completed nodes
        completed = state.get("completed_nodes", [])
        if "mutation_classifier" not in completed:
            completed.append("mutation_classifier")
        state["completed_nodes"] = completed
        
        logger.info(f"Classified {len(variants)} variants: {mutation_counts}")
        
    except Exception as e:
        logger.error(f"Error in mutation classification: {str(e)}")
        state["errors"] = state.get("errors", []) + [{
            "node": "mutation_classifier",
            "error": str(e),
            "timestamp": datetime.now().isoformat()
        }]
    
    return state 