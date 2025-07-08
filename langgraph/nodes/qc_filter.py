"""
Quality control filtering node.
Filters variants based on quality metrics.
"""
import logging
from datetime import datetime
from typing import Dict, Any, List

logger = logging.getLogger(__name__)

# QC Thresholds (from documentation)
QUAL_THRESHOLD = 30.0
DEPTH_THRESHOLD = 10
ALLELE_FREQ_THRESHOLD = 0.01


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Filter variants based on quality control metrics.
    
    Filters:
    - QUAL >= 30
    - Depth >= 10  
    - Allele frequency >= 0.01
    
    Updates state with:
    - filtered_variants: variants passing QC
    """
    logger.info("Starting QC filtering")
    state["current_node"] = "qc_filter"
    
    try:
        raw_variants = state["raw_variants"]
        filtered_variants = []
        
        for variant in raw_variants:
            # Apply QC filters
            passed_qc = True
            qc_failures = []
            
            # Quality score filter
            if variant.get("quality", 0) < QUAL_THRESHOLD:
                passed_qc = False
                qc_failures.append(f"Low quality: {variant.get('quality', 0)}")
            
            # Depth filter
            if variant.get("depth", 0) < DEPTH_THRESHOLD:
                passed_qc = False
                qc_failures.append(f"Low depth: {variant.get('depth', 0)}")
            
            # Allele frequency filter
            if variant.get("allele_freq", 0) < ALLELE_FREQ_THRESHOLD:
                passed_qc = False
                qc_failures.append(f"Low allele freq: {variant.get('allele_freq', 0)}")
            
            if passed_qc:
                filtered_variants.append(variant)
            else:
                # Log filtered variants for debugging
                logger.debug(f"Variant filtered: {variant['variant_id']} - {qc_failures}")
        
        # Update state
        state["filtered_variants"] = filtered_variants
        
        # Calculate filtering stats
        total_variants = len(raw_variants)
        passed_variants = len(filtered_variants)
        filter_rate = (total_variants - passed_variants) / total_variants * 100 if total_variants > 0 else 0
        
        logger.info(f"QC filtering complete: {passed_variants}/{total_variants} passed ({filter_rate:.1f}% filtered)")
        
        # Add filtering stats to metadata
        state["file_metadata"]["qc_stats"] = {
            "total_variants": total_variants,
            "passed_qc": passed_variants,
            "failed_qc": total_variants - passed_variants,
            "filter_rate_percent": filter_rate,
            "thresholds": {
                "qual": QUAL_THRESHOLD,
                "depth": DEPTH_THRESHOLD,
                "allele_freq": ALLELE_FREQ_THRESHOLD
            }
        }
        
        state["completed_nodes"].append("qc_filter")
        
    except Exception as e:
        logger.error(f"QC filtering failed: {str(e)}")
        state["errors"].append({
            "node": "qc_filter",
            "error": str(e),
            "timestamp": datetime.now()
        })
        state["pipeline_status"] = "failed"
    
    return state 