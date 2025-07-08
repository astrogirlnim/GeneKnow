"""
Formatter node.
Structures all results into JSON format for frontend consumption.
"""
import logging
from datetime import datetime
from typing import Dict, Any

logger = logging.getLogger(__name__)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Format all pipeline results into structured JSON.
    
    Updates state with:
    - structured_json: formatted data ready for frontend/report
    """
    logger.info("Starting result formatting")
    state["current_node"] = "formatter"
    
    try:
        # Extract high-risk findings
        high_risk_findings = []
        risk_threshold = state.get("risk_threshold_percentage", 50.0)
        
        for cancer_type, risk_score in state["risk_scores"].items():
            if risk_score >= risk_threshold:
                high_risk_findings.append({
                    "cancer_type": cancer_type,
                    "risk_percentage": risk_score,
                    "affected_genes": state["risk_genes"].get(cancer_type, [])
                })
        
        # Format variant summary with TCGA matches
        variant_summary = []
        for variant in state["filtered_variants"]:
            variant_info = {
                "gene": variant.get("gene", "Unknown"),
                "variant": variant.get("variant_id", "Unknown"),
                "consequence": variant.get("consequence", variant.get("variant_classification", "unknown")),
                "hgvs_c": variant.get("hgvs_c", ""),
                "hgvs_p": variant.get("hgvs_p", variant.get("protein_change", "")),
                "quality_metrics": {
                    "quality": variant.get("quality", variant.get("qual", 0)),
                    "depth": variant.get("depth", 0),
                    "allele_freq": variant.get("allele_freq", 0)
                },
                "tcga_matches": {}
            }
            
            # Add clinical significance if available (from MAF)
            if "clinical_significance" in variant:
                variant_info["clinical_significance"] = variant["clinical_significance"]
            
            # Add TCGA match info
            for cancer_type, matches in state["tcga_matches"].items():
                if variant.get("variant_id") in matches:
                    match_info = matches[variant["variant_id"]]
                    variant_info["tcga_matches"][cancer_type] = {
                        "frequency_percent": match_info["frequency"] * 100,
                        "patient_fraction": f"{match_info['sample_count']}/{match_info['total_samples']}"
                    }
            
            variant_summary.append(variant_info)
        
        # Build structured JSON
        # Calculate processing time if not set
        processing_time = None
        if state.get("pipeline_start_time"):
            processing_time = (datetime.now() - state["pipeline_start_time"]).total_seconds()
        
        # Build TCGA summary
        tcga_summary = {
            "cancer_types_analyzed": list(state.get("tcga_matches", {}).keys()),
            "cohort_sizes": state.get("tcga_cohort_sizes", {}),
            "variants_with_tcga_data": len([v for v in variant_summary if v["tcga_matches"]])
        }
        
        structured_json = {
            "report_metadata": {
                "pipeline_version": "1.0.0",
                "processing_time_seconds": processing_time,
                "generated_at": datetime.now().isoformat(),
                "language": state.get("language", "en")
            },
            "patient_data": {
                "file_type": state["file_type"],
                "file_metadata": state["file_metadata"]
            },
            "summary": {
                "total_variants_found": state["variant_count"],
                "variants_passed_qc": len(state["filtered_variants"]),
                "high_risk_findings": len(high_risk_findings)
            },
            "risk_assessment": {
                "scores": state["risk_scores"],
                "high_risk_findings": high_risk_findings,
                "risk_genes": state["risk_genes"]
            },
            "variant_details": variant_summary,
            "quality_control": state["file_metadata"].get("qc_stats", {}),
            "tcga_summary": tcga_summary,
            "warnings": state["warnings"]
        }
        
        state["structured_json"] = structured_json
        state["completed_nodes"].append("formatter")
        
        logger.info("Formatting complete")
        
    except Exception as e:
        logger.error(f"Formatting failed: {str(e)}")
        state["errors"].append({
            "node": "formatter",
            "error": str(e),
            "timestamp": datetime.now()
        })
        state["pipeline_status"] = "failed"
    
    return state 