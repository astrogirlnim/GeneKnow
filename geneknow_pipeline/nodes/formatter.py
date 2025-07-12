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
    # Note: Don't set current_node to avoid concurrent updates

    try:
        # Extract high-risk findings if risk scores are available
        high_risk_findings = []
        risk_threshold = state.get("risk_threshold_percentage", 50.0)

        risk_scores = state.get("risk_scores", {})
        risk_genes = state.get("risk_genes", {})

        # Only process risk findings if risk assessment was performed
        if risk_scores:
            for cancer_type, risk_score in risk_scores.items():
                if risk_score >= risk_threshold:
                    high_risk_findings.append(
                        {
                            "cancer_type": cancer_type,
                            "risk_percentage": risk_score,
                            "affected_genes": risk_genes.get(cancer_type, []),
                        }
                    )

        # Format variant summary with TCGA matches and transformations
        variant_summary = []

        # Use variant_details if available (from variant_transformer), otherwise use filtered_variants
        variants_to_format = state.get("variant_details") or state.get("filtered_variants", [])

        for variant in variants_to_format:
            variant_info = {
                "gene": variant.get("gene", "Unknown"),
                "variant": variant.get("variant_id", "Unknown"),
                "variant_id": variant.get("variant_id", "Unknown"),
                "consequence": variant.get("consequence", variant.get("variant_classification", "unknown")),
                "mutation_type": variant.get("mutation_type", "snv"),  # From mutation_classifier
                "hgvs_c": variant.get("hgvs_c", ""),
                "hgvs_p": variant.get("hgvs_p", variant.get("protein_change", "")),
                "protein_change": variant.get("protein_change", ""),  # From variant_transformer
                "functional_impact": variant.get("functional_impact", "Unknown"),
                "quality_metrics": {
                    "quality": variant.get("quality", variant.get("qual", 0)),
                    "depth": variant.get("depth", 0),
                    "allele_freq": variant.get("allele_freq", 0),
                },
                "tcga_matches": {},
            }

            # Add transformation data if available
            if "transformation" in variant:
                variant_info["transformation"] = variant["transformation"]

            # Add clinical significance if available (from MAF)
            if "clinical_significance" in variant:
                variant_info["clinical_significance"] = variant["clinical_significance"]

            # Add CADD scores if available
            if "cadd_phred" in variant:
                variant_info["cadd_scores"] = {
                    "phred": variant.get("cadd_phred"),
                    "raw": variant.get("cadd_raw"),
                    "risk_weight": variant.get("cadd_risk_weight"),
                }

            # Add TCGA match info
            tcga_matches = state.get("tcga_matches", {})
            for cancer_type, matches in tcga_matches.items():
                if variant.get("variant_id") in matches:
                    match_info = matches[variant["variant_id"]]
                    variant_info["tcga_matches"][cancer_type] = {
                        "frequency": match_info.get("frequency", 0),
                        "frequency_percent": match_info.get("frequency", 0) * 100,
                        "patient_fraction": f"{match_info['sample_count']}/{match_info['total_samples']}",
                    }

            variant_summary.append(variant_info)

        # Build structured JSON
        # Calculate processing time if not set
        processing_time = None
        if state.get("pipeline_start_time"):
            processing_time = (datetime.now() - state["pipeline_start_time"]).total_seconds()

        # Build TCGA summary
        tcga_matches = state.get("tcga_matches", {})
        tcga_summary = {
            "cancer_types_analyzed": list(tcga_matches.keys()),
            "cohort_sizes": state.get("tcga_cohort_sizes", {}),
            "variants_with_tcga_data": len([v for v in variant_summary if v["tcga_matches"]]),
        }

        # Include CADD statistics if available
        cadd_stats = state.get("cadd_stats", {})
        if cadd_stats and "variants_scored" in cadd_stats:
            cadd_summary = {
                "variants_scored": cadd_stats.get("variants_scored", 0),
                "mean_phred": cadd_stats.get("mean_phred", 0),
                "max_phred": cadd_stats.get("max_phred", 0),
                "high_impact_variants": cadd_stats.get("variants_gt20", 0),
                "cancer_gene_variants": cadd_stats.get("variants_in_cancer_genes", 0),
            }
        else:
            cadd_summary = None

        # Include SHAP validation results if available
        shap_validation = None
        if state.get("shap_validation_status"):
            shap_validation = {
                "status": state.get("shap_validation_status", "SKIPPED"),
                "reasons": state.get("shap_validation_reasons", []),
                "top_contributors": state.get("shap_top_contributors", []),
                "feature_importance": state.get("shap_feature_importance", {}),
                "details": state.get("shap_validation_details", {}),
            }

        # Build comprehensive summary including mutation types
        summary = {
            "total_variants_found": state["variant_count"],
            "variants_passed_qc": len(state["filtered_variants"]),
            "high_risk_findings": len(high_risk_findings),
        }

        # Add mutation type distribution if available
        if state.get("mutation_type_distribution"):
            summary["mutation_types"] = state["mutation_type_distribution"]

        structured_json = {
            "report_metadata": {
                "pipeline_version": "1.0.0",
                "processing_time_seconds": processing_time,
                "generated_at": datetime.now().isoformat(),
                "language": state.get("language", "en"),
            },
            "patient_data": {"file_type": state["file_type"], "file_metadata": state["file_metadata"]},
            "summary": summary,
            "risk_assessment": (
                {"scores": risk_scores, "high_risk_findings": high_risk_findings, "risk_genes": risk_genes}
                if risk_scores
                else None
            ),
            "variant_details": variant_summary,
            "quality_control": state["file_metadata"].get("qc_stats", {}),
            "tcga_summary": tcga_summary,
            "cadd_summary": cadd_summary,
            "shap_validation": shap_validation,
            "metrics": state.get("metrics", {}),
            "metrics_summary": state.get("metrics_summary", {}),
            "warnings": state["warnings"],
            # Add all new analysis results
            "mutation_signatures": state.get("mutational_signatures", []),
            "structural_variants": state.get("structural_variants", []),
            "copy_number_variants": state.get("copy_number_variants", []),
            "pathway_analysis": state.get("pathway_analysis"),
            "gene_interactions": state.get("gene_interactions", []),
            "gene_network_analysis": state.get("gene_network_analysis"),
            "survival_analysis": state.get("survival_analysis"),
            "clinical_recommendations": state.get("clinical_recommendations"),
        }

        # Log metrics status to help debug
        metrics_available = "metrics" in state and state["metrics"]
        metrics_summary_available = "metrics_summary" in state and state["metrics_summary"]
        logger.info(
            f"Formatting complete - Metrics available: {metrics_available}, Summary available: {metrics_summary_available}"
        )

        logger.info("Formatting complete")

        # Return only the keys this node updates
        return {"structured_json": structured_json}

    except Exception as e:
        logger.error(f"Formatting failed: {str(e)}")
        # Return error state updates
        return {
            "pipeline_status": "failed",
            "errors": [{"node": "formatter", "error": str(e), "timestamp": datetime.now()}],
        }
