"""
Report writer node.
Generates structured report sections for frontend consumption.
"""
import logging
from datetime import datetime
from typing import Dict, Any, List

logger = logging.getLogger(__name__)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate structured report sections.

    Updates state with:
    - report_sections: structured report data
    - report_markdown: full markdown (optional, for PDF generation)
    - pipeline_status: 'completed'
    """
    logger.info("Starting report generation")
    # Note: Don't set current_node to avoid concurrent updates

    try:
        # Check if structured_json exists, if not, create a basic one
        if "structured_json" not in state or not state["structured_json"]:
            logger.warning("No structured_json found, creating basic structure")
            # Create a basic structured_json from state data
            state["structured_json"] = {
                "report_metadata": {
                    "pipeline_version": "1.0.0",
                    "processing_time_seconds": (datetime.now() - state.get("pipeline_start_time", datetime.now())).total_seconds(),
                    "generated_at": datetime.now().isoformat(),
                    "language": state.get("language", "en")
                },
                "patient_data": {
                    "file_type": state.get("file_type", "unknown"),
                    "file_metadata": state.get("file_metadata", {})
                },
                "summary": {
                    "total_variants_found": state.get("variant_count", 0),
                    "variants_passed_qc": len(state.get("filtered_variants", [])),
                    "high_risk_findings": len([s for s in state.get("risk_scores", {}).values() if s >= 50])
                },
                "risk_assessment": {
                    "scores": state.get("risk_scores", {}),
                    "high_risk_findings": [],
                    "risk_genes": state.get("risk_genes", {})
                },
                "variant_details": [],
                "quality_control": {},
                "tcga_summary": {
                    "cancer_types_analyzed": list(state.get("tcga_matches", {}).keys()),
                    "cohort_sizes": state.get("tcga_cohort_sizes", {}),
                    "variants_with_tcga_data": 0
                },
                "warnings": state.get("warnings", [])
            }

            # Build high risk findings
            for cancer_type, score in state.get("risk_scores", {}).items():
                if score >= 50:
                    state["structured_json"]["risk_assessment"]["high_risk_findings"].append({
                        "cancer_type": cancer_type,
                        "risk_percentage": score,
                        "affected_genes": state.get("risk_genes", {}).get(cancer_type, [])
                    })

        structured_data = state["structured_json"]
        language = state.get("language", "en")
        include_technical = state.get("include_technical_details", True)

        # Build structured report sections
        report_sections = {
            "header": {
                "title": "Patient Genomic Risk Report",
                "generated_date": datetime.now().strftime('%Y-%m-%d'),
                "generated_time": datetime.now().strftime('%H:%M:%S'),
                "pipeline_version": structured_data['report_metadata']['pipeline_version'],
                "processing_time_seconds": structured_data['report_metadata'].get('processing_time_seconds', 0),
                "language": language
            },

            "summary": {
                "description": "This report details the genomic risk assessment based on provided genetic data. High-risk findings are highlighted.",
                "total_variants_found": structured_data['summary']['total_variants_found'],
                "variants_passed_qc": structured_data['summary']['variants_passed_qc'],
                "high_risk_findings_count": structured_data['summary']['high_risk_findings'],
                "patient_info": state.get("patient_data", {})
            },

            "risk_findings": _format_risk_findings(structured_data) if structured_data.get("risk_assessment") else [],

            "variant_table": _format_variant_table(structured_data),

            "quality_control": {
                "total_variants": structured_data['summary']['total_variants_found'],
                "passed_variants": structured_data['summary']['variants_passed_qc'],
                "filtered_variants": structured_data['summary']['total_variants_found'] - structured_data['summary']['variants_passed_qc'],
                "pass_rate_percent": _calculate_pass_rate(structured_data),
                "filter_criteria": {
                    "min_quality": 30,
                    "min_depth": 10,
                    "min_allele_frequency": 0.01
                }
            },

            "tcga_summary": _format_tcga_summary(structured_data),

            "cadd_summary": _format_cadd_summary(structured_data),

            "recommendations": _generate_recommendations(structured_data),

            "technical_details": _format_technical_details(state) if include_technical else None,

            "disclaimers": [
                "This report is for research purposes only and should not be used for clinical decision-making without consultation with a qualified healthcare provider.",
                "Genetic risk assessment is based on current scientific understanding and may change as new research emerges.",
                "Not all genetic variants are included in this analysis."
            ],

            "metadata": {
                "report_id": f"GR-{datetime.now().strftime('%Y%m%d-%H%M%S')}",
                "model_version": state.get("model_version", "1.0.0"),
                "warnings": [w.get("warning", str(w)) for w in state.get("warnings", []) if isinstance(w, dict)] + [str(w) for w in state.get("warnings", []) if isinstance(w, str)],
                "file_analyzed": state.get("file_path", "Unknown")
            }
        }

        # Also generate markdown version for PDF export
        report_markdown = _generate_markdown_from_sections(report_sections)

        # Debug logging
        logger.info(f"Report sections created with keys: {list(report_sections.keys())}")
        logger.info(f"Report sections header: {report_sections.get('header', {})}")

        logger.info("Report generation complete")

        # Return only the keys this node updates
        return {
            "report_sections": report_sections,
            "report_markdown": report_markdown,
            "pipeline_status": "completed"
        }

    except Exception as e:
        import traceback
        logger.error(f"Report generation failed: {str(e)}")
        logger.error(f"Full traceback: {traceback.format_exc()}")
        # Return error state updates
        return {
            "pipeline_status": "failed",
            "errors": [{
                "node": "report_writer",
                "error": str(e),
                "timestamp": datetime.now()
            }]
        }


def _format_risk_findings(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Format risk findings for easy frontend consumption."""
    findings = []

    risk_assessment = data.get("risk_assessment")
    if not risk_assessment:
        return findings

    for finding in risk_assessment.get("high_risk_findings", []):
        findings.append({
            "cancer_type": finding["cancer_type"],
            "risk_percentage": finding["risk_percentage"],
            "risk_level": _get_risk_level(finding["risk_percentage"]),
            "affected_genes": finding["affected_genes"],
            "gene_count": len(finding["affected_genes"]),
            "recommendation": _get_risk_recommendation(finding["cancer_type"], finding["risk_percentage"])
        })

    # Add low-risk findings
    all_scores = risk_assessment.get("scores", {})
    high_risk_types = [f["cancer_type"] for f in findings]

    for cancer_type, score in all_scores.items():
        if cancer_type not in high_risk_types:
            findings.append({
                "cancer_type": cancer_type,
                "risk_percentage": score,
                "risk_level": _get_risk_level(score),
                "affected_genes": risk_assessment.get("risk_genes", {}).get(cancer_type, []),
                "gene_count": len(risk_assessment.get("risk_genes", {}).get(cancer_type, [])),
                "recommendation": _get_risk_recommendation(cancer_type, score)
            })

    # Sort by risk percentage descending
    findings.sort(key=lambda x: x["risk_percentage"], reverse=True)

    return findings


def _format_variant_table(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Format variant details as structured table data."""
    table_data = []

    for variant in data["variant_details"][:10]:  # Limit to top 10
        # Find best TCGA match
        best_match = None
        best_frequency = 0

        for cancer_type, match_data in variant["tcga_matches"].items():
            if match_data and match_data["frequency_percent"] > best_frequency:
                best_match = {
                    "cancer_type": cancer_type,
                    "frequency": match_data["frequency_percent"],
                    "patient_fraction": match_data["patient_fraction"]
                }
                best_frequency = match_data["frequency_percent"]

        table_data.append({
            "gene": variant["gene"],
            "variant_id": variant["variant"],
            "consequence": variant.get("consequence", "unknown"),
            "hgvs_c": variant.get("hgvs_c", ""),
            "hgvs_p": variant.get("hgvs_p", ""),
            "quality_score": variant["quality_metrics"]["quality"],
            "depth": variant["quality_metrics"]["depth"],
            "allele_frequency": variant["quality_metrics"]["allele_freq"],
            "tcga_best_match": best_match,
            "cadd_score": variant.get("cadd_scores", {}).get("phred") if "cadd_scores" in variant else None,
            "cadd_interpretation": _interpret_cadd_score(variant.get("cadd_scores", {}).get("phred")) if "cadd_scores" in variant else None,
            "clinical_significance": _get_clinical_significance(variant["gene"])
        })

    return table_data


def _format_tcga_summary(data: Dict[str, Any]) -> Dict[str, Any]:
    """Format TCGA comparison summary."""
    tcga_summary = data.get("tcga_summary", {})

    return {
        "cancer_types_analyzed": tcga_summary.get("cancer_types_analyzed", []),
        "total_cohort_size": sum(tcga_summary.get("cohort_sizes", {}).values()),
        "cohort_sizes": tcga_summary.get("cohort_sizes", {}),
        "variants_with_matches": tcga_summary.get("variants_with_tcga_data", 0),
        "match_rate_percent": (tcga_summary.get("variants_with_tcga_data", 0) /
                              data['summary']['variants_passed_qc'] * 100) if data['summary']['variants_passed_qc'] > 0 else 0
    }


def _format_cadd_summary(data: Dict[str, Any]) -> Dict[str, Any]:
    """Format CADD summary."""
    cadd_summary = data.get("cadd_summary", {})

    if not cadd_summary:
        return {
            "enabled": False,
            "variants_scored": 0,
            "mean_phred_score": 0,
            "max_phred_score": 0,
            "high_impact_variants": 0,
            "cancer_gene_variants": 0,
            "description": "CADD scoring not available"
        }

    return {
        "enabled": True,
        "variants_scored": cadd_summary.get("variants_scored", 0),
        "mean_phred_score": round(cadd_summary.get("mean_phred", 0), 2),
        "max_phred_score": round(cadd_summary.get("max_phred", 0), 2),
        "high_impact_variants": cadd_summary.get("high_impact_variants", 0),
        "cancer_gene_variants": cadd_summary.get("cancer_gene_variants", 0),
        "description": "Combined Annotation Dependent Depletion (CADD) scores predict variant deleteriousness"
    }


def _generate_recommendations(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Generate actionable recommendations based on findings."""
    recommendations = []

    risk_assessment = data.get("risk_assessment")
    if not risk_assessment:
        # No risk assessment performed, provide general recommendations
        recommendations.append({
            "priority": "medium",
            "category": "general",
            "title": "Variant Analysis Complete",
            "description": "Genetic variants have been identified and annotated. Consider consulting with a healthcare provider to discuss these findings.",
            "genes_involved": []
        })
        return recommendations

    # Check high-risk findings
    for finding in risk_assessment.get("high_risk_findings", []):
        if finding["risk_percentage"] > 70:
            recommendations.append({
                "priority": "high",
                "category": "screening",
                "title": f"Enhanced {finding['cancer_type'].capitalize()} Cancer Screening",
                "description": f"Due to elevated risk ({finding['risk_percentage']}%), consider discussing enhanced screening protocols with your healthcare provider.",
                "genes_involved": finding["affected_genes"]
            })
        elif finding["risk_percentage"] > 50:
            recommendations.append({
                "priority": "medium",
                "category": "monitoring",
                "title": f"Regular {finding['cancer_type'].capitalize()} Cancer Monitoring",
                "description": f"Moderate risk detected ({finding['risk_percentage']}%). Regular monitoring is recommended.",
                "genes_involved": finding["affected_genes"]
            })

    # General recommendations
    if data['summary']['high_risk_findings'] > 0:
        recommendations.append({
            "priority": "high",
            "category": "consultation",
            "title": "Genetic Counseling Consultation",
            "description": "Consider scheduling a consultation with a genetic counselor to discuss these findings in detail.",
            "genes_involved": []
        })

    return recommendations


def _format_technical_details(state: Dict[str, Any]) -> Dict[str, Any]:
    """Format technical pipeline details."""
    return {
        "pipeline_nodes": state.get("completed_nodes", []),
        "file_metadata": state.get("file_metadata", {}),
        "alignment_stats": state.get("alignment_stats", {}),
        "variant_calling_method": "Mock DeepVariant (simulated)",
        "reference_genome": "GRCh38",
        "quality_thresholds": {
            "min_quality_score": 30,
            "min_read_depth": 10,
            "min_allele_frequency": 0.01
        },
        "model_details": {
            "type": state.get("model_type", "logistic_regression"),
            "version": state.get("model_version", "1.0.0"),
            "training_date": "2025-01-07"
        }
    }


def _calculate_pass_rate(data: Dict[str, Any]) -> float:
    """Calculate QC pass rate percentage."""
    total = data['summary']['total_variants_found']
    passed = data['summary']['variants_passed_qc']
    return (passed / total * 100) if total > 0 else 100.0


def _get_risk_level(percentage: float) -> str:
    """Categorize risk level based on percentage."""
    if percentage >= 70:
        return "high"
    elif percentage >= 30:
        return "moderate"
    else:
        return "low"


def _get_risk_recommendation(cancer_type: str, risk_percentage: float) -> str:
    """Generate brief recommendation based on risk level."""
    if risk_percentage >= 70:
        return f"Enhanced {cancer_type} cancer screening recommended"
    elif risk_percentage >= 30:
        return f"Regular {cancer_type} cancer monitoring advised"
    else:
        return "Standard screening guidelines apply"


def _get_clinical_significance(gene: str) -> str:
    """Determine clinical significance based on gene."""
    high_risk_genes = ["BRCA1", "BRCA2", "TP53", "APC", "MLH1", "MSH2"]
    moderate_risk_genes = ["KRAS", "PIK3CA", "PTEN", "JAK2", "FLT3"]

    if gene in high_risk_genes:
        return "pathogenic"
    elif gene in moderate_risk_genes:
        return "likely_pathogenic"
    else:
        return "uncertain_significance"


def _interpret_cadd_score(phred_score: float) -> str:
    """Interpret CADD PHRED score into human-readable impact."""
    if phred_score is None:
        return "Not available"
    elif phred_score >= 25:
        return "Highly pathogenic"
    elif phred_score >= 20:
        return "Pathogenic"
    elif phred_score >= 15:
        return "Damaging"
    elif phred_score >= 10:
        return "Uncertain"
    else:
        return "Benign"


def _generate_markdown_from_sections(sections: Dict[str, Any]) -> str:
    """Generate markdown report from structured sections (for PDF export)."""
    md = f"""# {sections['header']['title']}

**Generated**: {sections['header']['generated_date']} at {sections['header']['generated_time']}
**Pipeline Version**: {sections['header']['pipeline_version']}

## Summary

{sections['summary']['description']}

- Total variants analyzed: {sections['summary']['total_variants_found']}
- Variants passing quality control: {sections['summary']['variants_passed_qc']}
- High-risk findings: {sections['summary']['high_risk_findings_count']}

## Risk Assessment

"""

    # Add risk findings
    if 'risk_findings' in sections and sections['risk_findings']:
        for finding in sections['risk_findings']:
            level_emoji = "🔴" if finding['risk_level'] == 'high' else "🟡" if finding['risk_level'] == 'moderate' else "🟢"
            md += f"### {level_emoji} {finding['cancer_type'].capitalize()} Cancer\n"
            md += f"- **Risk Score**: {finding['risk_percentage']}%\n"
            md += f"- **Risk Level**: {finding['risk_level'].capitalize()}\n"
            md += f"- **Affected Genes**: {', '.join(finding['affected_genes']) if finding['affected_genes'] else 'None'}\n"
            md += f"- **Recommendation**: {finding['recommendation']}\n\n"
    else:
        md += "*No risk assessment performed*\n\n"

    # Add variant table
    md += "## Key Genetic Variants\n\n"
    md += "| Gene | Variant | Consequence | TCGA Match | Clinical Significance | CADD Score |\n"
    md += "|------|---------|-------------|------------|-----------------------|------------|\n"

    for variant in sections['variant_table'][:5]:
        tcga_text = "No match"
        if variant['tcga_best_match']:
            tcga_text = f"{variant['tcga_best_match']['frequency']:.1f}% in {variant['tcga_best_match']['cancer_type']}"

        cadd_text = "Score not available"
        if variant['cadd_score'] is not None:
            cadd_text = f"{variant['cadd_score']:.1f} ({variant['cadd_interpretation']})"

        md += f"| {variant['gene']} | {variant['variant_id'][:20]}... | {variant['consequence']} | {tcga_text} | {variant['clinical_significance']} | {cadd_text} |\n"

    # Add recommendations
    if sections['recommendations']:
        md += "\n## Recommendations\n\n"
        for rec in sections['recommendations']:
            priority_emoji = "❗" if rec['priority'] == 'high' else "⚠️"
            md += f"{priority_emoji} **{rec['title']}**\n"
            md += f"   {rec['description']}\n\n"

    # Add disclaimers
    md += "\n## Important Notice\n\n"
    for disclaimer in sections['disclaimers']:
        md += f"- {disclaimer}\n"

    md += f"\n---\n*Report ID: {sections['metadata']['report_id']}*\n"
    md += f"*Processing time: {sections['header']['processing_time_seconds']:.2f} seconds*\n"

    return md
