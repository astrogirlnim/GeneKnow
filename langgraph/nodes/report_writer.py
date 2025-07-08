"""
Report writer node.
Generates human-readable report using LLM (Llama 3.1).
"""
import logging
from datetime import datetime
from typing import Dict, Any

logger = logging.getLogger(__name__)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate final report using LLM.
    
    Updates state with:
    - report_markdown: generated markdown report
    - report_pdf_path: path to PDF (if generated)
    - pipeline_status: 'completed'
    """
    logger.info("Starting report generation")
    state["current_node"] = "report_writer"
    
    try:
        structured_data = state["structured_json"]
        language = state.get("language", "en")
        include_technical = state.get("include_technical_details", True)
        
        # ⚠️ MOCK IMPLEMENTATION - Replace with real LLM generation
        # TODO: Real implementation should:
        # 1. Load Llama 3.1 model (GGUF format)
        # 2. Build prompt with structured data
        # 3. Generate report in requested language
        # 4. Convert markdown to PDF
        
        # ⚠️ MOCK: Generate sample report
        high_risk = structured_data["risk_assessment"]["high_risk_findings"]
        
        # Build mock report based on template from documentation
        report_markdown = f"""# Patient Genomic Risk Report

## Summary
> This report details the genomic risk assessment based on provided genetic data. High-risk findings are highlighted.

**Report Date**: {datetime.now().strftime('%Y-%m-%d')}  
**Pipeline Version**: {structured_data['report_metadata']['pipeline_version']}

## Key Risk Findings

"""
        
        # Add high risk findings
        for finding in high_risk:
            cancer = finding["cancer_type"].capitalize()
            risk = finding["risk_percentage"]
            genes = ", ".join(finding["affected_genes"])
            report_markdown += f"- **{cancer} Cancer Risk: {risk}%**\n"
            report_markdown += f"  - Affected genes: {genes}\n"
        
        report_markdown += f"""
## Variant Match Summary

| Gene | Variant | TCGA Match | Clinical Note |
|------|---------|------------|----------------|
"""
        
        # Add variant details
        for variant in structured_data["variant_details"][:5]:  # Limit to top 5
            gene = variant["gene"]
            var_id = variant["variant"]
            
            # Get TCGA match for most relevant cancer
            tcga_match = "No match"
            for cancer_type, match_data in variant["tcga_matches"].items():
                if match_data:
                    tcga_match = f"{match_data['frequency_percent']:.1f}% of {cancer_type.upper()} patients"
                    break
            
            clinical_note = "Consider enhanced screening" if gene in ["BRCA1", "BRCA2", "TP53"] else "Monitor"
            
            report_markdown += f"| {gene} | {var_id} | {tcga_match} | {clinical_note} |\n"
        
        # Calculate QC pass rate
        total_variants = structured_data['summary']['total_variants_found']
        passed_variants = structured_data['summary']['variants_passed_qc']
        qc_pass_rate = (passed_variants / total_variants * 100) if total_variants > 0 else 100.0
        
        # Get processing time
        processing_time = structured_data['report_metadata'].get('processing_time_seconds', 0)
        
        report_markdown += f"""

## Quality Control Summary

- Total variants found: {total_variants}
- Variants passing QC: {passed_variants}
- QC pass rate: {qc_pass_rate:.1f}%

## Notes

⚠️ **MOCK REPORT** - This is a demonstration report using simulated data.

Generated with Llama 3.1 using local inference.
Processing time: {processing_time:.2f} seconds
"""
        
        # Update state
        state["report_markdown"] = report_markdown
        state["report_pdf_path"] = "reports/genomic_report_mock.pdf"  # ⚠️ MOCK PATH
        state["pipeline_status"] = "completed"
        state["completed_nodes"].append("report_writer")
        
        # ⚠️ MOCK WARNING
        state["warnings"].append({
            "node": "report_writer",
            "warning": "Using MOCK report generation - replace with real LLM implementation",
            "timestamp": datetime.now()
        })
        
        logger.info("Report generation complete")
        
    except Exception as e:
        logger.error(f"Report generation failed: {str(e)}")
        state["errors"].append({
            "node": "report_writer",
            "error": str(e),
            "timestamp": datetime.now()
        })
        state["pipeline_status"] = "failed"
    
    return state 