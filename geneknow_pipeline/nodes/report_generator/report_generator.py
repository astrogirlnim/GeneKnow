"""
Main report generator entry point.
Replaces the existing report_writer.process function with LLM-enhanced report generation.
"""

import logging
import os
from datetime import datetime
from typing import Dict, Any, Optional, Callable

from .config import load_config, ReportConfig, LLMBackend
from .model_interface import ModelInterface
from .prompt_builder import PromptBuilder
from .formatter import ReportFormatter

logger = logging.getLogger(__name__)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate enhanced genomic risk assessment reports.
    
    This function replaces the original report_writer.process and generates
    professional clinical reports using LLM enhancement when available,
    with intelligent fallback to template-based generation.
    
    Args:
        state: Pipeline state containing structured_json and other data
        
    Returns:
        Updated state with report generation results
    """
    logger.info("Starting enhanced report generation")
    
    try:
        # Load configuration
        config = load_config()
        logger.info(f"Using report config: backend={config.backend.value}, style={config.style.value}")
        
        # Get or create structured JSON data
        structured_json = state.get("structured_json")
        if not structured_json:
            logger.warning("No structured_json found, creating basic structure")
            structured_json = _create_basic_structured_json(state)
            state["structured_json"] = structured_json
        
        # Filter for high-risk findings only (>5% risk)
        filtered_data = _filter_high_risk_findings(structured_json, config.risk_threshold)
        
        # Initialize components
        model_interface = ModelInterface(config)
        prompt_builder = PromptBuilder(config.style)
        
        # Set up output directory
        output_dir = os.path.join(os.getcwd(), "reports")
        formatter = ReportFormatter(config, output_dir)
        
        # Generate report ID
        report_id = f"report_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        # Get backend info for logging
        backend_info = model_interface.get_backend_info()
        logger.info(f"Report generation backend: {backend_info}")
        
        # Generate report content
        llm_content = None
        stream_callback = state.get("stream_callback")
        
        if model_interface.is_available():
            logger.info("Generating LLM-enhanced report content")
            llm_content = _generate_llm_content(
                filtered_data, 
                model_interface, 
                prompt_builder, 
                stream_callback
            )
        else:
            logger.info("LLM not available, using template-based generation")
            if stream_callback:
                stream_callback({
                    "section": "fallback",
                    "content": f"Using template-based generation ({config.dev_mode_indicator})",
                    "total_length": 0
                })
        
        # Format and save report
        report_paths = formatter.format_report(
            data=filtered_data,
            llm_content=llm_content,
            report_id=report_id
        )
        
        logger.info(f"Report generation complete. Files: {list(report_paths.keys())}")
        
        # Update state with report information
        report_info = {
            "report_id": report_id,
            "report_paths": report_paths,
            "backend_used": backend_info["backend"],
            "model_used": backend_info.get("model"),
            "llm_enhanced": backend_info["available"],
            "generation_time": datetime.now().isoformat(),
            "high_risk_findings_count": _count_high_risk_findings(filtered_data, config.risk_threshold)
        }
        
        # Preserve original report_sections for dashboard compatibility
        original_report_sections = state.get("report_sections", {})
        
        # Return updated state
        return {
            "report_sections": original_report_sections,  # Keep original for dashboard
            "report_generator_info": report_info,  # New report generator info
            "enhanced_report_paths": report_paths,  # Paths to generated files
            "pipeline_status": "completed"
        }
        
    except Exception as e:
        logger.error(f"Report generation failed: {str(e)}")
        
        # Return error state but don't break the pipeline
        return {
            "report_sections": state.get("report_sections", {}),  # Preserve original
            "report_generator_info": {
                "error": str(e),
                "generation_time": datetime.now().isoformat(),
                "llm_enhanced": False
            },
            "enhanced_report_paths": {},
            "pipeline_status": "completed"  # Don't fail the entire pipeline
        }


def _create_basic_structured_json(state: Dict[str, Any]) -> Dict[str, Any]:
    """Create basic structured JSON from pipeline state."""
    
    return {
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
            "high_risk_findings": len([s for s in state.get("risk_scores", {}).values() if s >= 5.0])
        },
        "risk_assessment": {
            "scores": state.get("risk_scores", {}),
            "high_risk_findings": [],
            "risk_genes": state.get("risk_genes", {})
        },
        "variant_details": _format_variant_details(state.get("filtered_variants", [])),
        "quality_control": {},
        "tcga_summary": {
            "cancer_types_analyzed": list(state.get("tcga_matches", {}).keys()),
            "cohort_sizes": state.get("tcga_cohort_sizes", {}),
            "variants_with_tcga_data": 0
        },
        "cadd_summary": state.get("cadd_stats", {}),
        "warnings": state.get("warnings", [])
    }


def _format_variant_details(variants: list) -> list:
    """Format variant details for report generation."""
    
    formatted_variants = []
    
    for variant in variants[:20]:  # Limit to top 20 variants
        formatted_variant = {
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
            "tcga_matches": {},
            "cadd_scores": {}
        }
        
        # Add clinical significance if available
        if "clinical_significance" in variant:
            formatted_variant["clinical_significance"] = variant["clinical_significance"]
        
        # Add CADD scores if available
        if "cadd_phred" in variant:
            formatted_variant["cadd_scores"] = {
                "phred": variant.get("cadd_phred"),
                "raw": variant.get("cadd_raw")
            }
        
        formatted_variants.append(formatted_variant)
    
    return formatted_variants


def _filter_high_risk_findings(data: Dict[str, Any], risk_threshold: float) -> Dict[str, Any]:
    """Filter data to focus on high-risk findings above threshold."""
    
    filtered_data = data.copy()
    
    # Filter risk assessment for high-risk findings only
    risk_assessment = filtered_data.get("risk_assessment", {})
    if risk_assessment and "scores" in risk_assessment:
        scores = risk_assessment["scores"]
        risk_genes = risk_assessment.get("risk_genes", {})
        
        # Keep only high-risk cancers
        high_risk_scores = {
            cancer: score for cancer, score in scores.items() 
            if score > risk_threshold
        }
        
        high_risk_genes = {
            cancer: genes for cancer, genes in risk_genes.items()
            if cancer in high_risk_scores
        }
        
        # Update high_risk_findings
        high_risk_findings = []
        for cancer_type, score in high_risk_scores.items():
            high_risk_findings.append({
                "cancer_type": cancer_type,
                "risk_percentage": score,
                "affected_genes": high_risk_genes.get(cancer_type, [])
            })
        
        filtered_data["risk_assessment"] = {
            "scores": high_risk_scores,
            "high_risk_findings": high_risk_findings,
            "risk_genes": high_risk_genes
        }
    
    return filtered_data


def _count_high_risk_findings(data: Dict[str, Any], risk_threshold: float) -> int:
    """Count the number of high-risk findings."""
    
    risk_assessment = data.get("risk_assessment", {})
    if not risk_assessment or "scores" not in risk_assessment:
        return 0
    
    return sum(1 for score in risk_assessment["scores"].values() if score > risk_threshold)


def _generate_llm_content(data: Dict[str, Any], 
                         model_interface: ModelInterface,
                         prompt_builder: PromptBuilder,
                         stream_callback: Optional[Callable] = None) -> str:
    """Generate LLM-enhanced report content."""
    
    try:
        # Check if we have any high-risk findings to report on
        high_risk_count = _count_high_risk_findings(data, prompt_builder.risk_threshold)
        
        if high_risk_count == 0:
            # Generate low-risk report
            prompt = prompt_builder._build_low_risk_prompt()
        else:
            # Generate comprehensive report
            prompt = prompt_builder.build_full_report_prompt(data)
        
        logger.info(f"Generating LLM content for {high_risk_count} high-risk findings")
        
        # Generate content with optional streaming
        content = model_interface.generate(prompt, stream_callback)
        
        if not content:
            logger.warning("LLM generated empty content")
            return ""
        
        logger.info(f"Generated {len(content)} characters of LLM content")
        return content
        
    except Exception as e:
        logger.error(f"LLM content generation failed: {e}")
        return ""


def generate_report_standalone(json_file_path: str, 
                              config_path: Optional[str] = None,
                              output_dir: Optional[str] = None) -> Dict[str, str]:
    """
    Standalone function to generate reports from JSON files.
    Useful for testing and batch processing.
    
    Args:
        json_file_path: Path to structured JSON file
        config_path: Optional path to config file
        output_dir: Optional output directory
        
    Returns:
        Dictionary with paths to generated files
    """
    import json
    
    # Load JSON data
    with open(json_file_path, 'r') as f:
        data = json.load(f)
    
    # Load configuration
    config = load_config(config_path)
    
    # Set up components
    model_interface = ModelInterface(config)
    prompt_builder = PromptBuilder(config.style)
    
    if output_dir is None:
        output_dir = os.path.dirname(json_file_path)
    
    formatter = ReportFormatter(config, output_dir)
    
    # Generate report
    report_id = f"standalone_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    
    # Filter for high-risk findings
    filtered_data = _filter_high_risk_findings(data, config.risk_threshold)
    
    # Generate LLM content if available
    llm_content = None
    if model_interface.is_available():
        llm_content = _generate_llm_content(filtered_data, model_interface, prompt_builder)
    
    # Format and save
    report_paths = formatter.format_report(
        data=filtered_data,
        llm_content=llm_content,
        report_id=report_id
    )
    
    logger.info(f"Standalone report generated: {report_paths}")
    return report_paths


# Example usage for testing
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        json_file = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else None
        
        try:
            paths = generate_report_standalone(json_file, output_dir=output_dir)
            print(f"Report generated successfully:")
            for format_type, path in paths.items():
                print(f"  {format_type}: {path}")
        except Exception as e:
            print(f"Error generating report: {e}")
    else:
        print("Usage: python report_generator.py <json_file> [output_dir]") 