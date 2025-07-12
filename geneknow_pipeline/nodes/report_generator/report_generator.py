"""
Main report generator entry point.
Replaces the existing report_writer.process function with LLM-enhanced report generation.
"""

import logging
import os
import re
from datetime import datetime
from typing import Dict, Any, Optional, Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

from .config import load_config
from .model_interface import ModelInterface
from .prompt_builder import PromptBuilder
from .formatter import ReportFormatter

logger = logging.getLogger(__name__)


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate enhanced genomic risk assessment reports.

    This function replaces the original report_writer.process and generates
    professional reports using LLM enhancement when available,
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
        logger.info(
            f"Using report config: backend={config.backend.value}, style={config.style.value}"
        )

        # Get or create structured JSON data
        structured_json = state.get("structured_json")
        if not structured_json:
            logger.warning("No structured_json found, creating basic structure")
            structured_json = _create_basic_structured_json(state)
            state["structured_json"] = structured_json

        logger.info(
            f"Structured JSON keys: {list(structured_json.keys()) if structured_json else 'None'}"
        )

        # Filter for high-risk findings only (>5% risk)
        filtered_data = _filter_high_risk_findings(
            structured_json, config.risk_threshold
        )

        # Initialize components
        model_interface = ModelInterface(config)
        prompt_builder = PromptBuilder(config.style)

        # Set up output directory (not used for file storage, just for config)
        output_dir = os.path.join(os.getcwd(), "reports")
        formatter = ReportFormatter(config, output_dir)

        # Generate report ID
        report_id = f"report_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

        # Get backend info for logging
        backend_info = model_interface.get_backend_info()
        logger.info(f"Report generation backend: {backend_info}")

        # Create report info early (before generating content)
        report_info = {
            "report_id": report_id,
            "backend_used": backend_info["backend"],
            "model_used": backend_info.get("model"),
            "llm_enhanced": backend_info["available"],
            "generation_time": datetime.now().isoformat(),
            "high_risk_findings_count": _count_high_risk_findings(
                filtered_data, config.risk_threshold
            ),
        }

        # Generate report content
        llm_content = None
        stream_callback = state.get("stream_callback")

        if model_interface.is_available():
            logger.info("Generating LLM-enhanced report content")
            llm_content = _generate_llm_content(
                filtered_data, model_interface, prompt_builder, stream_callback
            )
        else:
            logger.info("LLM not available, using template-based generation")
            if stream_callback:
                stream_callback(
                    {
                        "section": "fallback",
                        "content": f"Using template-based generation ({config.dev_mode_indicator})",
                        "total_length": 0,
                    }
                )

        # Format report content (in-memory, no file storage for HIPAA compliance)
        report_content = formatter.format_report(
            data=filtered_data, llm_content=llm_content, report_id=report_id
        )

        logger.info(
            f"Report generation complete. Formats: {list(report_content.keys())}"
        )

        # Update report info with content metadata
        report_info["available_formats"] = list(report_content.keys())
        report_info["content_length"] = len(report_content.get("markdown", ""))

        # Create report_sections for dashboard compatibility
        dashboard_report_sections = _create_dashboard_report_sections(
            filtered_data, report_info
        )

        # Return updated state with in-memory content
        return {
            "report_sections": dashboard_report_sections,  # For dashboard display
            "report_generator_info": report_info,  # Report generator metadata
            "enhanced_report_content": report_content,  # In-memory report content
            "pipeline_status": "completed",
        }

    except Exception as e:
        logger.error(f"Report generation failed: {str(e)}")

        # Return error state but don't break the pipeline
        error_report_sections = {
            "error": {
                "title": "Report Generation Error",
                "content": f"Failed to generate enhanced report: {str(e)}",
                "severity": "high",
            }
        }

        return {
            "report_sections": error_report_sections,  # Error message for dashboard
            "report_generator_info": {
                "error": str(e),
                "generation_time": datetime.now().isoformat(),
                "llm_enhanced": False,
            },
            "enhanced_report_content": {},
            "pipeline_status": "completed",  # Don't fail the entire pipeline
        }


def _create_basic_structured_json(state: Dict[str, Any]) -> Dict[str, Any]:
    """Create basic structured JSON from pipeline state."""

    return {
        "report_metadata": {
            "pipeline_version": "1.0.0",
            "processing_time_seconds": (
                datetime.now() - state.get("pipeline_start_time", datetime.now())
            ).total_seconds(),
            "generated_at": datetime.now().isoformat(),
            "language": state.get("language", "en"),
        },
        "patient_data": {
            "file_type": state.get("file_type", "unknown"),
            "file_metadata": state.get("file_metadata", {}),
        },
        "summary": {
            "total_variants_found": state.get("variant_count", 0),
            "variants_passed_qc": len(state.get("filtered_variants", [])),
            "high_risk_findings": len(
                [s for s in state.get("risk_scores", {}).values() if s >= 5.0]
            ),
        },
        "risk_assessment": {
            "scores": state.get("risk_scores", {}),
            "high_risk_findings": [],
            "risk_genes": state.get("risk_genes", {}),
        },
        "variant_details": _format_variant_details(state.get("filtered_variants", [])),
        "quality_control": {},
        "tcga_summary": {
            "cancer_types_analyzed": list(state.get("tcga_matches", {}).keys()),
            "cohort_sizes": state.get("tcga_cohort_sizes", {}),
            "variants_with_tcga_data": 0,
        },
        "cadd_summary": state.get("cadd_stats", {}),
        "warnings": state.get("warnings", []),
    }


def _format_variant_details(variants: list) -> list:
    """Format variant details for report generation."""

    formatted_variants = []

    for variant in variants[:20]:  # Limit to top 20 variants
        formatted_variant = {
            "gene": variant.get("gene", "Unknown"),
            "variant": variant.get("variant_id", "Unknown"),
            "consequence": variant.get(
                "consequence", variant.get("variant_classification", "unknown")
            ),
            "hgvs_c": variant.get("hgvs_c", ""),
            "hgvs_p": variant.get("hgvs_p", variant.get("protein_change", "")),
            "quality_metrics": {
                "quality": variant.get("quality", variant.get("qual", 0)),
                "depth": variant.get("depth", 0),
                "allele_freq": variant.get("allele_freq", 0),
            },
            "tcga_matches": {},
            "cadd_scores": {},
        }

        # Add clinical significance if available
        if "clinical_significance" in variant:
            formatted_variant["clinical_significance"] = variant[
                "clinical_significance"
            ]

        # Add CADD scores if available
        if "cadd_phred" in variant:
            formatted_variant["cadd_scores"] = {
                "phred": variant.get("cadd_phred"),
                "raw": variant.get("cadd_raw"),
            }

        formatted_variants.append(formatted_variant)

    return formatted_variants


def _filter_high_risk_findings(
    data: Dict[str, Any], risk_threshold: float
) -> Dict[str, Any]:
    """Filter data to focus on high-risk findings above threshold."""

    filtered_data = data.copy()

    # Filter risk assessment for high-risk findings only
    risk_assessment = filtered_data.get("risk_assessment") or {}
    if risk_assessment and "scores" in risk_assessment:
        scores = risk_assessment["scores"]
        risk_genes = risk_assessment.get("risk_genes", {})

        # Keep only high-risk cancers
        high_risk_scores = {
            cancer: score for cancer, score in scores.items() if score > risk_threshold
        }

        high_risk_genes = {
            cancer: genes
            for cancer, genes in risk_genes.items()
            if cancer in high_risk_scores
        }

        # Update high_risk_findings
        high_risk_findings = []
        for cancer_type, score in high_risk_scores.items():
            high_risk_findings.append(
                {
                    "cancer_type": cancer_type,
                    "risk_percentage": score,
                    "affected_genes": high_risk_genes.get(cancer_type, []),
                }
            )

        filtered_data["risk_assessment"] = {
            "scores": high_risk_scores,
            "high_risk_findings": high_risk_findings,
            "risk_genes": high_risk_genes,
        }

    return filtered_data


def _count_high_risk_findings(data: Dict[str, Any], risk_threshold: float) -> int:
    """Count the number of high-risk findings."""

    risk_assessment = data.get("risk_assessment") or {}
    if not risk_assessment or "scores" not in risk_assessment:
        return 0

    return sum(
        1 for score in risk_assessment["scores"].values() if score > risk_threshold
    )


def _create_dashboard_report_sections(
    data: Dict[str, Any], report_info: Dict[str, Any]
) -> Dict[str, Any]:
    """Create report_sections structure for dashboard display."""

    logger.info(
        f"Creating dashboard sections with data keys: {list(data.keys()) if data else 'None'}"
    )
    logger.info(f"Report info available: {report_info is not None}")

    sections = {}

    # Overview section
    risk_assessment = data.get("risk_assessment") or {}
    high_risk_findings = risk_assessment.get("high_risk_findings", [])

    if high_risk_findings:
        # High-risk overview
        top_risk = max(high_risk_findings, key=lambda x: x["risk_percentage"])
        sections["overview"] = {
            "title": "Analysis Overview",
            "content": f"Genomic analysis identified {len(high_risk_findings)} cancer types with elevated risk. Highest risk: {top_risk['cancer_type']} at {top_risk['risk_percentage']:.1f}%.",
            "severity": "high" if top_risk["risk_percentage"] >= 50 else "medium",
        }
    else:
        # Low-risk overview
        sections["overview"] = {
            "title": "Analysis Overview",
            "content": "Genomic analysis completed successfully. All cancer risk levels are within normal baseline ranges.",
            "severity": "low",
        }

    # Key findings section
    if high_risk_findings:
        findings_content = []
        for finding in high_risk_findings[:3]:  # Top 3 findings
            findings_content.append(
                f"{finding['cancer_type'].title()}: {finding['risk_percentage']:.1f}% risk"
            )

        sections["key_findings"] = {
            "title": "Key Findings",
            "content": ". ".join(findings_content) + ".",
            "severity": (
                "high" if high_risk_findings[0]["risk_percentage"] >= 50 else "medium"
            ),
            "technical_details": f"Analysis included {data.get('summary', {}).get('total_variants_found', 0)} variants with {data.get('summary', {}).get('variants_passed_qc', 0)} passing quality control.",
        }

    # Variant analysis section
    variant_details = data.get("variant_details", [])
    if variant_details:
        key_genes = [v["gene"] for v in variant_details[:3]]
        sections["variant_analysis"] = {
            "title": "Variant Analysis",
            "content": f"Key genetic variants identified in {', '.join(key_genes)}. Detailed pathogenicity assessment completed.",
            "severity": "medium",
                            "technical_details": f"CADD scores and pathogenicity evaluated for {len(variant_details)} variants.",
        }

    # Report generation info
    if report_info:
        backend_used = report_info.get("backend_used", "template")
        if backend_used != "none":
            sections["ai_enhancement"] = {
                "title": "AI-Enhanced Analysis",
                "content": f"Report generated using {backend_used} LLM with model {report_info.get('model_used', 'unknown')}. Enhanced interpretations provided.",
                "severity": "low",
            }

    return sections


def _clean_llm_output(content: str, section_name: str) -> str:
    """Clean LLM output to remove instruction artifacts and prompt text."""
    
    # Remove common instruction phrases
    instruction_patterns = [
        r"^Here is the .* section content for.*?:\s*",
        r"^Here is the .* section content:\s*",
        r"^Generated content for.*?:\s*",
        r"^Content for the .* section:\s*",
        r"^The following is the .* section.*?:\s*",
        r"^Below is the .* section.*?:\s*",
        r"^\*\*.*Section.*\*\*\s*",
        r"^#+ .*Section.*\s*",
        r"^## .*\s*",  # Remove section headers that start with ##
        r"^# .*\s*",   # Remove section headers that start with #
        r"^OUTPUT:\s*",
        r"^CONTENT:\s*",
        r"^RESULT:\s*",
        r"^RESPONSE:\s*",
    ]
    
    # Apply cleaning patterns
    cleaned_content = content
    for pattern in instruction_patterns:
        cleaned_content = re.sub(pattern, "", cleaned_content, flags=re.IGNORECASE | re.MULTILINE)
    
    # Remove section name headers that might be at the beginning
    section_header_patterns = [
        rf"^{re.escape(section_name)}\s*:?\s*",
        rf"^\*\*{re.escape(section_name)}\*\*\s*:?\s*",
        rf"^#{1,6}\s*{re.escape(section_name)}\s*:?\s*",
    ]
    
    for pattern in section_header_patterns:
        cleaned_content = re.sub(pattern, "", cleaned_content, flags=re.IGNORECASE | re.MULTILINE)
    
    # Remove any leading/trailing whitespace or newlines
    cleaned_content = cleaned_content.strip()
    
    # If content is significantly reduced, log a warning
    if len(cleaned_content) < len(content) * 0.5:
        logger.warning(f"Significant content reduction after cleaning {section_name}: {len(content)} -> {len(cleaned_content)} characters")
    
    return cleaned_content


def _generate_single_section(section_name: str, 
                           data: Dict[str, Any],
                           model_interface: ModelInterface,
                           prompt_builder: PromptBuilder,
                           stream_callback: Optional[Callable] = None) -> tuple[str, str]:
    """Generate content for a single section. Returns (section_name, content)."""
    
    try:
        if stream_callback:
            stream_callback({
                "section": section_name.lower().replace(" ", "_"),
                "content": f"Generating {section_name} section...",
                "total_length": 0
            })
        
        # Get the appropriate prompt for this section
        if section_name == "Summary":
            prompt = prompt_builder.build_summary_section_prompt(data)
        elif section_name == "Key Variants":
            prompt = prompt_builder.build_key_variants_section_prompt(data)
        elif section_name == "Risk Summary":
            prompt = prompt_builder.build_risk_summary_section_prompt(data)
        elif section_name == "Interpretation":
            prompt = prompt_builder.build_clinical_interpretation_section_prompt(data)
        elif section_name == "Recommendations":
            prompt = prompt_builder.build_recommendations_section_prompt(data)
        else:
            logger.warning(f"Unknown section: {section_name}")
            return section_name, f"Section '{section_name}' not implemented"
        
        # Generate content for this section
        section_content = model_interface.generate(prompt, None)  # No streaming per section
        
        if section_content and section_content.strip():
            content = section_content.strip()
            # Clean up LLM output artifacts
            content = _clean_llm_output(content, section_name)
            logger.info(f"Generated {section_name} section: {len(content)} characters")
            return section_name, content
        else:
            logger.warning(f"Empty content generated for {section_name} section, using fallback")
            return section_name, _get_fallback_section_content(section_name, data)
            
    except Exception as e:
        logger.error(f"Failed to generate {section_name} section: {e}")
        return section_name, _get_fallback_section_content(section_name, data)


def _generate_llm_content(data: Dict[str, Any], 
                         model_interface: ModelInterface,
                         prompt_builder: PromptBuilder,
                         stream_callback: Optional[Callable] = None) -> str:
    """Generate LLM-enhanced report content with configurable parallel section generation."""
    
    try:
        # Check if we have any high-risk findings to report on
        high_risk_count = _count_high_risk_findings(data, prompt_builder.risk_threshold)
        
        # Get config from model_interface
        config = model_interface.config
        
        # Define sections to generate
        section_names = ["Summary", "Key Variants", "Risk Summary", "Interpretation", "Recommendations"]
        sections = {}
        
        # Record start time for performance monitoring
        start_time = time.time()
        
        if config.enable_parallel_generation and len(section_names) > 1:
            logger.info(f"Generating LLM content for {high_risk_count} high-risk findings using parallel section generation (max workers: {config.max_parallel_workers})")
            
            # Generate sections in parallel using ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=min(len(section_names), config.max_parallel_workers)) as executor:
                # Submit all section generation tasks
                future_to_section = {
                    executor.submit(_generate_single_section, section_name, data, model_interface, prompt_builder, stream_callback): section_name 
                    for section_name in section_names
                }
                
                # Collect results as they complete
                for future in as_completed(future_to_section):
                    section_name = future_to_section[future]
                    try:
                        returned_section_name, content = future.result()
                        sections[returned_section_name] = content
                        
                        if stream_callback:
                            stream_callback({
                                "section": returned_section_name.lower().replace(" ", "_"),
                                "content": f"Completed {returned_section_name} section",
                                "total_length": len(content)
                            })
                            
                    except Exception as e:
                        logger.error(f"Section generation failed for {section_name}: {e}")
                        sections[section_name] = _get_fallback_section_content(section_name, data)
        else:
            # Sequential generation (fallback or disabled parallel processing)
            logger.info(f"Generating LLM content for {high_risk_count} high-risk findings using sequential section generation")
            
            for section_name in section_names:
                try:
                    returned_section_name, content = _generate_single_section(
                        section_name, data, model_interface, prompt_builder, stream_callback
                    )
                    sections[returned_section_name] = content
                    
                    if stream_callback:
                        stream_callback({
                            "section": returned_section_name.lower().replace(" ", "_"),
                            "content": f"Completed {returned_section_name} section",
                            "total_length": len(content)
                        })
                        
                except Exception as e:
                    logger.error(f"Section generation failed for {section_name}: {e}")
                    sections[section_name] = _get_fallback_section_content(section_name, data)
        
        # Log performance improvement
        generation_time = time.time() - start_time
        generation_mode = "parallel" if config.enable_parallel_generation else "sequential"
        logger.info(f"{generation_mode.capitalize()} section generation completed in {generation_time:.2f} seconds")
        # Combine all sections into the final report
        final_content = _combine_sections_into_report(sections)

        if stream_callback:
            stream_callback(
                {
                    "section": "complete",
                    "content": "Report generation complete",
                    "total_length": len(final_content),
                }
            )

        logger.info(f"Generated complete report with {len(final_content)} characters")
        return final_content

    except Exception as e:
        logger.error(f"LLM content generation failed: {e}")
        return ""


def _get_fallback_section_content(section_name: str, data: Dict[str, Any]) -> str:
    """Generate fallback content for a section if LLM generation fails."""

    if section_name == "Summary":
        summary = data.get("summary", {})
        return f"This genomic risk assessment analyzed {summary.get('total_variants_found', 0)} genetic variants, with {summary.get('variants_passed_qc', 0)} variants passing quality control filters. The analysis completed successfully with results available for review."

    elif section_name == "Key Variants":
        variant_details = data.get("variant_details", [])
        if not variant_details:
            return "No key pathogenic variants associated with elevated cancer risk were identified in this analysis."
        else:
            return f"Analysis identified {len(variant_details)} variants for review. Detailed variant interpretation requires genetics expertise."

    elif section_name == "Risk Summary":
        risk_assessment = data.get("risk_assessment", {})
        scores = risk_assessment.get("scores", {}) if risk_assessment else {}

        if not scores:
            return "| Cancer Type | Risk (%) |\n|-------------|----------|\n| No data | N/A |\n\nRisk assessment data not available for this analysis."

        # Create simple table
        table_rows = []
        for cancer_type, score in sorted(
            scores.items(), key=lambda x: x[1], reverse=True
        ):
            table_rows.append(f"| {cancer_type.title()} | {score:.1f} |")

        table = "| Cancer Type | Risk (%) |\n|-------------|----------|\n" + "\n".join(
            table_rows
        )

        high_risk_count = sum(1 for score in scores.values() if score > 5.0)
        if high_risk_count > 0:
            table += "\n\nRisks above 5% are considered elevated and warrant further attention."
        else:
            table += "\n\nAll risks are within baseline population levels (<5%), indicating no elevated genetic predisposition was detected."

        return table

    elif section_name == "Interpretation":
        return "This genomic risk assessment utilized multiple computational approaches including population frequency analysis, TCGA tumor database comparison, pathogenicity prediction algorithms, polygenic risk scoring, and machine learning risk models. The analysis provides research-grade insights that should be interpreted within the context of current scientific understanding and guidelines. These findings should be correlated with family history, lifestyle factors, and presentation for comprehensive risk assessment."

    elif section_name == "Recommendations":
        risk_assessment = data.get("risk_assessment", {})
        scores = risk_assessment.get("scores", {}) if risk_assessment else {}
        high_risk_count = (
            sum(1 for score in scores.values() if score > 5.0) if scores else 0
        )

        if high_risk_count > 0:
            return "- Genetic counseling is recommended to discuss these findings and develop a personalized risk management plan\n- Enhanced screening protocols may be appropriate for elevated cancer risks\n- Consider consultation with oncology specialists for high-risk findings\n- Maintain regular follow-up care and health monitoring\n- Consult a qualified healthcare provider for personalized guidance."
        else:
            return "- Adhere to standard age-appropriate cancer screening protocols\n- Maintain healthy lifestyle including regular exercise and balanced diet\n- Continue routine preventive care and health maintenance\n- Consider genetic counseling if strong family history of cancer emerges\n- Consult a qualified healthcare provider for personalized guidance."

    return "Content not available."


def _combine_sections_into_report(sections: Dict[str, str]) -> str:
    """Combine individual sections into a complete report with consistent formatting."""

    report_parts = []

    # Summary section
    if "Summary" in sections:
        report_parts.append(f"**Summary**\n\n{sections['Summary']}")

    # Key Variants section
    if "Key Variants" in sections:
        report_parts.append(f"**Key Variants**\n\n{sections['Key Variants']}")

    # Risk Summary section
    if "Risk Summary" in sections:
        report_parts.append(f"**Risk Summary**\n\n{sections['Risk Summary']}")

    # Interpretation section
    if "Interpretation" in sections:
        report_parts.append(
            f"**Interpretation**\n\n{sections['Interpretation']}"
        )

    # Recommendations section
    if "Recommendations" in sections:
        report_parts.append(f"**Recommendations**\n\n{sections['Recommendations']}")

    # Join all sections with double line breaks
    return "\n\n".join(report_parts)


def generate_report_standalone(
    json_file_path: str,
    config_path: Optional[str] = None,
    output_dir: Optional[str] = None,
) -> Dict[str, str]:
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
    with open(json_file_path, "r") as f:
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
        llm_content = _generate_llm_content(
            filtered_data, model_interface, prompt_builder
        )

    # Format and save
    report_content = formatter.format_report(
        data=filtered_data, llm_content=llm_content, report_id=report_id
    )

    logger.info(f"Standalone report generated: {report_content}")
    return report_content


# Example usage for testing
if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        json_file = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else None

        try:
            content = generate_report_standalone(json_file, output_dir=output_dir)
            print("Report generated successfully:")
            print(f"  Content length: {len(content)}")
            print(f"  Formats available: {list(content.keys())}")
        except Exception as e:
            print(f"Error generating report: {e}")
    else:
        print("Usage: python report_generator.py <json_file> [output_dir]")
