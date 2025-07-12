"""
State definitions for the genomic pipeline.
This module defines the TypedDict schema for the pipeline state.
"""

from typing import TypedDict, Optional, List, Dict, Any
from typing_extensions import Annotated
from datetime import datetime


# Define a reducer for completed_nodes that merges lists
def merge_completed_nodes(existing: List[str], new: List[str]) -> List[str]:
    """Merge completed nodes lists, removing duplicates while preserving order."""
    # Convert to set to remove duplicates, then back to list
    combined = existing + new
    # Preserve order while removing duplicates
    seen = set()
    result = []
    for node in combined:
        if node not in seen:
            seen.add(node)
            result.append(node)
    return result


class GenomicState(TypedDict, total=False):
    """
    State schema for the genomic analysis pipeline.
    All fields are optional (total=False) to support partial updates.
    """

    # File information
    file_path: str
    file_type: Optional[str]
    file_metadata: Dict[str, Any]
    aligned_bam_path: Optional[str]

    # Variant data
    raw_variants: List[Dict[str, Any]]
    filtered_variants: List[Dict[str, Any]]
    variant_count: int

    # Analysis results from parallel nodes
    tcga_matches: Optional[Dict[str, Dict[str, Any]]]
    tcga_cohort_sizes: Dict[str, int]
    tcga_summary: Optional[Dict[str, Any]]

    cadd_enriched_variants: Optional[List[Dict[str, Any]]]
    cadd_stats: Optional[Dict[str, Any]]

    clinvar_annotations: Optional[Dict[str, Dict[str, Any]]]
    clinvar_stats: Optional[Dict[str, Any]]

    prs_results: Optional[Dict[str, Dict[str, Any]]]
    prs_summary: Optional[Dict[str, Any]]

    pathway_burden_results: Optional[Dict[str, Any]]
    pathway_burden_summary: Optional[Dict[str, Any]]
    pathway_enriched_variants: Optional[List[Dict[str, Any]]]

    # ML fusion fields
    ml_ready_variants: Optional[List[Dict[str, Any]]]
    ml_fusion_results: Optional[Dict[str, Any]]
    ml_fusion_model_instance: Optional[Any]
    ml_fusion_feature_matrix: Optional[Any]

    # SHAP validation fields
    shap_validation_status: Optional[str]
    shap_validation_reasons: Optional[List[str]]
    shap_top_contributors: Optional[List[Dict[str, Any]]]
    shap_feature_importance: Optional[Dict[str, Any]]
    shap_validation_details: Optional[Dict[str, Any]]

    # New fields for clinical view nodes
    classified_variants: Optional[List[Dict[str, Any]]]
    mutation_type_distribution: Optional[Dict[str, int]]

    mutational_signatures: Optional[Dict[str, Any]]
    mutational_signatures_summary: Optional[Dict[str, Any]]

    variant_details: Optional[List[Dict[str, Any]]]
    variant_transformation_summary: Optional[Dict[str, Any]]

    structural_variants: Optional[List[Dict[str, Any]]]
    structural_variant_summary: Optional[Dict[str, Any]]

    copy_number_variants: Optional[List[Dict[str, Any]]]
    cnv_summary: Optional[Dict[str, Any]]

    genomic_alterations: Optional[Dict[str, Any]]

    pathway_analysis: Optional[Dict[str, Any]]

    gene_interactions: Optional[Dict[str, Any]]
    gene_network_analysis: Optional[Dict[str, Any]]

    survival_analysis: Optional[Dict[str, Any]]

    clinical_recommendations: Optional[Dict[str, Any]]

    # Risk assessment
    risk_scores: Dict[str, float]  # {cancer_type: risk_percentage}
    risk_genes: Dict[str, List[str]]  # {cancer_type: [affected_genes]}
    risk_details: Dict[str, Any]
    ml_risk_assessment: Dict[str, Any]

    # Metrics and reporting
    metrics: Dict[str, Any]
    metrics_calculated: bool
    metrics_summary: Dict[str, Any]
    structured_json: Dict[str, Any]  # Formatted data for frontend
    report_markdown: Optional[str]  # LLM-generated report
    report_sections: Dict[str, Any]  # Structured report sections for frontend
    report_pdf_path: Optional[str]  # Final PDF location
    enhanced_report_content: Dict[str, str]  # In-memory report content (markdown, html, txt)
    report_generator_info: Dict[str, Any]  # Report generator metadata and info

    # Pipeline control
    pipeline_status: str
    current_node: Optional[str]
    completed_nodes: Annotated[List[str], merge_completed_nodes]  # Use Annotated for concurrent updates
    errors: List[Dict[str, Any]]
    warnings: List[Dict[str, Any]]

    # User preferences
    preferences: Optional[Dict[str, Any]]
    patient_data: Dict[str, Any]
    include_technical_details: bool
    language: str

    # Pipeline metadata
    pipeline_start_time: datetime
    pipeline_end_time: Optional[datetime]
    processing_time_seconds: Optional[float]

    # Internal control for merge logic
    _merge_call_count: Optional[int]
