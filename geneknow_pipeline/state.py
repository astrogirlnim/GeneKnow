"""
State schema for GeneKnow LangGraph pipeline.
Defines the data structure passed between nodes.
"""

from typing import TypedDict, List, Dict, Optional, Any, Annotated
from datetime import datetime
import operator


class GenomicState(TypedDict):
    """Main state object passed through the LangGraph pipeline."""

    # Input data
    file_path: str  # Read-only input field
    file_type: str  # 'fastq' or 'bam'
    file_metadata: Dict[str, Any]  # size, read count, etc.

    # Processing intermediates
    aligned_bam_path: Optional[str]  # Path to aligned BAM (if FASTQ input)
    vcf_path: Optional[str]  # Path to DeepVariant output

    # Variant data
    raw_variants: List[Dict[str, Any]]  # Unfiltered variants from VCF
    filtered_variants: List[Dict[str, Any]]  # Post-QC variants
    variant_count: int

    # TCGA matching results
    tcga_matches: Dict[str, Dict[str, Any]]  # {cancer_type: {variant_id: match_info}}
    tcga_cohort_sizes: Dict[str, int]  # {cancer_type: sample_count}
    tcga_summary: Optional[Dict[str, Any]]  # Summary statistics from TCGA mapping

    # CADD enrichment results
    cadd_enriched_variants: Optional[
        List[Dict[str, Any]]
    ]  # Variants with CADD annotations
    cadd_stats: Optional[Dict[str, Any]]  # CADD scoring statistics

    # ClinVar annotation results
    clinvar_annotations: Dict[str, Dict[str, Any]]  # {variant_id: clinvar_data}
    clinvar_stats: Optional[Dict[str, Any]]  # ClinVar annotation statistics
    clinvar_pathogenic_variants: Optional[List[Dict[str, Any]]]  # Pathogenic variants
    clinvar_likely_pathogenic_variants: Optional[
        List[Dict[str, Any]]
    ]  # Likely pathogenic variants
    clinvar_benign_variants: Optional[List[Dict[str, Any]]]  # Benign variants
    clinvar_vus_variants: Optional[
        List[Dict[str, Any]]
    ]  # Variants of uncertain significance

    # PRS (Polygenic Risk Score) results
    prs_results: Dict[str, Dict[str, Any]]  # {cancer_type: prs_data}
    prs_summary: Dict[str, Any]  # Overall PRS summary

    # ML Fusion results
    ml_ready_variants: Optional[List[Dict[str, Any]]]  # Variants prepared for ML fusion
    ml_fusion_results: Optional[Dict[str, Any]]  # ML fusion model predictions
    ml_fusion_model_instance: Optional[Any]  # ML fusion model for SHAP
    ml_fusion_feature_matrix: Optional[Any]  # Feature matrix for SHAP

    # SHAP validation results
    shap_validation_status: Optional[str]  # "PASS", "FLAG_FOR_REVIEW", "ERROR", "SKIPPED"
    shap_validation_reasons: Optional[List[str]]  # Reasons for flagging
    shap_top_contributors: Optional[List[Dict[str, Any]]]  # Top 3 features with names & contributions
    shap_feature_importance: Optional[Dict[str, float]]  # Full SHAP values for detailed analysis
    shap_validation_details: Optional[Dict[str, Any]]  # Detailed validation info

    # Risk assessment
    risk_scores: Dict[str, float]  # {cancer_type: risk_percentage}
    risk_genes: Dict[str, List[str]]  # {cancer_type: [affected_genes]}

    # Report generation
    structured_json: Dict[str, Any]  # Formatted data for frontend
    report_markdown: Optional[str]  # LLM-generated report
    report_sections: Dict[str, Any]  # Structured report sections for frontend
    report_pdf_path: Optional[str]  # Final PDF location
    enhanced_report_content: Dict[str, str]  # In-memory report content (markdown, html, txt)
    report_generator_info: Dict[str, Any]  # Report generator metadata and info

    # Pipeline metadata
    pipeline_start_time: datetime
    pipeline_end_time: Optional[datetime]
    processing_time_seconds: Optional[float]

    # Error tracking
    errors: Annotated[
        List[Dict[str, Any]], operator.add
    ]  # {node: str, error: str, timestamp: datetime}
    warnings: Annotated[List[Dict[str, Any]], operator.add]  # Non-fatal issues
    pipeline_status: str  # 'running', 'completed', 'failed', 'partial'

    # Node completion tracking
    completed_nodes: List[str]  # Remove operator.add to prevent duplicates
    current_node: str

    # User preferences (from UI)
    language: str  # 'en', 'hi', 'es'
    include_technical_details: bool
    risk_threshold_percentage: float  # For highlighting high-risk findings
