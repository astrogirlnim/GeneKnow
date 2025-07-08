"""
State schema for GeneKnow LangGraph pipeline.
Defines the data structure passed between nodes.
"""
from typing import TypedDict, List, Dict, Optional, Any
from datetime import datetime


class GenomicState(TypedDict):
    """Main state object passed through the LangGraph pipeline."""
    
    # Input data
    file_path: str
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
    
    # Risk assessment
    risk_scores: Dict[str, float]  # {cancer_type: risk_percentage}
    risk_genes: Dict[str, List[str]]  # {cancer_type: [affected_genes]}
    
    # Report generation
    structured_json: Dict[str, Any]  # Formatted data for frontend
    report_markdown: Optional[str]  # LLM-generated report
    report_sections: Dict[str, Any]  # Structured report sections for frontend
    report_pdf_path: Optional[str]  # Final PDF location
    
    # Pipeline metadata
    pipeline_start_time: datetime
    pipeline_end_time: Optional[datetime]
    processing_time_seconds: Optional[float]
    
    # Error tracking
    errors: List[Dict[str, Any]]  # {node: str, error: str, timestamp: datetime}
    warnings: List[Dict[str, Any]]  # Non-fatal issues
    pipeline_status: str  # 'running', 'completed', 'failed', 'partial'
    
    # Node completion tracking
    completed_nodes: List[str]
    current_node: str
    
    # User preferences (from UI)
    language: str  # 'en', 'hi', 'es'
    include_technical_details: bool
    risk_threshold_percentage: float  # For highlighting high-risk findings 