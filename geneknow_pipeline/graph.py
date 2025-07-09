"""
Main LangGraph pipeline for GeneKnow genomic risk assessment.
Orchestrates the flow from FASTQ/BAM/VCF/MAF input to final PDF report.
"""
from langgraph.graph import StateGraph, END
from datetime import datetime
import logging
import os

try:
    # Try relative imports first (when used as a package)
    from .state import GenomicState
    from .nodes import (
        file_input,
        preprocess,
        variant_calling,
        qc_filter,
        population_mapper,
        cadd_scoring,
        feature_vector_builder,
        risk_model,
        formatter,
        report_writer
    )
except ImportError:
    # Fall back to absolute imports (when running directly)
    from state import GenomicState
    from nodes import (
        file_input,
        preprocess,
        variant_calling,
        qc_filter,
        population_mapper,
        cadd_scoring,
        feature_vector_builder,
        risk_model,
        formatter,
        report_writer
    )

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def merge_variants(state: dict) -> dict:
    """
    Merge results from parallel variant processing.
    Combines variants from both variant_calling and qc_filter paths.
    """
    logger.info("Merging results from parallel variant processing")
    
    # Ensure we have filtered_variants for downstream nodes
    if not state.get("filtered_variants"):
        # If QC filter didn't run (no pre-existing variants), 
        # filter the newly called variants
        if state.get("raw_variants"):
            logger.info("Applying QC filters to variants")
            from nodes.qc_filter import apply_qc_filters
            state["filtered_variants"] = apply_qc_filters(state["raw_variants"])
            state["variant_count"] = len(state["filtered_variants"])
        else:
            logger.warning("No variants found after merge")
            state["filtered_variants"] = []
            state["variant_count"] = 0
    
    # Log merge results
    logger.info(f"Merge complete: {state['variant_count']} filtered variants ready for population frequency mapping")
    
    return state


def route_after_preprocess(state: dict) -> list:
    """
    Determine which nodes to run after preprocessing.
    Returns a list to enable parallel execution.
    """
    # Check if this is a MAF file that already has filtered variants
    if state.get("file_type") == "maf" and state.get("filtered_variants"):
        logger.info("MAF file detected with pre-processed variants, skipping to population mapping")
        return ["population_mapper"]
    
    nodes_to_run = []
    
    # Always run variant calling if we have a BAM file
    if state.get("aligned_bam_path"):
        nodes_to_run.append("variant_calling")
    
    # Always run QC filter if we have any raw variants
    if state.get("raw_variants"):
        nodes_to_run.append("qc_filter")
    
    # If neither condition is met, we still need variant calling
    if not nodes_to_run:
        nodes_to_run.append("variant_calling")
    
    logger.info(f"Routing to nodes: {nodes_to_run}")
    return nodes_to_run


def create_genomic_pipeline() -> StateGraph:
    """
    Creates the main LangGraph pipeline for genomic analysis.
    Now with parallel execution of variant calling and QC filtering.
    Supports FASTQ, BAM, VCF, and MAF input files.
    
    Returns:
        StateGraph configured with all nodes and edges
    """
    # Initialize the graph with our state schema
    workflow = StateGraph(GenomicState)
    
    # Add all nodes to the graph
    workflow.add_node("file_input", file_input.process)
    workflow.add_node("preprocess", preprocess.process)
    workflow.add_node("variant_calling", variant_calling.process)
    workflow.add_node("qc_filter", qc_filter.process)
    workflow.add_node("merge_parallel", merge_variants)
    workflow.add_node("population_mapper", population_mapper.process)
    workflow.add_node("cadd_scoring", cadd_scoring.process)
    workflow.add_node("feature_vector_builder", feature_vector_builder.process)
    workflow.add_node("risk_model", risk_model.process)
    workflow.add_node("formatter", formatter.process)
    workflow.add_node("report_writer", report_writer.process)
    
    # Define the entry point
    workflow.set_entry_point("file_input")
    
    # Linear flow up to preprocess
    workflow.add_edge("file_input", "preprocess")
    
    # After preprocess, use conditional edges for parallel execution
    # This will run the nodes returned by route_after_preprocess in parallel
    workflow.add_conditional_edges(
        "preprocess",
        route_after_preprocess,
    )
    
    # Both parallel paths converge at merge_parallel
    workflow.add_edge("variant_calling", "merge_parallel")
    workflow.add_edge("qc_filter", "merge_parallel")
    
    # MAF files skip directly to population mapping (no merge needed)
    # Continue with linear flow after merge
    workflow.add_edge("merge_parallel", "population_mapper")
    
    # New flow with CADD scoring
    workflow.add_edge("population_mapper", "cadd_scoring")
    workflow.add_edge("cadd_scoring", "feature_vector_builder")
    
    # Check for legacy risk model flag
    use_legacy_risk = os.environ.get("USE_LEGACY_RISK", "false").lower() == "true"
    
    if use_legacy_risk:
        # Legacy path: feature_vector_builder -> risk_model -> formatter
        workflow.add_edge("feature_vector_builder", "risk_model") 
        workflow.add_edge("risk_model", "formatter")
    else:
        # New path: feature_vector_builder -> formatter (bypass risk_model)
        workflow.add_edge("feature_vector_builder", "formatter")
    
    workflow.add_edge("formatter", "report_writer")
    workflow.add_edge("report_writer", END)
    
    # Compile the graph
    compiled = workflow.compile()
    
    return compiled


def run_pipeline(file_path: str, user_preferences: dict = None) -> dict:
    """
    Execute the genomic analysis pipeline.
    
    Args:
        file_path: Path to FASTQ, BAM, VCF, or MAF file
        user_preferences: Optional user settings (language, detail level, etc.)
    
    Returns:
        Final state dictionary with all results
    """
    # Initialize the pipeline
    pipeline = create_genomic_pipeline()
    
    # Prepare initial state
    initial_state = {
        "file_path": file_path,
        "file_type": None,
        "file_metadata": {},
        "raw_variants": [],
        "filtered_variants": [],
        "variant_count": 0,
        "tcga_matches": {},
        "cadd_enriched_variants": None,
        "cadd_stats": None,
        "risk_scores": {},
        "structured_json": {},
        "report_markdown": "",
        "report_sections": {},  # Add this line
        "pipeline_status": "in_progress",
        "current_node": None,
        "completed_nodes": [],
        "errors": [],
        "warnings": [],
        "preferences": user_preferences,
        "patient_data": user_preferences.get("patient_data", {}),
        "include_technical_details": user_preferences.get("include_technical", True),
        "language": user_preferences.get("language", "en"),
        "pipeline_start_time": datetime.now()
    }
    
    try:
        # Run the pipeline
        logger.info(f"Starting genomic pipeline for: {file_path}")
        final_state = pipeline.invoke(initial_state)
        
        # Calculate total processing time
        final_state["pipeline_end_time"] = datetime.now()
        final_state["processing_time_seconds"] = (
            final_state["pipeline_end_time"] - final_state["pipeline_start_time"]
        ).total_seconds()
        
        logger.info(f"Pipeline completed in {final_state['processing_time_seconds']:.2f} seconds")
        
        # Debug logging
        logger.info(f"Final state keys: {list(final_state.keys())}")
        logger.info(f"Report sections in final state: {final_state.get('report_sections', 'NOT FOUND')}")
        
        return final_state
        
    except Exception as e:
        logger.error(f"Pipeline failed with error: {str(e)}")
        initial_state["pipeline_status"] = "failed"
        initial_state["errors"].append({
            "node": initial_state["current_node"],
            "error": str(e),
            "timestamp": datetime.now()
        })
        return initial_state


if __name__ == "__main__":
    # Test the pipeline with a sample file
    # ⚠️ MOCK DATA - Replace with real file path for testing
    test_file = "test-data/sample.fastq"
    result = run_pipeline(test_file)
    print(f"Pipeline status: {result['pipeline_status']}")
    print(f"Completed nodes: {result['completed_nodes']}") 