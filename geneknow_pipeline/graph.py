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
        tcga_mapper,
        cadd_scoring,
        feature_vector_builder,
        prs_calculator,
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
        tcga_mapper,
        cadd_scoring,
        feature_vector_builder,
        prs_calculator,
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
    state["current_node"] = "merge_parallel"
    
    # Check what we got from the parallel nodes
    has_raw_variants = bool(state.get("raw_variants"))
    has_filtered_variants = bool(state.get("filtered_variants"))
    
    # Ensure we have filtered_variants for downstream nodes
    if not has_filtered_variants:
        # If QC filter didn't run (no pre-existing variants), 
        # filter the newly called variants
        if has_raw_variants:
            logger.info("Applying QC filters to variants")
            from nodes.qc_filter import apply_qc_filters
            state["filtered_variants"] = apply_qc_filters(state["raw_variants"])
            state["variant_count"] = len(state["filtered_variants"])
        else:
            logger.warning("No variants found after merge")
            state["filtered_variants"] = []
            state["variant_count"] = 0
    else:
        # Update variant count if we have filtered variants
        state["variant_count"] = len(state["filtered_variants"])
    
    # Track completion of both parallel nodes
    completed_nodes = state.get("completed_nodes", [])
    if has_raw_variants and "variant_calling" not in completed_nodes:
        completed_nodes.append("variant_calling")
    if has_filtered_variants and "qc_filter" not in completed_nodes:
        completed_nodes.append("qc_filter")
    state["completed_nodes"] = completed_nodes
    
    # Log merge results
    logger.info(f"Merge complete: {state['variant_count']} filtered variants ready for population frequency mapping")
    
    return state


def merge_tcga_cadd_results(state: dict) -> dict:
    """
    Merge results from parallel TCGA mapping, CADD scoring, and PRS calculation.
    Ensures all three enrichment steps are complete before proceeding.
    """
    logger.info("Merging results from parallel TCGA mapping, CADD scoring, and PRS calculation")
    
    # Check what results we have from parallel execution
    # Note: In LangGraph, parallel nodes may have already updated the state
    tcga_complete = bool(state.get("tcga_matches"))
    cadd_complete = bool(state.get("cadd_enriched_variants"))
    prs_complete = bool(state.get("prs_results"))
    
    # Also check if we got PRS results that haven't been merged yet
    if not prs_complete and state.get("prs_summary"):
        prs_complete = True
    
    if tcga_complete and cadd_complete and prs_complete:
        logger.info("All three parallel processes completed successfully")
        logger.info(f"TCGA matches: {len(state.get('tcga_matches', {}))}")
        logger.info(f"CADD enriched variants: {len(state.get('cadd_enriched_variants', []))}")
        logger.info(f"PRS results: {len(state.get('prs_results', {}))}")
        
        # Merge the enriched variants from CADD scoring with TCGA data
        # The CADD scoring node produces cadd_enriched_variants
        # The TCGA mapping node produces tcga_matches but doesn't modify filtered_variants
        # The PRS calculator produces prs_results and prs_summary
        cadd_variants = state.get("cadd_enriched_variants", state.get("filtered_variants", []))
        tcga_matches = state.get("tcga_matches", {})
        
        # Add TCGA annotations to each variant
        merged_variants = []
        for variant in cadd_variants:
            # Create a copy to avoid modifying the original
            enriched_variant = variant.copy()
            variant_id = variant.get("variant_id", f"{variant['chrom']}:{variant['pos']}")
            
            # Add TCGA cancer relevance data from matches
            tcga_cancer_relevance = 0.0
            tcga_best_match = None
            
            # Check all cancer types for this variant
            for cancer_type, matches in tcga_matches.items():
                if variant_id in matches:
                    match_data = matches[variant_id]
                    enrichment = match_data.get("enrichment_score", 1.0)
                    
                    # Calculate cancer relevance (normalized enrichment)
                    relevance = min(enrichment / 10.0, 1.0)
                    if relevance > tcga_cancer_relevance:
                        tcga_cancer_relevance = relevance
                        tcga_best_match = {
                            "cancer_type": cancer_type,
                            "frequency": match_data.get("tumor_frequency", 0.0),
                            "enrichment": enrichment,
                            "sample_count": match_data.get("sample_count", 0),
                            "total_samples": match_data.get("total_samples", 1000)
                        }
            
            # Add TCGA data to the variant
            enriched_variant["tcga_cancer_relevance"] = tcga_cancer_relevance
            enriched_variant["tcga_best_match"] = tcga_best_match
            
            merged_variants.append(enriched_variant)
        
        # Ensure filtered_variants is updated with the merged data
        state["filtered_variants"] = merged_variants
        
        logger.info(f"Merged {len(merged_variants)} variants with TCGA and CADD annotations")
        logger.info(f"PRS calculation completed for {len(state['prs_results'])} cancer types")
        
        # Track completion of all three parallel nodes
        state["completed_nodes"] = state.get("completed_nodes", []) + ["tcga_mapper", "cadd_scoring", "prs_calculator"]
        
    else:
        logger.warning(f"Incomplete parallel processing - TCGA: {tcga_complete}, CADD: {cadd_complete}, PRS: {prs_complete}")
        if not tcga_complete:
            state["warnings"].append("TCGA mapping did not complete successfully")
        if not cadd_complete:
            state["warnings"].append("CADD scoring did not complete successfully")
            # If CADD failed, use filtered_variants as-is
            if not state.get("cadd_enriched_variants"):
                state["cadd_enriched_variants"] = state.get("filtered_variants", [])
        if not prs_complete:
            state["warnings"].append("PRS calculation did not complete successfully")
            # Add empty PRS results so pipeline can continue
            state["prs_results"] = {}
            state["prs_summary"] = {"overall_confidence": "failed"}
        
        # Track partial completion
        completed_nodes = state.get("completed_nodes", [])
        if tcga_complete and "tcga_mapper" not in completed_nodes:
            completed_nodes.append("tcga_mapper")
        if cadd_complete and "cadd_scoring" not in completed_nodes:
            completed_nodes.append("cadd_scoring")
        if prs_complete and "prs_calculator" not in completed_nodes:
            completed_nodes.append("prs_calculator")
        state["completed_nodes"] = completed_nodes
    
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
    Features parallel execution at two stages:
    1. Variant calling and QC filtering run in parallel
    2. TCGA mapping, CADD scoring, and PRS calculation run in parallel
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
    workflow.add_node("tcga_mapper", tcga_mapper.process)
    workflow.add_node("cadd_scoring", cadd_scoring.process)
    workflow.add_node("merge_tcga_cadd", merge_tcga_cadd_results)
    workflow.add_node("feature_vector_builder", feature_vector_builder.process)
    workflow.add_node("prs_calculator", prs_calculator.process)
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
    
    # Connect parallel paths and direct MAF path to population mapping
    workflow.add_edge("merge_parallel", "population_mapper")
    
    # Parallel execution of TCGA mapping and CADD scoring
    workflow.add_edge("population_mapper", "tcga_mapper")
    workflow.add_edge("population_mapper", "cadd_scoring")
    workflow.add_edge("population_mapper", "prs_calculator") # Added PRS calculator to parallel path
    
    # Both parallel paths converge at merge_tcga_cadd
    workflow.add_edge("tcga_mapper", "merge_tcga_cadd")
    workflow.add_edge("cadd_scoring", "merge_tcga_cadd")
    workflow.add_edge("prs_calculator", "merge_tcga_cadd") # Added PRS calculator to merge
    
    # Continue with feature vector building after merge
    workflow.add_edge("merge_tcga_cadd", "feature_vector_builder")
    
    # Feature vector builder feeds directly into risk model
    workflow.add_edge("feature_vector_builder", "risk_model")
    workflow.add_edge("risk_model", "formatter")
    
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
        "tcga_cohort_sizes": {},  # Initialize TCGA cohort sizes
        "cadd_enriched_variants": None,
        "cadd_stats": None,
        "prs_results": {},  # Initialize PRS results
        "prs_summary": {},  # Initialize PRS summary
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