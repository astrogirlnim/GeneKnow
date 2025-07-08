"""
Main LangGraph pipeline for GeneKnow genomic risk assessment.
Orchestrates the flow from FASTQ/BAM input to final PDF report.
"""
from langgraph.graph import StateGraph, END
from datetime import datetime
import logging

try:
    # Try relative imports first (when used as a package)
    from .state import GenomicState
    from .nodes import (
        file_input,
        preprocess,
        variant_calling,
        qc_filter,
        tcga_mapper,
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
        tcga_mapper,
        risk_model,
        formatter,
        report_writer
    )

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def create_genomic_pipeline() -> StateGraph:
    """
    Creates the main LangGraph pipeline for genomic analysis.
    
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
    workflow.add_node("tcga_mapper", tcga_mapper.process)
    workflow.add_node("risk_model", risk_model.process)
    workflow.add_node("formatter", formatter.process)
    workflow.add_node("report_writer", report_writer.process)
    
    # Define the entry point
    workflow.set_entry_point("file_input")
    
    # Add edges (linear flow with conditional routing)
    workflow.add_edge("file_input", "preprocess")
    
    # Conditional edge: Skip alignment if input is already BAM
    workflow.add_conditional_edges(
        "preprocess",
        lambda state: "variant_calling" if state["aligned_bam_path"] else "FAILED",
        {
            "variant_calling": "variant_calling",
            "FAILED": END
        }
    )
    
    workflow.add_edge("variant_calling", "qc_filter")
    workflow.add_edge("qc_filter", "tcga_mapper")
    workflow.add_edge("tcga_mapper", "risk_model")
    workflow.add_edge("risk_model", "formatter")
    workflow.add_edge("formatter", "report_writer")
    workflow.add_edge("report_writer", END)
    
    # Compile the graph
    # Note: Removing checkpointer for now as it requires additional configuration
    # To enable checkpointing, you'd need to provide thread_id in the config
    compiled = workflow.compile()
    
    return compiled


def run_pipeline(file_path: str, user_preferences: dict = None) -> dict:
    """
    Execute the genomic analysis pipeline.
    
    Args:
        file_path: Path to FASTQ or BAM file
        user_preferences: Optional user settings (language, detail level, etc.)
    
    Returns:
        Final state dictionary with all results
    """
    # Initialize the pipeline
    pipeline = create_genomic_pipeline()
    
    # Prepare initial state
    initial_state = {
        "file_path": file_path,
        "file_type": "unknown",  # Will be determined by file_input node
        "file_metadata": {},
        "aligned_bam_path": None,
        "vcf_path": None,
        "raw_variants": [],
        "filtered_variants": [],
        "variant_count": 0,
        "tcga_matches": {},
        "risk_scores": {},
        "risk_genes": {},
        "structured_json": {},
        "report_markdown": None,
        "report_pdf_path": None,
        "pipeline_start_time": datetime.now(),
        "pipeline_end_time": None,
        "processing_time_seconds": None,
        "errors": [],
        "warnings": [],
        "pipeline_status": "running",
        "completed_nodes": [],
        "current_node": "file_input",
        "language": user_preferences.get("language", "en") if user_preferences else "en",
        "include_technical_details": user_preferences.get("include_technical_details", True) if user_preferences else True,
        "risk_threshold_percentage": user_preferences.get("risk_threshold_percentage", 50.0) if user_preferences else 50.0
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