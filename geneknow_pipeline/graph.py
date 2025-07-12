"""
Main LangGraph pipeline for GeneKnow genomic risk assessment.
Orchestrates the flow from FASTQ/BAM/VCF/MAF input to final PDF report.
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
        population_mapper,
        tcga_mapper,
        cadd_scoring,
        clinvar_annotator,
        feature_vector_builder,
        prs_calculator,
        pathway_burden,
        ml_fusion_node,
        risk_model,
        shap_validator,
        metrics_calculator,
        formatter,
        report_writer,
        # New nodes for clinical view data
        mutation_classifier,
        mutational_signatures,
        variant_transformer,
        structural_variant_detector,
        cnv_detector,
        pathway_analyzer,
        survival_analyzer,
        clinical_recommendations,
    )
    from .nodes.report_generator import process as report_generator_process
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
        clinvar_annotator,
        feature_vector_builder,
        prs_calculator,
        pathway_burden,
        ml_fusion_node,
        risk_model,
        shap_validator,
        metrics_calculator,
        formatter,
        report_writer,
        # New nodes for clinical view data
        mutation_classifier,
        mutational_signatures,
        variant_transformer,
        structural_variant_detector,
        cnv_detector,
        pathway_analyzer,
        survival_analyzer,
        clinical_recommendations,
    )
    from nodes.report_generator import process as report_generator_process

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


def merge_static_model_results(state: dict) -> dict:
    """
    Merge results from parallel static model execution.

    This function is called by LangGraph each time a parallel node completes.
    We need to accumulate results and only perform the full merge when all nodes are done.

    IMPORTANT: Due to LangGraph's execution model, we must always preserve existing
    state data to avoid losing results from parallel nodes.
    """
    # Track how many times this has been called
    call_count = state.get("_merge_call_count", 0) + 1
    logger.info(f"Merge static models called (call #{call_count})")

    # Check which nodes have completed by checking for their data
    # Note: None means not run yet, {} or other values mean the node has run
    has_tcga = state.get("tcga_matches") is not None
    has_cadd = state.get("cadd_stats") is not None
    has_clinvar = state.get("clinvar_annotations") is not None
    has_prs = state.get("prs_results") is not None
    has_pathway = state.get("pathway_burden_results") is not None

    completed_count = sum([has_tcga, has_cadd, has_clinvar, has_prs, has_pathway])
    logger.info(
        f"Parallel nodes completed: {completed_count}/5 (TCGA={has_tcga}, CADD={has_cadd}, ClinVar={has_clinvar}, PRS={has_prs}, Pathway={has_pathway})"
    )

    # If not all nodes have completed, preserve any existing data
    if completed_count < 5:
        logger.info(f"Waiting for {5 - completed_count} more parallel node(s) to complete")
        # CRITICAL: Return existing state data to preserve results from completed nodes
        result = {"_merge_call_count": call_count}

        # Preserve all existing parallel node results
        if state.get("tcga_matches") is not None:
            result["tcga_matches"] = state["tcga_matches"]
        if state.get("tcga_summary") is not None:
            result["tcga_summary"] = state["tcga_summary"]
        if state.get("cadd_stats") is not None:
            result["cadd_stats"] = state["cadd_stats"]
        if state.get("cadd_enriched_variants") is not None:
            result["cadd_enriched_variants"] = state["cadd_enriched_variants"]
        if state.get("clinvar_annotations") is not None:
            result["clinvar_annotations"] = state["clinvar_annotations"]
        if state.get("clinvar_stats") is not None:
            result["clinvar_stats"] = state["clinvar_stats"]
        # Preserve all ClinVar-specific keys
        for key in [
            "clinvar_pathogenic_variants",
            "clinvar_likely_pathogenic_variants",
            "clinvar_benign_variants",
            "clinvar_vus_variants",
            "clinvar_drug_response_variants",
            "clinvar_risk_factor_variants",
        ]:
            if state.get(key) is not None:
                result[key] = state[key]
        if state.get("prs_results") is not None:
            result["prs_results"] = state["prs_results"]
        if state.get("prs_summary") is not None:
            result["prs_summary"] = state["prs_summary"]
        if state.get("pathway_burden_results") is not None:
            result["pathway_burden_results"] = state["pathway_burden_results"]
        if state.get("pathway_burden_summary") is not None:
            result["pathway_burden_summary"] = state["pathway_burden_summary"]
        if state.get("pathway_enriched_variants") is not None:
            result["pathway_enriched_variants"] = state["pathway_enriched_variants"]

        return result

    # All nodes have completed - check if we've already done the full merge
    if "merge_static_models" in state.get("completed_nodes", []):
        logger.info("Full merge already completed, skipping")
        return {}  # Return empty dict to avoid duplicate processing

    # Perform the full merge
    logger.info("All 5 parallel nodes completed - performing full merge")

    # Log what we received
    logger.info(f"TCGA matches: {len(state.get('tcga_matches', {}))}")
    logger.info(f"CADD enriched variants: {len(state.get('cadd_enriched_variants', []))}")
    logger.info(f"ClinVar annotations: {len(state.get('clinvar_annotations', {}))}")
    logger.info(f"PRS results: {len(state.get('prs_results', {}))}")
    logger.info(f"Pathway burden results: {len(state.get('pathway_burden_results', {}))}")

    # Start with the most enriched version of variants
    # Priority: pathway_enriched_variants > cadd_enriched_variants > filtered_variants
    if state.get("pathway_enriched_variants"):
        base_variants = state["pathway_enriched_variants"]
        logger.info("Using pathway enriched variants as base")
    elif state.get("cadd_enriched_variants"):
        base_variants = state["cadd_enriched_variants"]
        logger.info("Using CADD enriched variants as base")
    else:
        base_variants = state.get("filtered_variants", [])
        logger.info("Using filtered variants as base")

    # Create a mapping of variant_id to variant for easier lookup
    variant_map = {}
    for variant in base_variants:
        variant_id = variant.get("variant_id", f"{variant['chrom']}:{variant['pos']}")
        variant_map[variant_id] = variant.copy()  # Make a copy to avoid modifying originals

    logger.info(f"Base variants for merging: {len(base_variants)}")
    logger.info(f"Variant map size: {len(variant_map)}")

    # If pathway burden returned modified variants, merge in the pathway_damage_assessment
    if "pathway_enriched_variants" in state:
        pathway_merged_count = 0
        for variant in state["pathway_enriched_variants"]:
            variant_id = variant.get("variant_id", f"{variant['chrom']}:{variant['pos']}")
            if variant_id in variant_map and "pathway_damage_assessment" in variant:
                variant_map[variant_id]["pathway_damage_assessment"] = variant["pathway_damage_assessment"]
                pathway_merged_count += 1
        logger.info(f"Merged pathway_damage_assessment for {pathway_merged_count} variants")

    # Check for pathway_damage_assessment presence after merge
    variants_with_pathway_damage = sum(1 for v in variant_map.values() if "pathway_damage_assessment" in v)
    logger.info(f"Total variants with pathway_damage_assessment: {variants_with_pathway_damage}")

    tcga_matches = state.get("tcga_matches", {})
    clinvar_annotations = state.get("clinvar_annotations", {})

    # Add TCGA and ClinVar annotations to each variant
    merged_variants = []
    for variant_id, variant in variant_map.items():
        # variant is already a copy from variant_map
        enriched_variant = variant

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
                        "total_samples": match_data.get("total_samples", 1000),
                    }

        # Add TCGA data to the variant
        enriched_variant["tcga_cancer_relevance"] = tcga_cancer_relevance
        enriched_variant["tcga_best_match"] = tcga_best_match

        # Add ClinVar data to the variant
        clinvar_data = clinvar_annotations.get(variant_id, {})
        if clinvar_data.get("found_in_clinvar"):
            enriched_variant["clinvar_clinical_significance"] = clinvar_data.get("clinical_significance")
            enriched_variant["clinvar_risk_score"] = clinvar_data.get("clinical_risk_score", 0.0)
            enriched_variant["clinvar_interpretation"] = clinvar_data.get("clinical_interpretation")
            enriched_variant["clinvar_confidence"] = clinvar_data.get("confidence")
            enriched_variant["clinvar_actionability"] = clinvar_data.get("actionability")
            enriched_variant["clinvar_recommendation"] = clinvar_data.get("recommendation")
            enriched_variant["clinvar_condition"] = clinvar_data.get("condition")
            enriched_variant["clinvar_cancer_related"] = clinvar_data.get("cancer_related", False)
        else:
            enriched_variant["clinvar_clinical_significance"] = None
            enriched_variant["clinvar_risk_score"] = 0.0
            enriched_variant["clinvar_interpretation"] = "no_clinical_evidence"
            enriched_variant["clinvar_confidence"] = "none"
            enriched_variant["clinvar_actionability"] = "none"
            enriched_variant["clinvar_recommendation"] = "No clinical evidence available"
            enriched_variant["clinvar_condition"] = None
            enriched_variant["clinvar_cancer_related"] = False

        # pathway_damage_assessment is already included from the variant_map merge above

        merged_variants.append(enriched_variant)

    # Update file_metadata with summary statistics from all parallel nodes
    file_metadata = state.get("file_metadata", {})

    # Add TCGA summary
    if "tcga_summary" in state:
        file_metadata["tcga_summary"] = state["tcga_summary"]

    # Add ClinVar summary
    if "clinvar_stats" in state:
        file_metadata["clinvar_summary"] = state["clinvar_stats"]

    # Add CADD summary
    if "cadd_stats" in state:
        file_metadata["cadd_summary"] = state["cadd_stats"]

    # Add pathway burden summary
    if "pathway_burden_summary" in state:
        file_metadata["pathway_burden_summary"] = state["pathway_burden_summary"]

    logger.info(f"Merged {len(merged_variants)} variants with TCGA, CADD, ClinVar, and pathway annotations")
    logger.info(f"PRS calculation completed for {len(state.get('prs_results', {}))} cancer types")

    # Track completion
    completed = state.get("completed_nodes", [])
    for node in [
        "tcga_mapper",
        "cadd_scoring",
        "clinvar_annotator",
        "prs_calculator",
        "pathway_burden",
        "merge_static_models",
    ]:
        if node not in completed:
            completed.append(node)

    logger.info("Merge static models completed successfully")
    logger.info(f"Updated completed nodes: {completed}")

    # Return the full merge results
    result = {"_merge_call_count": call_count}  # Return the call count update
    result.update(
        {
            "filtered_variants": merged_variants,
            "file_metadata": file_metadata,
            "completed_nodes": completed,
        }
    )

    return result


def merge_variant_analysis_results(state: dict) -> dict:
    """
    Merge results from structural variant and CNV detection nodes.
    """
    logger.info("Merging variant analysis results")
    
    # Check if we have results from both nodes
    structural_variants = state.get("structural_variants", [])
    copy_number_variants = state.get("copy_number_variants", [])
    
    # Create summary statistics
    variant_analysis_summary = {
        "structural_variants_found": len(structural_variants),
        "copy_number_variants_found": len(copy_number_variants),
        "total_genomic_alterations": len(structural_variants) + len(copy_number_variants)
    }
    
    # Combine into genomic alterations
    genomic_alterations = {
        "structural_variants": structural_variants,
        "copy_number_variants": copy_number_variants,
        "summary": variant_analysis_summary
    }
    
    # Track completion
    completed = state.get("completed_nodes", [])
    if "merge_variant_analysis" not in completed:
        completed.append("merge_variant_analysis")
    
    logger.info(f"Merged variant analysis: {len(structural_variants)} SVs, {len(copy_number_variants)} CNVs")
    
    return {
        "genomic_alterations": genomic_alterations,
        "completed_nodes": completed
    }


def route_after_preprocess(state: dict) -> list:
    """
    Determine which nodes to run after preprocessing.
    Returns a list to enable parallel execution.
    """
    # Check if this is a MAF file that already has filtered variants
    if state.get("file_type") == "ma" and state.get("filtered_variants"):
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


def check_merge_complete(state: dict) -> str:
    """
    Check if the merge_static_models has completed processing all parallel nodes.
    Returns either "continue" to proceed or "wait" to re-run merge.
    """
    # Check if all parallel nodes have completed by checking for their data
    has_tcga = state.get("tcga_matches") is not None
    has_cadd = state.get("cadd_stats") is not None
    has_clinvar = state.get("clinvar_annotations") is not None
    has_prs = state.get("prs_results") is not None
    has_pathway = state.get("pathway_burden_results") is not None

    # Also check if merge has marked itself complete
    merge_complete = "merge_static_models" in state.get("completed_nodes", [])

    if all([has_tcga, has_cadd, has_clinvar, has_prs, has_pathway]) or merge_complete:
        logger.info("All parallel nodes complete, proceeding to feature vector builder")
        return "continue"
    else:
        missing = []
        if not has_tcga:
            missing.append("tcga")
        if not has_cadd:
            missing.append("cadd")
        if not has_clinvar:
            missing.append("clinvar")
        if not has_prs:
            missing.append("prs")
        if not has_pathway:
            missing.append("pathway_burden")
        logger.info(f"Waiting for parallel nodes to complete. Missing: {missing}")
        return "wait"


def merge_variant_analysis_results(state: dict) -> dict:
    """
    Merge results from parallel structural variant and CNV detection.
    """
    logger.info("Merging structural variant and CNV results")
    state["current_node"] = "merge_variant_analysis"

    # Check if both nodes have completed
    has_svs = state.get("structural_variants") is not None
    has_cnvs = state.get("copy_number_variants") is not None

    # Track completion
    completed_nodes = state.get("completed_nodes", [])

    if has_svs and has_cnvs:
        logger.info("Both SV and CNV detection completed")

        # Combine all genomic alterations
        all_alterations = {
            "structural_variants": state.get("structural_variants", []),
            "copy_number_variants": state.get("copy_number_variants", []),
            "total_alterations": len(state.get("structural_variants", [])) + len(state.get("copy_number_variants", [])),
        }

        # Update completed nodes
        if "merge_variant_analysis" not in completed_nodes:
            completed_nodes.append("merge_variant_analysis")

        logger.info(f"Merged {all_alterations['total_alterations']} genomic alterations")

        return {"genomic_alterations": all_alterations, "completed_nodes": completed_nodes}
    else:
        logger.info(f"Waiting for nodes to complete (SV={has_svs}, CNV={has_cnvs})")
        return {}


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
    workflow.add_node("clinvar_annotator", clinvar_annotator.process)
    workflow.add_node("pathway_burden", pathway_burden.process)
    workflow.add_node("merge_static_models", merge_static_model_results)
    workflow.add_node("feature_vector_builder", feature_vector_builder.process)
    workflow.add_node("prs_calculator", prs_calculator.process)
    workflow.add_node("ml_fusion", ml_fusion_node.process)
    workflow.add_node("risk_model", risk_model.process)
    workflow.add_node("shap_validator", shap_validator.process)

    # Add new clinical view nodes
    workflow.add_node("mutation_classifier", mutation_classifier.process)
    workflow.add_node("mutational_signatures", mutational_signatures.process)
    workflow.add_node("variant_transformer", variant_transformer.process)
    workflow.add_node("structural_variant_detector", structural_variant_detector.process)
    workflow.add_node("cnv_detector", cnv_detector.process)
    # REMOVED: pathway_analyzer node - using pathway_burden results directly
    workflow.add_node("survival_analyzer", survival_analyzer.process)
    workflow.add_node("clinical_recommendations", clinical_recommendations.process)

    workflow.add_node("metrics_calculator", metrics_calculator.process)
    workflow.add_node("formatter", formatter.process)
    workflow.add_node("report_generator", report_generator_process)

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

    # Parallel execution of TCGA mapping, CADD scoring, ClinVar annotation, PRS calculation, and pathway burden
    workflow.add_edge("population_mapper", "tcga_mapper")
    workflow.add_edge("population_mapper", "cadd_scoring")
    workflow.add_edge("population_mapper", "clinvar_annotator")
    workflow.add_edge("population_mapper", "prs_calculator")
    workflow.add_edge("population_mapper", "pathway_burden")

    # All parallel paths converge at merge_static_models
    workflow.add_edge("tcga_mapper", "merge_static_models")
    workflow.add_edge("cadd_scoring", "merge_static_models")
    workflow.add_edge("clinvar_annotator", "merge_static_models")
    workflow.add_edge("prs_calculator", "merge_static_models")
    workflow.add_edge("pathway_burden", "merge_static_models")

    # Continue with feature vector building after merge
    workflow.add_edge("merge_static_models", "feature_vector_builder")

    # Feature vector builder -> ML fusion -> risk model -> SHAP validator
    workflow.add_edge("feature_vector_builder", "ml_fusion")
    workflow.add_edge("ml_fusion", "risk_model")
    workflow.add_edge("risk_model", "shap_validator")

    # After SHAP validator, run new clinical analysis nodes in sequence
    workflow.add_edge("shap_validator", "mutation_classifier")
    workflow.add_edge("mutation_classifier", "mutational_signatures")
    workflow.add_edge("mutational_signatures", "variant_transformer")

    # Run structural variant and CNV detection in parallel
    workflow.add_edge("variant_transformer", "structural_variant_detector")
    workflow.add_edge("variant_transformer", "cnv_detector")

    # Create a merge point for structural variant and CNV results
    workflow.add_node("merge_variant_analysis", merge_variant_analysis_results)
    workflow.add_edge("structural_variant_detector", "merge_variant_analysis")
    workflow.add_edge("cnv_detector", "merge_variant_analysis")

    # Continue with pathway analysis after merging
    workflow.add_edge("merge_variant_analysis", "survival_analyzer")
    workflow.add_edge("survival_analyzer", "clinical_recommendations")

    # Clinical recommendations -> metrics calculator -> formatter -> report writer
    workflow.add_edge("clinical_recommendations", "metrics_calculator")
    workflow.add_edge("metrics_calculator", "formatter")

    workflow.add_edge("formatter", "report_generator")
    workflow.add_edge("report_generator", END)

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

    # Clean preferences to remove non-serializable values (like callbacks)
    if user_preferences is None:
        user_preferences = {}

    # Create a copy and remove any non-serializable values
    clean_preferences = {}
    for key, value in user_preferences.items():
        # Skip function objects and other non-serializable types
        if not callable(value):
            clean_preferences[key] = value

    # Prepare initial state
    initial_state = {
        "file_path": file_path,
        "file_type": None,
        "file_metadata": {},
        "raw_variants": [],
        "filtered_variants": [],
        "variant_count": 0,
        "tcga_matches": None,  # Use None to detect when node hasn't run
        "tcga_cohort_sizes": {},  # Keep as {} for metadata
        "cadd_enriched_variants": None,
        "cadd_stats": None,
        "clinvar_annotations": None,  # Use None to detect when node hasn't run
        "clinvar_stats": None,
        "prs_results": None,  # Use None to detect when node hasn't run
        "prs_summary": None,
        "pathway_burden_results": None,  # Use None to detect when node hasn't run
        "pathway_burden_summary": None,
        "pathway_enriched_variants": None,
        "ml_ready_variants": None,  # Variants prepared for ML fusion
        "ml_fusion_results": None,  # ML fusion model predictions
        "ml_fusion_model_instance": None,  # ML fusion model for SHAP
        "ml_fusion_feature_matrix": None,  # Feature matrix for SHAP
        "shap_validation_status": None,  # SHAP validation status
        "shap_validation_reasons": None,  # Reasons for flagging
        "shap_top_contributors": None,  # Top contributing features
        "shap_feature_importance": None,  # Full SHAP values
        "shap_validation_details": None,  # Detailed validation info
        # NEW: Initialize state fields for new nodes
        "classified_variants": None,  # Mutation classifier
        "mutation_type_distribution": None,
        "mutational_signatures": None,  # Mutational signatures
        "mutational_signatures_summary": None,
        "variant_details": None,  # Variant transformer
        "variant_transformation_summary": None,
        "structural_variants": None,  # Structural variant detector
        "structural_variant_summary": None,
        "copy_number_variants": None,  # CNV detector
        "cnv_summary": None,
        "pathway_analysis": None,  # Pathway analyzer
        "survival_analysis": None,  # Survival analyzer
        "clinical_recommendations": None,  # Clinical recommendations
        "genomic_alterations": None,  # From merge_variant_analysis
        "risk_scores": {},
        "risk_details": {},  # Initialize risk details for metrics
        "ml_risk_assessment": {},  # Initialize ML risk assessment for metrics
        "metrics": {},  # Initialize metrics
        "metrics_calculated": False,  # Track if metrics have been calculated
        "metrics_summary": {},  # Initialize metrics summary
        "structured_json": {},
        "report_markdown": "",
        "report_sections": {},  # Add this line
        "enhanced_report_content": {},  # For in-memory report content (HIPAA compliant)
        "report_generator_info": {},  # For report generator
        "pipeline_status": "in_progress",
        "current_node": None,
        "completed_nodes": [],
        "errors": [],
        "warnings": [],
        "preferences": clean_preferences,  # Use cleaned preferences
        "patient_data": clean_preferences.get("patient_data", {}),
        "include_technical_details": clean_preferences.get("include_technical", True),
        "language": clean_preferences.get("language", "en"),
        "pipeline_start_time": datetime.now(),
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
        initial_state["errors"].append(
            {
                "node": initial_state["current_node"],
                "error": str(e),
                "timestamp": datetime.now(),
            }
        )
        return initial_state


if __name__ == "__main__":
    # Test the pipeline with a sample file
    # ⚠️ MOCK DATA - Replace with real file path for testing
    test_file = "test-data/sample.fastq"
    result = run_pipeline(test_file)
    print(f"Pipeline status: {result['pipeline_status']}")
    print(f"Completed nodes: {result['completed_nodes']}")
