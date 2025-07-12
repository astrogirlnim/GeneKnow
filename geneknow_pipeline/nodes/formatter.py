"""
Formatter node.
Structures all results into JSON format for frontend consumption.
"""

import logging
from datetime import datetime
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)


def _ensure_structural_variants_format(structural_variants: List[Dict]) -> List[Dict]:
    """Ensure all structural variants have the genes_affected field as an array"""
    formatted_svs = []
    
    for sv in structural_variants:
        # Create a copy to avoid modifying the original
        formatted_sv = sv.copy()
        
        # Ensure genes_affected is always an array
        if "genes_affected" not in formatted_sv:
            # Try to extract genes from other fields
            genes = []
            
            # Check various fields where genes might be stored
            if "gene" in formatted_sv:
                genes.append(formatted_sv["gene"])
            elif "genes" in formatted_sv:
                if isinstance(formatted_sv["genes"], list):
                    genes.extend(formatted_sv["genes"])
                else:
                    genes.append(formatted_sv["genes"])
            elif "fusion_name" in formatted_sv:
                # For fusions, try to extract gene names from fusion name
                fusion_name = formatted_sv["fusion_name"]
                if "-" in fusion_name:
                    genes.extend(fusion_name.split("-"))
                else:
                    genes.append(fusion_name)
            
            # If no genes found, use "Unknown"
            if not genes:
                genes = ["Unknown"]
            
            formatted_sv["genes_affected"] = genes
        elif not isinstance(formatted_sv["genes_affected"], list):
            # Ensure it's a list
            formatted_sv["genes_affected"] = [formatted_sv["genes_affected"]]
        
        formatted_svs.append(formatted_sv)
    
    return formatted_svs


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Format all pipeline results into structured JSON.

    Updates state with:
    - structured_json: formatted data ready for frontend/report
    """
    logger.info("Starting result formatting")
    # Note: Don't set current_node to avoid concurrent updates

    try:
        # Extract high-risk findings if risk scores are available
        high_risk_findings = []
        risk_threshold = state.get("risk_threshold_percentage", 50.0)

        risk_scores = state.get("risk_scores", {})
        risk_genes = state.get("risk_genes", {})

        # Only process risk findings if risk assessment was performed
        if risk_scores:
            for cancer_type, risk_score in risk_scores.items():
                if risk_score >= risk_threshold:
                    high_risk_findings.append(
                        {
                            "cancer_type": cancer_type,
                            "risk_percentage": risk_score,
                            "affected_genes": risk_genes.get(cancer_type, []),
                        }
                    )

        # Format variant summary with TCGA matches and transformations
        variant_summary = []

        # Use variant_details if available (from variant_transformer), otherwise use filtered_variants
        variants_to_format = state.get("variant_details") or state.get(
            "filtered_variants", []
        )

        for variant in variants_to_format:
            variant_info = {
                "gene": variant.get("gene", "Unknown"),
                "variant": variant.get("variant_id", "Unknown"),
                "variant_id": variant.get("variant_id", "Unknown"),
                "consequence": variant.get(
                    "consequence", variant.get("variant_classification", "unknown")
                ),
                "mutation_type": variant.get(
                    "mutation_type", "snv"
                ),  # From mutation_classifier
                "hgvs_c": variant.get("hgvs_c", ""),
                "hgvs_p": variant.get("hgvs_p", variant.get("protein_change", "")),
                "protein_change": variant.get(
                    "protein_change", ""
                ),  # From variant_transformer
                "functional_impact": variant.get("functional_impact", "Unknown"),
                "quality_metrics": {
                    "quality": variant.get("quality", variant.get("qual", 0)),
                    "depth": variant.get("depth", 0),
                    "allele_freq": variant.get("allele_freq", 0),
                },
                "tcga_matches": {},
            }

            # Add transformation data if available
            if "transformation" in variant:
                variant_info["transformation"] = variant["transformation"]

            # Add clinical significance if available (from MAF)
            if "clinical_significance" in variant:
                variant_info["clinical_significance"] = variant["clinical_significance"]

            # Add CADD scores if available
            if "cadd_phred" in variant:
                variant_info["cadd_scores"] = {
                    "phred": variant.get("cadd_phred"),
                    "raw": variant.get("cadd_raw"),
                    "risk_weight": variant.get("cadd_risk_weight"),
                }

            # Add TCGA match info
            tcga_matches = state.get("tcga_matches", {})
            for cancer_type, matches in tcga_matches.items():
                if variant.get("variant_id") in matches:
                    match_info = matches[variant["variant_id"]]
                    variant_info["tcga_matches"][cancer_type] = {
                        "frequency": match_info.get("frequency", 0),
                        "frequency_percent": match_info.get("frequency", 0) * 100,
                        "patient_fraction": f"{match_info['sample_count']}/{match_info['total_samples']}",
                    }

            variant_summary.append(variant_info)

        # Build structured JSON
        # Calculate processing time if not set
        processing_time = None
        if state.get("pipeline_start_time"):
            processing_time = (
                datetime.now() - state["pipeline_start_time"]
            ).total_seconds()

        # Build TCGA summary
        tcga_matches = state.get("tcga_matches", {})
        tcga_summary = {
            "cancer_types_analyzed": list(tcga_matches.keys()),
            "cohort_sizes": state.get("tcga_cohort_sizes", {}),
            "variants_with_tcga_data": len(
                [v for v in variant_summary if v["tcga_matches"]]
            ),
        }

        # Include CADD statistics if available
        cadd_stats = state.get("cadd_stats", {})
        if cadd_stats and "variants_scored" in cadd_stats:
            cadd_summary = {
                "variants_scored": cadd_stats.get("variants_scored", 0),
                "mean_phred": cadd_stats.get("mean_phred", 0),
                "max_phred": cadd_stats.get("max_phred", 0),
                "high_impact_variants": cadd_stats.get("variants_gt20", 0),
                "cancer_gene_variants": cadd_stats.get("variants_in_cancer_genes", 0),
            }
        else:
            cadd_summary = None

        # Include SHAP validation results if available
        shap_validation = None
        if state.get("shap_validation_status"):
            shap_validation = {
                "status": state.get("shap_validation_status", "SKIPPED"),
                "reasons": state.get("shap_validation_reasons", []),
                "top_contributors": state.get("shap_top_contributors", []),
                "feature_importance": state.get("shap_feature_importance", {}),
                "details": state.get("shap_validation_details", {}),
            }

        # Build comprehensive summary including mutation types
        # SIMPLIFIED: Use filtered_variants as the source of truth
        filtered_variants = state.get("filtered_variants", [])
        variant_count = len(filtered_variants)
        
        # Get mutation type distribution if available
        mutation_type_dist = state.get("mutation_type_distribution")
        
        # Add debug logging to track data flow
        print(f"\n=== FORMATTER DEBUG ===")
        print(f"State keys: {list(state.keys())}")
        print(f"filtered_variants length: {variant_count}")
        print(f"mutation_type_distribution: {mutation_type_dist}")
        print(f"variant_count from state: {state.get('variant_count', 0)}")
        
        # Use consistent logic: filtered_variants is the source of truth
        # Since QC filtering is disabled, total_variants_found equals variants_passed_qc
        total_variants_found = variant_count
        variants_passed_qc = variant_count
        
        print(f"Calculated total_variants_found: {total_variants_found}")
        print(f"Calculated variants_passed_qc: {variants_passed_qc}")
        print(f"=== END FORMATTER DEBUG ===\n")
        
        summary = {
            "total_variants_found": total_variants_found,
            "variants_passed_qc": variants_passed_qc,
            "high_risk_findings": len(high_risk_findings),
        }
        
        # Add mutation type distribution if available
        if mutation_type_dist:
            summary["mutation_types"] = mutation_type_dist

        # Process pathway analysis data from pathway burden results
        pathway_analysis = None
        pathway_burden_results = state.get("pathway_burden_results", {})
        pathway_burden_summary = state.get("pathway_burden_summary", {})

        # DEBUG: Write debug info to file
        with open("formatter_debug.log", "a") as f:
            f.write(f"\n=== FORMATTER DEBUG ===\n")
            f.write(
                f"pathway_burden_results available: {bool(pathway_burden_results)}\n"
            )
            f.write(
                f"pathway_burden_summary available: {bool(pathway_burden_summary)}\n"
            )
            if pathway_burden_results:
                f.write(f"Number of pathways: {len(pathway_burden_results)}\n")
                for name, result in pathway_burden_results.items():
                    burden_score = result.get("burden_score", 0)
                    f.write(f"  {name} = {burden_score}\n")
            f.write(f"About to start transformation...\n")

        # DEBUG: Log what we're working with
        logger.info(
            f"FORMATTER DEBUG: pathway_burden_results available: {bool(pathway_burden_results)}"
        )
        if pathway_burden_results:
            logger.info(
                f"FORMATTER DEBUG: Number of pathways: {len(pathway_burden_results)}"
            )
            for name, result in pathway_burden_results.items():
                burden_score = result.get("burden_score", 0)
                logger.info(f"FORMATTER DEBUG: {name} = {burden_score}")

        if pathway_burden_results:
            # Transform pathway burden results into the format expected by frontend
            disrupted_pathways = []
            cancer_pathway_associations = {}

            # Convert pathway burden results to disrupted pathways format
            for pathway_name, burden_result in pathway_burden_results.items():
                burden_score = burden_result.get("burden_score", 0)
                logger.info(
                    f"FORMATTER DEBUG: Processing {pathway_name} with burden_score {burden_score}"
                )

                if burden_score > 0.1:  # Only include pathways with significant burden
                    logger.info(
                        f"FORMATTER DEBUG: Adding {pathway_name} to disrupted_pathways"
                    )
                    # Create mutations list from damaging genes
                    mutations = []
                    if burden_result.get("damaging_genes"):
                        for gene in burden_result["damaging_genes"]:
                            mutations.append(
                                {
                                    "gene": gene,
                                    "type": "missense",  # Default type, could be enhanced
                                    "effect": f"Damaging variant in {gene}",
                                }
                            )

                    disrupted_pathways.append(
                        {
                            "name": pathway_name.replace("_", " ").title(),
                            "pathway_id": pathway_name,
                            "significance": round(
                                burden_score * 100, 1
                            ),  # Convert to percentage
                            "affected_genes": burden_result.get("damaging_genes", []),
                            "mutations": mutations,
                            "description": burden_result.get(
                                "description", f"{pathway_name} pathway"
                            ),
                            "genes_affected_ratio": f"{burden_result.get('genes_with_damaging', 0)}/{burden_result.get('genes_in_pathway', 0)}",
                        }
                    )
                else:
                    logger.info(
                        f"FORMATTER DEBUG: Skipping {pathway_name} (burden_score {burden_score} <= 0.1)"
                    )

            # Create cancer pathway associations based on high burden pathways
            high_burden_pathways = pathway_burden_summary.get(
                "high_burden_pathways", []
            )
            if high_burden_pathways:
                # Map pathways to cancer types based on common associations
                pathway_cancer_mapping = {
                    "oncogenes": ["lung", "colon", "breast"],
                    "tumor_suppressors": ["breast", "lung", "colon", "prostate"],
                    "dna_repair": ["breast", "colon"],
                    "chromatin_remodeling": ["blood", "lung"],
                    "ras_mapk": ["lung", "colon", "prostate"],
                    "cell_cycle": ["breast", "lung", "prostate"],
                    "apoptosis": ["breast", "lung", "colon"],
                    "mismatch_repair": ["colon"],
                    "wnt_signaling": ["colon"],
                    "pi3k_akt": ["breast", "prostate"],
                }

                for pathway in high_burden_pathways:
                    associated_cancers = pathway_cancer_mapping.get(pathway, [])
                    for cancer in associated_cancers:
                        if cancer not in cancer_pathway_associations:
                            cancer_pathway_associations[cancer] = []
                        cancer_pathway_associations[cancer].append(pathway)

            # Create pathway analysis structure
            pathway_analysis = {
                "disrupted_pathways": disrupted_pathways,
                "cancer_pathway_associations": cancer_pathway_associations,
                "pathway_interactions": [],  # Could be enhanced
                "clinical_recommendations": [],  # Could be enhanced
                "summary": {
                    "total_pathways_disrupted": len(disrupted_pathways),
                    "highly_disrupted_pathways": len(
                        [p for p in disrupted_pathways if p["significance"] > 50]
                    ),
                    "total_genes_affected": len(
                        set(
                            [
                                gene
                                for p in disrupted_pathways
                                for gene in p["affected_genes"]
                            ]
                        )
                    ),
                    "pathway_interaction_count": 0,
                    "overall_burden_score": pathway_burden_summary.get(
                        "overall_burden_score", 0
                    ),
                    "high_burden_pathways": high_burden_pathways,
                },
            }

            # DEBUG: Log the final transformation result
            logger.info(
                f"FORMATTER DEBUG: Created pathway_analysis with {len(disrupted_pathways)} disrupted pathways"
            )
            for pathway in disrupted_pathways:
                logger.info(
                    f"FORMATTER DEBUG: - {pathway['name']}: {pathway['significance']}%"
                )

            # DEBUG: Write transformation result to file
            with open("formatter_debug.log", "a") as f:
                f.write(
                    f"Transformation complete: {len(disrupted_pathways)} disrupted pathways\n"
                )
                for pathway in disrupted_pathways:
                    f.write(f"  - {pathway['name']}: {pathway['significance']}%\n")
                f.write(f"pathway_analysis created: {pathway_analysis is not None}\n")

        # SIMPLIFIED: Always use pathway_burden results since pathway_analyzer node is removed
        # The pathway_burden node already provides all the analysis we need
        final_pathway_analysis = (
            pathway_analysis  # This comes from pathway_burden transformation above
        )

        # DEBUG: Write final result to file
        with open("formatter_debug.log", "a") as f:
            f.write(f"final_pathway_analysis: {final_pathway_analysis is not None}\n")
            if final_pathway_analysis:
                f.write(
                    f"final disrupted_pathways: {len(final_pathway_analysis.get('disrupted_pathways', []))}\n"
                )
            f.write(f"=== END FORMATTER DEBUG ===\n")

        # COMPREHENSIVE DEBUG: Log the complete pathway_analysis object
        logger.info(f"=== COMPREHENSIVE PATHWAY_ANALYSIS DEBUG ===")
        logger.info(f"pathway_analysis is None: {pathway_analysis is None}")
        logger.info(f"final_pathway_analysis is None: {final_pathway_analysis is None}")

        if final_pathway_analysis:
            logger.info(
                f"final_pathway_analysis keys: {list(final_pathway_analysis.keys())}"
            )

            # Log disrupted_pathways in detail
            disrupted_pathways = final_pathway_analysis.get("disrupted_pathways", [])
            logger.info(f"disrupted_pathways count: {len(disrupted_pathways)}")
            for i, pathway in enumerate(disrupted_pathways):
                logger.info(f"  Pathway {i+1}:")
                logger.info(f"    name: {pathway.get('name', 'MISSING')}")
                logger.info(f"    pathway_id: {pathway.get('pathway_id', 'MISSING')}")
                logger.info(
                    f"    significance: {pathway.get('significance', 'MISSING')}"
                )
                logger.info(
                    f"    affected_genes: {pathway.get('affected_genes', 'MISSING')}"
                )
                logger.info(f"    mutations count: {len(pathway.get('mutations', []))}")
                logger.info(f"    description: {pathway.get('description', 'MISSING')}")
                logger.info(
                    f"    genes_affected_ratio: {pathway.get('genes_affected_ratio', 'MISSING')}"
                )

            # Log cancer pathway associations
            cancer_associations = final_pathway_analysis.get(
                "cancer_pathway_associations", {}
            )
            logger.info(
                f"cancer_pathway_associations count: {len(cancer_associations)}"
            )
            for cancer, pathways in cancer_associations.items():
                logger.info(f"  {cancer}: {pathways}")

            # Log summary
            summary = final_pathway_analysis.get("summary", {})
            logger.info(f"summary keys: {list(summary.keys())}")
            logger.info(
                f"  total_pathways_disrupted: {summary.get('total_pathways_disrupted', 'MISSING')}"
            )
            logger.info(
                f"  highly_disrupted_pathways: {summary.get('highly_disrupted_pathways', 'MISSING')}"
            )
            logger.info(
                f"  total_genes_affected: {summary.get('total_genes_affected', 'MISSING')}"
            )
            logger.info(
                f"  overall_burden_score: {summary.get('overall_burden_score', 'MISSING')}"
            )
            logger.info(
                f"  high_burden_pathways: {summary.get('high_burden_pathways', 'MISSING')}"
            )
        else:
            logger.info("final_pathway_analysis is None - no pathway data available")

        logger.info(f"=== END COMPREHENSIVE DEBUG ===")

        # Build the structured JSON output
        structured_json = {
            "report_metadata": {
                "pipeline_version": "1.0.0",
                "processing_time_seconds": processing_time,
                "generated_at": datetime.now().isoformat(),
                "language": state.get("language", "en"),
            },
            "patient_data": {
                "file_type": state["file_type"],
                "file_metadata": state["file_metadata"],
            },
            "summary": {
                "total_variants_found": total_variants_found,
                "variants_passed_qc": variants_passed_qc,
                "high_risk_findings": len(high_risk_findings),
                "mutation_types": mutation_type_dist,
            },
            "risk_assessment": (
                {
                    "scores": risk_scores,
                    "high_risk_findings": high_risk_findings,
                    "risk_genes": risk_genes,
                }
                if risk_scores
                else None
            ),
            "variant_details": variant_summary,
            "quality_control": state["file_metadata"].get("qc_stats", {}),
            "tcga_summary": tcga_summary,
            "cadd_summary": cadd_summary,
            "shap_validation": shap_validation,
            "metrics": state.get("metrics", {}),
            "metrics_summary": state.get("metrics_summary", {}),
            "warnings": state.get("warnings", []),
            # Add all new analysis results
            "mutation_signatures": state.get("mutational_signatures", []),
            "structural_variants": _ensure_structural_variants_format(state.get("structural_variants", [])),
            "copy_number_variants": state.get("copy_number_variants", []),
            "pathway_analysis": final_pathway_analysis,
            "pathway_summary": final_pathway_analysis.get("summary", {}) if final_pathway_analysis else {},
            "survival_analysis": state.get("survival_analysis"),
            "clinical_recommendations": state.get("clinical_recommendations"),
        }

        # Log metrics status to help debug
        metrics_available = "metrics" in state and state["metrics"]
        metrics_summary_available = (
            "metrics_summary" in state and state["metrics_summary"]
        )
        logger.info(
            f"Formatting complete - Metrics available: {metrics_available}, Summary available: {metrics_summary_available}"
        )

        logger.info("Formatting complete")

        # Return only the keys this node updates
        return {"structured_json": structured_json}

    except Exception as e:
        logger.error(f"Formatting failed: {str(e)}")
        # Return error state updates
        return {
            "pipeline_status": "failed",
            "errors": [
                {"node": "formatter", "error": str(e), "timestamp": datetime.now()}
            ],
        }
