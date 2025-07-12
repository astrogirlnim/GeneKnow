"""
LangGraph nodes for genomic pipeline.
Each node represents a discrete processing step.
"""

# Import all nodes for easy access from graph.py
try:
    # Package imports
    from . import file_input
    from . import preprocess
    from . import variant_calling
    from . import qc_filter
    from . import population_mapper
    from . import tcga_mapper
    from . import cadd_scoring
    from . import clinvar_annotator
    from . import feature_vector_builder
    from . import prs_calculator
    from . import pathway_burden
    from . import ml_fusion_node
    from . import risk_model
    from . import shap_validator
    from . import metrics_calculator
    from . import formatter
    from . import report_writer
    from . import maf_parser

    # New nodes for clinical view data
    from . import mutation_classifier
    from . import mutational_signatures
    from . import variant_transformer
    from . import structural_variant_detector
    from . import cnv_detector
    from . import pathway_analyzer
    from . import gene_interaction_network
    from . import survival_analyzer
    from . import clinical_recommendations
except ImportError:
    # Direct imports
    import file_input
    import preprocess
    import variant_calling
    import qc_filter
    import population_mapper
    import tcga_mapper
    import cadd_scoring
    import clinvar_annotator
    import feature_vector_builder
    import prs_calculator
    import pathway_burden
    import ml_fusion_node
    import risk_model
    import shap_validator
    import metrics_calculator
    import formatter
    import report_writer
    import maf_parser

    # New nodes for clinical view data
    import mutation_classifier
    import mutational_signatures
    import variant_transformer
    import structural_variant_detector
    import cnv_detector
    import pathway_analyzer
    import gene_interaction_network
    import survival_analyzer
    import clinical_recommendations

__all__ = [
    "file_input",
    "preprocess",
    "variant_calling",
    "qc_filter",
    "population_mapper",
    "tcga_mapper",
    "cadd_scoring",
    "clinvar_annotator",
    "feature_vector_builder",
    "prs_calculator",
    "pathway_burden",
    "ml_fusion_node",
    "risk_model",
    "shap_validator",
    "metrics_calculator",
    "formatter",
    "report_writer",
    "maf_parser",
    # New nodes
    "mutation_classifier",
    "mutational_signatures",
    "variant_transformer",
    "structural_variant_detector",
    "cnv_detector",
    "pathway_analyzer",
    "gene_interaction_network",
    "survival_analyzer",
    "clinical_recommendations",
]
