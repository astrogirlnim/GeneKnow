"""
Gene Interaction Network Node for GeneKnow pipeline.
Builds and analyzes gene interaction networks.
Provides offline gene interaction data for desktop operation.
"""

from typing import Dict, List
import logging
from datetime import datetime
from collections import defaultdict

logger = logging.getLogger(__name__)

# Known gene interactions (curated subset for offline use)
# In production, would query STRING or BioGRID databases
GENE_INTERACTIONS = {
    # DNA repair interactions
    ("BRCA1", "BRCA2"): {"type": "protein-protein", "confidence": 0.95, "pathway": "DNA_REPAIR"},
    ("BRCA1", "TP53"): {"type": "regulation", "confidence": 0.9, "pathway": "DNA_REPAIR"},
    ("BRCA1", "PALB2"): {"type": "protein-protein", "confidence": 0.9, "pathway": "DNA_REPAIR"},
    ("BRCA2", "PALB2"): {"type": "protein-protein", "confidence": 0.95, "pathway": "DNA_REPAIR"},
    ("ATM", "TP53"): {"type": "phosphorylation", "confidence": 0.95, "pathway": "DNA_REPAIR"},
    ("ATM", "CHEK2"): {"type": "phosphorylation", "confidence": 0.9, "pathway": "DNA_REPAIR"},
    ("MLH1", "MSH2"): {"type": "protein-protein", "confidence": 0.95, "pathway": "DNA_REPAIR"},
    ("MSH2", "MSH6"): {"type": "protein-protein", "confidence": 0.95, "pathway": "DNA_REPAIR"},
    # Cell cycle interactions
    ("TP53", "MDM2"): {"type": "inhibition", "confidence": 0.95, "pathway": "CELL_CYCLE"},
    ("TP53", "CDKN1A"): {"type": "transcription", "confidence": 0.9, "pathway": "CELL_CYCLE"},
    ("CDK4", "CCND1"): {"type": "protein-protein", "confidence": 0.95, "pathway": "CELL_CYCLE"},
    ("CDK6", "CCND1"): {"type": "protein-protein", "confidence": 0.9, "pathway": "CELL_CYCLE"},
    ("RB1", "E2F1"): {"type": "inhibition", "confidence": 0.95, "pathway": "CELL_CYCLE"},
    ("CDKN2A", "CDK4"): {"type": "inhibition", "confidence": 0.9, "pathway": "CELL_CYCLE"},
    # PI3K/AKT pathway
    ("PIK3CA", "AKT1"): {"type": "activation", "confidence": 0.95, "pathway": "PI3K_AKT"},
    ("PTEN", "AKT1"): {"type": "inhibition", "confidence": 0.95, "pathway": "PI3K_AKT"},
    ("AKT1", "MTOR"): {"type": "activation", "confidence": 0.9, "pathway": "PI3K_AKT"},
    ("TSC1", "TSC2"): {"type": "protein-protein", "confidence": 0.95, "pathway": "PI3K_AKT"},
    ("TSC2", "MTOR"): {"type": "inhibition", "confidence": 0.9, "pathway": "PI3K_AKT"},
    # RAS/MAPK pathway
    ("KRAS", "BRAF"): {"type": "activation", "confidence": 0.95, "pathway": "RAS_MAPK"},
    ("BRAF", "MAP2K1"): {"type": "phosphorylation", "confidence": 0.95, "pathway": "RAS_MAPK"},
    ("MAP2K1", "MAPK1"): {"type": "phosphorylation", "confidence": 0.95, "pathway": "RAS_MAPK"},
    ("NF1", "KRAS"): {"type": "inhibition", "confidence": 0.9, "pathway": "RAS_MAPK"},
    # Apoptosis
    ("TP53", "BAX"): {"type": "transcription", "confidence": 0.9, "pathway": "APOPTOSIS"},
    ("BCL2", "BAX"): {"type": "inhibition", "confidence": 0.95, "pathway": "APOPTOSIS"},
    ("CASP9", "CASP3"): {"type": "activation", "confidence": 0.95, "pathway": "APOPTOSIS"},
    # Hormone signaling
    ("ESR1", "FOXA1"): {"type": "cooperation", "confidence": 0.9, "pathway": "HORMONE_SIGNALING"},
    ("AR", "FOXA1"): {"type": "cooperation", "confidence": 0.85, "pathway": "HORMONE_SIGNALING"},
    ("ERBB2", "ERBB3"): {"type": "heterodimerization", "confidence": 0.95, "pathway": "HORMONE_SIGNALING"},
}

# Interaction type properties
INTERACTION_PROPERTIES = {
    "protein-protein": {"strength": 0.9, "directness": "direct"},
    "phosphorylation": {"strength": 0.95, "directness": "direct"},
    "transcription": {"strength": 0.8, "directness": "indirect"},
    "activation": {"strength": 0.85, "directness": "direct"},
    "inhibition": {"strength": 0.85, "directness": "direct"},
    "regulation": {"strength": 0.7, "directness": "indirect"},
    "cooperation": {"strength": 0.7, "directness": "indirect"},
    "heterodimerization": {"strength": 0.95, "directness": "direct"},
}


def build_interaction_network(genes: List[str]) -> List[Dict]:
    """Build interaction network for given genes"""
    interactions = []
    seen = set()

    # Create a set for faster lookup
    gene_set = set(genes)

    # Find all interactions involving our genes
    for (gene1, gene2), interaction_data in GENE_INTERACTIONS.items():
        # Check if both genes are in our list
        if gene1 in gene_set and gene2 in gene_set:
            # Create unique key to avoid duplicates
            key = tuple(sorted([gene1, gene2]))

            if key not in seen:
                seen.add(key)

                interaction = {
                    "gene1": gene1,
                    "gene2": gene2,
                    "interaction_type": interaction_data["type"],
                    "confidence": interaction_data["confidence"],
                    "pathway": interaction_data.get("pathway", "unknown"),
                    "evidence": "curated_database",
                }

                # Add interaction properties
                props = INTERACTION_PROPERTIES.get(interaction_data["type"], {})
                interaction.update(props)

                interactions.append(interaction)

        # Also check for indirect connections (one gene in our set)
        elif gene1 in gene_set or gene2 in gene_set:
            # For genes with many interactions, we might want to show key partners
            # This helps understand the network context
            if interaction_data["confidence"] >= 0.9:
                our_gene = gene1 if gene1 in gene_set else gene2
                partner_gene = gene2 if gene1 in gene_set else gene1

                # Only include key cancer genes as partners
                cancer_genes = {"TP53", "KRAS", "EGFR", "BRCA1", "BRCA2", "APC", "PTEN"}
                if partner_gene in cancer_genes:
                    key = tuple(sorted([our_gene, partner_gene]))
                    if key not in seen:
                        seen.add(key)
                        interaction = {
                            "gene1": our_gene,
                            "gene2": partner_gene,
                            "interaction_type": interaction_data["type"],
                            "confidence": interaction_data["confidence"] * 0.8,  # Lower confidence for indirect
                            "pathway": interaction_data.get("pathway", "unknown"),
                            "evidence": "extended_network",
                        }
                        interactions.append(interaction)

    return interactions


def analyze_network_properties(interactions: List[Dict], genes: List[str]) -> Dict:
    """Analyze properties of the gene network"""
    if not interactions:
        return {"network_density": 0, "hub_genes": [], "isolated_genes": genes, "pathway_crosstalk": []}

    # Count connections per gene
    gene_connections = defaultdict(int)
    pathway_connections = defaultdict(set)

    for interaction in interactions:
        gene_connections[interaction["gene1"]] += 1
        gene_connections[interaction["gene2"]] += 1

        # Track pathway connections
        pathway = interaction.get("pathway", "unknown")
        pathway_connections[pathway].add(interaction["gene1"])
        pathway_connections[pathway].add(interaction["gene2"])

    # Identify hub genes (highly connected)
    avg_connections = sum(gene_connections.values()) / len(gene_connections) if gene_connections else 0
    hub_genes = [gene for gene, count in gene_connections.items() if count > avg_connections * 1.5]

    # Identify isolated genes
    connected_genes = set(gene_connections.keys())
    isolated_genes = [g for g in genes if g not in connected_genes]

    # Calculate network density
    n_genes = len(genes)
    max_possible_edges = n_genes * (n_genes - 1) / 2
    actual_edges = len(interactions)
    network_density = actual_edges / max_possible_edges if max_possible_edges > 0 else 0

    # Analyze pathway crosstalk
    pathway_crosstalk = []
    pathways = list(pathway_connections.keys())
    for i in range(len(pathways)):
        for j in range(i + 1, len(pathways)):
            shared_genes = pathway_connections[pathways[i]] & pathway_connections[pathways[j]]
            if shared_genes:
                pathway_crosstalk.append(
                    {
                        "pathways": [pathways[i], pathways[j]],
                        "shared_genes": list(shared_genes),
                        "crosstalk_strength": len(shared_genes),
                    }
                )

    return {
        "network_density": round(network_density, 3),
        "hub_genes": hub_genes,
        "isolated_genes": isolated_genes,
        "pathway_crosstalk": pathway_crosstalk,
        "total_interactions": len(interactions),
        "genes_in_network": len(connected_genes),
    }


def identify_functional_modules(interactions: List[Dict]) -> List[Dict]:
    """Identify functional modules in the network"""
    # Group interactions by pathway
    pathway_modules = defaultdict(list)

    for interaction in interactions:
        pathway = interaction.get("pathway", "unknown")
        pathway_modules[pathway].append(interaction)

    # Create module descriptions
    modules = []
    for pathway, module_interactions in pathway_modules.items():
        if len(module_interactions) >= 2:  # At least 2 interactions
            # Get unique genes in module
            module_genes = set()
            for inter in module_interactions:
                module_genes.add(inter["gene1"])
                module_genes.add(inter["gene2"])

            modules.append(
                {
                    "module_name": pathway,
                    "genes": list(module_genes),
                    "interaction_count": len(module_interactions),
                    "module_type": "pathway_based",
                    "functional_annotation": get_pathway_function(pathway),
                }
            )

    return modules


def get_pathway_function(pathway: str) -> str:
    """Get functional description for pathway"""
    pathway_functions = {
        "DNA_REPAIR": "DNA damage response and repair mechanisms",
        "CELL_CYCLE": "Cell division control and checkpoints",
        "PI3K_AKT": "Growth factor signaling and cell survival",
        "RAS_MAPK": "Proliferation signaling cascade",
        "APOPTOSIS": "Programmed cell death regulation",
        "HORMONE_SIGNALING": "Hormone receptor signaling",
        "unknown": "Uncharacterized interactions",
    }
    return pathway_functions.get(pathway, "Unknown function")


def process(state: Dict) -> Dict:
    """Generate gene interaction data"""
    logger.info("Starting gene interaction network analysis")
    state["current_node"] = "gene_interaction_network"

    try:
        # Get affected genes from variants
        variants = state.get("variant_details", state.get("filtered_variants", []))
        affected_genes = list(set(v.get("gene") for v in variants if v.get("gene")))

        if not affected_genes:
            state["gene_interactions"] = []
            logger.info("No genes to analyze for interactions")
            return state

        # Limit to reasonable number of genes for network analysis
        if len(affected_genes) > 50:
            # Prioritize by variant count
            gene_counts = defaultdict(int)
            for v in variants:
                if v.get("gene"):
                    gene_counts[v["gene"]] += 1

            # Take top 50 genes
            sorted_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)
            affected_genes = [g[0] for g in sorted_genes[:50]]
            logger.info("Limited analysis to top 50 genes by variant count")

        # Build interaction network
        interactions = build_interaction_network(affected_genes)

        # Add pathway-specific interactions if available
        pathway_analysis = state.get("pathway_analysis", {})
        for pathway in pathway_analysis.get("disrupted_pathways", []):
            pathway_genes = pathway.get("affected_genes", [])
            if len(pathway_genes) >= 2:
                # Add any missing interactions for this pathway
                pathway_interactions = build_interaction_network(pathway_genes)
                for inter in pathway_interactions:
                    # Add pathway context
                    inter["pathway_context"] = pathway["name"]
                    inter["pathway_significance"] = pathway.get("significance", 0)
                interactions.extend(pathway_interactions)

        # Remove duplicates
        seen = set()
        unique_interactions = []
        for interaction in interactions:
            key = tuple(sorted([interaction["gene1"], interaction["gene2"]]))
            if key not in seen:
                seen.add(key)
                unique_interactions.append(interaction)

        # Analyze network properties
        network_analysis = analyze_network_properties(unique_interactions, affected_genes)

        # Identify functional modules
        modules = identify_functional_modules(unique_interactions)

        # Update state
        state["gene_interactions"] = unique_interactions
        state["gene_network_analysis"] = {
            "properties": network_analysis,
            "functional_modules": modules,
            "summary": {
                "total_genes": len(affected_genes),
                "total_interactions": len(unique_interactions),
                "network_density": network_analysis["network_density"],
                "hub_gene_count": len(network_analysis["hub_genes"]),
                "module_count": len(modules),
            },
        }

        # Add to completed nodes
        completed = state.get("completed_nodes", [])
        if "gene_interaction_network" not in completed:
            completed.append("gene_interaction_network")
        state["completed_nodes"] = completed

        logger.info(f"Found {len(unique_interactions)} gene interactions")

    except Exception as e:
        logger.error(f"Error in gene interaction network analysis: {str(e)}")
        state["errors"] = state.get("errors", []) + [
            {"node": "gene_interaction_network", "error": str(e), "timestamp": datetime.now().isoformat()}
        ]
        # Set empty results on error
        state["gene_interactions"] = []

    return state
