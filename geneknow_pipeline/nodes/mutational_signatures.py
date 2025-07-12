"""
Mutational Signatures Node for GeneKnow pipeline.
Analyzes mutational signatures using COSMIC signature patterns.
Currently provides simulated analysis for offline operation.
"""

from typing import Dict, List
import logging
from datetime import datetime
from collections import defaultdict

logger = logging.getLogger(__name__)


def load_cosmic_signatures() -> Dict:
    """Load COSMIC mutational signatures (simulated for offline use)"""
    # In production, load from COSMIC database
    # For now, return common cancer-associated signatures
    return {
        "SBS1": {
            "name": "Spontaneous deamination",
            "description": "C>T transitions at CpG sites",
            "etiology": "Aging-related mutations",
            "cancer_types": ["All cancers"],
            "pattern_weight": 0.3,
        },
        "SBS3": {
            "name": "BRCA1/BRCA2 deficiency",
            "description": "Homologous recombination deficiency",
            "etiology": "Inherited DNA repair defects",
            "cancer_types": ["Breast", "Ovarian", "Pancreatic"],
            "pattern_weight": 0.0,
        },
        "SBS4": {
            "name": "Tobacco smoking",
            "description": "C>A mutations from tobacco carcinogens",
            "etiology": "Tobacco smoke exposure",
            "cancer_types": ["Lung", "Head and neck"],
            "pattern_weight": 0.0,
        },
        "SBS7a": {
            "name": "UV exposure",
            "description": "C>T mutations at dipyrimidines",
            "etiology": "Ultraviolet light exposure",
            "cancer_types": ["Melanoma", "Skin"],
            "pattern_weight": 0.0,
        },
        "SBS13": {
            "name": "APOBEC activity",
            "description": "C>G and C>T mutations",
            "etiology": "APOBEC enzyme activity",
            "cancer_types": ["Breast", "Bladder", "Cervical"],
            "pattern_weight": 0.0,
        },
    }


def extract_mutation_context(variant: Dict) -> str:
    """Extract trinucleotide context for signature analysis"""
    # In production, would fetch from reference genome
    # For now, create a simple context based on mutation type
    ref = variant.get("re", "")
    alt = variant.get("alt", "")

    if len(ref) == 1 and len(alt) == 1:
        # SNV - create trinucleotide context
        return f"N[{ref}>{alt}]N"
    else:
        return "complex"


def analyze_mutation_patterns(variants: List[Dict]) -> Dict[str, float]:
    """Analyze mutation patterns in variants"""
    mutation_contexts = defaultdict(int)

    for variant in variants:
        if variant.get("mutation_type") == "snv":
            context = extract_mutation_context(variant)
            mutation_contexts[context] += 1

    # Calculate pattern weights based on mutation types
    pattern_weights = {}
    total_snvs = sum(mutation_contexts.values())

    if total_snvs > 0:
        # C>T transitions (aging, UV)
        c_to_t = sum(v for k, v in mutation_contexts.items() if "C>T" in k)
        pattern_weights["C>T"] = c_to_t / total_snvs

        # C>A transversions (smoking)
        c_to_a = sum(v for k, v in mutation_contexts.items() if "C>A" in k)
        pattern_weights["C>A"] = c_to_a / total_snvs

        # C>G mutations (APOBEC)
        c_to_g = sum(v for k, v in mutation_contexts.items() if "C>G" in k)
        pattern_weights["C>G"] = c_to_g / total_snvs

    return pattern_weights


def calculate_signature_contributions(
    pattern_weights: Dict[str, float], cosmic_signatures: Dict
) -> List[Dict]:
    """Calculate contribution of each signature based on mutation patterns"""
    results = []

    # Simple heuristic-based assignment
    # In production, would use NMF decomposition

    # Check for aging signature (C>T transitions)
    if pattern_weights.get("C>T", 0) > 0.3:
        results.append(
            {
                "signature": "SBS1",
                "name": cosmic_signatures["SBS1"]["name"],
                "contribution": min(pattern_weights.get("C>T", 0), 0.5),
                "description": cosmic_signatures["SBS1"]["description"],
                "etiology": cosmic_signatures["SBS1"]["etiology"],
            }
        )

    # Check for smoking signature (C>A transversions)
    if pattern_weights.get("C>A", 0) > 0.2:
        results.append(
            {
                "signature": "SBS4",
                "name": cosmic_signatures["SBS4"]["name"],
                "contribution": min(pattern_weights.get("C>A", 0) * 1.5, 0.4),
                "description": cosmic_signatures["SBS4"]["description"],
                "etiology": cosmic_signatures["SBS4"]["etiology"],
            }
        )

    # Check for APOBEC signature (C>G mutations)
    if pattern_weights.get("C>G", 0) > 0.15:
        results.append(
            {
                "signature": "SBS13",
                "name": cosmic_signatures["SBS13"]["name"],
                "contribution": min(pattern_weights.get("C>G", 0) * 2, 0.3),
                "description": cosmic_signatures["SBS13"]["description"],
                "etiology": cosmic_signatures["SBS13"]["etiology"],
            }
        )

    # Check for BRCA deficiency (based on variant genes)
    # This would be more sophisticated in production
    if len(results) == 0:
        # Default aging signature
        results.append(
            {
                "signature": "SBS1",
                "name": cosmic_signatures["SBS1"]["name"],
                "contribution": 0.2,
                "description": cosmic_signatures["SBS1"]["description"],
                "etiology": cosmic_signatures["SBS1"]["etiology"],
            }
        )

    # Normalize contributions
    total_contribution = sum(r["contribution"] for r in results)
    if total_contribution > 0:
        for r in results:
            r["contribution"] = round(r["contribution"] / total_contribution, 3)

    # Sort by contribution
    results.sort(key=lambda x: x["contribution"], reverse=True)

    return results


def process(state: Dict) -> Dict:
    """Analyze mutational signatures in variants"""
    logger.info("Starting mutational signature analysis")
    # Note: Don't set current_node to avoid concurrent updates

    try:
        # Use classified variants if available, otherwise use filtered variants
        variants = state.get("classified_variants", state.get("filtered_variants", []))

        # Load COSMIC signatures
        cosmic_signatures = load_cosmic_signatures()

        # Analyze mutation patterns
        pattern_weights = analyze_mutation_patterns(variants)

        # Calculate signature contributions
        signatures = calculate_signature_contributions(
            pattern_weights, cosmic_signatures
        )

        # Add special check for BRCA1/BRCA2 variants
        brca_genes = {"BRCA1", "BRCA2"}
        variant_genes = {v.get("gene") for v in variants if v.get("gene")}

        if variant_genes.intersection(brca_genes):
            # Add BRCA deficiency signature if not already present
            if not any(s["signature"] == "SBS3" for s in signatures):
                signatures.append(
                    {
                        "signature": "SBS3",
                        "name": cosmic_signatures["SBS3"]["name"],
                        "contribution": 0.25,
                        "description": cosmic_signatures["SBS3"]["description"],
                        "etiology": cosmic_signatures["SBS3"]["etiology"],
                    }
                )
                # Re-normalize
                total = sum(s["contribution"] for s in signatures)
                for s in signatures:
                    s["contribution"] = round(s["contribution"] / total, 3)

        logger.info(f"Identified {len(signatures)} mutational signatures")

        # Return only the fields this node updates
        return {
            "mutational_signatures": signatures,
            "mutational_signatures_summary": {
                "total_signatures": len(signatures),
                "dominant_signature": signatures[0] if signatures else None,
                "pattern_weights": pattern_weights,
            }
        }

    except Exception as e:
        logger.error(f"Error in mutational signature analysis: {str(e)}")
        # Return error state updates
        return {
            "mutational_signatures": [],
            "errors": [{
                "node": "mutational_signatures",
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            }]
        }
