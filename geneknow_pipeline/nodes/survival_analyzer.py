"""
Survival Analyzer Node for GeneKnow pipeline.
Performs survival analysis based on genomic features.
Provides simulated survival curves for offline operation.
"""

import numpy as np
from typing import Dict, List, Tuple
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

# Survival impact of known mutations (hazard ratios)
# HR > 1 = worse survival, HR < 1 = better survival
MUTATION_HAZARD_RATIOS = {
    "TP53": {"HR": 2.1, "confidence": 0.95, "cancer_types": ["breast", "lung", "colon"]},
    "KRAS": {"HR": 1.8, "confidence": 0.9, "cancer_types": ["lung", "colon", "pancreatic"]},
    "BRAF": {"HR": 1.7, "confidence": 0.85, "cancer_types": ["melanoma", "colon"]},
    "EGFR": {"HR": 0.6, "confidence": 0.9, "cancer_types": ["lung"]},  # Better with targeted therapy
    "BRCA1": {"HR": 1.4, "confidence": 0.85, "cancer_types": ["breast", "ovarian"]},
    "BRCA2": {"HR": 1.3, "confidence": 0.85, "cancer_types": ["breast", "ovarian"]},
    "PIK3CA": {"HR": 1.5, "confidence": 0.8, "cancer_types": ["breast", "colon"]},
    "PTEN": {"HR": 1.6, "confidence": 0.85, "cancer_types": ["prostate", "breast"]},
    "MYC": {"HR": 2.0, "confidence": 0.9, "cancer_types": ["blood", "breast"]},
    "BCL2": {"HR": 0.8, "confidence": 0.75, "cancer_types": ["blood"]},  # Better in some blood cancers
    "APC": {"HR": 1.4, "confidence": 0.85, "cancer_types": ["colon"]},
    "ERBB2": {"HR": 0.7, "confidence": 0.9, "cancer_types": ["breast"]},  # Better with HER2 therapy
    "ALK": {"HR": 0.5, "confidence": 0.85, "cancer_types": ["lung"]},  # Better with ALK inhibitors
    "ROS1": {"HR": 0.6, "confidence": 0.8, "cancer_types": ["lung"]},  # Better with targeted therapy
    "MET": {"HR": 1.7, "confidence": 0.8, "cancer_types": ["lung", "gastric"]},
}

# Pathway-based survival impact
PATHWAY_SURVIVAL_IMPACT = {
    "DNA_REPAIR": {"HR": 1.5, "description": "DNA repair deficiency"},
    "CELL_CYCLE": {"HR": 1.8, "description": "Cell cycle dysregulation"},
    "PI3K_AKT": {"HR": 1.6, "description": "PI3K/AKT activation"},
    "RAS_MAPK": {"HR": 1.7, "description": "RAS/MAPK activation"},
    "APOPTOSIS": {"HR": 1.4, "description": "Apoptosis evasion"},
    "IMMUNE_CHECKPOINT": {"HR": 0.7, "description": "Immunotherapy response"},
    "HORMONE_SIGNALING": {"HR": 0.8, "description": "Hormone therapy response"},
}

# Base survival rates by cancer type (5-year survival)
BASE_SURVIVAL_RATES = {
    "breast": 0.89,
    "colon": 0.64,
    "lung": 0.18,
    "prostate": 0.98,
    "blood": 0.60,
    "ovarian": 0.47,
    "pancreatic": 0.09,
    "melanoma": 0.92,
}


def calculate_combined_hazard_ratio(mutations: List[str], pathways: List[str]) -> float:
    """Calculate combined hazard ratio from multiple factors"""
    combined_hr = 1.0

    # Apply mutation-based HRs
    for mutation in mutations:
        if mutation in MUTATION_HAZARD_RATIOS:
            hr_data = MUTATION_HAZARD_RATIOS[mutation]
            # Weight by confidence
            weighted_hr = 1 + (hr_data["HR"] - 1) * hr_data["confidence"]
            combined_hr *= weighted_hr

    # Apply pathway-based HRs
    for pathway in pathways:
        if pathway in PATHWAY_SURVIVAL_IMPACT:
            pathway_hr = PATHWAY_SURVIVAL_IMPACT[pathway]["HR"]
            # Pathways have multiplicative effect but dampened
            combined_hr *= 1 + (pathway_hr - 1) * 0.5

    # Cap the combined HR to reasonable range
    return min(max(combined_hr, 0.1), 10.0)


def generate_survival_curve(base_survival: float, hazard_ratio: float, time_points: np.ndarray) -> np.ndarray:
    """Generate survival curve using exponential model"""
    # Convert 5-year survival to annual hazard rate
    if base_survival > 0:
        base_hazard = -np.log(base_survival) / 5.0
    else:
        base_hazard = 1.0

    # Apply hazard ratio
    adjusted_hazard = base_hazard * hazard_ratio

    # Generate survival probabilities
    survival_probs = np.exp(-adjusted_hazard * time_points)

    # Add some realistic variation
    noise = np.random.normal(0, 0.02, len(time_points))
    survival_probs = np.maximum(0, np.minimum(1, survival_probs + noise))

    # Ensure monotonic decrease
    for i in range(1, len(survival_probs)):
        survival_probs[i] = min(survival_probs[i], survival_probs[i - 1])

    return survival_probs


def calculate_median_survival(survival_curve: np.ndarray, time_points: np.ndarray) -> float:
    """Calculate median survival time from survival curve"""
    # Find where survival drops below 50%
    below_50 = np.where(survival_curve < 0.5)[0]

    if len(below_50) == 0:
        return time_points[-1]  # Median not reached

    idx = below_50[0]
    if idx == 0:
        return time_points[0]

    # Linear interpolation
    t1, t2 = time_points[idx - 1], time_points[idx]
    s1, s2 = survival_curve[idx - 1], survival_curve[idx]

    if s1 != s2:
        median_time = t1 + (0.5 - s1) * (t2 - t1) / (s2 - s1)
    else:
        median_time = (t1 + t2) / 2

    return median_time


def generate_confidence_bands(survival_curve: np.ndarray, sample_size: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """Generate confidence bands for survival curve"""
    # Standard error using Greenwood's formula (simplified)
    se = np.sqrt(survival_curve * (1 - survival_curve) / sample_size)

    # 95% confidence interval
    z_score = 1.96
    lower = np.maximum(0, survival_curve - z_score * se)
    upper = np.minimum(1, survival_curve + z_score * se)

    return lower, upper


def analyze_prognostic_factors(variants: List[Dict], pathways: List[Dict]) -> Dict:
    """Analyze prognostic factors from variants and pathways"""
    positive_factors = []
    negative_factors = []

    # Check variants
    seen_genes = set()
    for variant in variants:
        gene = variant.get("gene")
        if gene and gene not in seen_genes:
            seen_genes.add(gene)

            if gene in MUTATION_HAZARD_RATIOS:
                hr_data = MUTATION_HAZARD_RATIOS[gene]
                factor = {
                    "gene": gene,
                    "hazard_ratio": hr_data["HR"],
                    "impact": "positive" if hr_data["HR"] < 1 else "negative",
                    "confidence": hr_data["confidence"],
                }

                if hr_data["HR"] < 1:
                    positive_factors.append(factor)
                else:
                    negative_factors.append(factor)

    # Check pathways
    for pathway in pathways:
        pathway_id = pathway.get("pathway_id")
        if pathway_id in PATHWAY_SURVIVAL_IMPACT:
            impact_data = PATHWAY_SURVIVAL_IMPACT[pathway_id]
            factor = {
                "pathway": pathway.get("name", pathway_id),
                "hazard_ratio": impact_data["HR"],
                "impact": "positive" if impact_data["HR"] < 1 else "negative",
                "description": impact_data["description"],
            }

            if impact_data["HR"] < 1:
                positive_factors.append(factor)
            else:
                negative_factors.append(factor)

    return {
        "positive_prognostic_factors": positive_factors,
        "negative_prognostic_factors": negative_factors,
        "factor_summary": {
            "positive_count": len(positive_factors),
            "negative_count": len(negative_factors),
            "overall_prognosis": "favorable" if len(positive_factors) > len(negative_factors) else "unfavorable",
        },
    }


def process(state: Dict) -> Dict:
    """Perform survival analysis"""
    logger.info("Starting survival analysis")
    state["current_node"] = "survival_analyzer"

    try:
        # Get variants and pathways
        variants = state.get("variant_details", state.get("filtered_variants", []))
        
        # Get pathway burden results directly (not pathway_analysis which hasn't been created yet)
        pathway_burden_results = state.get("pathway_burden_results", {})
        pathway_burden_summary = state.get("pathway_burden_summary", {})
        
        # Transform pathway burden results into disrupted pathways format
        disrupted_pathways = []
        if pathway_burden_results:
            for pathway_name, burden_result in pathway_burden_results.items():
                burden_score = burden_result.get("burden_score", 0)
                
                if burden_score > 0.1:  # Only include pathways with significant burden
                    # Create mutations list from damaging genes
                    mutations = []
                    if burden_result.get("damaging_genes"):
                        for gene in burden_result["damaging_genes"]:
                            mutations.append({
                                "gene": gene,
                                "type": "missense",  # Default type, could be enhanced
                                "effect": f"Damaging variant in {gene}"
                            })
                    
                    disrupted_pathways.append({
                        "name": pathway_name.replace("_", " ").title(),
                        "pathway_id": pathway_name.upper(),  # Convert to uppercase for matching
                        "significance": round(burden_score * 100, 1),  # Convert to percentage
                        "affected_genes": burden_result.get("damaging_genes", []),
                        "mutations": mutations,
                        "description": burden_result.get("description", f"{pathway_name} pathway"),
                        "genes_affected_ratio": f"{burden_result.get('genes_with_damaging', 0)}/{burden_result.get('genes_in_pathway', 0)}"
                    })
        
        logger.info(f"Transformed {len(disrupted_pathways)} disrupted pathways from pathway burden results")

        # Get risk scores to determine relevant cancer types
        risk_scores = state.get("risk_scores", {})

        # Find top cancer types by risk
        top_cancers = sorted(risk_scores.items(), key=lambda x: x[1], reverse=True)[:3]

        # Extract relevant mutations
        mutations = list(set(v.get("gene") for v in variants if v.get("gene") in MUTATION_HAZARD_RATIOS))

        # Extract disrupted pathway IDs
        pathway_ids = [p.get("pathway_id") for p in disrupted_pathways if p.get("significance", 0) > 30]

        # Generate survival curves for each cancer type
        survival_curves = {}
        time_points = np.linspace(0, 10, 50)

        for cancer_type, risk_score in top_cancers:
            if cancer_type in BASE_SURVIVAL_RATES:
                # Calculate combined hazard ratio
                combined_hr = calculate_combined_hazard_ratio(mutations, pathway_ids)

                # Generate curves
                base_survival = BASE_SURVIVAL_RATES[cancer_type]

                # Population curve (no mutations)
                population_curve = generate_survival_curve(base_survival, 1.0, time_points)

                # Patient curve (with mutations)
                patient_curve = generate_survival_curve(base_survival, combined_hr, time_points)

                # Generate confidence bands
                lower_ci, upper_ci = generate_confidence_bands(patient_curve)

                # Calculate median survival
                population_median = calculate_median_survival(population_curve, time_points)
                patient_median = calculate_median_survival(patient_curve, time_points)

                survival_curves[cancer_type] = {
                    "time_points": time_points.tolist(),
                    "population_survival": population_curve.tolist(),
                    "patient_survival": patient_curve.tolist(),
                    "confidence_interval": {"lower": lower_ci.tolist(), "upper": upper_ci.tolist()},
                    "median_survival": {
                        "population": round(population_median, 1),
                        "patient": round(patient_median, 1),
                        "difference": round(population_median - patient_median, 1),
                    },
                    "hazard_ratio": round(combined_hr, 2),
                    "five_year_survival": {
                        "population": round(population_curve[25] * 100, 1),  # 5-year point
                        "patient": round(patient_curve[25] * 100, 1),
                    },
                }

        # Analyze prognostic factors
        prognostic_analysis = analyze_prognostic_factors(variants, disrupted_pathways)

        # Create survival analysis result
        survival_analysis = {
            "survival_curves": survival_curves,
            "prognostic_factors": prognostic_analysis,
            "clinical_interpretation": {
                "mutations_analyzed": len(mutations),
                "pathways_analyzed": len(pathway_ids),
                "recommendation": generate_survival_recommendation(survival_curves, prognostic_analysis),
            },
            "methodology_note": "Survival estimates based on genomic features and published hazard ratios",
        }
        
        # Add frontend-compatible format
        # Convert the first (highest risk) cancer type's survival curve to the expected format
        if survival_curves and len(top_cancers) > 0:
            highest_risk_cancer = top_cancers[0][0]
            if highest_risk_cancer in survival_curves:
                curve_data = survival_curves[highest_risk_cancer]
                
                # Create patient profile with estimated survival
                patient_profile = {
                    "estimated_survival": [],
                    "risk_category": "High Risk" if top_cancers[0][1] > 50 else "Moderate Risk"
                }
                
                # Create population average array
                population_average = []
                
                # Convert time points to age-based survival data
                base_age = 40  # Starting age for visualization
                for i, time_point in enumerate(curve_data["time_points"]):
                    age = base_age + int(time_point)
                    
                    # Patient survival data point
                    patient_profile["estimated_survival"].append({
                        "age": age,
                        "probability": curve_data["patient_survival"][i]
                    })
                    
                    # Population average data point
                    population_average.append({
                        "age": age,
                        "probability": curve_data["population_survival"][i]
                    })
                
                # Add frontend-compatible fields
                survival_analysis["patient_profile"] = patient_profile
                survival_analysis["population_average"] = population_average

        # Update state
        state["survival_analysis"] = survival_analysis

        # Add to completed nodes
        completed = state.get("completed_nodes", [])
        if "survival_analyzer" not in completed:
            completed.append("survival_analyzer")
        state["completed_nodes"] = completed

        logger.info(f"Generated survival curves for {len(survival_curves)} cancer types")

    except Exception as e:
        logger.error(f"Error in survival analysis: {str(e)}")
        state["errors"] = state.get("errors", []) + [
            {"node": "survival_analyzer", "error": str(e), "timestamp": datetime.now().isoformat()}
        ]
        # Set empty results on error
        state["survival_analysis"] = {}

    return state


def generate_survival_recommendation(survival_curves: Dict, prognostic_analysis: Dict) -> str:
    """Generate clinical recommendation based on survival analysis"""
    recommendations = []

    # Check if survival is significantly worse
    for cancer_type, curve_data in survival_curves.items():
        hr = curve_data.get("hazard_ratio", 1.0)
        if hr > 1.5:
            recommendations.append(
                f"Enhanced {cancer_type} cancer surveillance recommended due to " f"{hr:.1f}x increased risk"
            )
        elif hr < 0.7:
            recommendations.append(f"Favorable {cancer_type} cancer prognosis if targeted therapy available")

    # Check prognostic factors
    positive_count = prognostic_analysis["factor_summary"]["positive_count"]
    negative_count = prognostic_analysis["factor_summary"]["negative_count"]

    if negative_count > positive_count:
        recommendations.append(
            "Consider aggressive screening and prevention strategies due to "
            f"{negative_count} negative prognostic factors"
        )
    elif positive_count > 0:
        recommendations.append(
            f"Presence of {positive_count} positive prognostic factors may indicate "
            "better response to targeted therapies"
        )

    return "; ".join(recommendations) if recommendations else "Standard surveillance recommended"
