"""
Clinical Recommendations Node for GeneKnow pipeline.
Generates comprehensive clinical recommendations based on all analyses.
"""

from typing import Dict, List
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

# Evidence-based screening guidelines
SCREENING_GUIDELINES = {
    "breast": {
        "high_risk": {
            "start_age": 25,
            "frequency": "Annual MRI + mammography",
            "additional": "Consider prophylactic options",
        },
        "moderate_risk": {
            "start_age": 30,
            "frequency": "Annual mammography + breast ultrasound",
            "additional": "Clinical breast exam every 6 months",
        },
        "standard_risk": {
            "start_age": 40,
            "frequency": "Annual mammography",
            "additional": "Clinical breast exam annually",
        },
    },
    "colon": {
        "high_risk": {
            "start_age": 20,
            "frequency": "Colonoscopy every 1-2 years",
            "additional": "Consider genetic counseling for Lynch syndrome",
        },
        "moderate_risk": {
            "start_age": 40,
            "frequency": "Colonoscopy every 3-5 years",
            "additional": "Annual FIT testing between colonoscopies",
        },
        "standard_risk": {
            "start_age": 45,
            "frequency": "Colonoscopy every 10 years",
            "additional": "Annual FIT testing",
        },
    },
    "lung": {
        "high_risk": {
            "start_age": 40,
            "frequency": "Annual low-dose CT scan",
            "additional": "Smoking cessation program",
        },
        "moderate_risk": {
            "start_age": 50,
            "frequency": "Annual low-dose CT if smoking history",
            "additional": "Pulmonary function testing",
        },
        "standard_risk": {
            "start_age": 55,
            "frequency": "Discuss screening if smoking history",
            "additional": "Annual chest X-ray not recommended",
        },
    },
    "prostate": {
        "high_risk": {
            "start_age": 40,
            "frequency": "Annual PSA + DRE",
            "additional": "Consider MRI fusion biopsy if elevated PSA",
        },
        "moderate_risk": {
            "start_age": 45,
            "frequency": "Annual PSA testing",
            "additional": "DRE at physician discretion",
        },
        "standard_risk": {
            "start_age": 50,
            "frequency": "Discuss PSA screening",
            "additional": "Shared decision making",
        },
    },
    "ovarian": {
        "high_risk": {
            "start_age": 30,
            "frequency": "CA-125 + transvaginal ultrasound every 6 months",
            "additional": "Consider risk-reducing salpingo-oophorectomy by age 40",
        },
        "moderate_risk": {
            "start_age": 35,
            "frequency": "Annual CA-125 + ultrasound",
            "additional": "Genetic counseling recommended",
        },
        "standard_risk": {
            "start_age": "No routine screening",
            "frequency": "No routine screening",
            "additional": "Be aware of symptoms",
        },
    },
}

# Targeted therapy options by gene
TARGETED_THERAPIES = {
    "EGFR": {
        "drugs": ["Erlotinib", "Gefitinib", "Osimertinib"],
        "cancer_types": ["lung"],
        "description": "EGFR tyrosine kinase inhibitors",
    },
    "ALK": {
        "drugs": ["Crizotinib", "Alectinib", "Brigatinib"],
        "cancer_types": ["lung"],
        "description": "ALK inhibitors",
    },
    "BRAF": {
        "drugs": ["Vemurafenib", "Dabrafenib", "Encorafenib"],
        "cancer_types": ["melanoma", "colon"],
        "description": "BRAF inhibitors",
    },
    "ERBB2": {
        "drugs": ["Trastuzumab", "Pertuzumab", "T-DM1"],
        "cancer_types": ["breast", "gastric"],
        "description": "HER2-targeted therapies",
    },
    "BRCA1": {
        "drugs": ["Olaparib", "Rucaparib", "Niraparib"],
        "cancer_types": ["breast", "ovarian"],
        "description": "PARP inhibitors",
    },
    "BRCA2": {
        "drugs": ["Olaparib", "Rucaparib", "Niraparib"],
        "cancer_types": ["breast", "ovarian"],
        "description": "PARP inhibitors",
    },
    "PIK3CA": {
        "drugs": ["Alpelisib"],
        "cancer_types": ["breast"],
        "description": "PI3K inhibitor",
    },
}

# Prevention strategies
PREVENTION_STRATEGIES = {
    "lifestyle": [
        "Maintain healthy weight (BMI 18.5-24.9)",
        "Regular physical activity (150 min/week moderate intensity)",
        "Limit alcohol consumption",
        "Avoid tobacco products",
        "Healthy diet rich in fruits and vegetables",
        "Limit processed meat consumption",
        "Sun protection measures",
    ],
    "chemoprevention": {
        "breast": ["Tamoxifen", "Raloxifene", "Aromatase inhibitors"],
        "colon": ["Aspirin", "NSAIDs"],
        "prostate": ["5-alpha reductase inhibitors"],
    },
    "surgical": {
        "breast": ["Prophylactic mastectomy", "Prophylactic oophorectomy"],
        "colon": ["Prophylactic colectomy for FAP/Lynch"],
        "ovarian": ["Risk-reducing salpingo-oophorectomy"],
    },
}


def determine_risk_level(risk_score: float) -> str:
    """Determine risk level from risk score"""
    if risk_score >= 50:
        return "high_risk"
    elif risk_score >= 20:
        return "moderate_risk"
    else:
        return "standard_risk"


def generate_screening_recommendations(risk_scores: Dict) -> List[Dict]:
    """Generate personalized screening recommendations"""
    recommendations = []

    for cancer_type, risk_score in risk_scores.items():
        if cancer_type in SCREENING_GUIDELINES:
            risk_level = determine_risk_level(risk_score)
            guidelines = SCREENING_GUIDELINES[cancer_type][risk_level]

            recommendation = {
                "cancer_type": cancer_type.capitalize(),
                "risk_level": risk_level.replace("_", " ").title(),
                "risk_score": round(risk_score, 1),
                "screening_recommendation": {
                    "start_age": guidelines["start_age"],
                    "frequency": guidelines["frequency"],
                    "additional_measures": guidelines["additional"],
                },
                "priority": (
                    "high"
                    if risk_score >= 50
                    else "moderate" if risk_score >= 20 else "standard"
                ),
            }

            recommendations.append(recommendation)

    # Sort by priority and risk score
    recommendations.sort(key=lambda x: (-x["risk_score"]))

    return recommendations


def generate_therapy_recommendations(variants: List[Dict]) -> List[Dict]:
    """Generate targeted therapy recommendations based on variants"""
    therapy_recommendations = []
    seen_genes = set()

    for variant in variants:
        gene = variant.get("gene")

        if gene and gene not in seen_genes and gene in TARGETED_THERAPIES:
            seen_genes.add(gene)
            therapy_info = TARGETED_THERAPIES[gene]

            recommendation = {
                "gene": gene,
                "variant": f"{variant.get('chrom', '')}:{variant.get('pos', '')}",
                "therapy_options": therapy_info["drugs"],
                "applicable_cancers": therapy_info["cancer_types"],
                "description": therapy_info["description"],
                "clinical_significance": "actionable",
            }

            therapy_recommendations.append(recommendation)

    return therapy_recommendations


def generate_prevention_recommendations(
    risk_scores: Dict, pathways: List[Dict]
) -> Dict:
    """Generate prevention recommendations based on risk profile"""
    prevention_recs = {
        "lifestyle_modifications": [],
        "chemoprevention_options": [],
        "surgical_options": [],
    }

    # General lifestyle recommendations for everyone
    prevention_recs["lifestyle_modifications"] = PREVENTION_STRATEGIES["lifestyle"]

    # Check for high-risk cancers
    high_risk_cancers = [cancer for cancer, score in risk_scores.items() if score >= 50]

    # Add specific chemoprevention options
    for cancer in high_risk_cancers:
        if cancer in PREVENTION_STRATEGIES["chemoprevention"]:
            options = PREVENTION_STRATEGIES["chemoprevention"][cancer]
            for option in options:
                prevention_recs["chemoprevention_options"].append(
                    {
                        "medication": option,
                        "indication": f"{cancer.capitalize()} cancer prevention",
                        "discuss_with": "oncologist",
                    }
                )

    # Check for DNA repair deficiency
    dna_repair_disrupted = any(
        p.get("pathway_id") == "DNA_REPAIR"
        for p in pathways
        if p.get("significance", 0) > 50
    )

    if dna_repair_disrupted:
        # Add surgical prevention options for BRCA-like phenotype
        if "breast" in high_risk_cancers:
            prevention_recs["surgical_options"].append(
                {
                    "procedure": "Prophylactic mastectomy",
                    "indication": "High breast cancer risk with DNA repair deficiency",
                    "timing": "Consider after completing childbearing",
                }
            )

        if "ovarian" in risk_scores and risk_scores["ovarian"] >= 30:
            prevention_recs["surgical_options"].append(
                {
                    "procedure": "Risk-reducing salpingo-oophorectomy",
                    "indication": "Ovarian cancer risk with DNA repair deficiency",
                    "timing": "Recommend by age 40-45 or after childbearing",
                }
            )

    return prevention_recs


def generate_clinical_trials_recommendations(
    variants: List[Dict], pathways: List[Dict]
) -> List[str]:
    """Generate clinical trial recommendations based on genomic profile"""
    trials_recs = []

    # Check for specific targetable alterations
    targetable_genes = set(
        v.get("gene") for v in variants if v.get("gene") in TARGETED_THERAPIES
    )

    if targetable_genes:
        trials_recs.append(
            f"Consider clinical trials targeting: {', '.join(targetable_genes)}"
        )

    # Check for immunotherapy eligibility
    immune_pathway_disrupted = any(
        p.get("pathway_id") == "IMMUNE_CHECKPOINT" for p in pathways
    )

    if immune_pathway_disrupted:
        trials_recs.append("Consider immunotherapy clinical trials")

    # DNA repair deficiency trials
    dna_repair_disrupted = any(
        p.get("pathway_id") == "DNA_REPAIR"
        for p in pathways
        if p.get("significance", 0) > 30
    )

    if dna_repair_disrupted:
        trials_recs.append("Eligible for PARP inhibitor trials")
        trials_recs.append("Consider DNA damage response (DDR) targeted trials")

    return trials_recs


def generate_monitoring_plan(risk_scores: Dict, variants: List[Dict]) -> Dict:
    """Generate a monitoring plan based on risk profile"""
    monitoring_plan = {"biomarkers": [], "imaging": [], "clinical_assessments": []}

    # Add biomarker monitoring
    for cancer_type, risk_score in risk_scores.items():
        if risk_score >= 20:
            if cancer_type == "breast":
                monitoring_plan["biomarkers"].append(
                    {
                        "marker": "CA 15-3",
                        "frequency": "Every 6 months if high risk",
                        "indication": "Breast cancer monitoring",
                    }
                )
            elif cancer_type == "colon":
                monitoring_plan["biomarkers"].append(
                    {
                        "marker": "CEA",
                        "frequency": "Annual",
                        "indication": "Colorectal cancer monitoring",
                    }
                )
            elif cancer_type == "ovarian":
                monitoring_plan["biomarkers"].append(
                    {
                        "marker": "CA-125",
                        "frequency": "Every 6 months",
                        "indication": "Ovarian cancer monitoring",
                    }
                )
            elif cancer_type == "prostate":
                monitoring_plan["biomarkers"].append(
                    {
                        "marker": "PSA",
                        "frequency": "Annual",
                        "indication": "Prostate cancer monitoring",
                    }
                )

    # Add imaging recommendations
    high_risk_count = sum(1 for score in risk_scores.values() if score >= 50)

    if high_risk_count >= 2:
        monitoring_plan["imaging"].append(
            {
                "modality": "Whole body MRI",
                "frequency": "Consider annual",
                "indication": "Multiple high cancer risks",
            }
        )

    # Clinical assessments
    monitoring_plan["clinical_assessments"] = [
        {
            "assessment": "Genetic counseling",
            "frequency": "Initial consultation recommended",
            "indication": "Hereditary cancer risk assessment",
        },
        {
            "assessment": "Annual cancer risk assessment",
            "frequency": "Yearly",
            "indication": "Update risk profile and screening plan",
        },
    ]

    return monitoring_plan


def process(state: Dict) -> Dict:
    """Generate comprehensive clinical recommendations"""
    logger.info("Generating clinical recommendations")
    state["current_node"] = "clinical_recommendations"

    try:
        # Get all relevant data
        risk_scores = state.get("risk_scores", {})
        variants = state.get("variant_details", state.get("filtered_variants", []))

        # Get pathway burden results directly (not pathway_analysis which hasn't been created yet)
        pathway_burden_results = state.get("pathway_burden_results", {})
        # pathway_burden_summary = state.get("pathway_burden_summary", {})

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
                            "pathway_id": pathway_name.upper(),  # Convert to uppercase for matching
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

        # survival_analysis = state.get("survival_analysis", {})

        # Generate screening recommendations
        screening_recs = generate_screening_recommendations(risk_scores)

        # Generate therapy recommendations
        therapy_recs = generate_therapy_recommendations(variants)

        # Generate prevention recommendations
        prevention_recs = generate_prevention_recommendations(
            risk_scores, disrupted_pathways
        )

        # Generate clinical trial recommendations
        trials_recs = generate_clinical_trials_recommendations(
            variants, disrupted_pathways
        )

        # Generate monitoring plan
        monitoring_plan = generate_monitoring_plan(risk_scores, variants)

        # Create comprehensive clinical recommendations
        clinical_recommendations = {
            "screening_recommendations": screening_recs,
            "targeted_therapy_options": therapy_recs,
            "prevention_strategies": prevention_recs,
            "clinical_trial_eligibility": trials_recs,
            "monitoring_plan": monitoring_plan,
            "summary": {
                "high_risk_cancers": [c for c, s in risk_scores.items() if s >= 50],
                "actionable_variants": len(therapy_recs),
                "priority_actions": generate_priority_actions(
                    screening_recs, therapy_recs, prevention_recs
                ),
                "follow_up_timeline": generate_follow_up_timeline(screening_recs),
            },
        }

        # Update state
        state["clinical_recommendations"] = clinical_recommendations

        # Add to structured JSON for report
        if "structured_json" not in state:
            state["structured_json"] = {}

        state["structured_json"]["clinical_recommendations"] = clinical_recommendations

        # Add to completed nodes
        completed = state.get("completed_nodes", [])
        if "clinical_recommendations" not in completed:
            completed.append("clinical_recommendations")
        state["completed_nodes"] = completed

        logger.info(f"Generated {len(screening_recs)} screening recommendations")
        logger.info(f"Found {len(therapy_recs)} targeted therapy options")

    except Exception as e:
        logger.error(f"Error generating clinical recommendations: {str(e)}")
        state["errors"] = state.get("errors", []) + [
            {
                "node": "clinical_recommendations",
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            }
        ]
        # Set empty results on error
        state["clinical_recommendations"] = {}

    return state


def generate_priority_actions(
    screening_recs: List[Dict], therapy_recs: List[Dict], prevention_recs: Dict
) -> List[str]:
    """Generate list of priority actions"""
    actions = []

    # High priority screening
    high_priority_screening = [r for r in screening_recs if r["priority"] == "high"]
    if high_priority_screening:
        for rec in high_priority_screening[:2]:  # Top 2
            actions.append(
                f"Begin {rec['cancer_type']} cancer screening: {rec['screening_recommendation']['frequency']}"
            )

    # Actionable variants
    if therapy_recs:
        actions.append(
            f"Discuss targeted therapy options with oncologist ({len(therapy_recs)} actionable variants found)"
        )

    # Prevention
    if prevention_recs.get("surgical_options"):
        actions.append("Genetic counseling for risk-reducing surgery options")

    if prevention_recs.get("chemoprevention_options"):
        actions.append("Discuss chemoprevention options with physician")

    return actions[:5]  # Top 5 actions


def generate_follow_up_timeline(screening_recs: List[Dict]) -> Dict:
    """Generate follow-up timeline"""
    timeline = {
        "immediate": [],  # Within 1 month
        "short_term": [],  # 1-3 months
        "medium_term": [],  # 3-12 months
        "long_term": [],  # Annual
    }

    for rec in screening_recs:
        if rec["priority"] == "high":
            timeline["immediate"].append(
                f"Schedule {rec['cancer_type']} screening consultation"
            )
        elif rec["priority"] == "moderate":
            timeline["short_term"].append(
                f"Discuss {rec['cancer_type']} screening with physician"
            )

    timeline["immediate"].append("Genetic counseling consultation")
    timeline["medium_term"].append("Establish care team and monitoring plan")
    timeline["long_term"].append("Annual risk reassessment")

    return timeline
