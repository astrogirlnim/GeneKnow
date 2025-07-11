"""
Pathway Analyzer Node for GeneKnow pipeline.
Analyzes disruption of biological pathways by variants.
Extends the existing pathway burden analysis with detailed pathway impact scores.
"""
from typing import Dict, List, Set
import logging
from datetime import datetime
from collections import defaultdict

logger = logging.getLogger(__name__)

# Comprehensive pathway definitions
PATHWAYS = {
    "DNA_REPAIR": {
        "name": "DNA Repair",
        "genes": ["BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "RAD51", 
                  "PALB2", "RAD50", "RAD52", "XRCC1", "XRCC2", "XRCC3", "MLH1", 
                  "MSH2", "MSH6", "PMS2", "MUTYH", "POLE", "POLD1"],
        "description": "Maintains genomic stability",
        "cancer_relevance": 0.95
    },
    "CELL_CYCLE": {
        "name": "Cell Cycle Control",
        "genes": ["TP53", "RB1", "CDKN2A", "CCND1", "CDK4", "CDK6", "CDKN1A", 
                  "CDKN1B", "E2F1", "MDM2", "MDM4", "CCNE1", "CDC25A"],
        "description": "Regulates cell division",
        "cancer_relevance": 0.9
    },
    "APOPTOSIS": {
        "name": "Apoptosis",
        "genes": ["TP53", "BCL2", "BAX", "CASP3", "CASP9", "FAS", "FADD", 
                  "APAF1", "BCL2L1", "MCL1", "BID", "BIRC5"],
        "description": "Programmed cell death",
        "cancer_relevance": 0.85
    },
    "WNT_SIGNALING": {
        "name": "WNT Signaling",
        "genes": ["APC", "CTNNB1", "WNT1", "TCF7L2", "AXIN1", "AXIN2", "GSK3B", 
                  "DVL1", "DVL2", "DVL3", "LRP5", "LRP6"],
        "description": "Cell fate determination",
        "cancer_relevance": 0.8
    },
    "PI3K_AKT": {
        "name": "PI3K-AKT Signaling",
        "genes": ["PIK3CA", "PTEN", "AKT1", "AKT2", "MTOR", "TSC1", "TSC2", 
                  "PIK3CB", "PIK3R1", "RICTOR", "RAPTOR", "STK11"],
        "description": "Cell growth and survival",
        "cancer_relevance": 0.9
    },
    "RAS_MAPK": {
        "name": "RAS-MAPK Signaling",
        "genes": ["KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2", 
                  "MAPK1", "MAPK3", "NF1", "RASA1", "SOS1"],
        "description": "Cell proliferation signaling",
        "cancer_relevance": 0.95
    },
    "NOTCH_SIGNALING": {
        "name": "NOTCH Signaling",
        "genes": ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL1", "DLL3", "DLL4", 
                  "JAG1", "JAG2", "FBXW7", "MAML1", "RBPJ"],
        "description": "Cell differentiation",
        "cancer_relevance": 0.75
    },
    "HORMONE_SIGNALING": {
        "name": "Hormone Signaling",
        "genes": ["ESR1", "ESR2", "AR", "PGR", "ERBB2", "FOXA1", "GATA3", 
                  "NCOA3", "NCOR1", "NCOR2"],
        "description": "Hormone receptor pathways",
        "cancer_relevance": 0.85
    },
    "IMMUNE_CHECKPOINT": {
        "name": "Immune Checkpoint",
        "genes": ["PD1", "PDL1", "CTLA4", "LAG3", "TIM3", "TIGIT", "B7H3", 
                  "IDO1", "CD47", "SIRPA"],
        "description": "Immune regulation",
        "cancer_relevance": 0.8
    },
    "CHROMATIN_REMODELING": {
        "name": "Chromatin Remodeling",
        "genes": ["ARID1A", "SMARCA4", "SMARCB1", "EZH2", "SETD2", "KMT2D", 
                  "CREBBP", "EP300", "ASXL1", "DNMT3A", "TET2", "IDH1", "IDH2"],
        "description": "Epigenetic regulation",
        "cancer_relevance": 0.8
    }
}

# Cancer-specific pathway associations
CANCER_PATHWAYS = {
    "breast": ["DNA_REPAIR", "CELL_CYCLE", "PI3K_AKT", "HORMONE_SIGNALING", "APOPTOSIS"],
    "colon": ["WNT_SIGNALING", "APOPTOSIS", "DNA_REPAIR", "CELL_CYCLE", "RAS_MAPK"],
    "lung": ["CELL_CYCLE", "PI3K_AKT", "APOPTOSIS", "RAS_MAPK", "CHROMATIN_REMODELING"],
    "prostate": ["HORMONE_SIGNALING", "PI3K_AKT", "CELL_CYCLE", "DNA_REPAIR"],
    "blood": ["APOPTOSIS", "CELL_CYCLE", "NOTCH_SIGNALING", "CHROMATIN_REMODELING"],
    "ovarian": ["DNA_REPAIR", "CELL_CYCLE", "PI3K_AKT", "APOPTOSIS"],
    "pancreatic": ["RAS_MAPK", "CELL_CYCLE", "DNA_REPAIR", "CHROMATIN_REMODELING"]
}

# Variant impact weights for pathway scoring
IMPACT_WEIGHTS = {
    'frameshift': 1.0,
    'nonsense': 1.0,
    'stop_gained': 1.0,
    'stop_lost': 0.9,
    'splice_site': 0.9,
    'missense': 0.5,
    'inframe': 0.3,
    'synonymous': 0.1,
    'unknown': 0.3
}


def get_variant_impact_weight(variant: Dict) -> float:
    """Get impact weight for a variant based on its consequence"""
    consequence = variant.get('consequence', 'unknown').lower()
    
    # Check various possible consequence fields
    for impact_type, weight in IMPACT_WEIGHTS.items():
        if impact_type in consequence:
            return weight
    
    # Check functional impact field
    functional_impact = variant.get('functional_impact', '').lower()
    if 'nonsense' in functional_impact:
        return 1.0
    elif 'frameshift' in functional_impact:
        return 1.0
    elif 'missense' in functional_impact:
        return 0.5
    
    return 0.3  # Default weight


def calculate_pathway_impact(pathway: Dict, affected_genes: Set[str], 
                           variant_impacts: Dict[str, float]) -> float:
    """Calculate the impact score for a pathway"""
    pathway_genes = set(pathway['genes'])
    affected_in_pathway = affected_genes.intersection(pathway_genes)
    
    if not affected_in_pathway:
        return 0.0
    
    # Calculate base score (fraction of pathway affected)
    base_score = len(affected_in_pathway) / len(pathway_genes) * 100
    
    # Weight by variant impact
    total_weight = sum(variant_impacts.get(gene, 0.3) for gene in affected_in_pathway)
    avg_weight = total_weight / len(affected_in_pathway) if affected_in_pathway else 0
    
    # Apply cancer relevance
    cancer_relevance = pathway.get('cancer_relevance', 0.8)
    
    # Calculate final score
    weighted_score = base_score * avg_weight * cancer_relevance
    
    return min(weighted_score, 100)  # Cap at 100


def analyze_pathway_interactions(disrupted_pathways: List[Dict]) -> Dict:
    """Analyze interactions between disrupted pathways"""
    interactions = []
    
    # Define known pathway interactions
    pathway_interactions = {
        ('DNA_REPAIR', 'CELL_CYCLE'): 'DNA damage checkpoint activation',
        ('PI3K_AKT', 'APOPTOSIS'): 'Survival signaling vs cell death',
        ('RAS_MAPK', 'CELL_CYCLE'): 'Growth signal to proliferation',
        ('WNT_SIGNALING', 'CELL_CYCLE'): 'Developmental signals to proliferation',
        ('HORMONE_SIGNALING', 'PI3K_AKT'): 'Hormone-driven growth signaling'
    }
    
    # Check for interacting pathways
    disrupted_ids = {p['pathway_id'] for p in disrupted_pathways}
    
    for (p1, p2), description in pathway_interactions.items():
        if p1 in disrupted_ids and p2 in disrupted_ids:
            interactions.append({
                'pathways': [p1, p2],
                'interaction': description,
                'impact': 'Synergistic pathway disruption'
            })
    
    return {
        'pathway_interactions': interactions,
        'interaction_count': len(interactions)
    }


def generate_pathway_recommendations(disrupted_pathways: List[Dict], 
                                   cancer_associations: Dict) -> List[str]:
    """Generate clinical recommendations based on pathway disruptions"""
    recommendations = []
    
    # Check for specific pathway patterns
    pathway_ids = {p['pathway_id'] for p in disrupted_pathways}
    
    if 'DNA_REPAIR' in pathway_ids:
        recommendations.append("Consider PARP inhibitor sensitivity testing")
        recommendations.append("Enhanced screening for DNA repair deficiency syndrome")
    
    if 'PI3K_AKT' in pathway_ids:
        recommendations.append("Consider PI3K/AKT/mTOR pathway inhibitors")
        recommendations.append("Monitor for metabolic dysregulation")
    
    if 'RAS_MAPK' in pathway_ids:
        recommendations.append("Consider MEK inhibitor therapy options")
        recommendations.append("Test for RAS pathway activation markers")
    
    if 'IMMUNE_CHECKPOINT' in pathway_ids:
        recommendations.append("Evaluate for immunotherapy eligibility")
        recommendations.append("Consider PD-1/PD-L1 expression testing")
    
    if 'HORMONE_SIGNALING' in pathway_ids:
        recommendations.append("Consider hormone receptor testing")
        recommendations.append("Evaluate endocrine therapy options")
    
    # Cancer-specific recommendations
    for cancer_type, associated_pathways in cancer_associations.items():
        if len(associated_pathways) >= 3:
            recommendations.append(f"High pathway burden suggests increased {cancer_type} cancer risk")
    
    return list(set(recommendations))  # Remove duplicates


def process(state: Dict) -> Dict:
    """Analyze pathway disruptions"""
    logger.info("Starting pathway disruption analysis")
    state["current_node"] = "pathway_analyzer"
    
    try:
        # Get variants
        variants = state.get('variant_details', 
                           state.get('filtered_variants', []))
        
        # Extract affected genes and their impacts
        affected_genes = set()
        variant_impacts = defaultdict(float)
        gene_mutations = defaultdict(list)
        
        for variant in variants:
            gene = variant.get('gene')
            if gene:
                affected_genes.add(gene)
                # Track highest impact per gene
                weight = get_variant_impact_weight(variant)
                if weight > variant_impacts[gene]:
                    variant_impacts[gene] = weight
                
                # Track mutations per gene
                gene_mutations[gene].append({
                    'type': variant.get('consequence', 'unknown'),
                    'impact': variant.get('functional_impact', 'Unknown impact'),
                    'protein_change': variant.get('protein_change', '')
                })
        
        # Analyze each pathway
        disrupted_pathways = []
        
        for pathway_id, pathway_data in PATHWAYS.items():
            impact_score = calculate_pathway_impact(
                pathway_data, affected_genes, variant_impacts
            )
            
            if impact_score > 10:  # Significance threshold
                affected_in_pathway = affected_genes.intersection(pathway_data['genes'])
                
                # Get mutation details for affected genes
                mutations = []
                for gene in affected_in_pathway:
                    if gene in gene_mutations:
                        # Take the most severe mutation
                        gene_muts = gene_mutations[gene]
                        most_severe = max(gene_muts, 
                                        key=lambda m: IMPACT_WEIGHTS.get(
                                            m['type'].lower(), 0.3))
                        mutations.append({
                            'gene': gene,
                            'type': most_severe['type'],
                            'effect': most_severe['impact']
                        })
                
                disrupted_pathways.append({
                    'name': pathway_data['name'],
                    'pathway_id': pathway_id,
                    'significance': round(impact_score, 1),
                    'affected_genes': list(affected_in_pathway),
                    'mutations': mutations,
                    'description': pathway_data['description'],
                    'genes_affected_ratio': f"{len(affected_in_pathway)}/{len(pathway_data['genes'])}"
                })
        
        # Sort by significance
        disrupted_pathways.sort(key=lambda x: x['significance'], reverse=True)
        
        # Build cancer pathway associations
        cancer_associations = {}
        for cancer_type, pathway_ids in CANCER_PATHWAYS.items():
            associated_pathways = [p['pathway_id'] for p in disrupted_pathways 
                                 if p['pathway_id'] in pathway_ids]
            if associated_pathways:
                cancer_associations[cancer_type] = associated_pathways
        
        # Analyze pathway interactions
        interaction_analysis = analyze_pathway_interactions(disrupted_pathways)
        
        # Generate recommendations
        recommendations = generate_pathway_recommendations(
            disrupted_pathways, cancer_associations
        )
        
        # Create pathway analysis result
        pathway_analysis = {
            'disrupted_pathways': disrupted_pathways,
            'cancer_pathway_associations': cancer_associations,
            'pathway_interactions': interaction_analysis['pathway_interactions'],
            'clinical_recommendations': recommendations,
            'summary': {
                'total_pathways_disrupted': len(disrupted_pathways),
                'highly_disrupted_pathways': sum(1 for p in disrupted_pathways 
                                               if p['significance'] > 50),
                'total_genes_affected': len(affected_genes),
                'pathway_interaction_count': interaction_analysis['interaction_count']
            }
        }
        
        # Update state
        state['pathway_analysis'] = pathway_analysis
        
        # Add to completed nodes
        completed = state.get("completed_nodes", [])
        if "pathway_analyzer" not in completed:
            completed.append("pathway_analyzer")
        state["completed_nodes"] = completed
        
        logger.info(f"Analyzed {len(disrupted_pathways)} disrupted pathways")
        
    except Exception as e:
        logger.error(f"Error in pathway analysis: {str(e)}")
        state["errors"] = state.get("errors", []) + [{
            "node": "pathway_analyzer",
            "error": str(e),
            "timestamp": datetime.now().isoformat()
        }]
        # Set empty results on error
        state['pathway_analysis'] = {
            'disrupted_pathways': [],
            'cancer_pathway_associations': {}
        }
    
    return state 