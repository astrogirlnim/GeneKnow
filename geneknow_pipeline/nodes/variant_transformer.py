"""
Variant Transformer Node for GeneKnow pipeline.
Adds transformation details to variants including protein changes.
Uses hardcoded genetic tables for offline operation.
"""
from typing import Dict, List
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

# Genetic code table - this is universal and doesn't change
# This MUST be hardcoded for offline operation
CODON_TABLE = {
    # Phenylalanine
    'TTT': 'Phe', 'TTC': 'Phe',
    # Leucine
    'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    # Isoleucine
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile',
    # Methionine (Start)
    'ATG': 'Met',
    # Valine
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    # Serine
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser', 'AGT': 'Ser', 'AGC': 'Ser',
    # Proline
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    # Threonine
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    # Alanine
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    # Tyrosine
    'TAT': 'Tyr', 'TAC': 'Tyr',
    # Stop codons
    'TAA': '*', 'TAG': '*', 'TGA': '*',
    # Histidine
    'CAT': 'His', 'CAC': 'His',
    # Glutamine
    'CAA': 'Gln', 'CAG': 'Gln',
    # Asparagine
    'AAT': 'Asn', 'AAC': 'Asn',
    # Lysine
    'AAA': 'Lys', 'AAG': 'Lys',
    # Aspartic acid
    'GAT': 'Asp', 'GAC': 'Asp',
    # Glutamic acid
    'GAA': 'Glu', 'GAG': 'Glu',
    # Cysteine
    'TGT': 'Cys', 'TGC': 'Cys',
    # Tryptophan
    'TGG': 'Trp',
    # Arginine
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
    # Glycine
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

# Amino acid properties - standard biochemical properties
# Also necessary for offline functional impact assessment
AA_PROPERTIES = {
    'Phe': {'hydrophobic': True, 'size': 'large', 'aromatic': True, 'charge': 'neutral'},
    'Leu': {'hydrophobic': True, 'size': 'large', 'aromatic': False, 'charge': 'neutral'},
    'Ile': {'hydrophobic': True, 'size': 'medium', 'aromatic': False, 'charge': 'neutral'},
    'Met': {'hydrophobic': True, 'size': 'medium', 'aromatic': False, 'charge': 'neutral'},
    'Val': {'hydrophobic': True, 'size': 'small', 'aromatic': False, 'charge': 'neutral'},
    'Ser': {'hydrophobic': False, 'size': 'small', 'aromatic': False, 'charge': 'neutral', 'polar': True},
    'Pro': {'hydrophobic': False, 'size': 'small', 'aromatic': False, 'charge': 'neutral', 'special': 'helix_breaker'},
    'Thr': {'hydrophobic': False, 'size': 'small', 'aromatic': False, 'charge': 'neutral', 'polar': True},
    'Ala': {'hydrophobic': True, 'size': 'tiny', 'aromatic': False, 'charge': 'neutral'},
    'Tyr': {'hydrophobic': False, 'size': 'large', 'aromatic': True, 'charge': 'neutral', 'polar': True},
    'His': {'hydrophobic': False, 'size': 'medium', 'aromatic': True, 'charge': 'positive'},
    'Gln': {'hydrophobic': False, 'size': 'medium', 'aromatic': False, 'charge': 'neutral', 'polar': True},
    'Asn': {'hydrophobic': False, 'size': 'small', 'aromatic': False, 'charge': 'neutral', 'polar': True},
    'Lys': {'hydrophobic': False, 'size': 'large', 'aromatic': False, 'charge': 'positive'},
    'Asp': {'hydrophobic': False, 'size': 'small', 'aromatic': False, 'charge': 'negative'},
    'Glu': {'hydrophobic': False, 'size': 'medium', 'aromatic': False, 'charge': 'negative'},
    'Cys': {'hydrophobic': True, 'size': 'small', 'aromatic': False, 'charge': 'neutral', 'special': 'disulfide'},
    'Trp': {'hydrophobic': True, 'size': 'large', 'aromatic': True, 'charge': 'neutral'},
    'Arg': {'hydrophobic': False, 'size': 'large', 'aromatic': False, 'charge': 'positive'},
    'Gly': {'hydrophobic': False, 'size': 'tiny', 'aromatic': False, 'charge': 'neutral', 'special': 'flexible'},
    '*': {'hydrophobic': False, 'size': 'none', 'aromatic': False, 'charge': 'none', 'special': 'stop'}
}


def simulate_codon_context(variant: Dict) -> tuple:
    """
    Simulate codon context for a variant.
    In production, this would fetch from reference genome.
    For offline operation, we create plausible codons.
    """
    ref = variant.get('ref', '')
    alt = variant.get('alt', '')
    gene = variant.get('gene', '')
    
    # For known cancer genes, use common codons
    cancer_gene_codons = {
        'TP53': {'ref': 'CGT', 'alt': 'CAT'},  # R273H - common TP53 mutation
        'KRAS': {'ref': 'GGT', 'alt': 'GAT'},  # G12D - common KRAS mutation
        'BRAF': {'ref': 'GTG', 'alt': 'GAG'},  # V600E - common BRAF mutation
        'EGFR': {'ref': 'CTG', 'alt': 'CGG'},  # L858R - common EGFR mutation
        'PIK3CA': {'ref': 'CAT', 'alt': 'AAT'},  # H1047N - common PIK3CA mutation
    }
    
    if gene in cancer_gene_codons and len(ref) == 1 and len(alt) == 1:
        return cancer_gene_codons[gene]['ref'], cancer_gene_codons[gene]['alt']
    
    # For other variants, create plausible codons
    if len(ref) == 1 and len(alt) == 1:
        # Create a codon with the variant in the middle
        ref_codon = 'A' + ref + 'T'  # Simple default
        alt_codon = 'A' + alt + 'T'
        
        # Try to make it a valid codon
        if ref_codon in CODON_TABLE:
            return ref_codon, alt_codon
        else:
            # Use a common codon as template
            return 'ATT', 'A' + alt + 'T'
    
    return '', ''


def get_protein_change(variant: Dict) -> Dict:
    """Calculate protein-level changes from DNA variant"""
    ref_codon, alt_codon = simulate_codon_context(variant)
    
    if not ref_codon or not alt_codon:
        return {
            "original": "",
            "mutated": "",
            "amino_acid_change": "Unknown",
            "effect": "Unable to determine protein change"
        }
    
    ref_aa = CODON_TABLE.get(ref_codon, 'X')
    alt_aa = CODON_TABLE.get(alt_codon, 'X')
    
    # Determine functional impact
    impact = assess_functional_impact(ref_aa, alt_aa, variant.get('gene', ''))
    
    return {
        "original": ref_codon,
        "mutated": alt_codon,
        "amino_acid_change": f"{ref_aa}→{alt_aa}",
        "effect": impact
    }


def assess_functional_impact(ref_aa: str, alt_aa: str, gene: str) -> str:
    """Assess the functional impact of amino acid change"""
    # Synonymous mutation
    if ref_aa == alt_aa:
        return "Synonymous - no amino acid change"
    
    # Nonsense mutation
    if alt_aa == '*':
        return f"Nonsense mutation - premature stop codon in {gene}"
    
    # Loss of stop codon
    if ref_aa == '*':
        return f"Stop loss - extended protein in {gene}"
    
    # Get properties
    ref_props = AA_PROPERTIES.get(ref_aa, {})
    alt_props = AA_PROPERTIES.get(alt_aa, {})
    
    # Check for drastic changes
    impacts = []
    
    # Charge change
    if ref_props.get('charge') != alt_props.get('charge'):
        if ref_props.get('charge') == 'positive' and alt_props.get('charge') == 'negative':
            impacts.append("Charge reversal (positive to negative)")
        elif ref_props.get('charge') == 'negative' and alt_props.get('charge') == 'positive':
            impacts.append("Charge reversal (negative to positive)")
        else:
            impacts.append(f"Charge change ({ref_props.get('charge')} to {alt_props.get('charge')})")
    
    # Hydrophobicity change
    if ref_props.get('hydrophobic') != alt_props.get('hydrophobic'):
        if ref_props.get('hydrophobic'):
            impacts.append("Loss of hydrophobicity")
        else:
            impacts.append("Gain of hydrophobicity")
    
    # Size change
    size_order = {'tiny': 0, 'small': 1, 'medium': 2, 'large': 3}
    ref_size = size_order.get(ref_props.get('size', 'medium'), 2)
    alt_size = size_order.get(alt_props.get('size', 'medium'), 2)
    
    if abs(ref_size - alt_size) >= 2:
        impacts.append(f"Significant size change ({ref_props.get('size')} to {alt_props.get('size')})")
    
    # Special residues
    if ref_props.get('special') or alt_props.get('special'):
        if ref_props.get('special') == 'helix_breaker':
            impacts.append("Loss of proline - may affect protein structure")
        elif alt_props.get('special') == 'helix_breaker':
            impacts.append("Introduction of proline - may disrupt secondary structure")
        elif ref_props.get('special') == 'disulfide':
            impacts.append("Loss of cysteine - may disrupt disulfide bonds")
        elif alt_props.get('special') == 'disulfide':
            impacts.append("Introduction of cysteine - may form new disulfide bonds")
    
    # Known hotspot mutations in cancer genes
    cancer_hotspots = {
        'TP53': {'R273': 'DNA binding domain hotspot', 'R248': 'DNA binding domain hotspot'},
        'KRAS': {'G12': 'GTPase domain hotspot', 'G13': 'GTPase domain hotspot'},
        'BRAF': {'V600': 'Kinase domain activation hotspot'},
        'PIK3CA': {'H1047': 'Kinase domain hotspot', 'E545': 'Helical domain hotspot'}
    }
    
    if gene in cancer_hotspots:
        # Check if this might be a hotspot (simplified check)
        for hotspot, description in cancer_hotspots[gene].items():
            if ref_aa in hotspot or alt_aa in hotspot:
                impacts.append(f"Mutation at {description}")
                break
    
    # Generate final impact description
    if impacts:
        return " - ".join(impacts)
    elif ref_props.get('size') == alt_props.get('size') and ref_props.get('charge') == alt_props.get('charge'):
        return "Conservative substitution - minimal impact expected"
    else:
        return "Moderate change - potential functional impact"


def process(state: Dict) -> Dict:
    """Add transformation details to variants"""
    logger.info("Starting variant transformation")
    state["current_node"] = "variant_transformer"
    
    try:
        # Get variants - prefer classified variants if available
        variants = state.get('classified_variants', state.get('filtered_variants', []))
        
        # Process each variant
        transformed_variants = []
        for variant in variants:
            # Copy variant to avoid modifying original
            transformed_variant = variant.copy()
            
            # Add transformation details
            transformation = get_protein_change(transformed_variant)
            transformed_variant['transformation'] = transformation
            
            # Add functional impact assessment
            transformed_variant['functional_impact'] = transformation['effect']
            
            # Add protein change notation if available
            if transformation['amino_acid_change'] != 'Unknown':
                aa_change = transformation['amino_acid_change']
                transformed_variant['protein_change'] = f"p.{aa_change}"
            
            transformed_variants.append(transformed_variant)
        
        # Update state
        state['variant_details'] = transformed_variants
        
        # Add summary statistics
        impact_summary = {
            'total_variants': len(transformed_variants),
            'nonsense_mutations': sum(1 for v in transformed_variants 
                                    if 'nonsense' in v.get('functional_impact', '').lower()),
            'charge_changes': sum(1 for v in transformed_variants 
                                if 'charge' in v.get('functional_impact', '').lower()),
            'conservative_changes': sum(1 for v in transformed_variants 
                                      if 'conservative' in v.get('functional_impact', '').lower())
        }
        state['variant_transformation_summary'] = impact_summary
        
        # Add to completed nodes
        completed = state.get("completed_nodes", [])
        if "variant_transformer" not in completed:
            completed.append("variant_transformer")
        state["completed_nodes"] = completed
        
        logger.info(f"Transformed {len(transformed_variants)} variants")
        
    except Exception as e:
        logger.error(f"Error in variant transformation: {str(e)}")
        state["errors"] = state.get("errors", []) + [{
            "node": "variant_transformer",
            "error": str(e),
            "timestamp": datetime.now().isoformat()
        }]
        # Pass through variants without transformation on error
        state['variant_details'] = state.get('classified_variants', 
                                            state.get('filtered_variants', []))
    
    return state 