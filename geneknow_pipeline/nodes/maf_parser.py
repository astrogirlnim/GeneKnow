"""
MAF (Mutation Annotation Format) parser for TCGA data.
Converts MAF files into the same variant format used by VCF files.
"""
import os
import logging
import pandas as pd
from typing import Dict, Any, List
from datetime import datetime

logger = logging.getLogger(__name__)

# Clinical significance mapping based on variant classification
CLASSIFICATION_TO_SIGNIFICANCE = {
    'Missense_Mutation': 'Uncertain significance',
    'Nonsense_Mutation': 'Likely pathogenic',
    'Frame_Shift_Del': 'Pathogenic',
    'Frame_Shift_Ins': 'Pathogenic',
    'In_Frame_Del': 'Likely pathogenic',
    'In_Frame_Ins': 'Likely pathogenic',
    'Silent': 'Likely benign',
    'Splice_Site': 'Pathogenic',
    'Translation_Start_Site': 'Pathogenic',
    'Nonstop_Mutation': 'Pathogenic',
    '3\'UTR': 'Likely benign',
    '5\'UTR': 'Likely benign',
    'Intron': 'Benign',
    'IGR': 'Benign',
    'RNA': 'Uncertain significance',
    'lincRNA': 'Uncertain significance'
}


def parse_maf_file(file_path: str) -> List[Dict[str, Any]]:
    """
    Parse a MAF file and convert to standardized variant format.
    
    Args:
        file_path: Path to MAF file
        
    Returns:
        List of variant dictionaries
    """
    logger.info(f"Parsing MAF file: {file_path}")
    
    try:
        # Read MAF file (tab-separated)
        df = pd.read_csv(file_path, sep='\t', low_memory=False, comment='#')
        logger.info(f"Loaded MAF with {len(df)} variants")
        
        variants = []
        
        for idx, row in df.iterrows():
            try:
                # Get chromosome - handle different formats
                chrom = str(row.get('Chromosome', ''))
                if not chrom.startswith('chr'):
                    chrom = f"chr{chrom}"
                
                # Calculate allele frequency if possible
                t_depth = row.get('t_depth', 0)
                t_alt_count = row.get('t_alt_count', 0)
                if pd.notna(row.get('tumor_f')):
                    allele_freq = float(row.get('tumor_f'))
                elif t_depth > 0 and pd.notna(t_alt_count):
                    allele_freq = float(t_alt_count) / float(t_depth)
                else:
                    allele_freq = 0.5  # Default if no frequency data
                
                # Create variant in our standard format
                variant = {
                    'chrom': chrom,
                    'pos': int(row.get('Start_Position', 0)),
                    'ref': str(row.get('Reference_Allele', 'N')),
                    'alt': str(row.get('Tumor_Seq_Allele2', 'N')),
                    'gene': str(row.get('Hugo_Symbol', 'Unknown')),
                    'variant_id': f"{chrom}:{row.get('Start_Position')}:{row.get('Reference_Allele')}>{row.get('Tumor_Seq_Allele2')}",
                    
                    # MAF-specific fields
                    'variant_classification': row.get('Variant_Classification', ''),
                    'variant_type': row.get('Variant_Type', ''),
                    'consequence': row.get('Consequence', row.get('Variant_Classification', '')),
                    'protein_change': row.get('HGVSp', ''),
                    
                    # Quality metrics (MAF files don't have QUAL like VCF)
                    'qual': 100.0,  # Default high quality for MAF variants
                    'depth': int(t_depth) if pd.notna(t_depth) else 100,
                    'allele_freq': allele_freq,
                    
                    # Clinical interpretation
                    'clinical_significance': CLASSIFICATION_TO_SIGNIFICANCE.get(
                        row.get('Variant_Classification', ''), 
                        'Uncertain significance'
                    ),
                    
                    # Sample information
                    'sample_id': str(row.get('Tumor_Sample_Barcode', 'Unknown')),
                    'normal_sample_id': str(row.get('Matched_Norm_Sample_Barcode', '')),
                    
                    # Additional annotations if available
                    'dbsnp_rs': str(row.get('dbSNP_RS', '')),
                    'existing_variation': str(row.get('Existing_variation', ''))
                }
                
                # Only include variants that pass basic filters
                if variant['ref'] not in ['N', '-', ''] and variant['alt'] not in ['N', '-', '']:
                    variants.append(variant)
                    
            except Exception as e:
                logger.warning(f"Error parsing variant at row {idx}: {e}")
                continue
        
        logger.info(f"Successfully parsed {len(variants)} variants from MAF")
        return variants
        
    except Exception as e:
        logger.error(f"Failed to parse MAF file: {e}")
        raise


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process MAF file node for LangGraph pipeline.
    
    This node is used when the input file is a MAF file instead of VCF.
    It parses the MAF and converts it to our standard variant format.
    """
    logger.info("Processing MAF file")
    state["current_node"] = "maf_parser"
    
    try:
        file_path = state["file_path"]
        
        # Parse MAF file
        variants = parse_maf_file(file_path)
        
        # Update state with parsed variants
        state["variants"] = variants
        state["variant_count"] = len(variants)
        
        # MAF files are already filtered/processed, so we can use them as filtered_variants
        state["filtered_variants"] = variants
        
        # Add metadata
        state["file_metadata"]["variant_source"] = "MAF"
        state["file_metadata"]["total_variants"] = len(variants)
        
        # Calculate MAF-specific statistics
        if variants:
            df = pd.DataFrame(variants)
            state["file_metadata"]["maf_info"] = {
                "unique_genes": df['gene'].nunique(),
                "unique_samples": df['sample_id'].nunique() if 'sample_id' in df else 1,
                "variant_classifications": dict(df['variant_classification'].value_counts().head(10)),
                "top_mutated_genes": dict(df['gene'].value_counts().head(10))
            }
        
        state["completed_nodes"].append("maf_parser")
        logger.info(f"MAF processing complete: {len(variants)} variants")
        
    except Exception as e:
        logger.error(f"MAF processing failed: {str(e)}")
        state["errors"].append({
            "node": "maf_parser",
            "error": str(e),
            "timestamp": datetime.now()
        })
        state["pipeline_status"] = "failed"
    
    return state 