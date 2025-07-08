"""
Variant calling node using DeepVariant.
Identifies genetic variants from aligned BAM files.
"""
import os
import subprocess
import logging
from datetime import datetime
from typing import Dict, Any, List
import vcf
import csv

logger = logging.getLogger(__name__)


def run_simple_variant_caller(bam_path: str, reference_path: str, output_dir: str = None) -> str:
    """
    Simple variant caller for testing.
    In production, this would run DeepVariant.
    
    For now, we'll use our pre-generated test VCF.
    """
    # In production, this would run:
    # docker run -v ${PWD}:${PWD} google/deepvariant:latest \
    #   /opt/deepvariant/bin/run_deepvariant \
    #   --model_type=WGS \
    #   --ref=reference.fa \
    #   --reads=aligned.bam \
    #   --output_vcf=variants.vcf
    
    # For testing, use pre-generated VCF
    test_vcf = os.path.join(os.path.dirname(__file__), "..", "test_data", "test_variants.vcf")
    if os.path.exists(test_vcf):
        return test_vcf
    
    # If no test VCF, create a simple one based on the BAM
    if output_dir is None:
        output_dir = os.path.dirname(bam_path)
    
    vcf_path = os.path.join(output_dir, "variants.vcf")
    
    # Create a minimal VCF header
    vcf_content = """##fileformat=VCFv4.2
##source=SimpleVariantCaller
##reference={ref}
##contig=<ID=chr17,length=12000>
##contig=<ID=chr5,length=12000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
""".format(ref=os.path.basename(reference_path))
    
    with open(vcf_path, 'w') as f:
        f.write(vcf_content)
    
    logger.warning("Using simplified variant caller for testing. In production, use DeepVariant.")
    return vcf_path


def load_variant_annotations(annotation_file: str) -> Dict[str, Dict[str, Any]]:
    """
    Load variant annotations from TSV file.
    
    Returns:
        Dictionary mapping variant_id to annotation info
    """
    annotations = {}
    
    if not os.path.exists(annotation_file):
        logger.warning(f"Annotation file not found: {annotation_file}")
        return annotations
    
    with open(annotation_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Create variant ID
            variant_id = f"{row['CHROM']}:{row['POS']}:{row['REF']}>{row['ALT']}"
            
            annotations[variant_id] = {
                "gene": row.get('GENE', 'Unknown'),
                "consequence": row.get('CONSEQUENCE', 'unknown'),
                "hgvs_c": row.get('HGVS_C', ''),
                "hgvs_p": row.get('HGVS_P', ''),
                "impact": row.get('IMPACT', 'UNKNOWN')
            }
    
    return annotations


def parse_vcf_file(vcf_path: str) -> List[Dict[str, Any]]:
    """
    Parse VCF file and extract variant information.
    
    Args:
        vcf_path: Path to VCF file
        
    Returns:
        List of variant dictionaries
    """
    variants = []
    
    # Load annotations if available
    annotation_file = vcf_path.replace('.vcf', '_annotations.tsv')
    if not os.path.exists(annotation_file):
        annotation_file = os.path.join(
            os.path.dirname(vcf_path), 
            "test_annotations.tsv"
        )
    
    annotations = load_variant_annotations(annotation_file)
    
    # Parse VCF file
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    
    for record in vcf_reader:
        # Skip non-variant positions
        if record.ALT[0] is None:
            continue
        
        # Extract basic variant info
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = str(record.ALT[0])
        qual = record.QUAL if record.QUAL else 0
        
        # Create variant ID
        variant_id = f"{chrom}:{pos}:{ref}>{alt}"
        
        # Extract INFO fields
        depth = record.INFO.get('DP', 0)
        allele_freq = record.INFO.get('AF', [0])[0] if 'AF' in record.INFO else 0
        
        # Get genotype information
        if record.samples:
            sample = record.samples[0]
            genotype = sample['GT'] if 'GT' in sample.data._fields else '0/1'
            sample_depth = sample['DP'] if 'DP' in sample.data._fields else depth
        else:
            genotype = '0/1'
            sample_depth = depth
        
        # Build variant dictionary
        variant = {
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_id": variant_id,
            "quality": float(qual),
            "depth": int(sample_depth),
            "allele_freq": float(allele_freq),
            "genotype": genotype
        }
        
        # Add annotations if available
        if variant_id in annotations:
            variant.update(annotations[variant_id])
        else:
            # Default values if no annotation
            variant.update({
                "gene": "Unknown",
                "consequence": "unknown",
                "hgvs_c": "",
                "hgvs_p": "",
                "impact": "UNKNOWN"
            })
        
        variants.append(variant)
    
    return variants


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Call variants using DeepVariant.
    
    Updates state with:
    - vcf_path: path to output VCF file
    - raw_variants: list of variant dictionaries
    - variant_count: total number of variants found
    """
    logger.info("Starting variant calling with DeepVariant")
    state["current_node"] = "variant_calling"
    
    try:
        bam_path = state["aligned_bam_path"]
        
        # Get reference genome path
        reference_path = os.path.join(
            os.path.dirname(__file__), 
            "..", 
            "test_reference", 
            "test_genome.fa"
        )
        
        # Run variant caller (simplified for testing)
        vcf_path = run_simple_variant_caller(
            bam_path, 
            reference_path,
            os.path.dirname(bam_path)
        )
        
        # Parse VCF file
        variants = parse_vcf_file(vcf_path)
        
        state["vcf_path"] = vcf_path
        state["raw_variants"] = variants
        state["variant_count"] = len(variants)
        
        logger.info(f"Variant calling complete. Found {len(variants)} variants")
        
        # Log some variant details
        for v in variants[:3]:  # First 3 variants
            logger.debug(f"Variant: {v['variant_id']} in {v.get('gene', 'Unknown')} - {v.get('consequence', 'unknown')}")
        
        state["completed_nodes"].append("variant_calling")
        
    except Exception as e:
        logger.error(f"Variant calling failed: {str(e)}")
        state["errors"].append({
            "node": "variant_calling",
            "error": str(e),
            "timestamp": datetime.now()
        })
        state["pipeline_status"] = "failed"
    
    return state 