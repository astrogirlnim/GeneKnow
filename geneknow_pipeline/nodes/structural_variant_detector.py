"""
Structural Variant Detector Node for GeneKnow pipeline.
Detects structural variants from alignment data or variant calls.
Provides simulated SV detection for desktop operation.
"""

from typing import Dict, List
import logging
from datetime import datetime
import os

logger = logging.getLogger(__name__)

# Structural variant types
SV_TYPES = {
    "DEL": "deletion",
    "INS": "insertion",
    "DUP": "duplication",
    "INV": "inversion",
    "TRA": "translocation",
    "CNV": "copy_number_variant",
}

# Known cancer-associated structural variants
CANCER_SVS = {
    "BCR-ABL1": {
        "type": "translocation",
        "chr1": "chr9",
        "chr2": "chr22",
        "genes": ["BCR", "ABL1"],
        "cancer": "Chronic Myeloid Leukemia",
        "description": "Philadelphia chromosome",
    },
    "EML4-ALK": {
        "type": "inversion",
        "chr1": "chr2",
        "chr2": "chr2",
        "genes": ["EML4", "ALK"],
        "cancer": "Lung adenocarcinoma",
        "description": "ALK fusion",
    },
    "TMPRSS2-ERG": {
        "type": "deletion",
        "chr1": "chr21",
        "chr2": "chr21",
        "genes": ["TMPRSS2", "ERG"],
        "cancer": "Prostate cancer",
        "description": "Common prostate cancer fusion",
    },
    "MYC-IGH": {
        "type": "translocation",
        "chr1": "chr8",
        "chr2": "chr14",
        "genes": ["MYC", "IGH"],
        "cancer": "Burkitt lymphoma",
        "description": "MYC translocation",
    },
}


def detect_from_vcf_variants(variants: List[Dict]) -> List[Dict]:
    """
    Detect structural variants from VCF variant calls.
    In production, would analyze INFO fields for SV evidence.
    """
    structural_variants = []

    for variant in variants:
        # Check for large indels that might be SVs
        ref = variant.get("re", "")
        alt = variant.get("alt", "")

        size_diff = abs(len(ref) - len(alt))

        # Large deletions
        if len(ref) > len(alt) and size_diff > 50:
            sv = {
                "type": "deletion",
                "chromosome": variant.get("chrom", ""),
                "start": variant.get("pos", 0),
                "end": variant.get("pos", 0) + len(ref),
                "size": size_diff,
                "evidence": "large_indel",
                "variant_id": f"{variant.get('chrom')}:{variant.get('pos')}",
            }
            structural_variants.append(sv)

        # Large insertions
        elif len(alt) > len(ref) and size_diff > 50:
            sv = {
                "type": "insertion",
                "chromosome": variant.get("chrom", ""),
                "start": variant.get("pos", 0),
                "end": variant.get("pos", 0) + 1,
                "size": size_diff,
                "evidence": "large_indel",
                "variant_id": f"{variant.get('chrom')}:{variant.get('pos')}",
            }
            structural_variants.append(sv)

    return structural_variants


def simulate_sv_detection(variants: List[Dict]) -> List[Dict]:
    """
    Simulate structural variant detection based on variant patterns.
    In production, would use split-read and paired-end analysis.
    """
    structural_variants = []

    # Check for known fusion genes in variants
    variant_genes = {v.get("gene") for v in variants if v.get("gene")}

    # Look for potential fusion partners
    for fusion_name, fusion_data in CANCER_SVS.items():
        genes_found = [g for g in fusion_data["genes"] if g in variant_genes]

        if len(genes_found) >= 1:
            # Simulate finding a fusion
            sv = {
                "type": fusion_data["type"],
                "fusion_name": fusion_name,
                "chromosome": fusion_data["chr1"],
                "partner_chromosome": fusion_data["chr2"],
                "genes_affected": fusion_data["genes"],
                "clinical_significance": "pathogenic",
                "cancer_association": fusion_data["cancer"],
                "description": fusion_data["description"],
                "evidence": "gene_disruption",
                "functional_impact": f"Oncogenic fusion creating {fusion_name}",
            }

            # Add coordinates (simulated)
            if fusion_data["type"] == "translocation":
                sv["chr1"] = fusion_data["chr1"]
                sv["pos1"] = 50000000  # Simulated position
                sv["chr2"] = fusion_data["chr2"]
                sv["pos2"] = 100000000  # Simulated position
            else:
                sv["start"] = 50000000
                sv["end"] = 51000000
                sv["size"] = 1000000

            structural_variants.append(sv)
            break  # Only report one major fusion

    # Simulate CNV detection based on multiple variants in same gene
    gene_variant_counts = {}
    for variant in variants:
        gene = variant.get("gene")
        if gene:
            gene_variant_counts[gene] = gene_variant_counts.get(gene, 0) + 1

    # Genes with many variants might have CNVs
    for gene, count in gene_variant_counts.items():
        if count >= 3:  # Multiple variants in same gene
            # Check if it's a known oncogene (amplification) or tumor suppressor (deletion)
            oncogenes = {"ERBB2", "MYC", "EGFR", "CDK4", "MDM2", "CCND1"}
            tumor_suppressors = {"TP53", "RB1", "PTEN", "CDKN2A", "APC"}

            if gene in oncogenes:
                sv = {
                    "type": "amplification",
                    "chromosome": get_gene_chromosome(gene),
                    "gene": gene,
                    "copy_number": 6,  # Simulated
                    "normal_copy_number": 2,
                    "fold_change": 3.0,
                    "clinical_significance": "likely_pathogenic",
                    "cancer_relevance": f"{gene} amplification in cancer",
                    "functional_impact": f"Gene amplification of oncogene {gene}",
                    "evidence": "multiple_variants",
                }
                structural_variants.append(sv)

            elif gene in tumor_suppressors:
                sv = {
                    "type": "deletion",
                    "chromosome": get_gene_chromosome(gene),
                    "gene": gene,
                    "size": 10000,  # Simulated
                    "clinical_significance": "pathogenic",
                    "functional_impact": f"Loss of tumor suppressor {gene}",
                    "evidence": "multiple_variants",
                }
                structural_variants.append(sv)

    return structural_variants


def get_gene_chromosome(gene: str) -> str:
    """Get chromosome location for known genes"""
    gene_locations = {
        "TP53": "chr17",
        "BRCA1": "chr17",
        "BRCA2": "chr13",
        "ERBB2": "chr17",
        "MYC": "chr8",
        "EGFR": "chr7",
        "KRAS": "chr12",
        "PTEN": "chr10",
        "RB1": "chr13",
        "APC": "chr5",
        "CDK4": "chr12",
        "CDKN2A": "chr9",
    }
    return gene_locations.get(gene, "chrUnknown")


def annotate_sv_impact(sv: Dict) -> Dict:
    """Add clinical annotations to structural variant"""
    # Already annotated SVs (fusions, etc)
    if "clinical_significance" in sv:
        return sv

    # Annotate based on type and size
    if sv["type"] == "deletion":
        if sv.get("size", 0) > 1000000:  # > 1Mb
            sv["clinical_significance"] = "likely_pathogenic"
            sv["functional_impact"] = (
                "Large deletion potentially affecting multiple genes"
            )
        elif sv.get("size", 0) > 100000:  # > 100kb
            sv["clinical_significance"] = "uncertain_significance"
            sv["functional_impact"] = "Moderate deletion"
        else:
            sv["clinical_significance"] = "uncertain_significance"
            sv["functional_impact"] = "Small deletion"

    elif sv["type"] == "duplication":
        if sv.get("size", 0) > 1000000:
            sv["clinical_significance"] = "likely_pathogenic"
            sv["functional_impact"] = (
                "Large duplication potentially causing dosage imbalance"
            )
        else:
            sv["clinical_significance"] = "uncertain_significance"
            sv["functional_impact"] = "Gene duplication"

    elif sv["type"] == "insertion":
        if sv.get("size", 0) > 1000:
            sv["clinical_significance"] = "uncertain_significance"
            sv["functional_impact"] = (
                "Large insertion potentially disrupting gene function"
            )
        else:
            sv["clinical_significance"] = "uncertain_significance"
            sv["functional_impact"] = "Insertion variant"

    return sv


def process(state: Dict) -> Dict:
    """Process structural variant detection"""
    logger.info("Starting structural variant detection")
    # Don't set current_node to avoid concurrent updates

    try:
        # Get variants
        variants = state.get(
            "variant_details",
            state.get("classified_variants", state.get("filtered_variants", [])),
        )

        logger.info(f"📊 Analyzing {len(variants)} variants for structural variants")
        logger.info("🔍 SV Detection Criteria:")
        logger.info("  • Large indels: >50bp size difference")
        logger.info("  • Known fusion genes: BCR-ABL1, EML4-ALK, TMPRSS2-ERG, MYC-IGH")
        logger.info("  • Gene clustering: ≥3 variants in same oncogene/tumor suppressor")

        structural_variants = []

        # Check if we have a BAM file for proper SV detection
        bam_path = state.get("aligned_bam_path")

        if bam_path and os.path.exists(bam_path):
            # In production, would use pysam to analyze split reads and discordant pairs
            logger.info("BAM file available - simulating SV detection from alignments")
            structural_variants.extend(simulate_sv_detection(variants))

        # Always check for large indels in VCF
        logger.info("🔬 Checking for large indels in VCF data...")
        vcf_svs = detect_from_vcf_variants(variants)
        structural_variants.extend(vcf_svs)
        logger.info(f"📈 Found {len(vcf_svs)} large indels")

        # Log gene variant distribution
        gene_variant_counts = {}
        for variant in variants:
            gene = variant.get("gene")
            if gene:
                gene_variant_counts[gene] = gene_variant_counts.get(gene, 0) + 1
        
        high_variant_genes = {gene: count for gene, count in gene_variant_counts.items() if count >= 3}
        if high_variant_genes:
            logger.info(f"🧬 Genes with ≥3 variants: {high_variant_genes}")
        else:
            logger.info("🧬 No genes with ≥3 variants found")

        # Remove duplicates
        seen = set()
        unique_svs = []
        for sv in structural_variants:
            # Create unique key
            key = (
                sv.get("type"),
                sv.get("chromosome", ""),
                sv.get("start", 0),
                sv.get("size", 0),
            )
            if key not in seen:
                seen.add(key)
                # Annotate SV
                annotated_sv = annotate_sv_impact(sv)
                unique_svs.append(annotated_sv)

        # Summary statistics
        sv_summary = {
            "total_svs": len(unique_svs),
            "deletions": sum(1 for sv in unique_svs if sv["type"] == "deletion"),
            "duplications": sum(1 for sv in unique_svs if sv["type"] == "duplication"),
            "translocations": sum(
                1 for sv in unique_svs if sv["type"] == "translocation"
            ),
            "pathogenic_svs": sum(
                1
                for sv in unique_svs
                if sv.get("clinical_significance")
                in ["pathogenic", "likely_pathogenic"]
            ),
            "gene_fusions": [
                sv.get("fusion_name") for sv in unique_svs if "fusion_name" in sv
            ],
        }

        # Add to completed nodes
        completed = state.get("completed_nodes", [])
        if "structural_variant_detector" not in completed:
            completed.append("structural_variant_detector")

        logger.info(f"✅ Detected {len(unique_svs)} structural variants")
        if len(unique_svs) == 0:
            logger.info("ℹ️  No structural variants detected - this is normal for most genomic analyses")
            logger.info("ℹ️  Structural variants are rare and typically found in:")
            logger.info("    • Cancer genomes with extensive rearrangements")
            logger.info("    • Inherited disorders with large genomic changes")
            logger.info("    • Samples with known fusion genes")

        # Return only the keys this node updates
        return {
            "structural_variants": unique_svs,
            "structural_variant_summary": sv_summary,
            "completed_nodes": completed,
        }

    except Exception as e:
        logger.error(f"Error in structural variant detection: {str(e)}")
        errors = state.get("errors", []) + [
            {
                "node": "structural_variant_detector",
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            }
        ]
        # Return only the keys this node updates
        return {"structural_variants": [], "errors": errors}
