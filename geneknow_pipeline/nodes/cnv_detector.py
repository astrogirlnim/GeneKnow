"""
CNV Detector Node for GeneKnow pipeline.
Detects copy number variants using depth-based analysis.
Provides simulated CNV detection for desktop operation.
"""

import numpy as np
from typing import Dict, List
import logging
from datetime import datetime
from collections import defaultdict

logger = logging.getLogger(__name__)

# Known cancer-associated CNVs
CANCER_CNVS = {
    "ERBB2": {
        "cancer_type": "breast",
        "expected_cn": "amplification",
        "description": "HER2 amplification in breast cancer",
        "therapeutic_relevance": "Trastuzumab (Herceptin) target",
    },
    "MYC": {
        "cancer_type": "multiple",
        "expected_cn": "amplification",
        "description": "MYC amplification in multiple cancers",
        "therapeutic_relevance": "Poor prognosis marker",
    },
    "CDKN2A": {
        "cancer_type": "melanoma",
        "expected_cn": "deletion",
        "description": "CDKN2A deletion in melanoma",
        "therapeutic_relevance": "CDK4/6 inhibitor sensitivity",
    },
    "PTEN": {
        "cancer_type": "multiple",
        "expected_cn": "deletion",
        "description": "PTEN loss in multiple cancers",
        "therapeutic_relevance": "PI3K/AKT pathway activation",
    },
    "AR": {
        "cancer_type": "prostate",
        "expected_cn": "amplification",
        "description": "AR amplification in prostate cancer",
        "therapeutic_relevance": "Androgen deprivation therapy resistance",
    },
    "EGFR": {
        "cancer_type": "glioblastoma",
        "expected_cn": "amplification",
        "description": "EGFR amplification in GBM",
        "therapeutic_relevance": "EGFR inhibitor target",
    },
    "MDM2": {
        "cancer_type": "sarcoma",
        "expected_cn": "amplification",
        "description": "MDM2 amplification in liposarcoma",
        "therapeutic_relevance": "MDM2 inhibitor target",
    },
    "CCND1": {
        "cancer_type": "breast",
        "expected_cn": "amplification",
        "description": "Cyclin D1 amplification",
        "therapeutic_relevance": "CDK4/6 inhibitor sensitivity",
    },
}


def simulate_depth_profile(gene: str, has_cnv: bool = False) -> np.ndarray:
    """
    Simulate read depth profile for a gene region.
    In production, would calculate from BAM file pileup.
    """
    # Baseline depth (diploid)
    baseline_depth = 100
    noise_level = 10

    if not has_cnv:
        # Normal diploid - add some noise
        depths = np.random.normal(baseline_depth, noise_level, 100)
    else:
        # Simulate CNV
        cnv_info = CANCER_CNVS.get(gene, {})

        if cnv_info.get("expected_cn") == "amplification":
            # Amplification - 3-8 copies
            copy_number = np.random.choice([3, 4, 5, 6, 8])
            expected_depth = baseline_depth * (copy_number / 2.0)
            depths = np.random.normal(expected_depth, noise_level, 100)

        elif cnv_info.get("expected_cn") == "deletion":
            # Deletion - 0 or 1 copy
            copy_number = np.random.choice([0, 1])
            expected_depth = baseline_depth * (copy_number / 2.0)
            depths = np.random.normal(max(expected_depth, 10), noise_level, 100)

        else:
            # Unknown - random CNV
            copy_number = np.random.choice([1, 3, 4])
            expected_depth = baseline_depth * (copy_number / 2.0)
            depths = np.random.normal(expected_depth, noise_level, 100)

    # Ensure positive depths
    depths = np.maximum(depths, 0)

    return depths


def calculate_copy_number(depths: np.ndarray, baseline_depth: float) -> int:
    """Calculate copy number from depth profile"""
    mean_depth = np.mean(depths)

    # Calculate copy number estimate
    cn_estimate = round(2 * mean_depth / baseline_depth) if baseline_depth > 0 else 2

    # Ensure valid copy number
    return max(0, min(cn_estimate, 10))


def segment_cnv_regions(gene_cnvs: Dict[str, Dict]) -> List[Dict]:
    """
    Segment and merge adjacent CNV regions.
    In production, would use CBS or similar algorithm.
    """
    cnvs = []

    # Group by chromosome
    chr_cnvs = defaultdict(list)
    for gene, cnv_data in gene_cnvs.items():
        if cnv_data["copy_number"] != 2:  # Not diploid
            chr_cnvs[cnv_data["chromosome"]].append(cnv_data)

    # Process each chromosome
    for chrom, cnv_list in chr_cnvs.items():
        # Sort by position (simulated)
        cnv_list.sort(key=lambda x: x.get("start", 0))

        # For now, report each CNV separately
        # In production, would merge adjacent CNVs
        cnvs.extend(cnv_list)

    return cnvs


def get_gene_info(gene: str) -> Dict:
    """Get chromosome and position info for gene"""
    # Simulated gene locations
    gene_info = {
        "ERBB2": {"chr": "chr17", "start": 37844000, "end": 37886000},
        "MYC": {"chr": "chr8", "start": 128748000, "end": 128753000},
        "CDKN2A": {"chr": "chr9", "start": 21967000, "end": 21995000},
        "TP53": {"chr": "chr17", "start": 7565000, "end": 7590000},
        "PTEN": {"chr": "chr10", "start": 89623000, "end": 89728000},
        "AR": {"chr": "chrX", "start": 66763000, "end": 66950000},
        "EGFR": {"chr": "chr7", "start": 55086000, "end": 55324000},
        "MDM2": {"chr": "chr12", "start": 69201000, "end": 69239000},
        "CCND1": {"chr": "chr11", "start": 69455000, "end": 69469000},
        "BRCA1": {"chr": "chr17", "start": 41196000, "end": 41277000},
        "BRCA2": {"chr": "chr13", "start": 32889000, "end": 32973000},
    }

    default = {"chr": "chrUnknown", "start": 0, "end": 1000}
    return gene_info.get(gene, default)


def annotate_cnv(cnv: Dict) -> Dict:
    """Add clinical annotations to CNV"""
    gene = cnv.get("gene", "")
    copy_number = cnv.get("copy_number", 2)

    # Check if it's a known cancer CNV
    if gene in CANCER_CNVS:
        cnv_info = CANCER_CNVS[gene]

        # Check if CNV matches expected pattern
        if (cnv_info["expected_cn"] == "amplification" and copy_number >= 4) or (
            cnv_info["expected_cn"] == "deletion" and copy_number <= 1
        ):

            cnv["cancer_relevance"] = cnv_info["description"]
            cnv["therapeutic_relevance"] = cnv_info.get("therapeutic_relevance", "")
            cnv["cancer_type"] = cnv_info["cancer_type"]

    # Assess clinical significance based on copy number
    if copy_number == 0:
        cnv["clinical_significance"] = "pathogenic"
        cnv["type"] = "homozygous_deletion"
    elif copy_number == 1:
        cnv["clinical_significance"] = "likely_pathogenic"
        cnv["type"] = "heterozygous_deletion"
    elif copy_number >= 4:
        cnv["clinical_significance"] = "likely_pathogenic"
        cnv["type"] = "amplification"
    elif copy_number == 3:
        cnv["clinical_significance"] = "uncertain_significance"
        cnv["type"] = "duplication"
    else:
        cnv["clinical_significance"] = "benign"
        cnv["type"] = "normal"

    return cnv


def process(state: Dict) -> Dict:
    """Detect and annotate CNVs"""
    logger.info("Starting CNV detection")
    # Don't set current_node to avoid concurrent updates

    try:
        # Get variants and genes
        variants = state.get("variant_details", state.get("filtered_variants", []))

        # Count variants per gene (proxy for CNV detection)
        gene_variant_counts = defaultdict(int)
        for variant in variants:
            gene = variant.get("gene")
            if gene:
                gene_variant_counts[gene] += 1

        # Simulate CNV detection for genes
        gene_cnvs = {}

        # Check known cancer genes
        cancer_genes = set(CANCER_CNVS.keys())
        genes_to_check = set(gene_variant_counts.keys()) | cancer_genes

        for gene in genes_to_check:
            # Decide if gene has CNV based on variant count and known associations
            variant_count = gene_variant_counts.get(gene, 0)

            # Higher variant count might indicate CNV
            has_cnv = (
                variant_count >= 3
                or (gene in CANCER_CNVS and variant_count >= 1)
                or (gene in CANCER_CNVS and np.random.random() < 0.3)
            )  # 30% chance for known genes

            if has_cnv or gene in ["ERBB2", "MYC"]:  # Always check key oncogenes
                # Simulate depth profile
                depths = simulate_depth_profile(gene, has_cnv)

                # Calculate copy number
                baseline = 100  # Simulated baseline depth
                copy_number = calculate_copy_number(depths, baseline)

                if copy_number != 2:  # Not diploid
                    gene_info = get_gene_info(gene)

                    cnv = {
                        "gene": gene,
                        "chromosome": gene_info["chr"],
                        "start": gene_info["start"],
                        "end": gene_info["end"],
                        "size": gene_info["end"] - gene_info["start"],
                        "copy_number": copy_number,
                        "normal_copy_number": 2,
                        "fold_change": copy_number / 2.0,
                        "depth_ratio": np.mean(depths) / baseline,
                        "genes_affected": [gene],
                        "evidence": "depth_analysis",
                    }

                    gene_cnvs[gene] = cnv

        # Segment and merge CNVs
        cnvs = segment_cnv_regions(gene_cnvs)

        # Annotate CNVs
        annotated_cnvs = []
        for cnv in cnvs:
            annotated_cnv = annotate_cnv(cnv)
            annotated_cnvs.append(annotated_cnv)

        # Sort by clinical significance
        significance_order = {"pathogenic": 0, "likely_pathogenic": 1, "uncertain_significance": 2, "benign": 3}
        annotated_cnvs.sort(
            key=lambda x: significance_order.get(x.get("clinical_significance", "uncertain_significance"), 2)
        )

        # Summary statistics
        cnv_summary = {
            "total_cnvs": len(annotated_cnvs),
            "amplifications": sum(1 for cnv in annotated_cnvs if cnv.get("type") == "amplification"),
            "deletions": sum(
                1 for cnv in annotated_cnvs if cnv.get("type") in ["heterozygous_deletion", "homozygous_deletion"]
            ),
            "pathogenic_cnvs": sum(
                1 for cnv in annotated_cnvs if cnv.get("clinical_significance") in ["pathogenic", "likely_pathogenic"]
            ),
            "cancer_associated": sum(1 for cnv in annotated_cnvs if "cancer_relevance" in cnv),
            "therapeutic_targets": [cnv["gene"] for cnv in annotated_cnvs if cnv.get("therapeutic_relevance")],
        }

        # Add to completed nodes
        completed = state.get("completed_nodes", [])
        if "cnv_detector" not in completed:
            completed.append("cnv_detector")

        logger.info(f"Detected {len(annotated_cnvs)} CNVs")

        # Return only the keys this node updates
        return {"copy_number_variants": annotated_cnvs, "cnv_summary": cnv_summary, "completed_nodes": completed}

    except Exception as e:
        logger.error(f"Error in CNV detection: {str(e)}")
        errors = state.get("errors", []) + [
            {"node": "cnv_detector", "error": str(e), "timestamp": datetime.now().isoformat()}
        ]
        # Return only the keys this node updates
        return {"copy_number_variants": [], "errors": errors}
