"""
Preprocessing node for alignment.
Converts FASTQ to aligned BAM using BWA-MEM2 or BWA.
"""

import os
import subprocess
import logging
from datetime import datetime
from typing import Dict, Any
import tempfile

logger = logging.getLogger(__name__)

# Configuration
REFERENCE_GENOME = os.path.join(os.path.dirname(__file__), "..", "test_reference", "test_genome.fa")
# For production, use real reference: "/path/to/hg38.fa"


def check_bwa_installed():
    """Check if BWA or BWA-MEM2 is installed."""
    try:
        # Try BWA-MEM2 first (faster)
        result = subprocess.run(["bwa-mem2", "version"], capture_output=True, text=True)
        if result.returncode == 0:
            return "bwa-mem2"
    except FileNotFoundError:
        pass

    try:
        # Fall back to regular BWA
        result = subprocess.run(["bwa"], capture_output=True, text=True)
        if result.returncode == 1:  # BWA returns 1 when run without args
            return "bwa"
    except FileNotFoundError:
        pass

    return None


def align_fastq(fastq_path: str, output_dir: str = None) -> tuple[str, Dict[str, Any]]:
    """
    Align FASTQ file to reference genome using BWA.

    Args:
        fastq_path: Path to input FASTQ file
        output_dir: Directory for output files (uses temp if None)

    Returns:
        Tuple of (aligned BAM path, alignment stats)
    """
    bwa_cmd = check_bwa_installed()
    if not bwa_cmd:
        raise RuntimeError("Neither BWA nor BWA-MEM2 is installed")

    # Check reference genome exists
    if not os.path.exists(REFERENCE_GENOME):
        raise FileNotFoundError(f"Reference genome not found: {REFERENCE_GENOME}")

    # Check reference index exists
    if not os.path.exists(f"{REFERENCE_GENOME}.bwt"):
        raise FileNotFoundError(f"Reference genome not indexed. Run: bwa index {REFERENCE_GENOME}")

    # Prepare output directory
    if output_dir is None:
        output_dir = tempfile.gettempdir()

    # Generate output filename
    base_name = os.path.basename(fastq_path).replace(".fastq.gz", "").replace(".fastq", "")
    sam_path = os.path.join(output_dir, f"{base_name}_aligned.sam")
    bam_path = os.path.join(output_dir, f"{base_name}_aligned.bam")
    sorted_bam_path = os.path.join(output_dir, f"{base_name}_aligned_sorted.bam")

    try:
        # Run BWA alignment
        logger.info(f"Running {bwa_cmd} alignment...")
        align_cmd = [bwa_cmd, "mem", "-t", "4", REFERENCE_GENOME, fastq_path]

        with open(sam_path, "w") as sam_file:
            process = subprocess.Popen(align_cmd, stdout=sam_file, stderr=subprocess.PIPE)
            _, stderr = process.communicate()

            if process.returncode != 0:
                raise RuntimeError(f"BWA alignment failed: {stderr.decode()}")

        # Convert SAM to BAM
        logger.info("Converting SAM to BAM...")
        subprocess.run(["samtools", "view", "-bS", sam_path, "-o", bam_path], check=True)

        # Sort BAM
        logger.info("Sorting BAM file...")
        subprocess.run(["samtools", "sort", bam_path, "-o", sorted_bam_path], check=True)

        # Index BAM
        logger.info("Indexing BAM file...")
        subprocess.run(["samtools", "index", sorted_bam_path], check=True)

        # Get alignment statistics
        logger.info("Calculating alignment statistics...")
        stats_output = subprocess.run(
            ["samtools", "flagstat", sorted_bam_path], capture_output=True, text=True, check=True
        ).stdout

        # Parse statistics
        stats = parse_samtools_flagstat(stats_output)

        # Clean up intermediate files
        os.remove(sam_path)
        os.remove(bam_path)

        return sorted_bam_path, stats

    except subprocess.CalledProcessError as e:
        logger.error(f"Alignment pipeline failed: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during alignment: {str(e)}")
        raise


def parse_samtools_flagstat(flagstat_output: str) -> Dict[str, Any]:
    """Parse samtools flagstat output into statistics dictionary."""
    stats = {}
    lines = flagstat_output.strip().split("\n")

    for line in lines:
        if "in total" in line:
            stats["total_reads"] = int(line.split()[0])
        elif "mapped (" in line and "primary mapped" not in line:
            stats["mapped_reads"] = int(line.split()[0])
            # Extract percentage
            pct = line.split("(")[1].split("%")[0]
            stats["mapping_rate"] = float(pct) / 100
        elif "properly paired" in line:
            stats["properly_paired"] = int(line.split()[0])

    # Calculate additional stats
    stats["unmapped_reads"] = stats.get("total_reads", 0) - stats.get("mapped_reads", 0)
    stats["reference_genome"] = os.path.basename(REFERENCE_GENOME).replace(".fa", "")

    return stats


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Preprocess the input file.
    - For FASTQ: Run BWA alignment to create BAM
    - For BAM: Validate and pass through
    - For VCF: Load variants directly
    - For MAF: Pass to MAF parser

    Updates state with:
    - aligned_bam_path: Path to aligned BAM file (for FASTQ/BAM inputs)
    - raw_variants: Loaded variants (for VCF inputs)
    - alignment_stats: Alignment statistics
    """
    logger.info("Starting preprocessing node")
    state["current_node"] = "preprocess"

    try:
        file_type = state["file_type"]
        file_path = state["file_path"]

        # Check if input is MAF - pass to MAF parser
        if file_type == "maf":
            logger.info("Input is MAF, passing to MAF parser")
            # Import and call MAF parser
            from . import maf_parser

            # MAF parser will handle its own state updates
            return maf_parser.process(state)

        # Check if input is VCF - load variants directly
        if file_type == "vcf":
            logger.info("Input is VCF, loading variants directly")
            import vcf as pyvcf

            variants = []
            with open(file_path, "r") as vcf_file:
                vcf_reader = pyvcf.Reader(vcf_file)
                for record in vcf_reader:
                    variant_data = {
                        "chrom": record.CHROM,
                        "pos": record.POS,
                        "re": record.REF,
                        "alt": str(record.ALT[0]) if record.ALT else ".",
                        "qual": record.QUAL or 0,
                        "quality": record.QUAL or 0,  # Map to field name expected by QC filter
                        "filter": record.FILTER or [],
                        "info": dict(record.INFO) if record.INFO else {},
                        "variant_id": f"{record.CHROM}:{record.POS}:{record.REF}>{record.ALT[0] if record.ALT else '.'}",
                    }

                    # Add sample-specific data if available
                    if record.samples:
                        sample = record.samples[0]
                        if hasattr(sample.data, "DP"):
                            variant_data["depth"] = sample.data.DP
                        if hasattr(sample.data, 'AF'):
                            # Handle AF as list (take first value) or single value
                            af_value = sample.data.AF
                            if isinstance(af_value, list):
                                variant_data["allele_freq"] = af_value[0] if af_value else 0.0
                            else:
                                variant_data["allele_freq"] = af_value or 0.0

                    variants.append(variant_data)

            logger.info(f"Loaded {len(variants)} variants from VCF")

            # Return only the fields this node updates
            # This is important for LangGraph parallel execution
            return {"raw_variants": variants, "variant_count": len(variants)}

        # For FASTQ files, run alignment
        if file_type == "fastq":
            logger.info("Aligning FASTQ file with BWA")

            # Use the same directory as input file for output
            output_dir = os.path.dirname(os.path.abspath(file_path))

            # Perform alignment
            aligned_bam_path, alignment_stats = align_fastq(file_path, output_dir)

            state["aligned_bam_path"] = aligned_bam_path

            # Update metadata with alignment stats
            state["file_metadata"]["alignment_stats"] = alignment_stats

            logger.info(
                f"Alignment complete: {alignment_stats['mapped_reads']}/{alignment_stats['total_reads']} reads mapped ({alignment_stats['mapping_rate']*100:.1f}%)"
            )

        elif file_type == "bam":
            logger.info("BAM file detected, validating alignment")

            # For BAM files, validate and pass through
            # In production, should check:
            # 1. BAM has proper header with reference info
            # 2. File is sorted and indexed
            # 3. Read groups are properly formatted

            # Basic validation
            try:
                result = subprocess.run(["samtools", "quickcheck", file_path], capture_output=True)
                if result.returncode != 0:
                    raise ValueError("BAM file validation failed")
            except subprocess.CalledProcessError:
                raise ValueError("Invalid BAM file")

            state["aligned_bam_path"] = file_path

        else:
            raise ValueError(f"Unexpected file type in preprocess: {file_type}")

        state["completed_nodes"].append("preprocess")
        logger.info(f"Preprocessing complete. BAM path: {state['aligned_bam_path']}")

        # For FASTQ/BAM files, return only the fields we updated
        return {"aligned_bam_path": state["aligned_bam_path"], "file_metadata": state["file_metadata"]}

    except Exception as e:
        logger.error(f"Preprocessing failed: {str(e)}")
        state["errors"].append({"node": "preprocess", "error": str(e), "timestamp": datetime.now()})
        state["pipeline_status"] = "failed"

        # Return error state
        return {"errors": state["errors"], "pipeline_status": state["pipeline_status"]}

    # This should never be reached
    return {}
