"""
File input validation node.
Validates uploaded FASTQ/BAM/VCF/MAF files and extracts metadata.
"""

import os
import gzip
import logging
from datetime import datetime
from typing import Dict, Any
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pysam

logger = logging.getLogger(__name__)


def validate_fastq(file_path: str) -> Dict[str, Any]:
    """
    Validate FASTQ file and extract metadata.

    Args:
        file_path: Path to FASTQ file

    Returns:
        Dictionary with file metadata
    """
    metadata = {
        "format": "FASTQ",
        "compression": "gzip" if file_path.endswith(".gz") else "none",
        "file_size_mb": os.path.getsize(file_path) / (1024 * 1024),
    }

    # Sample the file to get read statistics
    read_count = 0
    read_lengths = []
    quality_scores = []

    try:
        # Open file with appropriate handler
        if file_path.endswith(".gz"):
            handle = gzip.open(file_path, "rt")
        else:
            handle = open(file_path, "r")

        # Use FastqGeneralIterator for efficient parsing
        for title, seq, qual in FastqGeneralIterator(handle):
            read_count += 1
            read_lengths.append(len(seq))

            # Sample first 1000 reads for statistics
            if read_count <= 1000:
                # Calculate average quality score for this read
                qual_scores = [ord(c) - 33 for c in qual]  # Assuming Phred+33
                quality_scores.extend(qual_scores)

            # For large files, estimate from first 10k reads
            if read_count >= 10000:
                break

        handle.close()

        # Calculate statistics
        if read_lengths:
            metadata["read_length"] = int(sum(read_lengths) / len(read_lengths))
            metadata["read_length_min"] = min(read_lengths)
            metadata["read_length_max"] = max(read_lengths)
        else:
            raise ValueError("No reads found in FASTQ file")

        # Estimate total reads if we sampled
        if read_count >= 10000:
            # Estimate based on file size
            sample_size = handle.tell() if hasattr(handle, "tell") else 0
            if sample_size > 0:
                total_size = os.path.getsize(file_path)
                metadata["estimated_read_count"] = int((total_size / sample_size) * read_count)
            else:
                metadata["estimated_read_count"] = read_count * 10  # Rough estimate
        else:
            metadata["estimated_read_count"] = read_count

        # Quality encoding detection
        if quality_scores:
            min_qual = min(quality_scores)
            max_qual = max(quality_scores)
            avg_qual = sum(quality_scores) / len(quality_scores)

            metadata["quality_encoding"] = "Phred+33"  # Modern standard
            metadata["quality_score_min"] = min_qual
            metadata["quality_score_max"] = max_qual
            metadata["quality_score_avg"] = round(avg_qual, 1)

        # Check if paired-end (simple heuristic: /1 or /2 in read names)
        # Note: This is a simplified check
        metadata["paired_end"] = False  # Would need paired file to confirm

    except Exception as e:
        raise ValueError(f"Failed to parse FASTQ file: {str(e)}")

    return metadata


def validate_bam(file_path: str) -> Dict[str, Any]:
    """
    Validate BAM file and extract metadata.

    Args:
        file_path: Path to BAM file

    Returns:
        Dictionary with file metadata
    """
    metadata = {
        "format": "BAM" if file_path.endswith(".bam") else "SAM",
        "file_size_mb": os.path.getsize(file_path) / (1024 * 1024),
    }

    try:
        # Open BAM file
        bamfile = pysam.AlignmentFile(file_path, "rb")

        # Extract header information
        header = bamfile.header
        if "SQ" in header:
            # Get reference information
            ref_sequences = header["SQ"]
            metadata["reference_sequences"] = len(ref_sequences)

            # Try to identify reference genome
            if ref_sequences:
                first_seq = ref_sequences[0]["SN"]
                if first_seq.startswith("chr"):
                    metadata["reference_genome"] = "hg38"  # Assumption based on naming
                else:
                    metadata["reference_genome"] = "unknown"

        # Get alignment statistics
        total_reads = 0
        mapped_reads = 0
        properly_paired = 0

        # Sample first 100k reads for statistics
        for read in bamfile.fetch(until_eof=True):
            total_reads += 1

            if not read.is_unmapped:
                mapped_reads += 1

            if read.is_proper_pair:
                properly_paired += 1

            if total_reads >= 100000:
                break

        # Estimate total if we only sampled
        if total_reads >= 100000:
            # Use index if available for better estimate
            if os.path.exists(file_path + ".bai"):
                try:
                    idx_stats = bamfile.get_index_statistics()
                    total_estimate = sum(stat.total for stat in idx_stats)
                    if total_estimate > 0:
                        metadata["total_reads"] = total_estimate
                        metadata["mapped_reads"] = sum(stat.mapped for stat in idx_stats)
                    else:
                        # Rough estimate
                        metadata["total_reads"] = total_reads * 10
                        metadata["mapped_reads"] = mapped_reads * 10
                except Exception:
                    metadata["total_reads"] = total_reads * 10
                    metadata["mapped_reads"] = mapped_reads * 10
            else:
                metadata["total_reads"] = total_reads * 10
                metadata["mapped_reads"] = mapped_reads * 10
                metadata["has_index"] = False
        else:
            metadata["total_reads"] = total_reads
            metadata["mapped_reads"] = mapped_reads
            metadata["has_index"] = os.path.exists(file_path + ".bai")

        # Calculate mapping rate
        if metadata["total_reads"] > 0:
            metadata["mapping_rate"] = round(metadata["mapped_reads"] / metadata["total_reads"], 3)
        else:
            metadata["mapping_rate"] = 0.0

        # Check for read groups
        if "RG" in header:
            metadata["read_groups"] = len(header["RG"])

        bamfile.close()

    except Exception as e:
        raise ValueError(f"Failed to parse BAM file: {str(e)}")

    return metadata


def validate_maf(file_path: str) -> Dict[str, Any]:
    """
    Validate MAF file and extract metadata.

    Args:
        file_path: Path to MAF file

    Returns:
        Dictionary with file metadata
    """
    metadata = {
        "format": "MAF",
        "compression": "gzip" if file_path.endswith(".gz") else "none",
        "file_size_mb": os.path.getsize(file_path) / (1024 * 1024),
    }

    try:
        # Open file with appropriate handler
        if file_path.endswith(".gz"):
            f = gzip.open(file_path, "rt")
        else:
            f = open(file_path, "r")

        with f:
            # Skip comment lines
            header_line = None
            line_count = 0
            for line in f:
                line_count += 1
                if not line.startswith("#"):
                    header_line = line.strip()
                    break
                if line_count > 100:  # Safety check
                    break

            if not header_line:
                raise ValueError("No header found in MAF file")

            # Check for required MAF columns
            columns = header_line.split("\t")
            required_columns = ["Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]

            missing = [col for col in required_columns if col not in columns]
            if missing:
                raise ValueError(f"Missing required MAF columns: {missing}")

            metadata["column_count"] = len(columns)
            metadata["has_required_columns"] = True

            # Count total lines (variants)
            variant_count = 0
            for line in f:
                if not line.startswith("#") and line.strip():
                    variant_count += 1

            metadata["variant_count"] = variant_count

    except Exception as e:
        raise ValueError(f"Failed to validate MAF file: {str(e)}")

    return metadata


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate input file and extract metadata.

    Updates state with:
    - file_type: 'fastq', 'bam', 'vc', or 'maf'
    - file_metadata: size, format validation, read count estimate
    - errors: any validation failures
    """
    logger.info(f"Processing file input: {state['file_path']}")
    state["current_node"] = "file_input"

    try:
        file_path = state["file_path"]

        # Check file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        # Check file is readable
        if not os.access(file_path, os.R_OK):
            raise PermissionError(f"File is not readable: {file_path}")

        # Determine file type from extension
        file_lower = file_path.lower()
        if file_lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            file_type = "fastq"
        elif file_lower.endswith((".bam", ".sam")):
            file_type = "bam"
        elif file_lower.endswith((".vc", ".vcf.gz")):
            file_type = "vc"
        elif file_lower.endswith((".ma", ".maf.gz")):
            file_type = "ma"
        else:
            raise ValueError(f"Unsupported file type: {file_path}")

        # Validate based on file type
        if file_type == "fastq":
            metadata = validate_fastq(file_path)
        elif file_type == "bam":
            metadata = validate_bam(file_path)
        elif file_type == "vcf":
            # Basic VCF metadata
            metadata = {"format": "VCF", "file_size_mb": os.path.getsize(file_path) / (1024 * 1024)}
        elif file_type == "ma":
            metadata = validate_maf(file_path)
        else:
            raise ValueError(f"Unsupported file type: {file_path}")

        # Update state
        state["file_type"] = file_type
        state["file_metadata"] = metadata
        state["completed_nodes"].append("file_input")

        logger.info(f"File validation successful: {file_type} format detected")
        logger.info(f"File metadata: {metadata}")

    except Exception as e:
        logger.error(f"File validation failed: {str(e)}")
        state["errors"].append({"node": "file_input", "error": str(e), "timestamp": datetime.now()})
        state["pipeline_status"] = "failed"

    return state
