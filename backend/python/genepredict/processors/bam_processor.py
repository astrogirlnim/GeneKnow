"""
BAM File Processor

This module handles processing of BAM (Binary Alignment Map) files
for genomic alignment analysis.
"""

from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
from loguru import logger

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False
    logger.warning("pysam not available, BAM processing disabled")

from ..models.base import GenomicData


class BAMProcessor:
    """
    Processor for BAM (Binary Alignment Map) files
    """
    
    def __init__(self):
        """Initialize BAM processor"""
        self.supported_extensions = ['.bam']
        self.sample_id = None
        
        logger.info("🧬 BAM Processor initialized")
        if not HAS_PYSAM:
            logger.warning("⚠️ pysam not available - BAM processing disabled")
    
    def can_process(self, file_path: Path) -> bool:
        """Check if file can be processed by this processor"""
        return file_path.suffix.lower() in self.supported_extensions and HAS_PYSAM
    
    def process(self, file_path: Path, sample_id: Optional[str] = None) -> GenomicData:
        """
        Process BAM file and extract alignment statistics
        
        Args:
            file_path: Path to BAM file
            sample_id: Optional sample identifier
            
        Returns:
            GenomicData object containing alignment statistics
        """
        logger.info(f"📄 Processing BAM file: {file_path}")
        
        if not HAS_PYSAM:
            raise ImportError("pysam is required for BAM processing")
        
        if not file_path.exists():
            raise FileNotFoundError(f"BAM file not found: {file_path}")
        
        # Determine sample ID
        self.sample_id = sample_id or file_path.stem.replace('.bam', '')
        
        # Process BAM file
        alignment_stats = self._process_bam(file_path)
        
        # Create metadata
        metadata = {
            'processor': 'BAMProcessor',
            'file_size': file_path.stat().st_size,
            'sample_id': self.sample_id,
            'alignment_stats': alignment_stats
        }
        
        logger.info(f"✅ BAM processing complete")
        
        return GenomicData(
            variants=[],  # BAM files don't contain variants directly
            metadata=metadata,
            file_path=str(file_path),
            file_type='bam',
            sample_id=self.sample_id
        )
    
    def _process_bam(self, file_path: Path) -> Dict[str, Any]:
        """Process BAM file and extract alignment statistics"""
        logger.debug("🔧 Processing BAM alignment data")
        
        try:
            with pysam.AlignmentFile(str(file_path), "rb") as bamfile:
                # Get basic statistics
                stats = {
                    'total_reads': 0,
                    'mapped_reads': 0,
                    'unmapped_reads': 0,
                    'duplicate_reads': 0,
                    'quality_failed': 0,
                    'paired_reads': 0,
                    'properly_paired': 0
                }
                
                # Sample first 100k reads for statistics
                for i, read in enumerate(bamfile.fetch()):
                    if i >= 100000:  # Sample limit
                        break
                    
                    stats['total_reads'] += 1
                    
                    if read.is_unmapped:
                        stats['unmapped_reads'] += 1
                    else:
                        stats['mapped_reads'] += 1
                    
                    if read.is_duplicate:
                        stats['duplicate_reads'] += 1
                    
                    if read.is_qcfail:
                        stats['quality_failed'] += 1
                    
                    if read.is_paired:
                        stats['paired_reads'] += 1
                        
                        if read.is_proper_pair:
                            stats['properly_paired'] += 1
                
                # Calculate percentages
                total = stats['total_reads']
                if total > 0:
                    stats['mapping_rate'] = stats['mapped_reads'] / total
                    stats['duplication_rate'] = stats['duplicate_reads'] / total
                    stats['quality_fail_rate'] = stats['quality_failed'] / total
                
                logger.debug(f"📊 BAM statistics: {stats}")
                return stats
                
        except Exception as e:
            logger.error(f"❌ Error processing BAM file: {str(e)}")
            raise
    
    def validate_file(self, file_path: Path) -> Tuple[bool, List[str]]:
        """Validate BAM file format"""
        logger.debug(f"🔍 Validating BAM file: {file_path}")
        
        errors = []
        
        if not HAS_PYSAM:
            errors.append("pysam library not available")
            return False, errors
        
        try:
            if not file_path.exists():
                errors.append(f"File does not exist: {file_path}")
                return False, errors
            
            # Try to open BAM file
            with pysam.AlignmentFile(str(file_path), "rb") as bamfile:
                # Check if file has header
                if not bamfile.header:
                    errors.append("BAM file missing header")
                
                # Try to read first record
                try:
                    next(bamfile.fetch())
                except StopIteration:
                    # Empty file is valid
                    pass
                
        except Exception as e:
            errors.append(f"Error reading BAM file: {str(e)}")
        
        is_valid = len(errors) == 0
        
        if is_valid:
            logger.debug("✅ BAM file validation passed")
        else:
            logger.warning(f"⚠️ BAM file validation failed: {errors}")
        
        return is_valid, errors 