"""
FASTQ File Processor

This module handles processing of FASTQ files for sequencing read analysis.
"""

import gzip
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
from loguru import logger

try:
    from Bio import SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    logger.warning("BioPython not available, using basic FASTQ parsing")

from ..models.base import GenomicData


class FASTQProcessor:
    """
    Processor for FASTQ files
    """
    
    def __init__(self):
        """Initialize FASTQ processor"""
        self.supported_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
        self.sample_id = None
        
        logger.info("🧬 FASTQ Processor initialized")
        if not HAS_BIOPYTHON:
            logger.warning("⚠️ BioPython not available - using basic parsing")
    
    def can_process(self, file_path: Path) -> bool:
        """Check if file can be processed by this processor"""
        return any(str(file_path).endswith(ext) for ext in self.supported_extensions)
    
    def process(self, file_path: Path, sample_id: Optional[str] = None) -> GenomicData:
        """
        Process FASTQ file and extract read statistics
        
        Args:
            file_path: Path to FASTQ file
            sample_id: Optional sample identifier
            
        Returns:
            GenomicData object containing read statistics
        """
        logger.info(f"📄 Processing FASTQ file: {file_path}")
        
        if not file_path.exists():
            raise FileNotFoundError(f"FASTQ file not found: {file_path}")
        
        # Determine sample ID
        self.sample_id = sample_id or file_path.stem.replace('.fastq', '').replace('.fq', '')
        
        # Process FASTQ file
        read_stats = self._process_fastq(file_path)
        
        # Create metadata
        metadata = {
            'processor': 'FASTQProcessor',
            'file_size': file_path.stat().st_size,
            'sample_id': self.sample_id,
            'read_stats': read_stats
        }
        
        logger.info(f"✅ FASTQ processing complete: {read_stats['total_reads']} reads analyzed")
        
        return GenomicData(
            variants=[],  # FASTQ files don't contain variants directly
            metadata=metadata,
            file_path=str(file_path),
            file_type='fastq',
            sample_id=self.sample_id
        )
    
    def _process_fastq(self, file_path: Path) -> Dict[str, Any]:
        """Process FASTQ file and extract read statistics"""
        logger.debug("🔧 Processing FASTQ read data")
        
        try:
            if HAS_BIOPYTHON:
                return self._process_with_biopython(file_path)
            else:
                return self._process_basic(file_path)
                
        except Exception as e:
            logger.error(f"❌ Error processing FASTQ file: {str(e)}")
            raise
    
    def _process_with_biopython(self, file_path: Path) -> Dict[str, Any]:
        """Process FASTQ using BioPython"""
        logger.debug("🔧 Using BioPython for FASTQ processing")
        
        stats = {
            'total_reads': 0,
            'total_bases': 0,
            'avg_read_length': 0,
            'min_read_length': float('inf'),
            'max_read_length': 0,
            'avg_quality': 0,
            'quality_distribution': {},
            'base_composition': {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
        }
        
        # Open file (handle gzip compression)
        if str(file_path).endswith('.gz'):
            file_handle = gzip.open(file_path, 'rt')
        else:
            file_handle = open(file_path, 'r')
        
        try:
            total_quality = 0
            read_count = 0
            
            # Sample first 10k reads for statistics
            for i, record in enumerate(SeqIO.parse(file_handle, "fastq")):
                if i >= 10000:  # Sample limit
                    break
                
                read_count += 1
                read_length = len(record.seq)
                read_quality = record.letter_annotations.get("phred_quality", [])
                
                # Update basic stats
                stats['total_reads'] += 1
                stats['total_bases'] += read_length
                stats['min_read_length'] = min(stats['min_read_length'], read_length)
                stats['max_read_length'] = max(stats['max_read_length'], read_length)
                
                # Quality statistics
                if read_quality:
                    avg_read_quality = sum(read_quality) / len(read_quality)
                    total_quality += avg_read_quality
                
                # Base composition
                sequence = str(record.seq).upper()
                for base in sequence:
                    if base in stats['base_composition']:
                        stats['base_composition'][base] += 1
                    else:
                        stats['base_composition']['N'] += 1
            
            # Calculate averages
            if read_count > 0:
                stats['avg_read_length'] = stats['total_bases'] / read_count
                stats['avg_quality'] = total_quality / read_count
            
            # Fix infinite min_read_length
            if stats['min_read_length'] == float('inf'):
                stats['min_read_length'] = 0
            
        finally:
            file_handle.close()
        
        logger.debug(f"📊 FASTQ statistics: {stats}")
        return stats
    
    def _process_basic(self, file_path: Path) -> Dict[str, Any]:
        """Process FASTQ using basic text parsing"""
        logger.debug("🔧 Using basic text parsing for FASTQ processing")
        
        stats = {
            'total_reads': 0,
            'total_bases': 0,
            'avg_read_length': 0,
            'min_read_length': float('inf'),
            'max_read_length': 0,
            'base_composition': {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
        }
        
        # Open file (handle gzip compression)
        if str(file_path).endswith('.gz'):
            file_handle = gzip.open(file_path, 'rt')
        else:
            file_handle = open(file_path, 'r')
        
        try:
            line_count = 0
            
            # Sample first 40k lines (10k reads)
            for i, line in enumerate(file_handle):
                if i >= 40000:  # Sample limit (4 lines per read)
                    break
                
                line_count += 1
                
                # Every 4th line starting from 2nd is the sequence
                if line_count % 4 == 2:
                    sequence = line.strip().upper()
                    read_length = len(sequence)
                    
                    # Update statistics
                    stats['total_reads'] += 1
                    stats['total_bases'] += read_length
                    stats['min_read_length'] = min(stats['min_read_length'], read_length)
                    stats['max_read_length'] = max(stats['max_read_length'], read_length)
                    
                    # Base composition
                    for base in sequence:
                        if base in stats['base_composition']:
                            stats['base_composition'][base] += 1
                        else:
                            stats['base_composition']['N'] += 1
            
            # Calculate averages
            if stats['total_reads'] > 0:
                stats['avg_read_length'] = stats['total_bases'] / stats['total_reads']
            
            # Fix infinite min_read_length
            if stats['min_read_length'] == float('inf'):
                stats['min_read_length'] = 0
        
        finally:
            file_handle.close()
        
        logger.debug(f"📊 FASTQ statistics: {stats}")
        return stats
    
    def validate_file(self, file_path: Path) -> Tuple[bool, List[str]]:
        """Validate FASTQ file format"""
        logger.debug(f"🔍 Validating FASTQ file: {file_path}")
        
        errors = []
        
        try:
            if not file_path.exists():
                errors.append(f"File does not exist: {file_path}")
                return False, errors
            
            # Check file extension
            if not self.can_process(file_path):
                errors.append(f"Invalid file extension: {file_path.suffix}")
            
            # Open file and check basic FASTQ format
            if str(file_path).endswith('.gz'):
                file_handle = gzip.open(file_path, 'rt')
            else:
                file_handle = open(file_path, 'r')
            
            try:
                # Check first few records
                line_count = 0
                for i, line in enumerate(file_handle):
                    if i >= 20:  # Check first 5 records (4 lines each)
                        break
                    
                    line_count += 1
                    line = line.strip()
                    
                    # Check FASTQ format
                    if line_count % 4 == 1:  # Header line
                        if not line.startswith('@'):
                            errors.append(f"Invalid FASTQ header on line {line_count}")
                    elif line_count % 4 == 3:  # Quality header line
                        if not line.startswith('+'):
                            errors.append(f"Invalid FASTQ quality header on line {line_count}")
                
                # Check if we have complete records
                if line_count % 4 != 0:
                    errors.append("Incomplete FASTQ record at end of file")
                
            finally:
                file_handle.close()
                
        except Exception as e:
            errors.append(f"Error reading FASTQ file: {str(e)}")
        
        is_valid = len(errors) == 0
        
        if is_valid:
            logger.debug("✅ FASTQ file validation passed")
        else:
            logger.warning(f"⚠️ FASTQ file validation failed: {errors}")
        
        return is_valid, errors 