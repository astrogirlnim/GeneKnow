"""
VCF File Processor

This module handles processing of VCF (Variant Call Format) files
for genomic variant analysis.
"""

import os
import re
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
import pandas as pd
from loguru import logger

try:
    import cyvcf2
    HAS_CYVCF2 = True
except ImportError:
    HAS_CYVCF2 = False
    logger.warning("cyvcf2 not available, using basic VCF parsing")

from ..models.base import GenomicData


class VCFProcessor:
    """
    Processor for VCF (Variant Call Format) files
    """
    
    def __init__(self):
        """Initialize VCF processor"""
        self.supported_extensions = ['.vcf', '.vcf.gz']
        self.sample_id = None
        
        logger.info("🧬 VCF Processor initialized")
        if not HAS_CYVCF2:
            logger.warning("⚠️ cyvcf2 not available - using basic parsing")
    
    def can_process(self, file_path: Path) -> bool:
        """
        Check if file can be processed by this processor
        
        Args:
            file_path: Path to the file
            
        Returns:
            True if file can be processed
        """
        return file_path.suffix.lower() in self.supported_extensions
    
    def process(self, file_path: Path, sample_id: Optional[str] = None) -> GenomicData:
        """
        Process VCF file and extract variant information
        
        Args:
            file_path: Path to VCF file
            sample_id: Optional sample identifier
            
        Returns:
            GenomicData object containing parsed variants
        """
        logger.info(f"📄 Processing VCF file: {file_path}")
        
        if not file_path.exists():
            raise FileNotFoundError(f"VCF file not found: {file_path}")
        
        # Determine sample ID
        if sample_id:
            self.sample_id = sample_id
        else:
            self.sample_id = file_path.stem.replace('.vcf', '')
        
        # Process based on available libraries
        if HAS_CYVCF2:
            variants = self._process_with_cyvcf2(file_path)
        else:
            variants = self._process_basic(file_path)
        
        # Create metadata
        metadata = {
            'processor': 'VCFProcessor',
            'file_size': file_path.stat().st_size,
            'variant_count': len(variants),
            'sample_id': self.sample_id
        }
        
        logger.info(f"✅ VCF processing complete: {len(variants)} variants found")
        
        return GenomicData(
            variants=variants,
            metadata=metadata,
            file_path=str(file_path),
            file_type='vcf',
            sample_id=self.sample_id
        )
    
    def _process_with_cyvcf2(self, file_path: Path) -> List[Dict[str, Any]]:
        """
        Process VCF file using cyvcf2 library
        
        Args:
            file_path: Path to VCF file
            
        Returns:
            List of variant dictionaries
        """
        logger.debug("🔧 Using cyvcf2 for VCF processing")
        
        variants = []
        
        try:
            vcf = cyvcf2.VCF(str(file_path))
            
            for i, variant in enumerate(vcf):
                # Extract basic variant information
                variant_dict = {
                    'chrom': variant.CHROM,
                    'pos': variant.POS,
                    'id': variant.ID or f"var_{i}",
                    'ref': variant.REF,
                    'alt': variant.ALT,
                    'qual': variant.QUAL,
                    'filter': variant.FILTER,
                    'info': dict(variant.INFO)
                }
                
                # Extract genotype information
                if variant.genotypes:
                    genotypes = []
                    for sample_gt in variant.genotypes:
                        genotypes.append({
                            'gt': sample_gt[:-1],  # Remove phase info
                            'phased': sample_gt[-1]
                        })
                    variant_dict['genotypes'] = genotypes
                
                # Extract specific annotations
                variant_dict.update(self._extract_annotations(variant_dict))
                
                variants.append(variant_dict)
                
                # Log progress for large files
                if i % 10000 == 0 and i > 0:
                    logger.debug(f"📊 Processed {i} variants...")
            
            vcf.close()
            
        except Exception as e:
            logger.error(f"❌ Error processing VCF with cyvcf2: {str(e)}")
            raise
        
        return variants
    
    def _process_basic(self, file_path: Path) -> List[Dict[str, Any]]:
        """
        Process VCF file using basic text parsing
        
        Args:
            file_path: Path to VCF file
            
        Returns:
            List of variant dictionaries
        """
        logger.debug("🔧 Using basic text parsing for VCF processing")
        
        variants = []
        
        try:
            with open(file_path, 'r') as f:
                for line_num, line in enumerate(f):
                    line = line.strip()
                    
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    # Parse variant line
                    fields = line.split('\t')
                    if len(fields) < 8:
                        continue
                    
                    # Extract basic fields
                    variant_dict = {
                        'chrom': fields[0],
                        'pos': int(fields[1]),
                        'id': fields[2] if fields[2] != '.' else f"var_{line_num}",
                        'ref': fields[3],
                        'alt': fields[4].split(','),
                        'qual': float(fields[5]) if fields[5] != '.' else None,
                        'filter': fields[6],
                        'info': self._parse_info_field(fields[7])
                    }
                    
                    # Extract genotype information if available
                    if len(fields) >= 10:
                        format_fields = fields[8].split(':')
                        genotype_data = fields[9].split(':')
                        
                        genotype_dict = {}
                        for i, fmt in enumerate(format_fields):
                            if i < len(genotype_data):
                                genotype_dict[fmt] = genotype_data[i]
                        
                        variant_dict['genotypes'] = [genotype_dict]
                    
                    # Extract annotations
                    variant_dict.update(self._extract_annotations(variant_dict))
                    
                    variants.append(variant_dict)
                    
                    # Log progress for large files
                    if line_num % 10000 == 0 and line_num > 0:
                        logger.debug(f"📊 Processed {line_num} lines...")
        
        except Exception as e:
            logger.error(f"❌ Error processing VCF with basic parsing: {str(e)}")
            raise
        
        return variants
    
    def _parse_info_field(self, info_str: str) -> Dict[str, Any]:
        """
        Parse INFO field from VCF file
        
        Args:
            info_str: INFO field string
            
        Returns:
            Dictionary of INFO field values
        """
        info_dict = {}
        
        if info_str == '.':
            return info_dict
        
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                # Try to convert to appropriate type
                try:
                    if ',' in value:
                        # Multiple values
                        info_dict[key] = value.split(',')
                    elif value.replace('.', '').isdigit():
                        # Numeric value
                        info_dict[key] = float(value) if '.' in value else int(value)
                    else:
                        # String value
                        info_dict[key] = value
                except ValueError:
                    info_dict[key] = value
            else:
                # Flag field
                info_dict[item] = True
        
        return info_dict
    
    def _extract_annotations(self, variant_dict: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract and standardize variant annotations
        
        Args:
            variant_dict: Basic variant information
            
        Returns:
            Dictionary of additional annotations
        """
        annotations = {}
        info = variant_dict.get('info', {})
        
        # Extract gene information
        gene_fields = ['GENE', 'SYMBOL', 'Gene', 'gene']
        for field in gene_fields:
            if field in info:
                annotations['gene'] = info[field]
                break
        
        # Extract clinical significance
        clinical_fields = ['CLNSIG', 'clinical_significance', 'Clinical_significance']
        for field in clinical_fields:
            if field in info:
                annotations['clinical_significance'] = info[field]
                break
        
        # Extract impact/consequence
        impact_fields = ['IMPACT', 'Consequence', 'EFFECT', 'effect']
        for field in impact_fields:
            if field in info:
                annotations['impact'] = info[field]
                break
        
        # Extract allele frequency
        af_fields = ['AF', 'allele_frequency', 'gnomAD_AF']
        for field in af_fields:
            if field in info:
                annotations['allele_frequency'] = info[field]
                break
        
        # Extract HGVS notation
        hgvs_fields = ['HGVS_C', 'HGVSc', 'hgvs_c']
        for field in hgvs_fields:
            if field in info:
                annotations['hgvs_c'] = info[field]
                break
        
        # Extract quality metrics
        if 'DP' in info:
            annotations['depth'] = info['DP']
        
        if 'GQ' in info:
            annotations['genotype_quality'] = info['GQ']
        
        # Add variant type based on REF/ALT
        ref = variant_dict.get('ref', '')
        alt = variant_dict.get('alt', [])
        
        if isinstance(alt, list) and len(alt) > 0:
            alt_allele = alt[0]
            
            if len(ref) == 1 and len(alt_allele) == 1:
                annotations['variant_type'] = 'SNV'
            elif len(ref) > len(alt_allele):
                annotations['variant_type'] = 'deletion'
            elif len(ref) < len(alt_allele):
                annotations['variant_type'] = 'insertion'
            else:
                annotations['variant_type'] = 'complex'
        
        return annotations
    
    def validate_file(self, file_path: Path) -> Tuple[bool, List[str]]:
        """
        Validate VCF file format
        
        Args:
            file_path: Path to VCF file
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        logger.debug(f"🔍 Validating VCF file: {file_path}")
        
        errors = []
        
        try:
            if not file_path.exists():
                errors.append(f"File does not exist: {file_path}")
                return False, errors
            
            # Check file extension
            if not self.can_process(file_path):
                errors.append(f"Invalid file extension: {file_path.suffix}")
            
            # Check if file is readable
            with open(file_path, 'r') as f:
                # Check for VCF header
                first_line = f.readline().strip()
                if not first_line.startswith('##fileformat=VCF'):
                    errors.append("Missing or invalid VCF header")
                
                # Check for column header
                found_header = False
                for _ in range(100):  # Check first 100 lines
                    line = f.readline().strip()
                    if not line:
                        break
                    if line.startswith('#CHROM'):
                        found_header = True
                        break
                
                if not found_header:
                    errors.append("Missing column header line")
        
        except Exception as e:
            errors.append(f"Error reading file: {str(e)}")
        
        is_valid = len(errors) == 0
        
        if is_valid:
            logger.debug("✅ VCF file validation passed")
        else:
            logger.warning(f"⚠️ VCF file validation failed: {errors}")
        
        return is_valid, errors 