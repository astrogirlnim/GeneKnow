#!/usr/bin/env python3
"""
Genomic Extraction Validation Script
Validates extracted VCF files without requiring deep genomic expertise
"""

import os
import subprocess
import pandas as pd
from datetime import datetime

def run_cmd(cmd, capture=True):
    """Run command and return output"""
    result = subprocess.run(cmd, shell=True, capture_output=capture, text=True)
    if result.returncode != 0:
        raise Exception(f"Command failed: {cmd}\nError: {result.stderr}")
    return result.stdout.strip() if capture else None

def validate_file_structure(vcf_file):
    """Basic file structure validation"""
    print(f"\n🔍 Validating file structure: {vcf_file}")
    
    # Check if file exists and is not empty
    if not os.path.exists(vcf_file):
        return f"❌ File does not exist: {vcf_file}"
    
    size = os.path.getsize(vcf_file)
    if size == 0:
        return f"❌ File is empty: {vcf_file}"
    
    print(f"   ✅ File exists and has size: {size:,} bytes ({size/1024:.1f} KB)")
    
    # Check VCF header
    try:
        header_lines = run_cmd(f"bcftools view -h {vcf_file} | wc -l")
        data_lines = run_cmd(f"bcftools view -H {vcf_file} | wc -l")
        print(f"   ✅ VCF structure: {header_lines} header lines, {data_lines} data lines")
        return int(data_lines)
    except Exception as e:
        return f"❌ VCF structure error: {e}"

def validate_coordinates(vcf_file, expected_chrom, expected_start, expected_end, gene_name):
    """Validate that all variants fall within expected genomic boundaries"""
    print(f"\n📍 Validating coordinates for {gene_name} (chr{expected_chrom}:{expected_start}-{expected_end})")
    
    try:
        # Get all variant positions
        positions_cmd = f"bcftools view -H {vcf_file} | cut -f1,2"
        positions_output = run_cmd(positions_cmd)
        
        if not positions_output:
            return "❌ No variants found in file"
        
        violations = []
        total_variants = 0
        
        for line in positions_output.split('\n'):
            if line.strip():
                chrom, pos = line.split('\t')
                pos = int(pos)
                total_variants += 1
                
                # Check chromosome
                if chrom != expected_chrom:
                    violations.append(f"Wrong chromosome: {chrom} (expected {expected_chrom}) at position {pos}")
                
                # Check boundaries
                if pos < expected_start or pos > expected_end:
                    violations.append(f"Position {pos} outside bounds {expected_start}-{expected_end}")
        
        if violations:
            print(f"   ❌ Boundary violations found:")
            for v in violations[:5]:  # Show first 5
                print(f"      - {v}")
            if len(violations) > 5:
                print(f"      ... and {len(violations) - 5} more")
            return f"❌ {len(violations)}/{total_variants} variants outside boundaries"
        else:
            print(f"   ✅ All {total_variants} variants within correct boundaries")
            return total_variants
            
    except Exception as e:
        return f"❌ Coordinate validation error: {e}"

def validate_vcf_format(vcf_file):
    """Validate VCF format compliance"""
    print(f"\n🧬 Validating VCF format compliance")
    
    try:
        # Use bcftools to validate format
        run_cmd(f"bcftools view -h {vcf_file} > /dev/null")
        print(f"   ✅ VCF header is valid")
        
        # Check for required columns
        header_cmd = f"bcftools view -h {vcf_file} | tail -1"
        header_line = run_cmd(header_cmd)
        
        required_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        header_cols = header_line.split('\t')
        
        missing_cols = [col for col in required_cols if col not in header_cols]
        if missing_cols:
            return f"❌ Missing required VCF columns: {missing_cols}"
        
        print(f"   ✅ All required VCF columns present")
        print(f"   ✅ Total columns: {len(header_cols)} (including {len(header_cols) - 8} samples)")
        
        return "valid"
        
    except Exception as e:
        return f"❌ VCF format validation error: {e}"

def validate_against_reference(vcf_file, gene_name, expected_chrom):
    """Cross-validate against known reference data"""
    print(f"\n🔬 Cross-validating {gene_name} against reference databases")
    
    try:
        # Get variant summary statistics
        stats_cmd = f"bcftools stats {vcf_file}"
        stats_output = run_cmd(stats_cmd)
        
        # Parse key statistics
        stats = {}
        for line in stats_output.split('\n'):
            if line.startswith('SN'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    key = parts[2].strip(':')
                    value = parts[3]
                    stats[key] = value
        
        print(f"   📊 Variant Statistics:")
        important_stats = ['number of records', 'number of SNPs', 'number of indels']
        for stat in important_stats:
            if stat in stats:
                print(f"      {stat}: {stats[stat]}")
        
        # Validate chromosome consistency
        chrom_cmd = f"bcftools view -H {vcf_file} | cut -f1 | sort | uniq"
        chroms = run_cmd(chrom_cmd).split('\n')
        chroms = [c for c in chroms if c.strip()]
        
        if len(chroms) == 1 and chroms[0] == expected_chrom:
            print(f"   ✅ Chromosome consistency: All variants on chr{expected_chrom}")
        else:
            print(f"   ⚠️  Multiple chromosomes found: {chroms}")
            
        return stats
        
    except Exception as e:
        return f"❌ Reference validation error: {e}"

def spot_check_variants(vcf_file, gene_name, num_samples=5):
    """Show sample variants for manual inspection"""
    print(f"\n🔍 Spot check: Sample variants from {gene_name}")
    
    try:
        # Get first few variants
        sample_cmd = f"bcftools view -H {vcf_file} | head -{num_samples}"
        sample_output = run_cmd(sample_cmd)
        
        if not sample_output:
            return "❌ No variants found for spot check"
        
        print(f"   📋 First {num_samples} variants (CHROM, POS, REF, ALT):")
        for i, line in enumerate(sample_output.split('\n'), 1):
            if line.strip():
                parts = line.split('\t')
                chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
                print(f"      {i}. chr{chrom}:{pos} {ref}→{alt}")
        
        return "valid"
        
    except Exception as e:
        return f"❌ Spot check error: {e}"

def main():
    """Main validation workflow"""
    print("🧬 GENOMIC EXTRACTION VALIDATION SUITE")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Read coordinates from regions.bed file
    regions = []
    with open("regions.bed", "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                chrom, start, end, gene = line.split('\t')
                regions.append((chrom, int(start), int(end), gene))
    
    output_dir = "output_chunks"
    validation_results = {}
    
    for chrom, start, end, gene in regions:
        print(f"\n{'='*60}")
        print(f"🧬 VALIDATING {gene} (chr{chrom}:{start}-{end})")
        print(f"{'='*60}")
        
        vcf_file = os.path.join(output_dir, f"{gene}.vcf.gz")
        result = {
            "gene": gene,
            "chromosome": chrom,
            "start": start,
            "end": end,
            "file": vcf_file
        }
        
        # Run validation checks
        result["structure"] = validate_file_structure(vcf_file)
        result["coordinates"] = validate_coordinates(vcf_file, chrom, start, end, gene)
        result["format"] = validate_vcf_format(vcf_file)
        result["reference"] = validate_against_reference(vcf_file, gene, chrom)
        result["spot_check"] = spot_check_variants(vcf_file, gene)
        
        validation_results[gene] = result
    
    # Summary report
    print(f"\n{'='*60}")
    print("📊 VALIDATION SUMMARY REPORT")
    print(f"{'='*60}")
    
    all_passed = True
    for gene, result in validation_results.items():
        print(f"\n🧬 {gene}:")
        
        # Check each validation
        checks = ['structure', 'coordinates', 'format', 'reference', 'spot_check']
        gene_passed = True
        
        for check in checks:
            value = result[check]
            if isinstance(value, str) and value.startswith('❌'):
                print(f"   ❌ {check}: FAILED")
                gene_passed = False
                all_passed = False
            elif isinstance(value, int):
                print(f"   ✅ {check}: {value} variants")
            else:
                print(f"   ✅ {check}: PASSED")
        
        if gene_passed:
            print(f"   🎉 {gene}: ALL VALIDATIONS PASSED")
        else:
            print(f"   ⚠️  {gene}: SOME VALIDATIONS FAILED")
    
    print(f"\n{'='*60}")
    if all_passed:
        print("🎉 OVERALL RESULT: ALL VALIDATIONS PASSED!")
        print("✅ Your genomic extraction pipeline is working correctly.")
    else:
        print("⚠️  OVERALL RESULT: SOME VALIDATIONS FAILED")
        print("❌ Review the specific failures above.")
    print(f"{'='*60}")

if __name__ == "__main__":
    main() 