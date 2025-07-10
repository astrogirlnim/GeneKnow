#!/usr/bin/env python3
"""
Quick setup and test script for TCGA Frequency Match implementation.
Creates the TCGA database and tests the integration.
"""

import os
import sys
import logging
from pathlib import Path

# Add current directory to path for imports
sys.path.append(str(Path(__file__).parent))

def setup_tcga_database():
    """Create the TCGA database."""
    print("🔧 Setting up TCGA database...")
    
    try:
        from create_tcga_database import create_tcga_database
        db_path = create_tcga_database()
        print(f"✅ TCGA database created: {db_path}")
        return True
    except Exception as e:
        print(f"❌ Failed to create TCGA database: {e}")
        return False

def test_tcga_mapper():
    """Test the TCGA mapper node directly."""
    print("\n🧪 Testing TCGA mapper...")
    
    try:
        from nodes.tcga_mapper import process
        
        # Create test state with sample variants
        test_state = {
            "filtered_variants": [
                {
                    "chrom": "17",
                    "pos": 43044295,
                    "ref": "G",
                    "alt": "A", 
                    "gene": "BRCA1",
                    "variant_id": "17:43044295:G>A"
                },
                {
                    "chrom": "17", 
                    "pos": 7674221,
                    "ref": "G",
                    "alt": "A",
                    "gene": "TP53",
                    "variant_id": "17:7674221:G>A"
                },
                {
                    "chrom": "12",
                    "pos": 25245350,
                    "ref": "G", 
                    "alt": "A",
                    "gene": "KRAS",
                    "variant_id": "12:25245350:G>A"
                }
            ],
            "file_metadata": {},
            "completed_nodes": [],
            "errors": [],
            "warnings": []
        }
        
        # Run TCGA mapping
        result_state = process(test_state)
        
        # Check results
        if "tcga_matches" in result_state:
            tcga_matches = result_state["tcga_matches"]
            print(f"✅ TCGA mapping successful!")
            print(f"   Cancer types analyzed: {list(tcga_matches.keys())}")
            
            # Count matches
            total_matches = sum(len(matches) for matches in tcga_matches.values())
            print(f"   Total variant matches: {total_matches}")
            
            # Show sample matches
            if total_matches > 0:
                print("\n📊 Sample matches:")
                for cancer_type, matches in tcga_matches.items():
                    if matches:
                        for variant_id, match_data in list(matches.items())[:2]:  # Show first 2
                            enrichment = match_data.get("enrichment_score", 1.0)
                            tumor_freq = match_data.get("tumor_frequency", 0.0)
                            print(f"   {variant_id} in {cancer_type}: {tumor_freq:.1%} frequency, {enrichment:.1f}x enrichment")
            
            return True
        else:
            print("❌ No TCGA matches found in result")
            return False
            
    except Exception as e:
        print(f"❌ TCGA mapper test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_full_pipeline():
    """Test a small variant through the full pipeline."""
    print("\n🚀 Testing full pipeline integration...")
    
    try:
        from graph import run_pipeline
        import tempfile
        
        # Create a minimal test VCF
        vcf_content = """##fileformat=VCFv4.2
##source=TestVCF
##reference=GRCh38
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
17	43044295	.	G	A	60	PASS	DP=30;AF=0.5	GT	0/1
17	7674221	.	G	A	50	PASS	DP=25;AF=0.4	GT	0/1
"""
        
        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(vcf_content)
            test_vcf_path = f.name
        
        try:
            # Run pipeline
            result = run_pipeline(test_vcf_path, {
                "patient_data": {"age": 45, "sex": "F"},
                "language": "en"
            })
            
            # Check results
            if result.get("pipeline_status") != "failed":
                print("✅ Full pipeline test successful!")
                print(f"   Nodes completed: {result.get('completed_nodes', [])}")
                print(f"   Risk scores: {result.get('risk_scores', {})}")
                
                # Check TCGA integration
                tcga_summary = result.get("file_metadata", {}).get("tcga_summary", {})
                if tcga_summary:
                    print(f"   TCGA matches: {tcga_summary.get('variants_matched', 0)}")
                    print(f"   Match rate: {tcga_summary.get('match_rate', 0)*100:.1f}%")
                
                return True
            else:
                print(f"❌ Pipeline failed: {result.get('errors', [])}")
                return False
                
        finally:
            # Clean up temp file
            os.unlink(test_vcf_path)
            
    except Exception as e:
        print(f"❌ Full pipeline test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all setup and tests."""
    print("🧬 TCGA Frequency Match Setup & Test")
    print("=" * 50)
    
    # Configure logging
    logging.basicConfig(level=logging.WARNING)  # Reduce noise during testing
    
    success_count = 0
    total_tests = 3
    
    # 1. Setup database
    if setup_tcga_database():
        success_count += 1
    
    # 2. Test TCGA mapper directly
    if test_tcga_mapper():
        success_count += 1
    
    # 3. Test full pipeline integration
    if test_full_pipeline():
        success_count += 1
    
    # Summary
    print(f"\n📋 Test Summary: {success_count}/{total_tests} passed")
    
    if success_count == total_tests:
        print("🎉 All tests passed! TCGA Frequency Match is ready to use.")
        print("\nNext steps:")
        print("1. The TCGA database has been created")
        print("2. The tcga_mapper node is integrated into your pipeline")
        print("3. You can now run your full pipeline with TCGA enrichment analysis")
    else:
        print("⚠️  Some tests failed. Check the errors above.")
        
    return success_count == total_tests

if __name__ == "__main__":
    main() 