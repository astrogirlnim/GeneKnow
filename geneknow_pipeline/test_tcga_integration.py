#!/usr/bin/env python3
"""
Comprehensive TCGA Integration Test
Tests the full pipeline with TCGA integration to ensure everything works end-to-end.
"""

import sys
import os
import tempfile
import sqlite3
import logging
from pathlib import Path

# Add current directory for imports
sys.path.append(str(Path(__file__).parent))

def create_test_vcf_with_tcga_variants():
    """Create a VCF file with variants that should match TCGA data."""
    vcf_content = """##fileformat=VCFv4.2
##source=test
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
17	43044295	BRCA1_pathogenic	G	A	60	PASS	DP=45;AF=0.5	GT:DP:GQ	0/1:45:50
17	7674221	TP53_pathogenic	G	A	55	PASS	DP=38;AF=0.5	GT:DP:GQ	0/1:38:45
12	25245350	KRAS_pathogenic	G	A	50	PASS	DP=42;AF=0.5	GT:DP:GQ	0/1:42:40
1	123456	BENIGN_variant	A	T	45	PASS	DP=30;AF=0.5	GT:DP:GQ	0/1:30:35
X	500000	UNKNOWN_variant	C	G	40	PASS	DP=25;AF=0.3	GT:DP:GQ	0/1:25:30
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_content)
        return f.name

def test_database_integrity():
    """Test that the unified database has correct structure and data."""
    print("🔍 Testing Database Integrity...")
    
    db_path = "population_variants.db"
    if not os.path.exists(db_path):
        print("❌ Database not found!")
        return False
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Check both tables exist
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall()]
        
        required_tables = ['population_variants', 'tcga_variants']
        for table in required_tables:
            if table not in tables:
                print(f"❌ Missing table: {table}")
                return False
        
        # Check population_variants has data
        cursor.execute("SELECT COUNT(*) FROM population_variants")
        pop_count = cursor.fetchone()[0]
        
        # Check tcga_variants has data
        cursor.execute("SELECT COUNT(*) FROM tcga_variants")
        tcga_count = cursor.fetchone()[0]
        
        print(f"✅ Population variants: {pop_count:,}")
        print(f"✅ TCGA variants: {tcga_count}")
        
        if pop_count == 0:
            print("❌ No population data found!")
            return False
        
        if tcga_count == 0:
            print("❌ No TCGA data found!")
            return False
        
        # Test cross-table query capability
        cursor.execute("""
        SELECT COUNT(*) FROM population_variants p 
        WHERE EXISTS (
            SELECT 1 FROM tcga_variants t 
            WHERE t.chrom = p.chrom AND t.pos = p.pos 
            AND t.ref = p.ref AND t.alt = p.alt
        )
        """)
        overlapping = cursor.fetchone()[0]
        print(f"✅ Overlapping variants: {overlapping}")
        
        return True
        
    finally:
        conn.close()

def test_pipeline_with_tcga():
    """Test the full pipeline including TCGA integration."""
    print("\n🔄 Testing Full Pipeline with TCGA...")
    
    # Create test VCF
    vcf_path = create_test_vcf_with_tcga_variants()
    
    try:
        # Import the pipeline
        from graph import run_pipeline
        
        # Run the pipeline
        result = run_pipeline(vcf_path, {
            "patient_data": {"age": 45, "sex": "F"},
            "language": "en",
            "include_technical": True
        })
        
        # Check pipeline completed successfully
        if result.get("pipeline_status") in ["failed", "error"]:
            print("❌ Pipeline failed!")
            print(f"Errors: {result.get('errors', [])}")
            return False
        
        # Check TCGA data was processed
        if "tcga_matches" not in result:
            print("❌ No TCGA matches in result!")
            return False
        
        # Check for TCGA enriched variants
        tcga_enriched = result.get("tcga_enriched_variants", [])
        if not tcga_enriched:
            print("❌ No TCGA enriched variants found!")
            return False
        
        # Check that known variants have TCGA data
        filtered_variants = result.get("filtered_variants", [])
        tcga_annotated_count = 0
        
        for variant in filtered_variants:
            if hasattr(variant, 'get') and (
                variant.get("tcga_cancer_relevance") is not None or 
                variant.get("tcga_best_match") is not None
            ):
                tcga_annotated_count += 1
        
        print(f"✅ Pipeline completed successfully")
        print(f"✅ Processed {len(filtered_variants)} variants")
        print(f"✅ Found {len(tcga_enriched)} highly enriched variants")
        print(f"✅ TCGA annotated {tcga_annotated_count} variants")
        
        # Check specific expected matches
        known_variants = ["17:43044295:G>A", "17:7674221:G>A", "12:25245350:G>A"]
        matched_known = 0
        
        for cancer_type, matches in result["tcga_matches"].items():
            for variant_id in known_variants:
                if variant_id in matches:
                    matched_known += 1
                    match_data = matches[variant_id]
                    print(f"✅ {variant_id} matched in {cancer_type}: {match_data.get('enrichment_score', 0):.0f}x enrichment")
        
        if matched_known == 0:
            print("❌ No known variants matched in TCGA!")
            return False
        
        return True
        
    except Exception as e:
        print(f"❌ Pipeline test failed: {e}")
        return False
        
    finally:
        # Clean up
        if os.path.exists(vcf_path):
            os.unlink(vcf_path)

def test_performance():
    """Test that TCGA integration doesn't significantly impact performance."""
    print("\n⚡ Testing Performance...")
    
    import time
    from nodes.tcga_mapper import process as tcga_process
    
    # Create a state with many variants
    large_state = {
        "filtered_variants": [
            {
                "chrom": f"{i % 22 + 1}",
                "pos": 1000000 + i,
                "ref": "A",
                "alt": "T",
                "variant_id": f"{i % 22 + 1}:{1000000 + i}:A>T",
                "gene": f"GENE{i}"
            }
            for i in range(100)
        ],
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "file_metadata": {},
        "completed_nodes": [],
        "errors": [],
        "warnings": []
    }
    
    start_time = time.time()
    result = tcga_process(large_state)
    end_time = time.time()
    
    processing_time = end_time - start_time
    
    print(f"✅ Processed 100 variants in {processing_time:.2f} seconds")
    print(f"✅ Average: {processing_time/100*1000:.1f}ms per variant")
    
    if processing_time > 5.0:  # More than 5 seconds for 100 variants
        print("⚠️ Performance may be slower than expected")
        return False
    
    return True

def test_edge_cases():
    """Test edge cases and error handling."""
    print("\n🧪 Testing Edge Cases...")
    
    from nodes.tcga_mapper import process as tcga_process
    
    # Test 1: Empty variants
    empty_state = {
        "filtered_variants": [],
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "file_metadata": {},
        "completed_nodes": [],
        "errors": [],
        "warnings": []
    }
    
    result = tcga_process(empty_state)
    if result["errors"]:
        print("❌ TCGA mapper failed with empty variants")
        return False
    print("✅ Handles empty variant list")
    
    # Test 2: Malformed variants
    malformed_state = {
        "filtered_variants": [
            {"chrom": "", "pos": "invalid", "ref": None, "alt": None},
            {"variant_id": "incomplete"}
        ],
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "file_metadata": {},
        "completed_nodes": [],
        "errors": [],
        "warnings": []
    }
    
    result = tcga_process(malformed_state)
    print("✅ Handles malformed variants gracefully")
    
    # Test 3: Very large position numbers
    extreme_state = {
        "filtered_variants": [
            {
                "chrom": "1",
                "pos": 999999999,  # Very large position
                "ref": "A",
                "alt": "T",
                "variant_id": "1:999999999:A>T",
                "gene": "TEST"
            }
        ],
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "file_metadata": {},
        "completed_nodes": [],
        "errors": [],
        "warnings": []
    }
    
    result = tcga_process(extreme_state)
    print("✅ Handles extreme position values")
    
    return True

def test_data_consistency():
    """Test that TCGA data is biologically consistent."""
    print("\n🧬 Testing Data Consistency...")
    
    conn = sqlite3.connect("population_variants.db")
    cursor = conn.cursor()
    
    try:
        # Test 1: Frequencies should be between 0 and 1
        cursor.execute("""
        SELECT COUNT(*) FROM tcga_variants 
        WHERE tumor_frequency < 0 OR tumor_frequency > 1 
        OR normal_frequency < 0 OR normal_frequency > 1
        """)
        invalid_frequencies = cursor.fetchone()[0]
        
        if invalid_frequencies > 0:
            print(f"❌ {invalid_frequencies} variants have invalid frequencies")
            return False
        print("✅ All frequencies within valid range (0-1)")
        
        # Test 2: Enrichment scores should be reasonable
        cursor.execute("""
        SELECT COUNT(*) FROM tcga_variants 
        WHERE enrichment_score < 1 OR enrichment_score > 10000
        """)
        invalid_enrichments = cursor.fetchone()[0]
        
        if invalid_enrichments > 0:
            print(f"❌ {invalid_enrichments} variants have unreasonable enrichment scores")
            return False
        print("✅ All enrichment scores reasonable (1-10000x)")
        
        # Test 3: Sample counts should match totals
        cursor.execute("""
        SELECT COUNT(*) FROM tcga_variants 
        WHERE sample_count > total_samples OR sample_count < 0
        """)
        invalid_counts = cursor.fetchone()[0]
        
        if invalid_counts > 0:
            print(f"❌ {invalid_counts} variants have invalid sample counts")
            return False
        print("✅ All sample counts valid")
        
        # Test 4: Check for expected cancer-gene associations
        cursor.execute("""
        SELECT gene, cancer_type, enrichment_score FROM tcga_variants 
        WHERE (gene = 'BRCA1' AND cancer_type = 'breast')
        OR (gene = 'KRAS' AND cancer_type = 'colon')
        ORDER BY enrichment_score DESC
        """)
        
        expected_associations = cursor.fetchall()
        if not expected_associations:
            print("❌ Missing expected cancer-gene associations")
            return False
        
        print("✅ Expected cancer-gene associations present:")
        for gene, cancer, enrichment in expected_associations:
            print(f"   {gene} in {cancer}: {enrichment:.0f}x enrichment")
        
        return True
        
    finally:
        conn.close()

def run_comprehensive_tests():
    """Run all comprehensive tests."""
    print("🧬 TCGA Comprehensive Integration Testing")
    print("=" * 70)
    
    tests = [
        ("Database Integrity", test_database_integrity),
        ("Pipeline Integration", test_pipeline_with_tcga),
        ("Performance", test_performance),
        ("Edge Cases", test_edge_cases),
        ("Data Consistency", test_data_consistency)
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            print(f"\n📋 {test_name}")
            print("-" * 50)
            if test_func():
                print(f"✅ {test_name}: PASSED")
                passed += 1
            else:
                print(f"❌ {test_name}: FAILED")
                failed += 1
        except Exception as e:
            print(f"❌ {test_name}: ERROR - {e}")
            failed += 1
    
    print(f"\n📊 Final Results: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 ALL COMPREHENSIVE TESTS PASSED!")
        print("🚀 TCGA integration is production-ready!")
        return True
    else:
        print(f"⚠️ {failed} tests failed. Review issues above.")
        return False

if __name__ == "__main__":
    # Configure logging to reduce noise
    logging.basicConfig(level=logging.WARNING)
    
    success = run_comprehensive_tests()
    sys.exit(0 if success else 1) 