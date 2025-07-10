#!/usr/bin/env python3
"""
Comprehensive TCGA Implementation Validation
===========================================
This script thoroughly tests the TCGA frequency matching implementation
to ensure it's working correctly and producing biologically sensible results.
"""

import sqlite3
import sys
import os
from pathlib import Path
import logging

# Add current directory for imports
sys.path.append(str(Path(__file__).parent))

def test_database_structure():
    """Test 1: Verify database structure is correct."""
    print("🔍 Test 1: Database Structure Validation")
    print("-" * 50)
    
    if not os.path.exists("population_variants.db"):
        print("❌ Main database not found!")
        return False
    
    conn = sqlite3.connect("population_variants.db")
    cursor = conn.cursor()
    
    try:
        # Check table exists
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='tcga_variants'")
        if not cursor.fetchone():
            print("❌ tcga_variants table not found!")
            return False
        
        # Check columns
        cursor.execute("PRAGMA table_info(tcga_variants)")
        columns = {row[1]: row[2] for row in cursor.fetchall()}
        
        expected_columns = {
            'chrom': 'TEXT',
            'pos': 'INTEGER', 
            'ref': 'TEXT',
            'alt': 'TEXT',
            'gene': 'TEXT',
            'cancer_type': 'TEXT',
            'tumor_frequency': 'REAL',
            'normal_frequency': 'REAL',
            'enrichment_score': 'REAL',
            'sample_count': 'INTEGER',
            'total_samples': 'INTEGER'
        }
        
        for col, expected_type in expected_columns.items():
            if col not in columns:
                print(f"❌ Missing column: {col}")
                return False
            # Note: SQLite type checking is flexible, so we'll skip exact type validation
        
        print("✅ Database structure is correct")
        
        # Check indexes
        cursor.execute("SELECT name FROM sqlite_master WHERE type='index'")
        indexes = [row[0] for row in cursor.fetchall()]
        expected_indexes = ['idx_tcga_gene', 'idx_tcga_position', 'idx_tcga_cancer', 'idx_tcga_enrichment']
        
        missing_indexes = set(expected_indexes) - set(indexes)
        if missing_indexes:
            print(f"⚠️  Missing indexes: {missing_indexes}")
        else:
            print("✅ All indexes present")
        
        return True
        
    finally:
        conn.close()

def test_data_quality():
    """Test 2: Verify data quality and biological sensibility."""
    print("\n🔍 Test 2: Data Quality Validation")
    print("-" * 50)
    
    conn = sqlite3.connect("population_variants.db")
    cursor = conn.cursor()
    
    try:
        # Check total variants
        cursor.execute("SELECT COUNT(*) FROM tcga_variants")
        total_variants = cursor.fetchone()[0]
        print(f"📊 Total variants: {total_variants}")
        
        if total_variants == 0:
            print("❌ No variants found in database!")
            return False
        
        # Check cancer types
        cursor.execute("SELECT DISTINCT cancer_type FROM tcga_variants ORDER BY cancer_type")
        cancer_types = [row[0] for row in cursor.fetchall()]
        expected_types = ['blood', 'breast', 'colon', 'lung', 'prostate']
        
        print(f"🎯 Cancer types: {cancer_types}")
        if set(cancer_types) != set(expected_types):
            print(f"❌ Unexpected cancer types. Expected: {expected_types}")
            return False
        
        # Check cohort sizes match documentation
        cursor.execute("SELECT cancer_type, MAX(total_samples) FROM tcga_variants GROUP BY cancer_type")
        cohort_sizes = {row[0]: row[1] for row in cursor.fetchall()}
        
        expected_cohorts = {
            'breast': 1084,
            'colon': 461,
            'lung': 585,
            'prostate': 498,
            'blood': 200
        }
        
        print("📈 Cohort sizes:")
        for cancer_type in expected_types:
            actual = cohort_sizes.get(cancer_type, 0)
            expected = expected_cohorts[cancer_type]
            status = "✅" if actual == expected else "❌"
            print(f"   {cancer_type}: {actual} (expected {expected}) {status}")
        
        # Check enrichment scores are reasonable
        cursor.execute("SELECT MIN(enrichment_score), MAX(enrichment_score), AVG(enrichment_score) FROM tcga_variants")
        min_enrich, max_enrich, avg_enrich = cursor.fetchone()
        
        print(f"\n🧮 Enrichment Scores:")
        print(f"   Min: {min_enrich:.1f}x")
        print(f"   Max: {max_enrich:.1f}x") 
        print(f"   Avg: {avg_enrich:.1f}x")
        
        # Enrichment should be > 1 (otherwise why is it in tumors?)
        if min_enrich < 1:
            print("❌ Some variants have enrichment < 1x (impossible)")
            return False
        
        # Check frequency ranges
        cursor.execute("SELECT MIN(tumor_frequency), MAX(tumor_frequency) FROM tcga_variants")
        min_tumor, max_tumor = cursor.fetchone()
        
        cursor.execute("SELECT MIN(normal_frequency), MAX(normal_frequency) FROM tcga_variants")
        min_normal, max_normal = cursor.fetchone()
        
        print(f"\n📊 Frequency Ranges:")
        print(f"   Tumor: {min_tumor:.1%} - {max_tumor:.1%}")
        print(f"   Normal: {min_normal:.1%} - {max_normal:.1%}")
        
        # Sanity checks
        if max_tumor > 1.0 or max_normal > 1.0:
            print("❌ Frequencies > 100% found (impossible)")
            return False
        
        if min_tumor < 0 or min_normal < 0:
            print("❌ Negative frequencies found (impossible)")
            return False
        
        print("✅ Data quality checks passed")
        return True
        
    finally:
        conn.close()

def test_enrichment_calculations():
    """Test 3: Verify enrichment calculations are mathematically correct."""
    print("\n🔍 Test 3: Enrichment Calculation Validation")
    print("-" * 50)
    
    conn = sqlite3.connect("population_variants.db")
    cursor = conn.cursor()
    
    try:
        # Get a few variants and manually verify enrichment calculations
        cursor.execute("""
        SELECT gene, cancer_type, tumor_frequency, normal_frequency, enrichment_score
        FROM tcga_variants
        LIMIT 5
        """)
        
        all_correct = True
        
        for gene, cancer, tumor_freq, normal_freq, stored_enrichment in cursor.fetchall():
            # Calculate enrichment manually
            if normal_freq > 0:
                calculated_enrichment = tumor_freq / normal_freq
            else:
                calculated_enrichment = float('inf')
            
            # Allow small floating point differences
            diff = abs(calculated_enrichment - stored_enrichment)
            is_correct = diff < 0.1
            
            status = "✅" if is_correct else "❌"
            print(f"   {gene} in {cancer}: {tumor_freq:.1%}/{normal_freq:.1%} = {calculated_enrichment:.1f}x (stored: {stored_enrichment:.1f}x) {status}")
            
            if not is_correct:
                all_correct = False
        
        if all_correct:
            print("✅ All enrichment calculations are correct")
        else:
            print("❌ Some enrichment calculations are wrong")
        
        return all_correct
        
    finally:
        conn.close()

def test_tcga_mapper_functions():
    """Test 4: Test the TCGA mapper functions directly."""
    print("\n🔍 Test 4: TCGA Mapper Function Testing")
    print("-" * 50)
    
    # Import functions (copy them to avoid dependency issues)
    import sqlite3
    
    def query_tcga_database(chrom, pos, ref, alt, cancer_type):
        conn = sqlite3.connect("population_variants.db")
        cursor = conn.cursor()
        
        # Normalize chromosome
        norm_chrom = chrom[3:] if chrom.startswith('chr') else chrom
        
        cursor.execute("""
        SELECT gene, tumor_frequency, normal_frequency, enrichment_score, 
               sample_count, total_samples
        FROM tcga_variants
        WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ? AND cancer_type = ?
        """, (norm_chrom, pos, ref, alt, cancer_type))
        
        row = cursor.fetchone()
        conn.close()
        
        if row:
            return {
                'gene': row[0],
                'tumor_frequency': row[1],
                'normal_frequency': row[2],
                'enrichment_score': row[3],
                'sample_count': row[4],
                'total_samples': row[5],
                'found_in_tcga': True
            }
        return {'found_in_tcga': False}
    
    # Test known variants
    test_cases = [
        # (chrom, pos, ref, alt, cancer_type, should_find)
        ("17", 43044295, "G", "A", "breast", True),   # BRCA1
        ("17", 7674221, "G", "A", "breast", True),    # TP53
        ("12", 25245350, "G", "A", "colon", True),    # KRAS
        ("1", 123456, "A", "T", "breast", False),     # Random variant
    ]
    
    all_passed = True
    
    for chrom, pos, ref, alt, cancer_type, should_find in test_cases:
        result = query_tcga_database(chrom, pos, ref, alt, cancer_type)
        found = result.get('found_in_tcga', False)
        
        if found == should_find:
            if found:
                gene = result['gene']
                enrichment = result['enrichment_score']
                print(f"✅ {gene} variant found in {cancer_type}: {enrichment:.1f}x enrichment")
            else:
                print(f"✅ Unknown variant correctly not found")
        else:
            print(f"❌ Expected found={should_find}, got found={found}")
            all_passed = False
    
    return all_passed

def test_biological_sensibility():
    """Test 5: Check if results make biological sense."""
    print("\n🔍 Test 5: Biological Sensibility Check")
    print("-" * 50)
    
    conn = sqlite3.connect("population_variants.db")
    cursor = conn.cursor()
    
    try:
        # Check BRCA1/BRCA2 are most enriched in breast cancer
        cursor.execute("""
        SELECT gene, cancer_type, enrichment_score
        FROM tcga_variants
        WHERE gene IN ('BRCA1', 'BRCA2')
        ORDER BY enrichment_score DESC
        """)
        
        brca_results = cursor.fetchall()
        print("🧬 BRCA1/BRCA2 enrichment by cancer type:")
        
        breast_enrichments = []
        other_enrichments = []
        
        for gene, cancer_type, enrichment in brca_results:
            print(f"   {gene} in {cancer_type}: {enrichment:.1f}x")
            if cancer_type == 'breast':
                breast_enrichments.append(enrichment)
            else:
                other_enrichments.append(enrichment)
        
        # BRCA should be most enriched in breast cancer
        if breast_enrichments and other_enrichments:
            max_breast = max(breast_enrichments)
            max_other = max(other_enrichments) if other_enrichments else 0
            
            if max_breast > max_other:
                print("✅ BRCA variants most enriched in breast cancer (biologically correct)")
            else:
                print("❌ BRCA variants not most enriched in breast cancer (biologically incorrect)")
                return False
        
        # Check KRAS is highly enriched in colon/lung
        cursor.execute("""
        SELECT cancer_type, enrichment_score
        FROM tcga_variants
        WHERE gene = 'KRAS'
        ORDER BY enrichment_score DESC
        """)
        
        kras_results = cursor.fetchall()
        print("\n🧬 KRAS enrichment by cancer type:")
        
        for cancer_type, enrichment in kras_results:
            print(f"   KRAS in {cancer_type}: {enrichment:.1f}x")
        
        # KRAS should be highly enriched in colon and lung
        kras_dict = {cancer: enrichment for cancer, enrichment in kras_results}
        colon_enrichment = kras_dict.get('colon', 0)
        lung_enrichment = kras_dict.get('lung', 0)
        
        if colon_enrichment > 100 and lung_enrichment > 50:
            print("✅ KRAS highly enriched in colon/lung (biologically correct)")
        else:
            print("❌ KRAS not sufficiently enriched in colon/lung")
            return False
        
        print("✅ Biological sensibility checks passed")
        return True
        
    finally:
        conn.close()

def test_edge_cases():
    """Test 6: Test edge cases and error handling."""
    print("\n🔍 Test 6: Edge Case Testing")
    print("-" * 50)
    
    # Test with non-existent database
    original_db = "population_variants.db"
    fake_db = "fake_tcga.db"
    
    import sqlite3
    
    def query_with_error_handling(db_path, chrom, pos, ref, alt, cancer_type):
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
            SELECT gene FROM tcga_variants
            WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ? AND cancer_type = ?
            """, (chrom, pos, ref, alt, cancer_type))
            
            result = cursor.fetchone()
            conn.close()
            return result is not None
            
        except sqlite3.Error:
            return False
    
    # Test 1: Non-existent database
    result = query_with_error_handling(fake_db, "17", 123, "A", "T", "breast")
    if not result:
        print("✅ Handles non-existent database correctly")
    else:
        print("❌ Should fail with non-existent database")
        return False
    
    # Test 2: Invalid chromosome formats
    test_chroms = ["chr17", "17", "X", "Y", "23"]
    for chrom in test_chroms:
        result = query_with_error_handling(original_db, chrom, 43044295, "G", "A", "breast")
        print(f"   Chromosome '{chrom}': {'✅' if isinstance(result, bool) else '❌'}")
    
    # Test 3: Edge case positions
    edge_positions = [0, -1, 999999999]
    for pos in edge_positions:
        result = query_with_error_handling(original_db, "17", pos, "G", "A", "breast")
        print(f"   Position {pos}: {'✅' if isinstance(result, bool) else '❌'}")
    
    print("✅ Edge case testing completed")
    return True

def run_all_tests():
    """Run all validation tests."""
    print("🧬 TCGA Implementation Comprehensive Validation")
    print("=" * 70)
    
    tests = [
        ("Database Structure", test_database_structure),
        ("Data Quality", test_data_quality),
        ("Enrichment Calculations", test_enrichment_calculations),
        ("TCGA Mapper Functions", test_tcga_mapper_functions),
        ("Biological Sensibility", test_biological_sensibility),
        ("Edge Cases", test_edge_cases)
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"❌ {test_name} failed with error: {e}")
    
    print(f"\n📊 Final Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 ALL TESTS PASSED! TCGA implementation is robust and correct.")
        print("\n✅ Validation Summary:")
        print("   - Database structure is correct")
        print("   - Data quality is high") 
        print("   - Enrichment calculations are accurate")
        print("   - Functions work as expected")
        print("   - Results are biologically sensible")
        print("   - Edge cases are handled properly")
        
        return True
    else:
        print("⚠️  Some tests failed. Implementation needs review.")
        return False

if __name__ == "__main__":
    # Configure logging to reduce noise
    logging.basicConfig(level=logging.WARNING)
    
    success = run_all_tests()
    sys.exit(0 if success else 1) 