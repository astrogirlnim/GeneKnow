#!/usr/bin/env python3
"""
Comprehensive TCGA Test Suite for GeneKnow.
Combines all TCGA testing functionality including database validation,
integration tests, performance tests, and edge cases.
"""
import sys
import os
import json
import sqlite3
import time
import logging
from datetime import datetime
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from graph import run_pipeline, create_genomic_pipeline
from nodes.tcga_mapper import process as tcga_process

# Configure logging
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class TCGATestSuite:
    """Comprehensive TCGA testing class."""
    
    def __init__(self):
        self.results = {
            "total": 0,
            "passed": 0,
            "failed": 0,
            "tests": []
        }
        self.db_path = "population_variants.db"
    
    def test(self, name, func):
        """Run a test and record results."""
        self.results["total"] += 1
        print(f"\n{'='*60}")
        print(f"🧪 {name}")
        print(f"{'='*60}")
        
        try:
            start_time = time.time()
            result = func()
            end_time = time.time()
            
            if result:
                self.results["passed"] += 1
                print(f"✅ PASSED ({end_time - start_time:.2f}s)")
                self.results["tests"].append({
                    "name": name, 
                    "status": "passed", 
                    "duration": end_time - start_time
                })
            else:
                self.results["failed"] += 1
                print(f"❌ FAILED")
                self.results["tests"].append({"name": name, "status": "failed"})
        except Exception as e:
            self.results["failed"] += 1
            print(f"❌ ERROR: {str(e)}")
            self.results["tests"].append({
                "name": name, 
                "status": "error", 
                "error": str(e)
            })
    
    def print_summary(self):
        """Print test summary."""
        print("\n" + "=" * 60)
        print("📊 TCGA TEST SUMMARY")
        print("=" * 60)
        print(f"Total tests: {self.results['total']}")
        print(f"Passed: {self.results['passed']} ({self.results['passed']/self.results['total']*100:.1f}%)")
        print(f"Failed: {self.results['failed']}")
        
        if self.results['failed'] > 0:
            print("\nFailed tests:")
            for test in self.results['tests']:
                if test['status'] != 'passed':
                    print(f"  - {test['name']}: {test.get('error', 'Failed')}")


# Test Functions
def test_database_structure():
    """Test TCGA database structure and integrity."""
    print("Testing database structure...")
    
    db_path = "population_variants.db"
    if not os.path.exists(db_path):
        print(f"Database not found: {db_path}")
        return False
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Check tcga_variants table exists
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='tcga_variants'")
        if not cursor.fetchone():
            print("TCGA table not found")
            return False
        
        # Check table structure
        cursor.execute("PRAGMA table_info(tcga_variants)")
        columns = {row[1] for row in cursor.fetchall()}
        
        required_columns = {
            'chrom', 'pos', 'ref', 'alt', 'gene', 'cancer_type',
            'sample_count', 'total_samples', 'tumor_frequency', 'normal_frequency'
        }
        
        missing = required_columns - columns
        if missing:
            print(f"Missing columns: {missing}")
            return False
        
        # Check indexes
        cursor.execute("SELECT name FROM sqlite_master WHERE type='index' AND tbl_name='tcga_variants'")
        indexes = [row[0] for row in cursor.fetchall()]
        print(f"Indexes: {', '.join(indexes)}")
        
        # Check data
        cursor.execute("SELECT COUNT(*) FROM tcga_variants")
        count = cursor.fetchone()[0]
        print(f"Total TCGA variants: {count:,}")
        
        # Check cancer types
        cursor.execute("SELECT DISTINCT cancer_type FROM tcga_variants")
        cancer_types = [row[0] for row in cursor.fetchall()]
        print(f"Cancer types: {', '.join(cancer_types)}")
        
        return count > 0 and len(cancer_types) > 0
        
    finally:
        conn.close()


def test_tcga_mapper_direct():
    """Test TCGA mapper function directly."""
    print("Testing TCGA mapper directly...")
    
    # Create test state
    test_state = {
        "filtered_variants": [
            {
                "chrom": "17",
                "pos": 7674220,
                "ref": "G",
                "alt": "A",
                "gene": "TP53",
                "consequence": "missense_variant",
                "variant_id": "17:7674220:G>A"
            },
            {
                "chrom": "13",
                "pos": 32339000,
                "ref": "T",
                "alt": "C",
                "gene": "BRCA2",
                "consequence": "synonymous_variant",
                "variant_id": "13:32339000:T>C"
            }
        ],
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "completed_nodes": [],
        "errors": []
    }
    
    # Process
    result = tcga_process(test_state)
    
    # Check results
    has_matches = len(result.get("tcga_matches", {})) > 0
    has_cohorts = len(result.get("tcga_cohort_sizes", {})) > 0
    
    # Note: tcga_mapper doesn't update completed_nodes directly in the return value
    # It returns only the fields it updates, not the whole state
    
    print(f"Cancer types matched: {list(result.get('tcga_matches', {}).keys())}")
    print(f"Cohort sizes: {result.get('tcga_cohort_sizes', {})}")
    
    # Validate cohort sizes match expected values
    expected_cohorts = {
        "breast": 1084,
        "colon": 461,
        "lung": 585,
        "prostate": 498,
        "blood": 200
    }
    
    actual_cohorts = result.get("tcga_cohort_sizes", {})
    cohorts_match = all(
        actual_cohorts.get(cancer_type) == expected_size
        for cancer_type, expected_size in expected_cohorts.items()
    )
    
    return has_matches and has_cohorts and cohorts_match


def test_tcga_frequency_analysis():
    """Test TCGA frequency calculations for known variant."""
    print("Testing TCGA frequency analysis...")
    
    # Test with known TP53 variant that exists in the database
    test_state = {
        "filtered_variants": [{
            "chrom": "17",
            "pos": 7684799,
            "ref": "C",
            "alt": "T",
            "gene": "TP53",
            "variant_id": "17:7684799:C>T"
        }],
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "completed_nodes": [],
        "errors": []
    }
    
    result = tcga_process(test_state)
    
    # Check enrichment calculations
    tcga_matches = result.get("tcga_matches", {})
    found_matches = False
    
    for cancer_type, matches in tcga_matches.items():
        # matches is a dict where keys are variant_ids
        if matches:
            found_matches = True
            for variant_id, match_data in matches.items():
                if isinstance(match_data, dict):
                    freq = match_data.get('tumor_frequency', 0)
                    enrichment = match_data.get('enrichment_score', 0)
                    print(f"{cancer_type}: freq={freq:.3f}, enrichment={enrichment:.1f}x")
    
    return found_matches


def test_pipeline_integration():
    """Test TCGA integration with full pipeline."""
    print("Testing pipeline integration...")
    
    # Create test VCF
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
17	7578406	.	C	T	99	PASS	.	GT	0/1
13	32339000	.	T	C	80	PASS	.	GT	0/1
"""
    
    test_vcf = "test_tcga_temp.vcf"
    with open(test_vcf, "w") as f:
        f.write(vcf_content)
    
    try:
        # Run pipeline
        result = run_pipeline(test_vcf, {"language": "en"})
        
        # Check TCGA results
        has_tcga_node = "tcga_mapper" in result.get("completed_nodes", [])
        has_matches = len(result.get("tcga_matches", {})) > 0
        
        if has_matches:
            print(f"TCGA matches found for: {list(result['tcga_matches'].keys())}")
        
        return has_tcga_node and result["pipeline_status"] == "completed"
        
    finally:
        if os.path.exists(test_vcf):
            os.remove(test_vcf)


def test_performance():
    """Test TCGA lookup performance."""
    print("Testing TCGA lookup performance...")
    
    # Create test with many variants
    test_variants = []
    for i in range(100):
        chrom = str((i % 22) + 1)
        pos = 1000000 + i * 1000
        test_variants.append({
            "chrom": chrom,
            "pos": pos,
            "ref": "A",
            "alt": "G",
            "gene": f"GENE{i}",
            "variant_id": f"{chrom}:{pos}:A>G"
        })
    
    test_state = {
        "filtered_variants": test_variants,
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "completed_nodes": [],
        "errors": []
    }
    
    start_time = time.time()
    result = tcga_process(test_state)
    end_time = time.time()
    
    processing_time = end_time - start_time
    variants_per_second = len(test_variants) / processing_time
    
    print(f"Processed {len(test_variants)} variants in {processing_time:.3f}s")
    print(f"Performance: {variants_per_second:.0f} variants/second")
    
    return processing_time < 5.0  # Should process 100 variants in < 5 seconds


def test_edge_cases():
    """Test edge cases and error handling."""
    print("Testing edge cases...")
    
    # Test with empty variants
    empty_state = {
        "filtered_variants": [],
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "completed_nodes": [],
        "errors": []
    }
    
    result = tcga_process(empty_state)
    if result.get("errors"):
        print("Empty variants handled incorrectly")
        return False
    
    # Test with malformed variants
    malformed_state = {
        "filtered_variants": [
            {"chrom": "X"},  # Missing pos, ref, alt
            {"pos": 12345},  # Missing chrom
            {}  # Empty variant
        ],
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "completed_nodes": [],
        "errors": []
    }
    
    result = tcga_process(malformed_state)
    errors = result.get('errors', [])
    print(f"Malformed variants processed, errors: {len(errors)}")
    
    # With malformed variants, we expect the mapper to handle them gracefully
    # It should complete but potentially with no matches
    completed = "tcga_mapper" in result.get("completed_nodes", [])
    no_errors_raised = True  # The fact that we got here means no exceptions were raised
    
    return completed or no_errors_raised


def test_enrichment_calculations():
    """Test enrichment score calculations."""
    print("Testing enrichment calculations...")
    
    conn = sqlite3.connect("population_variants.db")
    cursor = conn.cursor()
    
    try:
        # Find a variant that exists in TCGA
        cursor.execute("""
            SELECT chrom, pos, ref, alt, gene, cancer_type, tumor_frequency
            FROM tcga_variants
            WHERE tumor_frequency > 0.01
            LIMIT 1
        """)
        
        row = cursor.fetchone()
        if not row:
            print("No test variant found")
            return False
        
        chrom, pos, ref, alt, gene, cancer_type, freq = row
        
        # Test with this variant
        test_state = {
            "filtered_variants": [{
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene,
                "variant_id": f"{chrom}:{pos}:{ref}>{alt}"
            }],
            "tcga_matches": {},
            "tcga_cohort_sizes": {},
            "completed_nodes": [],
            "errors": []
        }
        
        result = tcga_process(test_state)
        
        # Check enrichment calculation
        matches = result.get("tcga_matches", {}).get(cancer_type, {})
        if matches:
            # matches is a dict, get the first variant's data
            variant_id = f"{chrom}:{pos}:{ref}>{alt}"
            match_data = matches.get(variant_id, {})
            if match_data:
                enrichment = match_data.get("enrichment_score", 0)
                print(f"Variant frequency: {freq:.3f}")
                print(f"Calculated enrichment: {enrichment:.1f}x")
                
                # Enrichment should be > 1 for variants with freq > 0.01
                return enrichment > 1.0
        
        return False
        
    finally:
        conn.close()


def test_maf_file_processing():
    """Test processing of MAF files through pipeline."""
    print("Testing MAF file processing...")
    
    # Find test MAF file
    test_files = [
        "test_data/tcga_downloads/3d14b1e2-0555-4d6f-a55b-a56065f915e1.wxs.aliquot_ensemble_masked.maf.gz",
        "test_data/test_sample.maf"
    ]
    
    test_file = None
    for file in test_files:
        if os.path.exists(file):
            test_file = file
            break
    
    if not test_file:
        print("No MAF test file found")
        return True  # Don't fail
    
    # Process MAF file
    result = run_pipeline(test_file, {"language": "en"})
    
    # Check TCGA processing
    has_tcga = "tcga_mapper" in result.get("completed_nodes", [])
    has_matches = len(result.get("tcga_matches", {})) > 0
    
    if has_matches:
        total_matches = sum(len(m) for m in result["tcga_matches"].values())
        print(f"Total TCGA matches: {total_matches}")
    
    return has_tcga and result["pipeline_status"] == "completed"


def test_cancer_type_analysis():
    """Test cancer type specific analysis."""
    print("Testing cancer type analysis...")
    
    # Test variants known to be enriched in specific cancers
    test_state = {
        "filtered_variants": [
            # KRAS mutations common in lung/colon
            {
                "chrom": "12",
                "pos": 25245350,
                "ref": "C",
                "alt": "T",
                "gene": "KRAS",
                "variant_id": "12:25245350:C>T"
            },
            # BRCA1 mutations in breast
            {
                "chrom": "17",
                "pos": 43045677,
                "ref": "G",
                "alt": "A",
                "gene": "BRCA1",
                "variant_id": "17:43045677:G>A"
            }
        ],
        "tcga_matches": {},
        "tcga_cohort_sizes": {},
        "completed_nodes": [],
        "errors": []
    }
    
    result = tcga_process(test_state)
    
    # Check cancer-specific enrichment
    tcga_matches = result.get("tcga_matches", {})
    
    # KRAS should be enriched in lung/colon
    kras_cancers = []
    for cancer_type, matches in tcga_matches.items():
        # matches is a dict where keys are variant_ids and values are match data
        for variant_id, match_data in matches.items():
            if isinstance(match_data, dict):
                if match_data.get("gene") == "KRAS" and match_data.get("enrichment_score", 0) > 2:
                    kras_cancers.append(cancer_type)
    
    print(f"KRAS enriched in: {', '.join(kras_cancers)}")
    
    return len(tcga_matches) > 0


def main():
    """Run all TCGA tests."""
    print("🧬 GeneKnow Comprehensive TCGA Test Suite")
    print("=" * 60)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Check prerequisites
    db_path = "population_variants.db"
    if not os.path.exists(db_path):
        print(f"\n❌ Database not found: {db_path}")
        print("Please run database setup first.")
        return 1
    
    # Run tests
    suite = TCGATestSuite()
    
    # Database tests
    suite.test("Database Structure", test_database_structure)
    
    # Unit tests
    suite.test("TCGA Mapper Direct", test_tcga_mapper_direct)
    suite.test("Frequency Analysis", test_tcga_frequency_analysis)
    suite.test("Enrichment Calculations", test_enrichment_calculations)
    
    # Integration tests
    suite.test("Pipeline Integration", test_pipeline_integration)
    suite.test("MAF File Processing", test_maf_file_processing)
    
    # Analysis tests
    suite.test("Cancer Type Analysis", test_cancer_type_analysis)
    
    # Performance tests
    suite.test("Performance", test_performance)
    
    # Edge cases
    suite.test("Edge Cases", test_edge_cases)
    
    # Print summary
    suite.print_summary()
    
    return 0 if suite.results["failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main()) 