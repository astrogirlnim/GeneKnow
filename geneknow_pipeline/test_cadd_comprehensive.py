#!/usr/bin/env python3
"""
Comprehensive test suite for CADD scoring node implementation.
Tests all aspects: unit, integration, performance, memory, edge cases.
"""
import os
import sys
import json
import sqlite3
import tempfile
import unittest
import time
import psutil
import subprocess
from datetime import datetime
from unittest.mock import patch, MagicMock
import logging

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import CADD module
from nodes.cadd_scoring import (
    normalize_chromosome,
    lookup_cadd_score,
    calculate_risk_weight,
    process,
    create_job_record,
    update_job_status,
    query_remote_cadd,
    CADD_PHRED_THRESHOLDS
)

# Set up logging for tests
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestCADDUnit(unittest.TestCase):
    """Enhanced unit tests for CADD scoring functionality."""
    
    def setUp(self):
        """Set up test environment with comprehensive test data."""
        # Create temporary database
        self.temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        self.temp_db_path = self.temp_db.name
        self.temp_db.close()
        
        # Create comprehensive test database
        self._create_test_database()
        
        # Store original DB path
        self.original_db_path = os.environ.get("CADD_DB_PATH")
        os.environ["CADD_DB_PATH"] = self.temp_db_path
        
    def tearDown(self):
        """Clean up test environment."""
        if os.path.exists(self.temp_db_path):
            os.unlink(self.temp_db_path)
            
        if self.original_db_path:
            os.environ["CADD_DB_PATH"] = self.original_db_path
        elif "CADD_DB_PATH" in os.environ:
            del os.environ["CADD_DB_PATH"]
            
    def _create_test_database(self):
        """Create comprehensive test database with various scenarios."""
        conn = sqlite3.connect(self.temp_db_path)
        cursor = conn.cursor()
        
        # Create tables
        cursor.execute("""
            CREATE TABLE population_variants (
                chrom TEXT,
                pos INTEGER,
                ref TEXT,
                alt TEXT,
                gene TEXT,
                gnomad_af REAL DEFAULT 0.0,
                clinical_significance TEXT,
                is_pathogenic INTEGER DEFAULT 0,
                consequence TEXT,
                review_status TEXT,
                PRIMARY KEY (chrom, pos, ref, alt)
            )
        """)
        
        cursor.execute("""
            CREATE TABLE cadd_scores (
                chrom TEXT NOT NULL,
                pos INTEGER NOT NULL,
                ref TEXT NOT NULL,
                alt TEXT NOT NULL,
                raw_score REAL,
                phred_score REAL,
                job_id TEXT,
                source TEXT DEFAULT 'local',
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (chrom, pos, ref, alt)
            )
        """)
        
        cursor.execute("""
            CREATE TABLE cadd_jobs (
                job_id TEXT PRIMARY KEY,
                job_type TEXT NOT NULL,
                started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                completed_at TIMESTAMP,
                variant_count INTEGER DEFAULT 0,
                status TEXT DEFAULT 'running',
                metadata TEXT
            )
        """)
        
        # Insert comprehensive test data
        test_cadd_scores = [
            # Cancer genes with various scores
            ('17', 43044295, 'A', 'T', 2.345, 24.7, 'test_job_001'),  # BRCA1 - pathogenic
            ('17', 43045719, 'G', 'A', 3.456, 34.1, 'test_job_001'),  # BRCA1 - highly pathogenic
            ('13', 32315474, 'G', 'A', 1.234, 15.3, 'test_job_001'),  # BRCA2 - uncertain
            ('7', 55019017, 'T', 'G', 2.890, 28.5, 'test_job_001'),   # EGFR - pathogenic
            ('17', 7578406, 'C', 'T', 2.123, 22.8, 'test_job_001'),   # TP53 - pathogenic
            
            # Benign variants
            ('5', 112839936, 'C', 'T', 0.456, 8.2, 'test_job_001'),   # APC - benign
            ('1', 100000, 'G', 'A', 0.123, 3.5, 'test_job_001'),      # Non-cancer gene - benign
            
            # Edge cases
            ('X', 153296777, 'C', 'T', 1.789, 18.9, 'test_job_001'),  # X chromosome
            ('Y', 2786989, 'G', 'A', 1.456, 16.2, 'test_job_001'),    # Y chromosome
            ('MT', 8697, 'G', 'A', 0.789, 9.8, 'test_job_001'),       # Mitochondrial
            
            # High-impact variants
            ('2', 212578380, 'G', 'A', 3.890, 38.9, 'test_job_001'),  # ERBB4 - very high
            ('3', 178952085, 'A', 'G', 4.123, 41.2, 'test_job_001'),  # PIK3CA - extremely high
        ]
        
        cursor.executemany(
            "INSERT INTO cadd_scores (chrom, pos, ref, alt, raw_score, phred_score, job_id) VALUES (?, ?, ?, ?, ?, ?, ?)",
            test_cadd_scores
        )
        
        conn.commit()
        conn.close()
        
    def test_chromosome_normalization_comprehensive(self):
        """Test chromosome normalization with all possible formats."""
        test_cases = [
            # Standard chromosomes
            ("1", "1"), ("chr1", "1"),
            ("22", "22"), ("chr22", "22"),
            # Sex chromosomes
            ("X", "X"), ("chrX", "X"),
            ("Y", "Y"), ("chrY", "Y"),
            # Mitochondrial
            ("MT", "MT"), ("chrMT", "MT"),
            ("M", "M"), ("chrM", "M"),
            # Already normalized
            ("1", "1"), ("X", "X"),
        ]
        
        for input_chr, expected in test_cases:
            with self.subTest(chromosome=input_chr):
                result = normalize_chromosome(input_chr)
                self.assertEqual(result, expected)
                
    def test_risk_weight_calculation_comprehensive(self):
        """Test risk weight calculation across full PHRED range."""
        # Test boundary values
        test_cases = [
            # Benign range (< 10)
            (0.0, 0.1), (5.0, 0.1), (9.9, 0.1),
            # Uncertain range (10-15)
            (10.0, 0.1), (12.5, 0.2), (15.0, 0.3),
            # Damaging range (15-20)
            (15.1, 0.306), (17.5, 0.45), (19.9, 0.594),
            # Pathogenic range (20-25)
            (20.0, 0.6), (22.5, 0.7), (24.9, 0.796),
            # Highly pathogenic (> 25)
            (25.0, 0.8), (30.0, 0.9), (35.0, 1.0),
            (40.0, 1.0), (50.0, 1.0),  # Should cap at 1.0
        ]
        
        for phred, expected in test_cases:
            with self.subTest(phred=phred):
                result = calculate_risk_weight(phred)
                self.assertAlmostEqual(result, expected, places=2,
                    msg=f"PHRED {phred} should give risk weight ~{expected}")


class TestCADDPerformance(unittest.TestCase):
    """Performance tests for CADD scoring."""
    
    def setUp(self):
        """Set up performance test environment."""
        self.temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        self.temp_db_path = self.temp_db.name
        self.temp_db.close()
        self._create_large_test_database()
        
    def tearDown(self):
        """Clean up test environment."""
        if os.path.exists(self.temp_db_path):
            os.unlink(self.temp_db_path)
            
    def _create_large_test_database(self):
        """Create large test database for performance testing."""
        conn = sqlite3.connect(self.temp_db_path)
        cursor = conn.cursor()
        
        # Create tables with indexes
        cursor.execute("""
            CREATE TABLE cadd_scores (
                chrom TEXT NOT NULL,
                pos INTEGER NOT NULL,
                ref TEXT NOT NULL,
                alt TEXT NOT NULL,
                raw_score REAL,
                phred_score REAL,
                job_id TEXT,
                PRIMARY KEY (chrom, pos, ref, alt)
            )
        """)
        
        # Create index for fast lookups
        cursor.execute("""
            CREATE INDEX idx_cadd_lookup ON cadd_scores(chrom, pos, ref, alt)
        """)
        
        # Insert 100k test variants
        logger.info("Creating large test database with 100k variants...")
        test_data = []
        for i in range(100000):
            chrom = str((i % 22) + 1)
            pos = 1000000 + i * 100
            ref = 'ACGT'[i % 4]
            alt = 'ACGT'[(i + 1) % 4]
            raw = (i % 100) / 20.0
            phred = (i % 40) + 5
            test_data.append((chrom, pos, ref, alt, raw, phred, 'perf_test'))
            
            if len(test_data) >= 10000:
                cursor.executemany(
                    "INSERT INTO cadd_scores VALUES (?, ?, ?, ?, ?, ?, ?)",
                    test_data
                )
                test_data = []
                
        conn.commit()
        conn.close()
        logger.info("Large test database created")
        
    def test_lookup_performance(self):
        """Test that lookups meet performance requirements (<200ms per 10k)."""
        conn = sqlite3.connect(self.temp_db_path)
        
        # Prepare 10k lookup queries
        queries = []
        for i in range(10000):
            chrom = str((i % 22) + 1)
            pos = 1000000 + i * 100
            ref = 'ACGT'[i % 4]
            alt = 'ACGT'[(i + 1) % 4]
            queries.append((chrom, pos, ref, alt))
            
        # Time the lookups
        start_time = time.time()
        successful_lookups = 0
        
        for chrom, pos, ref, alt in queries:
            result = lookup_cadd_score(chrom, pos, ref, alt, conn)
            if result:
                successful_lookups += 1
                
        end_time = time.time()
        elapsed_ms = (end_time - start_time) * 1000
        
        conn.close()
        
        # Log results
        logger.info(f"Performance test: {successful_lookups}/10000 lookups in {elapsed_ms:.1f}ms")
        logger.info(f"Average: {elapsed_ms/10000:.3f}ms per lookup")
        
        # Assert performance requirement
        self.assertLess(elapsed_ms, 200,
            f"10k lookups took {elapsed_ms:.1f}ms, requirement is <200ms")
        self.assertEqual(successful_lookups, 10000,
            "All lookups should succeed in performance test")


class TestCADDMemory(unittest.TestCase):
    """Memory usage tests for CADD scoring."""
    
    def test_memory_usage(self):
        """Test that SQLite page cache stays under 200MB."""
        # Create a moderate test database
        temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        temp_db_path = temp_db.name
        temp_db.close()
        
        try:
            # Create database with 50k variants
            conn = sqlite3.connect(temp_db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
                CREATE TABLE cadd_scores (
                    chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
                    raw_score REAL, phred_score REAL, job_id TEXT,
                    PRIMARY KEY (chrom, pos, ref, alt)
                )
            """)
            
            # Insert test data
            for i in range(50000):
                cursor.execute(
                    "INSERT INTO cadd_scores VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (str((i % 22) + 1), 1000000 + i, 'A', 'G', 1.0, 15.0, 'mem_test')
                )
                
            conn.commit()
            
            # Get process memory before operations
            process = psutil.Process()
            mem_before = process.memory_info().rss / 1024 / 1024  # MB
            
            # Perform many lookups to load cache
            for i in range(10000):
                cursor.execute(
                    "SELECT * FROM cadd_scores WHERE chrom=? AND pos=? AND ref=? AND alt=?",
                    (str((i % 22) + 1), 1000000 + i, 'A', 'G')
                )
                cursor.fetchone()
                
            # Get memory after operations
            mem_after = process.memory_info().rss / 1024 / 1024  # MB
            mem_increase = mem_after - mem_before
            
            conn.close()
            
            logger.info(f"Memory test: {mem_before:.1f}MB -> {mem_after:.1f}MB "
                       f"(+{mem_increase:.1f}MB)")
            
            # Assert memory requirement
            self.assertLess(mem_increase, 200,
                f"Memory increase {mem_increase:.1f}MB exceeds 200MB limit")
                
        finally:
            if os.path.exists(temp_db_path):
                os.unlink(temp_db_path)


class TestCADDEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""
    
    def test_missing_database(self):
        """Test behavior when database is missing."""
        # Use non-existent database path
        import nodes.cadd_scoring as cadd_module
        original_path = cadd_module.POP_DB_PATH
        cadd_module.POP_DB_PATH = "/non/existent/database.db"
        
        try:
            state = {
                "filtered_variants": [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}],
                "completed_nodes": [],
                "errors": []
            }
            
            result = process(state)
            
            # Should not fail, just pass through
            self.assertIn("cadd_scoring", result["completed_nodes"])
            self.assertIn("error", result["cadd_stats"])
            
        finally:
            cadd_module.POP_DB_PATH = original_path
            
    def test_malformed_variants(self):
        """Test handling of malformed variant data."""
        state = {
            "filtered_variants": [
                {"chrom": "1"},  # Missing pos, ref, alt
                {"pos": 100, "ref": "A", "alt": "G"},  # Missing chrom
                {"chrom": "X", "pos": "not_a_number", "ref": "A", "alt": "G"},  # Invalid pos
                {},  # Empty variant
            ],
            "completed_nodes": [],
            "errors": []
        }
        
        # Should handle gracefully
        result = process(state)
        self.assertIn("cadd_scoring", result["completed_nodes"])
        
    @patch('nodes.cadd_scoring.requests')
    def test_remote_fallback(self, mock_requests):
        """Test remote CADD lookup fallback."""
        # Mock remote response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "1\t12345\tA\tG\t1.234\t15.6\n"
        mock_requests.get.return_value = mock_response
        
        # Enable remote lookups
        os.environ["USE_REMOTE_CADD"] = "true"
        
        conn = sqlite3.connect(":memory:")
        result = query_remote_cadd("1", 12345, "A", "G", "test_job", conn)
        
        # Currently returns None (not implemented), but structure is there
        self.assertIsNone(result)  # Will change when implemented


class TestCADDIntegration(unittest.TestCase):
    """Integration tests with full pipeline."""
    
    def test_pipeline_integration(self):
        """Test CADD scoring in full pipeline context."""
        from graph import create_genomic_pipeline
        
        # Create test state
        test_state = {
            "file_path": "test.maf",
            "file_type": "maf",
            "filtered_variants": [
                {
                    "variant_id": "17:43044295:A>T",
                    "chrom": "17",
                    "pos": 43044295,
                    "ref": "A",
                    "alt": "T",
                    "gene": "BRCA1",
                    "consequence": "missense_variant",
                    "quality": 100,
                    "depth": 50,
                    "allele_freq": 0.5
                }
            ],
            "variant_count": 1,
            "completed_nodes": ["file_input", "preprocess", "qc_filter", "population_mapper"],
            "errors": [],
            "warnings": [],
            "tcga_matches": {},
            "risk_genes": {"breast": ["BRCA1", "BRCA2"]},
            "pipeline_start_time": datetime.now()
        }
        
        # Run CADD scoring node directly
        from nodes import cadd_scoring
        result = cadd_scoring.process(test_state)
        
        # Verify integration
        self.assertIn("cadd_enriched_variants", result)
        self.assertIn("cadd_stats", result)
        self.assertIn("cadd_scoring", result["completed_nodes"])
        
        # Check that enriched variants are properly formatted
        if result["cadd_enriched_variants"]:
            variant = result["cadd_enriched_variants"][0]
            if "cadd_phred" in variant:
                self.assertIn("cadd_raw", variant)
                self.assertIn("cadd_risk_weight", variant)
                self.assertIsInstance(variant["cadd_phred"], (int, float))


def run_all_tests():
    """Run all test suites and generate report."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestCADDUnit))
    suite.addTests(loader.loadTestsFromTestCase(TestCADDPerformance))
    suite.addTests(loader.loadTestsFromTestCase(TestCADDMemory))
    suite.addTests(loader.loadTestsFromTestCase(TestCADDEdgeCases))
    suite.addTests(loader.loadTestsFromTestCase(TestCADDIntegration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Generate summary
    print("\n" + "=" * 70)
    print("CADD SCORING TEST SUMMARY")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {((result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100):.1f}%")
    print("=" * 70)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1) 