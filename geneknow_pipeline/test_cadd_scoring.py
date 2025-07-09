#!/usr/bin/env python3
"""
Unit tests for the CADD scoring node.
Tests lookup functionality, risk weight calculation, and pipeline integration.
"""
import os
import sys
import json
import sqlite3
import tempfile
import unittest
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import directly
from nodes.cadd_scoring import (
    normalize_chromosome,
    lookup_cadd_score,
    calculate_risk_weight,
    process
)


class TestCADDScoring(unittest.TestCase):
    """Test suite for CADD scoring functionality."""
    
    def setUp(self):
        """Set up test environment."""
        # Create temporary database for testing
        self.temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        self.temp_db_path = self.temp_db.name
        self.temp_db.close()
        
        # Create test database
        conn = sqlite3.connect(self.temp_db_path)
        cursor = conn.cursor()
        
        # Create schema
        cursor.execute("""
            CREATE TABLE cadd (
                chrom TEXT NOT NULL,
                pos INTEGER NOT NULL,
                ref TEXT NOT NULL,
                alt TEXT NOT NULL,
                raw_score REAL,
                phred_score REAL,
                PRIMARY KEY (chrom, pos, ref, alt)
            )
        """)
        
        # Insert test data
        test_variants = [
            ('chr17', 43044295, 'A', 'T', 2.345, 24.7),  # BRCA1 - pathogenic
            ('chr17', 43045719, 'G', 'A', 3.456, 34.1),  # BRCA1 - highly pathogenic
            ('chr13', 32315474, 'G', 'A', 1.234, 15.3),  # BRCA2 - uncertain
            ('chr5', 112839936, 'C', 'T', 0.456, 8.2),   # APC - benign
            ('chr7', 55019017, 'T', 'G', 2.890, 28.5),   # EGFR - pathogenic
        ]
        
        cursor.executemany(
            "INSERT INTO cadd VALUES (?, ?, ?, ?, ?, ?)",
            test_variants
        )
        
        conn.commit()
        conn.close()
        
        # Set environment variable to use test database
        os.environ['CADD_DB_PATH'] = self.temp_db_path
        os.environ['USE_REMOTE_CADD'] = 'false'
        
    def tearDown(self):
        """Clean up test environment."""
        # Remove temporary database
        if os.path.exists(self.temp_db_path):
            os.unlink(self.temp_db_path)
    
    def test_normalize_chromosome(self):
        """Test chromosome normalization."""
        self.assertEqual(normalize_chromosome('1'), 'chr1')
        self.assertEqual(normalize_chromosome('chr1'), 'chr1')
        self.assertEqual(normalize_chromosome('X'), 'chrX')
        self.assertEqual(normalize_chromosome('chrMT'), 'chrMT')
    
    def test_lookup_cadd_score(self):
        """Test CADD score lookup functionality."""
        # Test successful lookup
        result = lookup_cadd_score('chr17', 43044295, 'A', 'T')
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result['raw'], 2.345, places=3)
        self.assertAlmostEqual(result['phred'], 24.7, places=1)
        
        # Test lookup with unnormalized chromosome
        result = lookup_cadd_score('17', 43044295, 'A', 'T')
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result['phred'], 24.7, places=1)
        
        # Test lookup miss
        result = lookup_cadd_score('chr1', 12345, 'G', 'C')
        self.assertIsNone(result)
    
    def test_calculate_risk_weight(self):
        """Test risk weight calculation from PHRED scores."""
        # Test benign (< 10)
        self.assertAlmostEqual(calculate_risk_weight(5.0), 0.1, places=2)
        self.assertAlmostEqual(calculate_risk_weight(9.9), 0.1, places=2)
        
        # Test uncertain (10-15)
        self.assertAlmostEqual(calculate_risk_weight(10.0), 0.1, places=2)
        self.assertAlmostEqual(calculate_risk_weight(12.5), 0.2, places=2)
        self.assertAlmostEqual(calculate_risk_weight(15.0), 0.3, places=2)
        
        # Test damaging (15-20)
        self.assertAlmostEqual(calculate_risk_weight(17.5), 0.45, places=2)
        self.assertAlmostEqual(calculate_risk_weight(20.0), 0.6, places=2)
        
        # Test pathogenic (20-25)
        self.assertAlmostEqual(calculate_risk_weight(22.5), 0.7, places=2)
        self.assertAlmostEqual(calculate_risk_weight(25.0), 0.8, places=2)
        
        # Test highly pathogenic (> 25)
        self.assertAlmostEqual(calculate_risk_weight(30.0), 0.9, places=2)
        self.assertAlmostEqual(calculate_risk_weight(50.0), 1.0, places=2)
    
    def test_process_node(self):
        """Test the main process function."""
        # Create test state
        state = {
            "current_node": None,
            "completed_nodes": [],
            "errors": [],
            "filtered_variants": [
                {
                    "variant_id": "chr17:43044295:A>T",
                    "chrom": "chr17",
                    "pos": 43044295,
                    "ref": "A",
                    "alt": "T",
                    "gene": "BRCA1"
                },
                {
                    "variant_id": "chr5:112839936:C>T",
                    "chrom": "chr5",
                    "pos": 112839936,
                    "ref": "C",
                    "alt": "T",
                    "gene": "APC"
                },
                {
                    "variant_id": "chr1:12345:G>C",
                    "chrom": "chr1",
                    "pos": 12345,
                    "ref": "G",
                    "alt": "C",
                    "gene": "UNKNOWN"
                }
            ]
        }
        
        # Process the state
        result = process(state)
        
        # Check that node was marked as complete
        self.assertIn("cadd_scoring", result["completed_nodes"])
        
        # Check enriched variants
        self.assertIn("cadd_enriched_variants", result)
        self.assertEqual(len(result["cadd_enriched_variants"]), 3)
        
        # Check first variant (should have CADD score)
        var1 = result["cadd_enriched_variants"][0]
        self.assertIn("cadd_phred", var1)
        self.assertAlmostEqual(var1["cadd_phred"], 24.7, places=1)
        self.assertIn("cadd_risk_weight", var1)
        self.assertTrue(0.6 <= var1["cadd_risk_weight"] <= 0.8)
        
        # Check second variant (benign)
        var2 = result["cadd_enriched_variants"][1]
        self.assertAlmostEqual(var2["cadd_phred"], 8.2, places=1)
        self.assertAlmostEqual(var2["cadd_risk_weight"], 0.1, places=2)
        
        # Check third variant (no CADD score)
        var3 = result["cadd_enriched_variants"][2]
        self.assertIsNone(var3["cadd_phred"])
        self.assertEqual(var3["cadd_risk_weight"], 0.2)  # Default
        
        # Check statistics
        self.assertIn("cadd_stats", result)
        stats = result["cadd_stats"]
        self.assertEqual(stats["total_variants"], 3)
        self.assertEqual(stats["variants_scored"], 2)
        self.assertEqual(stats["lookup_missing"], 1)
        self.assertEqual(stats["variants_gt20"], 1)
        self.assertAlmostEqual(stats["mean_phred"], (24.7 + 8.2) / 2, places=1)
        self.assertAlmostEqual(stats["max_phred"], 24.7, places=1)
    
    def test_error_handling(self):
        """Test error handling in the process function."""
        # Set invalid database path
        os.environ['CADD_DB_PATH'] = '/invalid/path/to/database.db'
        
        state = {
            "current_node": None,
            "completed_nodes": [],
            "errors": [],
            "filtered_variants": [
                {
                    "variant_id": "chr17:43044295:A>T",
                    "chrom": "chr17",
                    "pos": 43044295,
                    "ref": "A",
                    "alt": "T"
                }
            ]
        }
        
        # Process should not crash
        result = process(state)
        
        # Should still mark as complete (graceful failure)
        self.assertIn("cadd_scoring", result["completed_nodes"])
        
        # Should have default risk weight
        var = result["filtered_variants"][0]
        self.assertEqual(var.get("cadd_risk_weight", 0.2), 0.2)


def run_tests():
    """Run all tests and print results."""
    # Create test suite
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCADDScoring)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "="*70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success: {result.wasSuccessful()}")
    print("="*70)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1) 