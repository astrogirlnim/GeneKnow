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
    process,
    create_job_record,
    update_job_status
)


class TestCADDScoring(unittest.TestCase):
    """Test suite for CADD scoring functionality."""
    
    def setUp(self):
        """Set up test environment."""
        # Create temporary database for testing
        self.temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        self.temp_db_path = self.temp_db.name
        self.temp_db.close()
        
        # Create test database with both population_variants and cadd_scores tables
        conn = sqlite3.connect(self.temp_db_path)
        cursor = conn.cursor()
        
        # Create population_variants table (simplified schema)
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
        
        # Insert some test population variants
        test_pop_variants = [
            ('17', 43044295, 'A', 'T', 'BRCA1', 0.001, 'Pathogenic', 1, 'missense_variant', 'reviewed'),
            ('17', 43045719, 'G', 'A', 'BRCA1', 0.0001, 'Pathogenic', 1, 'nonsense_variant', 'reviewed'),
            ('13', 32315474, 'G', 'A', 'BRCA2', 0.005, 'Uncertain', 0, 'missense_variant', 'reviewed'),
            ('5', 112839936, 'C', 'T', 'APC', 0.01, 'Benign', 0, 'synonymous_variant', 'reviewed'),
        ]
        cursor.executemany(
            "INSERT INTO population_variants VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            test_pop_variants
        )
        
        # Create cadd_scores table
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
                PRIMARY KEY (chrom, pos, ref, alt),
                FOREIGN KEY (chrom, pos, ref, alt) 
                    REFERENCES population_variants(chrom, pos, ref, alt)
            )
        """)
        
        # Create cadd_jobs table
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
        
        # Insert test CADD data
        test_cadd_scores = [
            ('17', 43044295, 'A', 'T', 2.345, 24.7, 'test_job_001'),  # BRCA1 - pathogenic
            ('17', 43045719, 'G', 'A', 3.456, 34.1, 'test_job_001'),  # BRCA1 - highly pathogenic
            ('13', 32315474, 'G', 'A', 1.234, 15.3, 'test_job_001'),  # BRCA2 - uncertain
            ('5', 112839936, 'C', 'T', 0.456, 8.2, 'test_job_001'),   # APC - benign
            ('7', 55019017, 'T', 'G', 2.890, 28.5, 'test_job_001'),   # EGFR - pathogenic
        ]
        
        cursor.executemany(
            "INSERT INTO cadd_scores (chrom, pos, ref, alt, raw_score, phred_score, job_id) VALUES (?, ?, ?, ?, ?, ?, ?)",
            test_cadd_scores
        )
        
        conn.commit()
        conn.close()
        
        # Set environment variable to use test database
        self.original_db_path = os.environ.get("CADD_DB_PATH")
        os.environ["CADD_DB_PATH"] = self.temp_db_path
        
    def tearDown(self):
        """Clean up test environment."""
        # Remove temporary database
        if os.path.exists(self.temp_db_path):
            os.unlink(self.temp_db_path)
            
        # Restore original environment variable
        if self.original_db_path:
            os.environ["CADD_DB_PATH"] = self.original_db_path
        elif "CADD_DB_PATH" in os.environ:
            del os.environ["CADD_DB_PATH"]
            
    def test_normalize_chromosome(self):
        """Test chromosome normalization."""
        self.assertEqual(normalize_chromosome("chr1"), "1")
        self.assertEqual(normalize_chromosome("chrX"), "X")
        self.assertEqual(normalize_chromosome("1"), "1")
        self.assertEqual(normalize_chromosome("X"), "X")
        
    def test_lookup_cadd_score(self):
        """Test CADD score lookup from database."""
        conn = sqlite3.connect(self.temp_db_path)
        
        # Test existing variant - note we use normalized chromosome
        result = lookup_cadd_score("chr17", 43044295, "A", "T", conn)
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result["raw"], 2.345, places=3)
        self.assertAlmostEqual(result["phred"], 24.7, places=1)
        self.assertEqual(result["job_id"], "test_job_001")
        
        # Test non-existing variant
        result = lookup_cadd_score("chr1", 1000000, "G", "A", conn)
        self.assertIsNone(result)
        
        conn.close()
        
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
        self.assertAlmostEqual(calculate_risk_weight(40.0), 1.0, places=2)
        
    def test_job_tracking(self):
        """Test job creation and update functions."""
        conn = sqlite3.connect(self.temp_db_path)
        
        # Create job
        job_id = create_job_record(conn, "unit_test")
        self.assertTrue(job_id.startswith("cadd_"))
        
        # Check job was created
        cursor = conn.cursor()
        cursor.execute("SELECT job_type, status FROM cadd_jobs WHERE job_id = ?", (job_id,))
        row = cursor.fetchone()
        self.assertEqual(row[0], "unit_test")
        self.assertEqual(row[1], "running")
        
        # Update job status
        update_job_status(conn, job_id, "completed", 10)
        
        # Check job was updated
        cursor.execute("SELECT status, variant_count FROM cadd_jobs WHERE job_id = ?", (job_id,))
        row = cursor.fetchone()
        self.assertEqual(row[0], "completed")
        self.assertEqual(row[1], 10)
        
        conn.close()
        
    def test_process_function(self):
        """Test the main process function."""
        # Override the database path for this test
        import nodes.cadd_scoring as cadd_module
        original_path = cadd_module.POP_DB_PATH
        cadd_module.POP_DB_PATH = self.temp_db_path
        
        # Create test state
        state = {
            "filtered_variants": [
                {
                    "variant_id": "chr17:43044295:A>T",
                    "chrom": "chr17",
                    "pos": 43044295,
                    "ref": "A",
                    "alt": "T",
                    "gene": "BRCA1",
                    "risk_weight": 0.5
                },
                {
                    "variant_id": "chr1:1000000:G>A",
                    "chrom": "chr1",
                    "pos": 1000000,
                    "ref": "G",
                    "alt": "A",
                    "gene": "TEST1"
                }
            ],
            "benign_variants": [],
            "pathogenic_variants": ["chr17:43044295:A>T"],
            "risk_genes": {"breast": ["BRCA1", "BRCA2"]},
            "completed_nodes": [],
            "errors": []
        }
        
        # Process variants
        result = process(state)
        
        # Check results
        self.assertIn("cadd_enriched_variants", result)
        self.assertIn("cadd_stats", result)
        self.assertIn("cadd_scoring", result["completed_nodes"])
        
        # Check enriched variants
        enriched = result["cadd_enriched_variants"]
        self.assertEqual(len(enriched), 2)
        
        # First variant should have CADD scores
        var1 = enriched[0]
        self.assertIn("cadd_phred", var1)
        self.assertAlmostEqual(var1["cadd_phred"], 24.7, places=1)
        self.assertIn("cadd_risk_weight", var1)
        self.assertGreater(var1["risk_weight"], 0.5)  # Should be max of original and CADD
        
        # Second variant should not have CADD scores
        var2 = enriched[1]
        self.assertNotIn("cadd_phred", var2)
        
        # Check statistics
        stats = result["cadd_stats"]
        self.assertEqual(stats["total_variants"], 2)
        self.assertEqual(stats["variants_scored"], 1)
        self.assertEqual(stats["lookup_missing"], 1)
        self.assertEqual(stats["variants_in_cancer_genes"], 1)
        self.assertAlmostEqual(stats["mean_phred"], 24.7, places=1)
        self.assertAlmostEqual(stats["max_phred"], 24.7, places=1)
        self.assertEqual(stats["variants_gt20"], 1)
        
        # Restore original path
        cadd_module.POP_DB_PATH = original_path


if __name__ == "__main__":
    unittest.main() 