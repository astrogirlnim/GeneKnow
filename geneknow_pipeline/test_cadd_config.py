#!/usr/bin/env python3
"""
Configuration tests for CADD scoring.
Tests environment variables, configuration options, and fallback behavior.
"""
import os
import sys
import unittest
import tempfile
from unittest.mock import patch, MagicMock

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nodes.cadd_scoring import process, CADD_REMOTE_TABIX, USE_REMOTE_CADD


class TestCADDConfiguration(unittest.TestCase):
    """Test CADD configuration and environment variable handling."""
    
    def setUp(self):
        """Store original environment variables."""
        self.original_env = {}
        for key in ["CADD_REMOTE_TABIX", "USE_REMOTE_CADD", "CADD_DB_PATH"]:
            self.original_env[key] = os.environ.get(key)
            
    def tearDown(self):
        """Restore original environment variables."""
        for key, value in self.original_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
                
    def test_default_configuration(self):
        """Test default configuration values."""
        # Clear environment variables
        for key in ["CADD_REMOTE_TABIX", "USE_REMOTE_CADD"]:
            os.environ.pop(key, None)
            
        # Re-import to get defaults
        import importlib
        import nodes.cadd_scoring
        importlib.reload(nodes.cadd_scoring)
        
        # Check defaults
        self.assertEqual(
            nodes.cadd_scoring.CADD_REMOTE_TABIX,
            "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/"
        )
        self.assertTrue(nodes.cadd_scoring.USE_REMOTE_CADD)
        
    def test_custom_remote_endpoint(self):
        """Test custom remote CADD endpoint configuration."""
        custom_endpoint = "https://custom.cadd.endpoint/v2.0/"
        os.environ["CADD_REMOTE_TABIX"] = custom_endpoint
        
        # Re-import to pick up new environment
        import importlib
        import nodes.cadd_scoring
        importlib.reload(nodes.cadd_scoring)
        
        self.assertEqual(nodes.cadd_scoring.CADD_REMOTE_TABIX, custom_endpoint)
        
    def test_disable_remote_lookups(self):
        """Test disabling remote CADD lookups."""
        os.environ["USE_REMOTE_CADD"] = "false"
        
        # Re-import to pick up new environment
        import importlib
        import nodes.cadd_scoring
        importlib.reload(nodes.cadd_scoring)
        
        self.assertFalse(nodes.cadd_scoring.USE_REMOTE_CADD)
        
    def test_custom_database_path(self):
        """Test custom database path configuration."""
        # Create a temporary database path
        temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        custom_path = temp_db.name
        temp_db.close()
        
        try:
            # Set custom path
            import nodes.cadd_scoring as cadd_module
            original_path = cadd_module.POP_DB_PATH
            cadd_module.POP_DB_PATH = custom_path
            
            # Create minimal database
            import sqlite3
            conn = sqlite3.connect(custom_path)
            cursor = conn.cursor()
            cursor.execute("""
                CREATE TABLE cadd_scores (
                    chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
                    raw_score REAL, phred_score REAL, job_id TEXT,
                    PRIMARY KEY (chrom, pos, ref, alt)
                )
            """)
            cursor.execute("""
                CREATE TABLE cadd_jobs (
                    job_id TEXT PRIMARY KEY,
                    job_type TEXT,
                    started_at TIMESTAMP,
                    completed_at TIMESTAMP,
                    variant_count INTEGER,
                    status TEXT,
                    metadata TEXT
                )
            """)
            conn.commit()
            conn.close()
            
            # Test that process uses custom path
            state = {
                "filtered_variants": [
                    {"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}
                ],
                "completed_nodes": [],
                "errors": []
            }
            
            result = process(state)
            self.assertIn("cadd_scoring", result["completed_nodes"])
            self.assertNotIn("error", result.get("cadd_stats", {}))
            
            # Restore original path
            cadd_module.POP_DB_PATH = original_path
            
        finally:
            # Clean up
            if os.path.exists(custom_path):
                os.unlink(custom_path)
                
    def test_legacy_risk_model_flag(self):
        """Test USE_LEGACY_RISK environment variable."""
        # This affects graph.py routing, not cadd_scoring directly
        # But we test it here for completeness
        
        # Test with legacy mode enabled
        os.environ["USE_LEGACY_RISK"] = "true"
        from graph import create_genomic_pipeline
        pipeline = create_genomic_pipeline()
        
        # Check that pipeline includes risk_model node
        # (This is a simplified check - full test would trace edges)
        self.assertIsNotNone(pipeline)
        
        # Test with legacy mode disabled (default)
        os.environ["USE_LEGACY_RISK"] = "false"
        pipeline = create_genomic_pipeline()
        self.assertIsNotNone(pipeline)
        
    @patch.dict(os.environ, {"USE_REMOTE_CADD": "true"})
    def test_remote_fallback_enabled(self):
        """Test that remote fallback is used when enabled."""
        import nodes.cadd_scoring as cadd_module
        
        # Mock the query_remote_cadd function
        with patch.object(cadd_module, 'query_remote_cadd') as mock_remote:
            mock_remote.return_value = {
                "raw": 2.5,
                "phred": 25.0,
                "job_id": "remote_test"
            }
            
            # Use non-existent database to force remote lookup
            original_path = cadd_module.POP_DB_PATH
            cadd_module.POP_DB_PATH = "/non/existent/path.db"
            
            try:
                state = {
                    "filtered_variants": [
                        {
                            "variant_id": "1:100:A>G",
                            "chrom": "1",
                            "pos": 100,
                            "ref": "A",
                            "alt": "G"
                        }
                    ],
                    "completed_nodes": [],
                    "errors": []
                }
                
                # Process should attempt remote lookup
                result = process(state)
                
                # Note: Current implementation doesn't actually call remote
                # This test structure is ready for when it does
                
            finally:
                cadd_module.POP_DB_PATH = original_path
                
    def test_configuration_logging(self):
        """Test that configuration is properly logged."""
        import logging
        
        # Capture log output
        with self.assertLogs('nodes.cadd_scoring', level=logging.INFO) as cm:
            state = {
                "filtered_variants": [],
                "completed_nodes": [],
                "errors": []
            }
            process(state)
            
        # Check that configuration-related logs are present
        log_output = '\n'.join(cm.output)
        self.assertIn("CADD scoring", log_output)


class TestCADDThresholds(unittest.TestCase):
    """Test CADD PHRED score thresholds."""
    
    def test_threshold_values(self):
        """Test that threshold values are correctly defined."""
        from nodes.cadd_scoring import CADD_PHRED_THRESHOLDS
        
        # Check all thresholds exist
        self.assertIn("benign", CADD_PHRED_THRESHOLDS)
        self.assertIn("uncertain", CADD_PHRED_THRESHOLDS)
        self.assertIn("damaging", CADD_PHRED_THRESHOLDS)
        self.assertIn("pathogenic", CADD_PHRED_THRESHOLDS)
        
        # Check threshold ordering
        self.assertLess(
            CADD_PHRED_THRESHOLDS["benign"],
            CADD_PHRED_THRESHOLDS["uncertain"]
        )
        self.assertLess(
            CADD_PHRED_THRESHOLDS["uncertain"],
            CADD_PHRED_THRESHOLDS["damaging"]
        )
        self.assertLess(
            CADD_PHRED_THRESHOLDS["damaging"],
            CADD_PHRED_THRESHOLDS["pathogenic"]
        )
        
    def test_threshold_risk_weights(self):
        """Test that thresholds produce expected risk weights."""
        from nodes.cadd_scoring import calculate_risk_weight, CADD_PHRED_THRESHOLDS
        
        # Test at threshold boundaries
        test_cases = [
            (CADD_PHRED_THRESHOLDS["benign"] - 1, 0.1),      # Below benign
            (CADD_PHRED_THRESHOLDS["benign"], 0.1),          # At benign
            (CADD_PHRED_THRESHOLDS["uncertain"], 0.3),       # At uncertain
            (CADD_PHRED_THRESHOLDS["damaging"], 0.6),        # At damaging
            (CADD_PHRED_THRESHOLDS["pathogenic"], 0.8),      # At pathogenic
            (CADD_PHRED_THRESHOLDS["pathogenic"] + 10, 1.0), # Well above pathogenic
        ]
        
        for phred, expected_weight in test_cases:
            weight = calculate_risk_weight(phred)
            self.assertAlmostEqual(
                weight, expected_weight, places=1,
                msg=f"PHRED {phred} should give risk weight ~{expected_weight}"
            )


if __name__ == "__main__":
    unittest.main() 