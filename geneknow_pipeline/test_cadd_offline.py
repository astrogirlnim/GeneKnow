#!/usr/bin/env python3
"""
Test suite for offline CADD scoring implementation.
Tests the local CADD-like scoring algorithm.
"""
import unittest
import sys
import os
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import CADD module
from nodes.cadd_scoring import (
    calculate_risk_weight,
    compute_cadd_score,
    process
)


class TestOfflineCADDScoring(unittest.TestCase):
    """Test offline CADD scoring functionality."""
    
    def test_compute_cadd_score_frameshift(self):
        """Test scoring of frameshift variants."""
        variant = {
            "consequence": "frameshift_variant",
            "gene": "TP53",
            "allele_frequency": 0.001,
            "quality": 100,
            "depth": 50
        }
        
        result = compute_cadd_score(variant)
        self.assertIn("raw", result)
        self.assertIn("phred", result)
        # Frameshift in TP53 should have high score
        self.assertGreater(result["phred"], 25.0)
        
    def test_compute_cadd_score_synonymous(self):
        """Test scoring of synonymous variants."""
        variant = {
            "consequence": "synonymous_variant",
            "gene": "UNKNOWN",
            "allele_frequency": 0.1,
            "quality": 100,
            "depth": 30
        }
        
        result = compute_cadd_score(variant)
        # Synonymous should have low score
        self.assertLess(result["phred"], 10.0)
        
    def test_compute_cadd_score_clinical_significance(self):
        """Test that clinical significance affects scores."""
        # Pathogenic variant
        variant_path = {
            "consequence": "missense_variant",
            "gene": "BRCA1",
            "clinical_significance": "Pathogenic",
            "allele_frequency": 0.01,
            "quality": 100,
            "depth": 40
        }
        
        # Benign variant
        variant_benign = {
            "consequence": "missense_variant",
            "gene": "BRCA1",
            "clinical_significance": "Benign",
            "allele_frequency": 0.01,
            "quality": 100,
            "depth": 40
        }
        
        result_path = compute_cadd_score(variant_path)
        result_benign = compute_cadd_score(variant_benign)
        
        # Pathogenic should score higher than benign
        self.assertGreater(result_path["phred"], 20.0)
        self.assertLessEqual(result_benign["phred"], 10.0)
        
    def test_calculate_risk_weight(self):
        """Test risk weight calculation from PHRED scores."""
        # Test different PHRED score ranges
        self.assertAlmostEqual(calculate_risk_weight(5.0), 0.1, places=2)
        self.assertAlmostEqual(calculate_risk_weight(12.5), 0.2, places=2)
        self.assertAlmostEqual(calculate_risk_weight(17.5), 0.45, places=2)
        self.assertAlmostEqual(calculate_risk_weight(22.5), 0.7, places=2)
        self.assertAlmostEqual(calculate_risk_weight(30.0), 0.9, places=2)
        
    def test_process_function(self):
        """Test the main process function."""
        state = {
            "filtered_variants": [
                {
                    "chrom": "17",
                    "pos": 7577121,
                    "ref": "G",
                    "alt": "A",
                    "gene": "TP53",
                    "consequence": "missense_variant",
                    "allele_frequency": 0.001,
                    "quality": 100,
                    "depth": 50
                },
                {
                    "chrom": "13",
                    "pos": 32914437,
                    "ref": "T",
                    "alt": "C",
                    "gene": "BRCA2",
                    "consequence": "synonymous_variant",
                    "allele_frequency": 0.05,
                    "quality": 100,
                    "depth": 40
                }
            ],
            "risk_genes": {
                "breast": ["TP53", "BRCA2"]
            },
            "completed_nodes": [],
            "errors": []
        }
        
        # Process variants
        result_state = process(state)
        
        # Check results
        self.assertIn("cadd_enriched_variants", result_state)
        self.assertIn("cadd_stats", result_state)
        self.assertEqual(len(result_state["cadd_enriched_variants"]), 2)
        
        # Check stats
        stats = result_state["cadd_stats"]
        self.assertEqual(stats["variants_scored"], 2)
        self.assertEqual(stats["scoring_method"], "offline_algorithm")
        self.assertGreater(stats["mean_phred"], 0)
        self.assertGreater(stats["max_phred"], 0)
        
        # Check that variants have CADD scores
        for variant in result_state["cadd_enriched_variants"]:
            self.assertIn("cadd_phred", variant)
            self.assertIn("cadd_raw", variant)
            self.assertIn("cadd_risk_weight", variant)
            self.assertEqual(variant["cadd_source"], "computed_offline")
            
    def test_rare_variant_scoring(self):
        """Test that rare variants score higher."""
        # Common variant
        common_variant = {
            "consequence": "missense_variant",
            "gene": "KRAS",
            "allele_frequency": 0.5,
            "quality": 100,
            "depth": 50
        }
        
        # Rare variant (same otherwise)
        rare_variant = {
            "consequence": "missense_variant",
            "gene": "KRAS",
            "allele_frequency": 0.0001,
            "quality": 100,
            "depth": 50
        }
        
        common_score = compute_cadd_score(common_variant)
        rare_score = compute_cadd_score(rare_variant)
        
        # Rare variant should score higher
        self.assertGreater(rare_score["phred"], common_score["phred"])
        
    def test_cancer_gene_detection(self):
        """Test that cancer genes are properly tracked."""
        state = {
            "filtered_variants": [
                {
                    "chrom": "12",
                    "pos": 25398285,
                    "ref": "C",
                    "alt": "T",
                    "gene": "KRAS",
                    "consequence": "missense_variant",
                    "allele_frequency": 0.01,
                    "quality": 100,
                    "depth": 60
                },
                {
                    "chrom": "1",
                    "pos": 1000000,
                    "ref": "A",
                    "alt": "G",
                    "gene": "UNKNOWN_GENE",
                    "consequence": "missense_variant",
                    "allele_frequency": 0.01,
                    "quality": 100,
                    "depth": 60
                }
            ],
            "risk_genes": {
                "lung": ["KRAS", "EGFR"],
                "colon": ["KRAS", "APC"]
            },
            "completed_nodes": [],
            "errors": []
        }
        
        result_state = process(state)
        stats = result_state["cadd_stats"]
        
        # Should detect KRAS as a cancer gene
        self.assertEqual(stats["variants_in_cancer_genes"], 1)


if __name__ == "__main__":
    unittest.main() 