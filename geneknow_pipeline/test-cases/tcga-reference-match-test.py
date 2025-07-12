#!/usr/bin/env python3
"""
TCGA Reference Match Test Suite
Tests the TCGA frequency matching functionality to ensure variants are correctly
matched against tumor frequency data.
"""

import sys
import logging
from pathlib import Path

# Add the parent directory to Python path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


# Import the TCGA mapper directly to avoid dependency issues
def import_tcga_mapper():
    """Import TCGA mapper without triggering other node dependencies."""
    import importlib.util

    # Direct import of tcga_mapper module
    tcga_mapper_path = Path(__file__).parent.parent / "nodes" / "tcga_mapper.py"
    spec = importlib.util.spec_from_file_location("tcga_mapper", tcga_mapper_path)
    tcga_mapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tcga_mapper)

    return tcga_mapper.process


tcga_mapper_process = import_tcga_mapper()


def run_tcga_frequency_analysis(state):
    """
    Wrapper function to run TCGA frequency analysis and format results
    for the test expectations.

    This function calls our TCGA mapper and translates the output format
    to match what the tests expect (tcga_frequency/tcga_patient_count fields).
    """
    # Ensure state has required structure
    if "tcga_matches" not in state:
        state["tcga_matches"] = {}
    if "tcga_cohort_sizes" not in state:
        state["tcga_cohort_sizes"] = {}
    if "file_metadata" not in state:
        state["file_metadata"] = {}
    if "completed_nodes" not in state:
        state["completed_nodes"] = []
    if "errors" not in state:
        state["errors"] = []
    if "warnings" not in state:
        state["warnings"] = []

    # Run the TCGA mapper
    updated_state = tcga_mapper_process(state)

    # Translate TCGA mapper output to test-expected format
    for variant in updated_state["filtered_variants"]:
        variant_id = variant.get("variant_id")

        # Find the best TCGA match for this variant
        best_frequency = 0
        best_patient_count = 0

        for cancer_type, matches in updated_state.get("tcga_matches", {}).items():
            if variant_id in matches:
                match_data = matches[variant_id]
                tumor_freq = match_data.get("tumor_frequency", 0)
                sample_count = match_data.get("sample_count", 0)

                # Keep the highest frequency found across all cancer types
                if tumor_freq > best_frequency:
                    best_frequency = tumor_freq
                    best_patient_count = sample_count

        # Add the expected fields to the variant
        variant["tcga_frequency"] = best_frequency
        variant["tcga_patient_count"] = best_patient_count

    return updated_state


def test_tcga_frequency_analysis_known_variant():
    """Test TCGA frequency analysis for a known variant that should be in the dataset."""
    # Use a real variant from our comprehensive TCGA dataset
    state = {
        "filtered_variants": [
            {
                "chrom": "17",
                "pos": 43111347,
                "re": "T",
                "alt": "G",
                "variant_id": "17:43111347:T>G",
                "gene": "BRCA1",
            }
        ]
    }

    # Call the function that performs TCGA lookup and modifies the variant
    updated_state = run_tcga_frequency_analysis(state)

    variant = updated_state["filtered_variants"][0]

    # Assert that a TCGA frequency value was added
    assert (
        "tcga_frequency" in variant or "tcga_patient_count" in variant
    ), "TCGA frequency data not added to variant"

    # Optional: assert non-zero if you know the test variant exists
    freq = variant.get("tcga_frequency", 0) or variant.get("tcga_patient_count", 0)
    assert freq > 0, "Expected known variant to have non-zero TCGA frequency"


def test_tcga_frequency_analysis_unknown_variant():
    """Test TCGA frequency analysis for an unknown variant."""
    state = {
        "filtered_variants": [
            {
                "chrom": "1",
                "pos": 12345678,
                "re": "T",
                "alt": "C",
                "variant_id": "1:12345678:T>C",
                "gene": "GENE_UNKNOWN",
            }
        ]
    }

    updated_state = run_tcga_frequency_analysis(state)

    variant = updated_state["filtered_variants"][0]
    assert "tcga_frequency" in variant or "tcga_patient_count" in variant

    # If the variant is unknown, frequency should be zero
    freq = variant.get("tcga_frequency", 0) or variant.get("tcga_patient_count", 0)
    assert freq == 0, "Unknown variant should have zero TCGA frequency"


def test_tcga_multiple_variants():
    """Test TCGA frequency analysis with multiple variants."""
    state = {
        "filtered_variants": [
            {
                "chrom": "17",
                "pos": 43111347,
                "re": "T",
                "alt": "G",
                "variant_id": "17:43111347:T>G",
                "gene": "BRCA1",
            },
            {
                "chrom": "17",
                "pos": 7675897,
                "re": "G",
                "alt": "A",
                "variant_id": "17:7675897:G>A",
                "gene": "TP53",
            },
            {
                "chrom": "12",
                "pos": 25210833,
                "re": "G",
                "alt": "A",
                "variant_id": "12:25210833:G>A",
                "gene": "KRAS",
            },
        ]
    }

    updated_state = run_tcga_frequency_analysis(state)

    # Check that all variants got frequency data
    for i, variant in enumerate(updated_state["filtered_variants"]):
        print(variant)
        assert "tcga_frequency" in variant, f"Variant {i} missing tcga_frequency"
        assert (
            "tcga_patient_count" in variant
        ), f"Variant {i} missing tcga_patient_count"


def test_tcga_state_structure():
    """Test that TCGA analysis creates the expected state structure."""
    state = {
        "filtered_variants": [
            {
                "chrom": "17",
                "pos": 43111347,
                "re": "T",
                "alt": "G",
                "variant_id": "17:43111347:T>G",
                "gene": "BRCA1",
            }
        ]
    }

    updated_state = run_tcga_frequency_analysis(state)

    # Check that the expected state keys were added
    assert "tcga_matches" in updated_state, "tcga_matches not added to state"
    assert "tcga_cohort_sizes" in updated_state, "tcga_cohort_sizes not added to state"

    # Check that cohort sizes are reasonable
    cohort_sizes = updated_state["tcga_cohort_sizes"]
    assert isinstance(cohort_sizes, dict), "tcga_cohort_sizes should be a dictionary"

    # Check for expected cancer types
    expected_cancer_types = ["breast", "colon", "lung", "prostate", "blood"]
    for cancer_type in expected_cancer_types:
        assert cancer_type in cohort_sizes, f"Missing cancer type: {cancer_type}"
        assert cohort_sizes[cancer_type] > 0, f"Invalid cohort size for {cancer_type}"


def run_all_tests():
    """Run all test functions and report results."""
    test_functions = [
        test_tcga_frequency_analysis_known_variant,
        test_tcga_frequency_analysis_unknown_variant,
        test_tcga_multiple_variants,
        test_tcga_state_structure,
    ]

    print("🧪 Running TCGA Reference Match Test Suite")
    print("=" * 60)

    passed = 0
    failed = 0

    for test_func in test_functions:
        try:
            print(f"Running {test_func.__name__}...", end=" ")
            test_func()
            print("✅ PASSED")
            passed += 1
        except Exception as e:
            print(f"❌ FAILED: {e}")
            failed += 1

    print(f"\n📊 Test Results: {passed} passed, {failed} failed")

    if failed == 0:
        print("🎉 All tests passed!")
        return True
    else:
        print(f"⚠️ {failed} tests failed.")
        return False


if __name__ == "__main__":
    # Configure logging to reduce noise during testing
    logging.basicConfig(level=logging.WARNING)

    # Check if unified database exists
    db_path = Path(__file__).parent.parent / "population_variants.db"
    if not db_path.exists():
        print("❌ Main database not found!")
        print(f"Expected location: {db_path}")
        print("Run database setup first.")
        sys.exit(1)

    # Check if TCGA table exists in the database
    import sqlite3

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='tcga_variants'"
    )
    if not cursor.fetchone():
        print("❌ TCGA table not found in database!")
        print("Run 'python setup_tcga.py' first to create the TCGA table.")
        sys.exit(1)
    conn.close()

    # Run the tests
    success = run_all_tests()
    sys.exit(0 if success else 1)
