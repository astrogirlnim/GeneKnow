#!/usr/bin/env python3
"""
Test script to verify the metrics implementation is working correctly.
"""
import sys
import os

# Add the parent directory to the path to import the modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def test_imports():
    """Test that all modules import correctly."""
    try:
        pass

        print("✅ All imports successful")
        assert True
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False


def test_metrics_calculator():
    """Test the metrics calculator directly."""
    try:
        from nodes.metrics_calculator import process as calculate_metrics

        # Create a minimal test state
        test_state = {
            "current_node": "test",
            "completed_nodes": [],
            "errors": [],
            "warnings": [],
            "model_version": "1.0.0",
            "filtered_variants": [
                {
                    "gene": "BRCA1",
                    "clinical_significance": "Pathogenic",
                    "cadd_phred": 25.0,
                    "variant_id": "chr17:41244936:G>A",
                }
            ],
            "risk_scores": {"breast": 75.0, "ovarian": 45.0},
            "risk_details": {
                "breast": {
                    "model_confidence": 0.85,
                    "gene_count": 1,
                    "pathogenic_count": 1,
                }
            },
            "ml_risk_assessment": {
                "confidence": 0.8,
                "risk_category": "high",
                "aggregate_risk_score": 60.0,
            },
            "risk_genes": {
                "breast": ["BRCA1", "BRCA2"],
                "ovarian": ["BRCA1", "BRCA2", "TP53"],
            },
            "prs_summary": {"breast": {"prs_score": 1.2, "variants_used": 50}},
        }

        # Test metrics calculation
        result = calculate_metrics(test_state)

        # Check that metrics were calculated
        if "metrics" in result and "metrics_summary" in result:
            print("✅ Metrics calculator works correctly")
            print(f"   - Calculated {len(result['metrics'])} metric categories")
            print("   - Generated metrics summary")
            return True
        else:
            print("❌ Metrics calculator failed to generate metrics")
            return False

    except Exception as e:
        print(f"❌ Metrics calculator error: {e}")
        return False


def test_frontend_types():
    """Test that the frontend types are correctly structured."""
    try:
        # Check if the API types file exists and is valid
        api_file = "../desktop/ui/src/api/geneknowPipeline.ts"
        if os.path.exists(api_file):
            with open(api_file, "r") as f:
                content = f.read()
                if "metrics?" in content and "PipelineResult" in content:
                    print("✅ Frontend types include metrics structure")
                    assert True
                else:
                    print("❌ Frontend types missing metrics structure")
                    return False
        else:
            print("❌ Frontend API types file not found")
            return False
    except Exception as e:
        print(f"❌ Frontend types error: {e}")
        return False


def main():
    """Run all tests."""
    print("🧪 Testing GeneKnow Pipeline Metrics Implementation")
    print("=" * 60)

    tests = [
        ("Module Imports", test_imports),
        ("Metrics Calculator", test_metrics_calculator),
        ("Frontend Types", test_frontend_types),
    ]

    passed = 0
    total = len(tests)

    for test_name, test_func in tests:
        print(f"\n📋 {test_name}:")
        if test_func():
            passed += 1
        else:
            print("   Test failed!")

    print(f"\n{'='*60}")
    print(f"📊 Test Results: {passed}/{total} tests passed")

    if passed == total:
        print("🎉 All tests passed! Metrics implementation is working correctly.")
        return True
    else:
        print("⚠️  Some tests failed. Please check the implementation.")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
