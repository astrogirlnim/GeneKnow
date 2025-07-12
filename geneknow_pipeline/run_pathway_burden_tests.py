#!/usr/bin/env python3
"""
Comprehensive test runner for the Pathway Burden Model.
Executes all types of tests: unit, integration, and performance.
"""

import os
import sys
import subprocess
import time
from datetime import datetime

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))


def run_command(cmd, description):
    """Run a command and capture its output."""
    print(f"\n{'='*60}")
    print(f"🔧 {description}")
    print(f"{'='*60}")
    print(f"Command: {cmd}")
    print(f"Started: {datetime.now().strftime('%H:%M:%S')}")

    start_time = time.time()
    try:
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=300
        )  # 5 minute timeout

        end_time = time.time()
        duration = end_time - start_time

        print(f"Duration: {duration:.2f} seconds")
        print(f"Return code: {result.returncode}")

        if result.stdout:
            print("\n📋 STDOUT:")
            print(result.stdout)

        if result.stderr:
            print("\n⚠️  STDERR:")
            print(result.stderr)

        success = result.returncode == 0
        print(f"\n{'✅ SUCCESS' if success else '❌ FAILED'}")

        return success, duration, result.stdout, result.stderr

    except subprocess.TimeoutExpired:
        print("\n⏰ TIMEOUT: Test exceeded 5 minutes")
        return False, 300, "", "Timeout"
    except Exception as e:
        print(f"\n❌ ERROR: {str(e)}")
        return False, 0, "", str(e)


def check_prerequisites():
    """Check if all required components are available."""

    print("🔍 Checking Prerequisites")
    print("=" * 60)

    checks = []

    # Check if pathway_burden.py exists
    pathway_burden_file = "nodes/pathway_burden.py"
    if os.path.exists(pathway_burden_file):
        print(f"✅ {pathway_burden_file} exists")
        checks.append(True)
    else:
        print(f"❌ {pathway_burden_file} missing")
        checks.append(False)

    # Check if test file exists
    test_file = "test_pathway_burden.py"
    if os.path.exists(test_file):
        print(f"✅ {test_file} exists")
        checks.append(True)
    else:
        print(f"❌ {test_file} missing")
        checks.append(False)

    # Check if graph.py has pathway_burden
    try:
        with open("graph.py", "r") as f:
            content = f.read()
            if "pathway_burden" in content:
                print("✅ graph.py includes pathway_burden")
                checks.append(True)
            else:
                print("❌ graph.py missing pathway_burden")
                checks.append(False)
    except Exception:
        print("❌ graph.py not readable")
        checks.append(False)

    # Check if __init__.py has pathway_burden
    try:
        with open("nodes/__init__.py", "r") as f:
            content = f.read()
            if "pathway_burden" in content:
                print("✅ nodes/__init__.py includes pathway_burden")
                checks.append(True)
            else:
                print("❌ nodes/__init__.py missing pathway_burden")
                checks.append(False)
    except:
        print("❌ nodes/__init__.py not readable")
        checks.append(False)

    all_good = all(checks)
    print(
        f"\n{'✅ All prerequisites met' if all_good else '❌ Some prerequisites missing'}"
    )

    return all_good


def run_all_tests():
    """Run all pathway burden tests."""

    print("🧬 Pathway Burden Model - Complete Test Suite")
    print("=" * 70)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Check prerequisites
    if not check_prerequisites():
        print("\n❌ Prerequisites not met. Cannot run tests.")
        return False

    # Test results
    results = {}
    total_duration = 0

    # Run comprehensive test
    success, duration, stdout, stderr = run_command(
        "python test_pathway_burden.py", "Comprehensive Pathway Burden Tests"
    )
    results["pathway_burden_tests"] = success
    total_duration += duration

    # Quick Pipeline Test (if available)
    if os.path.exists("test_pipeline_comprehensive.py"):
        success, duration, stdout, stderr = run_command(
            "python test_pipeline_comprehensive.py", "Comprehensive Pipeline Test"
        )
        results["pipeline_test"] = success
        total_duration += duration

    # Import Test
    success, duration, stdout, stderr = run_command(
        "python -c \"from nodes.pathway_burden import process, CANCER_PATHWAYS; print('Import successful'); print(f'Pathways: {len(CANCER_PATHWAYS)}')\"",
        "Import Test - Module Loading",
    )
    results["import_test"] = success
    total_duration += duration

    # Syntax Check
    success, duration, stdout, stderr = run_command(
        "python -m py_compile nodes/pathway_burden.py",
        "Syntax Check - Code Compilation",
    )
    results["syntax_check"] = success
    total_duration += duration

    # Generate summary
    print("\n" + "=" * 70)
    print("🏆 TEST RESULTS SUMMARY")
    print("=" * 70)

    passed = sum(1 for success in results.values() if success)
    total = len(results)

    print(f"📊 Overall: {passed}/{total} tests passed")
    print(f"⏱️  Total duration: {total_duration:.2f} seconds")

    for test_name, success in results.items():
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"  {test_name}: {status}")

    # Overall success
    overall_success = all(results.values())

    if overall_success:
        print("\n🎉 ALL TESTS PASSED!")
        print("The Pathway Burden Model is fully functional and ready for use.")
    else:
        print("\n❌ Some tests failed.")
        print("Please check the detailed output above for specific issues.")

    # Generate test report
    report_file = "pathway_burden_test_report.txt"
    with open(report_file, "w") as f:
        f.write("Pathway Burden Model Test Report\n")
        f.write("=" * 40 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Tests: {total}\n")
        f.write(f"Passed: {passed}\n")
        f.write(f"Failed: {total - passed}\n")
        f.write(f"Duration: {total_duration:.2f} seconds\n")
        f.write(f"Overall: {'PASS' if overall_success else 'FAIL'}\n\n")

        f.write("Detailed Results:\n")
        for test_name, success in results.items():
            f.write(f"  {test_name}: {'PASS' if success else 'FAIL'}\n")

    print(f"\n📄 Test report saved to: {report_file}")

    return overall_success


def quick_test():
    """Run a quick smoke test."""

    print("🚀 Quick Smoke Test")
    print("=" * 40)

    try:
        # Test import
        from nodes.pathway_burden import process, CANCER_PATHWAYS, is_damaging_variant

        print("✅ Module import: SUCCESS")
        print(f"  - Pathways defined: {len(CANCER_PATHWAYS)}")
        print("  - Functions available: process, is_damaging_variant, etc.")

        # Test basic functionality
        test_variant = {
            "variant_id": "test:123:A>G",
            "gene": "BRCA1",
            "cadd_phred": 25,
            "clinical_significance": "Pathogenic",
            "allele_frequency": 0.001,
            "consequence": "missense_variant",
        }

        damage_result = is_damaging_variant(test_variant)
        print("✅ Damage assessment: SUCCESS")
        print(f"  - Damage score: {damage_result['damage_score']:.3f}")
        print(f"  - Is damaging: {damage_result['is_damaging']}")

        # Test with minimal state
        test_state = {"filtered_variants": [test_variant]}
        result = process(test_state)

        print("✅ Node processing: SUCCESS")
        print(f"  - Pathways analyzed: {len(result.get('pathway_burden_results', {}))}")
        print(
            f"  - Overall burden score: {result.get('pathway_burden_summary', {}).get('overall_burden_score', 0):.3f}"
        )

        print("\n🎉 Quick test PASSED!")
        return True

    except Exception as e:
        print(f"\n❌ Quick test FAILED: {str(e)}")
        import traceback

        traceback.print_exc()
        return False


def main():
    """Main function."""

    if len(sys.argv) > 1:
        if sys.argv[1] == "quick":
            success = quick_test()
        elif sys.argv[1] == "prereq":
            success = check_prerequisites()
        else:
            print("Usage: python run_pathway_burden_tests.py [quick|prereq]")
            print("  quick  - Run quick smoke test")
            print("  prereq - Check prerequisites only")
            print("  (no args) - Run full test suite")
            return
    else:
        success = run_all_tests()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
