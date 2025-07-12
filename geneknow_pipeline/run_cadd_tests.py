#!/usr/bin/env python3
"""
CADD Test Runner - Executes all CADD tests and generates comprehensive reports.
Includes unit tests, integration tests, performance benchmarks, and validation.
"""
import os
import sys
import json
import time
import subprocess
from datetime import datetime
import logging
import traceback

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CADDTestRunner:
    """Comprehensive test runner for CADD scoring implementation."""

    def __init__(self):
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "tests": {},
            "summary": {
                "total_tests": 0,
                "passed": 0,
                "failed": 0,
                "errors": 0,
                "performance_metrics": {},
                "coverage": {}
            }
        }

    def run_all_tests(self):
        """Execute all CADD test suites."""
        print("=" * 80)
        print("🧬 CADD SCORING COMPREHENSIVE TEST SUITE")
        print("=" * 80)
        print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print()

        # Now we only need to run the comprehensive test file
        self.run_test_file("Comprehensive CADD Tests", "test_cadd_comprehensive.py")

        # Run performance benchmarks
        self.run_performance_benchmarks()

        # Run database validation
        self.validate_database_schema()

        # Generate report
        self.generate_report()

    def run_test_file(self, test_name, test_file):
        """Run a specific test file and capture results."""
        print(f"\n{'='*60}")
        print(f"Running {test_name}: {test_file}")
        print(f"{'='*60}")

        start_time = time.time()

        try:
            # Run the test file
            result = subprocess.run(
                [sys.executable, test_file],
                capture_output=True,
                text=True,
                cwd=os.path.dirname(os.path.abspath(__file__))
            )

            end_time = time.time()
            duration = end_time - start_time

            # Parse results
            success = result.returncode == 0
            output = result.stdout
            errors = result.stderr

            # Store results
            self.results["tests"][test_name] = {
                "file": test_file,
                "success": success,
                "duration": duration,
                "return_code": result.returncode,
                "output": output,
                "errors": errors
            }

            # Update summary
            self.results["summary"]["total_tests"] += 1
            if success:
                self.results["summary"]["passed"] += 1
                print(f"✅ {test_name} PASSED ({duration:.2f}s)")
            else:
                self.results["summary"]["failed"] += 1
                print(f"❌ {test_name} FAILED ({duration:.2f}s)")

            # Print key output lines
            if not success and errors:
                print("\nErrors:")
                print(errors[:500])  # First 500 chars of error

        except Exception as e:
            self.results["tests"][test_name] = {
                "file": test_file,
                "success": False,
                "error": str(e),
                "traceback": traceback.format_exc()
            }
            self.results["summary"]["errors"] += 1
            print(f"❌ {test_name} ERROR: {e}")

    def run_performance_benchmarks(self):
        """Run specific performance benchmarks."""
        print(f"\n{'='*60}")
        print("Running Performance Benchmarks")
        print(f"{'='*60}")

        try:
            # Import and run performance tests directly
            from test_cadd_comprehensive import TestCADDPerformance
            import unittest

            # Create test suite
            suite = unittest.TestLoader().loadTestsFromTestCase(TestCADDPerformance)
            runner = unittest.TextTestRunner(verbosity=2)

            start_time = time.time()
            result = runner.run(suite)
            duration = time.time() - start_time

            # Store performance metrics
            self.results["summary"]["performance_metrics"] = {
                "duration": duration,
                "tests_run": result.testsRun,
                "success": result.wasSuccessful(),
                "details": "See test_cadd_comprehensive.py output for detailed metrics"
            }

            if result.wasSuccessful():
                print(f"✅ Performance benchmarks PASSED ({duration:.2f}s)")
            else:
                print(f"❌ Performance benchmarks FAILED ({duration:.2f}s)")

        except Exception as e:
            print(f"❌ Performance benchmark error: {e}")
            self.results["summary"]["performance_metrics"]["error"] = str(e)

    def validate_database_schema(self):
        """Validate the database schema and data integrity."""
        print(f"\n{'='*60}")
        print("Validating Database Schema")
        print(f"{'='*60}")

        db_path = os.path.join(os.path.dirname(__file__), "population_variants.db")

        if not os.path.exists(db_path):
            print(f"⚠️  Database not found at: {db_path}")
            self.results["summary"]["database_validation"] = {
                "exists": False,
                "path": db_path
            }
            return

        try:
            import sqlite3
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()

            # Check tables exist
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
            tables = [row[0] for row in cursor.fetchall()]

            required_tables = ["cadd_scores", "cadd_jobs"]
            missing_tables = [t for t in required_tables if t not in tables]

            # Check CADD scores count
            if "cadd_scores" in tables:
                cursor.execute("SELECT COUNT(*) FROM cadd_scores")
                cadd_count = cursor.fetchone()[0]
            else:
                cadd_count = 0

            # Check indexes
            cursor.execute("SELECT name FROM sqlite_master WHERE type='index'")
            indexes = [row[0] for row in cursor.fetchall()]

            conn.close()

            # Store validation results
            validation = {
                "exists": True,
                "path": db_path,
                "size_mb": os.path.getsize(db_path) / 1024 / 1024,
                "tables": tables,
                "missing_tables": missing_tables,
                "cadd_scores_count": cadd_count,
                "indexes": indexes,
                "valid": len(missing_tables) == 0
            }

            self.results["summary"]["database_validation"] = validation

            if validation["valid"]:
                print("✅ Database validation PASSED")
                print(f"   - Size: {validation['size_mb']:.1f} MB")
                print(f"   - CADD scores: {validation['cadd_scores_count']:,}")
                print(f"   - Tables: {', '.join(tables)}")
            else:
                print("❌ Database validation FAILED")
                print(f"   - Missing tables: {', '.join(missing_tables)}")

        except Exception as e:
            print(f"❌ Database validation error: {e}")
            self.results["summary"]["database_validation"] = {
                "error": str(e)
            }

    def generate_report(self):
        """Generate comprehensive test report."""
        print(f"\n{'='*80}")
        print("TEST SUMMARY REPORT")
        print(f"{'='*80}")

        summary = self.results["summary"]

        # Overall results
        print("\nOverall Results:")
        print(f"  Total Test Suites: {summary['total_tests']}")
        print(f"  ✅ Passed: {summary['passed']}")
        print(f"  ❌ Failed: {summary['failed']}")
        print(f"  ⚠️  Errors: {summary['errors']}")

        if summary['total_tests'] > 0:
            success_rate = (summary['passed'] / summary['total_tests']) * 100
            print(f"  Success Rate: {success_rate:.1f}%")

        # Performance metrics
        if "performance_metrics" in summary and summary["performance_metrics"]:
            print("\nPerformance Metrics:")
            perf = summary["performance_metrics"]
            if "error" not in perf:
                print(f"  Duration: {perf.get('duration', 0):.2f}s")
                print(f"  Tests Run: {perf.get('tests_run', 0)}")
                print(f"  Success: {'✅' if perf.get('success') else '❌'}")
            else:
                print(f"  ❌ Error: {perf['error']}")

        # Database validation
        if "database_validation" in summary:
            print("\nDatabase Validation:")
            db_val = summary["database_validation"]
            if db_val.get("exists"):
                print(f"  Database Size: {db_val.get('size_mb', 0):.1f} MB")
                print(f"  CADD Scores: {db_val.get('cadd_scores_count', 0):,}")
                print(f"  Valid: {'✅' if db_val.get('valid') else '❌'}")
            else:
                print("  ❌ Database not found")

        # Individual test results
        print("\nIndividual Test Results:")
        for test_name, result in self.results["tests"].items():
            status = "✅" if result.get("success") else "❌"
            duration = result.get("duration", 0)
            print(f"  {status} {test_name}: {duration:.2f}s")

        # Save detailed report
        report_file = f"cadd_test_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        print(f"\n📄 Detailed report saved to: {report_file}")

        # Check if all tests passed
        all_passed = (
            summary["failed"] == 0 and
            summary["errors"] == 0 and
            summary.get("performance_metrics", {}).get("success", True) and
            summary.get("database_validation", {}).get("valid", True)
        )

        print(f"\n{'='*80}")
        if all_passed:
            print("✅ ALL TESTS PASSED! 🎉")
        else:
            print("❌ SOME TESTS FAILED - Please review the report")
        print(f"{'='*80}")

        return all_passed


def main():
    """Main entry point."""
    runner = CADDTestRunner()
    success = runner.run_all_tests()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
