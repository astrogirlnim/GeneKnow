#!/usr/bin/env python3
"""
Run all offline CADD scoring tests.
"""
import unittest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import test modules
from test_cadd_offline import TestOfflineCADDScoring
from test_cadd_comprehensive import (
    TestCADDComprehensive,
    TestCADDPerformance,
    TestCADDMemory,
    TestCADDEdgeCases
)

if __name__ == "__main__":
    # Create test suite
    suite = unittest.TestSuite()
    
    # Add offline CADD tests
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestOfflineCADDScoring))
    
    # Add comprehensive tests (if they don't reference old functions)
    try:
        suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCADDComprehensive))
        suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCADDPerformance))
        suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCADDMemory))
        suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCADDEdgeCases))
    except Exception as e:
        print(f"Note: Some comprehensive tests skipped due to: {e}")
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Exit with appropriate code
    sys.exit(0 if result.wasSuccessful() else 1) 