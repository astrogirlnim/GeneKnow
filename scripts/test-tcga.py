#!/usr/bin/env python3
"""
Test script for TCGA integration
Tests connectivity and basic functionality of the TCGA processor
"""

import asyncio
import sys
import os
from pathlib import Path

# Add backend/python to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'backend', 'python'))

try:
    from genepredict.processors.tcga_processor import TCGAProcessor
except ImportError as e:
    print(f"❌ Failed to import TCGA processor: {e}")
    print("Make sure you're in the project root and Python dependencies are installed")
    sys.exit(1)

def main():
    """Main test function."""
    print("🧬 Testing TCGA Integration for GenePredict")
    print("=" * 50)
    
    # Test data directory creation
    data_dir = Path("./data/tcga")
    print(f"📁 Creating data directory: {data_dir}")
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Test TCGA processor initialization
    print("🏗️  Initializing TCGA processor...")
    try:
        processor = TCGAProcessor(data_dir)
        print("✅ TCGA processor initialized successfully")
    except Exception as e:
        print(f"❌ Failed to initialize TCGA processor: {e}")
        return 1
    
    # Run async tests
    print("🌐 Running async connectivity tests...")
    try:
        asyncio.run(run_async_tests(processor))
        print("✅ All TCGA tests passed!")
    except Exception as e:
        print(f"❌ TCGA tests failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

async def run_async_tests(processor):
    """Run async test functions."""
    
    # Test 1: Get cancer types
    print("   Testing cancer types retrieval...")
    cancer_types = await processor.get_available_cancer_types()
    print(f"   Found {len(cancer_types)} cancer types: {cancer_types[:5]}...")
    
    # Test 2: Test clinical data processing (small sample)
    print("   Testing clinical data processing...")
    clinical_df = await processor.process_clinical_outcomes("BRCA")
    print(f"   Clinical data shape: {clinical_df.shape}")
    if not clinical_df.empty:
        print(f"   Columns: {list(clinical_df.columns)}")
        print(f"   Sample data:")
        print(clinical_df.head(2))
    
    # Test 3: Test population frequencies (mock data)
    print("   Testing population frequencies...")
    test_variants = [
        {"chromosome": "17", "position": 41234000},
        {"chromosome": "13", "position": 32400000}
    ]
    frequencies = await processor.get_population_frequencies(test_variants)
    print(f"   Got frequencies for {len(frequencies)} variants")
    for var_id, freq in frequencies.items():
        print(f"     {var_id}: {freq:.4f}")
    
    print("   All async tests completed successfully!")

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code) 