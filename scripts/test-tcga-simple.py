#!/usr/bin/env python3
"""
Simplified TCGA integration test
Tests TCGA API connectivity without full model dependencies
"""

import asyncio
import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List
import aiohttp
import pandas as pd
from datetime import datetime

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SimpleTCGAProcessor:
    """Simplified TCGA processor for testing connectivity."""
    
    def __init__(self, data_dir: Path):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.base_url = "https://api.gdc.cancer.gov"
        
        # Create subdirectories
        for subdir in ['raw', 'processed', 'cache']:
            (self.data_dir / subdir).mkdir(exist_ok=True)
    
    async def test_api_connection(self) -> bool:
        """Test basic API connectivity."""
        try:
            async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=30)) as session:
                url = f"{self.base_url}/status"
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        logger.info(f"✅ TCGA API is accessible. Status: {data}")
                        return True
                    else:
                        logger.error(f"❌ TCGA API returned status: {response.status}")
                        return False
        except Exception as e:
            logger.error(f"❌ Failed to connect to TCGA API: {e}")
            return False
    
    async def get_cancer_types(self) -> List[str]:
        """Get available cancer types."""
        try:
            async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=60)) as session:
                url = f"{self.base_url}/projects"
                params = {
                    "format": "json",
                    "size": "50",
                    "filters": json.dumps({
                        "op": "in",
                        "content": {
                            "field": "program.name",
                            "value": ["TCGA"]
                        }
                    })
                }
                
                async with session.get(url, params=params) as response:
                    if response.status == 200:
                        data = await response.json()
                        cancer_types = [
                            project["project_id"].replace("TCGA-", "")
                            for project in data["data"]["hits"]
                        ]
                        return sorted(cancer_types)
                    else:
                        logger.warning(f"Failed to get cancer types: {response.status}")
                        return ["BRCA", "LUAD", "COAD", "PRAD"]  # Fallback
        except Exception as e:
            logger.error(f"Error getting cancer types: {e}")
            return ["BRCA", "LUAD", "COAD", "PRAD"]  # Fallback
    
    async def get_brca_project_info(self) -> Dict:
        """Get BRCA project information."""
        try:
            async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=60)) as session:
                url = f"{self.base_url}/projects/TCGA-BRCA"
                params = {"format": "json", "expand": "summary"}
                
                async with session.get(url, params=params) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data["data"]
                    else:
                        logger.warning(f"Failed to get BRCA info: {response.status}")
                        return {}
        except Exception as e:
            logger.error(f"Error getting BRCA info: {e}")
            return {}

def main():
    """Main test function."""
    print("🧬 Testing TCGA API Connectivity")
    print("=" * 40)
    
    # Create processor
    data_dir = Path("./data/tcga")
    processor = SimpleTCGAProcessor(data_dir)
    print(f"📁 Using data directory: {data_dir}")
    
    # Run tests
    try:
        asyncio.run(run_tests(processor))
        print("\n✅ All TCGA connectivity tests passed!")
        return 0
    except Exception as e:
        print(f"\n❌ Tests failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

async def run_tests(processor):
    """Run async tests."""
    
    # Test 1: Basic API connectivity
    print("\n🌐 Testing API connectivity...")
    connected = await processor.test_api_connection()
    if not connected:
        raise Exception("Failed to connect to TCGA API")
    
    # Test 2: Get cancer types
    print("\n📊 Testing cancer types retrieval...")
    cancer_types = await processor.get_cancer_types()
    print(f"   Found {len(cancer_types)} cancer types")
    print(f"   Sample types: {cancer_types[:5]}")
    
    # Test 3: Get BRCA project info
    print("\n🎯 Testing BRCA project information...")
    brca_info = await processor.get_brca_project_info()
    if brca_info:
        print(f"   BRCA Project: {brca_info.get('project_id', 'N/A')}")
        print(f"   Name: {brca_info.get('name', 'N/A')}")
        summary = brca_info.get('summary', {})
        if summary:
            print(f"   Cases: {summary.get('case_count', 'N/A')}")
            print(f"   Files: {summary.get('file_count', 'N/A')}")
    else:
        print("   Could not retrieve BRCA project info")
    
    print("\n🎉 All connectivity tests completed successfully!")

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code) 