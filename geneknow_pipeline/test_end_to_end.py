#!/usr/bin/env python3
"""
Test end-to-end flow with the new report generator.
This simulates what the API server does.
"""

import sys
import json
import tempfile
from pathlib import Path

# Add the parent directory to the path so we can import the module
sys.path.insert(0, str(Path(__file__).parent))

from graph import run_pipeline
from datetime import datetime

def create_test_vcf():
    """Create a minimal test VCF file."""
    vcf_content = """##fileformat=VCFv4.2
##source=GeneKnowTest
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr12,length=133275309>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr17	41246481	.	C	T	99	PASS	GENE=BRCA1;CONSEQUENCE=frameshift
chr17	7577120	.	G	A	85	PASS	GENE=TP53;CONSEQUENCE=missense
chr12	25398284	.	C	T	70	PASS	GENE=KRAS;CONSEQUENCE=missense
"""
    
    # Create temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_content)
        return f.name

def test_end_to_end():
    """Test the complete pipeline flow."""
    print("=== End-to-End Pipeline Test ===")
    
    # Create test file
    test_file = create_test_vcf()
    print(f"Created test VCF: {test_file}")
    
    # User preferences
    user_preferences = {
        "language": "en",
        "include_technical": True,
        "patient_data": {
            "case_id": "E2E_TEST_001"
        }
    }
    
    try:
        print("\n=== Running Full Pipeline ===")
        start_time = datetime.now()
        
        # Run the pipeline (this should now use our new report generator)
        result = run_pipeline(test_file, user_preferences)
        
        end_time = datetime.now()
        processing_time = (end_time - start_time).total_seconds()
        
        print(f"✅ Pipeline completed in {processing_time:.2f} seconds")
        print(f"Pipeline status: {result.get('pipeline_status', 'unknown')}")
        
        # Check what we got back
        print(f"\n=== Pipeline Results ===")
        print(f"Result keys: {list(result.keys())}")
        
        # Check report_sections (for dashboard)
        if "report_sections" in result:
            print("✅ report_sections available for dashboard")
            sections = result["report_sections"]
            print(f"  Sections: {list(sections.keys())}")
            
            for key, section in sections.items():
                if isinstance(section, dict):
                    title = section.get('title', 'No title')
                    severity = section.get('severity', 'none')
                    content_preview = section.get('content', '')[:100] + "..." if len(section.get('content', '')) > 100 else section.get('content', '')
                    print(f"  ✅ {key}: {title} ({severity}) - {content_preview}")
        else:
            print("❌ No report_sections found")
        
        # Check enhanced report paths
        if "enhanced_report_paths" in result:
            print("✅ enhanced_report_paths available")
            paths = result["enhanced_report_paths"]
            for fmt, path in paths.items():
                if Path(path).exists():
                    size = Path(path).stat().st_size
                    print(f"  ✅ {fmt}: {path} ({size} bytes)")
                else:
                    print(f"  ❌ {fmt}: {path} (file not found)")
        else:
            print("❌ No enhanced_report_paths found")
        
        # Check report generator info
        if "report_generator_info" in result:
            print("✅ report_generator_info available")
            info = result["report_generator_info"]
            print(f"  Backend: {info.get('backend_used', 'unknown')}")
            print(f"  Model: {info.get('model_used', 'none')}")
            print(f"  LLM Enhanced: {info.get('llm_enhanced', False)}")
            print(f"  High-risk findings: {info.get('high_risk_findings_count', 0)}")
        else:
            print("❌ No report_generator_info found")
        
        # Check risk scores (should be preserved)
        if "risk_scores" in result:
            print("✅ risk_scores preserved")
            scores = result["risk_scores"]
            high_risk_scores = {k: v for k, v in scores.items() if v > 5.0}
            if high_risk_scores:
                print(f"  High-risk cancers: {high_risk_scores}")
            else:
                print("  No high-risk cancers detected")
        else:
            print("❌ No risk_scores found")
        
        # Simulate what the API would return
        print(f"\n=== API Response Simulation ===")
        api_response = {
            "job_id": "test_job_001",
            "status": "completed",
            "file_type": result.get("file_type", "vcf"),
            "processing_time_seconds": processing_time,
            "variant_count": result.get("variant_count", 0),
            "risk_scores": result.get("risk_scores", {}),
            "report_sections": result.get("report_sections", {}),
            "enhanced_report_paths": result.get("enhanced_report_paths", {}),
            "report_generator_info": result.get("report_generator_info", {}),
            "structured_json": result.get("structured_json", {}),
            "metrics": result.get("metrics", {}),
            "pipeline_status": result.get("pipeline_status", "unknown")
        }
        
        print("API response structure ready for frontend!")
        print(f"Response size: {len(json.dumps(api_response, default=str))} characters")
        
        return True
        
    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    finally:
        # Clean up test file
        try:
            Path(test_file).unlink()
            print(f"\nCleaned up test file: {test_file}")
        except:
            pass

if __name__ == "__main__":
    success = test_end_to_end()
    sys.exit(0 if success else 1) 