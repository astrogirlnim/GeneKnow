#!/usr/bin/env python3
"""
Test pipeline integration with the new report generator.
"""

import sys
from pathlib import Path

# Add the parent directory to the path so we can import the module
sys.path.insert(0, str(Path(__file__).parent))

from graph import create_genomic_pipeline
from datetime import datetime

def test_pipeline_integration():
    """Test that the pipeline works with the new report generator."""
    
    print("=== Testing Pipeline Integration ===")
    
    # Create the pipeline
    try:
        pipeline = create_genomic_pipeline()
        print("✅ Pipeline created successfully")
    except Exception as e:
        print(f"❌ Failed to create pipeline: {e}")
        return
    
    # Create a minimal test state (simulating what would come from earlier nodes)
    test_state = {
        "file_path": "test_sample.vcf",
        "file_type": "vcf",
        "variant_count": 1250,
        "risk_scores": {
            "breast": 82.5,
            "ovarian": 45.2,
            "colon": 15.3,
            "lung": 8.7
        },
        "structured_json": {
            "case_id": "TEST_INTEGRATION_001",
            "processing_time_seconds": 45.2,
            "variant_count": 1250,
            "risk_assessment": {
                "scores": {
                    "breast": 82.5,
                    "ovarian": 45.2,
                    "colon": 15.3,
                    "lung": 8.7
                },
                "high_risk_findings": [
                    {"cancer_type": "breast", "risk_percentage": 82.5, "affected_genes": ["BRCA1"]},
                    {"cancer_type": "ovarian", "risk_percentage": 45.2, "affected_genes": ["BRCA1"]},
                    {"cancer_type": "colon", "risk_percentage": 15.3, "affected_genes": ["TP53"]},
                    {"cancer_type": "lung", "risk_percentage": 8.7, "affected_genes": ["TP53"]}
                ]
            },
            "summary": {
                "total_variants_found": 1250,
                "variants_passed_qc": 1250,
                "high_risk_findings": 4
            },
            "variant_details": [
                {
                    "gene": "BRCA1",
                    "variant": "chr17:41246747",
                    "consequence": "frameshift_variant",
                    "quality_metrics": {"quality": 99, "depth": 45, "allele_freq": 0.673},
                    "cadd_scores": {"phred": 35.0}
                }
            ]
        },
        "pipeline_start_time": datetime.now(),
        "pipeline_status": "in_progress",
        "completed_nodes": ["file_input", "preprocess", "variant_calling", "qc_filter", 
                           "merge_parallel", "population_mapper", "tcga_mapper", "cadd_scoring",
                           "clinvar_annotator", "prs_calculator", "pathway_burden", 
                           "merge_static_models", "feature_vector_builder", "ml_fusion",
                           "risk_model", "metrics_calculator", "formatter"]
    }
    
    # Test just the report generator node
    print("\n=== Testing Report Generator Node ===")
    try:
        from nodes.report_generator import process as report_generator_process
        result = report_generator_process(test_state)
        
        print("✅ Report generator executed successfully")
        print(f"Result keys: {list(result.keys())}")
        
        # Check expected outputs
        if "report_sections" in result:
            print("✅ report_sections created for dashboard")
            sections = result["report_sections"]
            print(f"  Sections: {list(sections.keys())}")
            
            # Check structure
            for key, section in sections.items():
                if isinstance(section, dict) and "title" in section and "content" in section:
                    print(f"  ✅ {key}: {section['title']} (severity: {section.get('severity', 'none')})")
                else:
                    print(f"  ⚠️ {key}: Invalid section structure")
        else:
            print("❌ No report_sections in result")
        
        if "enhanced_report_paths" in result:
            print("✅ enhanced_report_paths created")
            paths = result["enhanced_report_paths"]
            for fmt, path in paths.items():
                print(f"  {fmt}: {path}")
        else:
            print("❌ No enhanced_report_paths in result")
        
        if "report_generator_info" in result:
            print("✅ report_generator_info created")
            info = result["report_generator_info"]
            print(f"  Backend: {info.get('backend_used', 'unknown')}")
            print(f"  LLM Enhanced: {info.get('llm_enhanced', False)}")
            print(f"  High-risk findings: {info.get('high_risk_findings_count', 0)}")
        else:
            print("❌ No report_generator_info in result")
            
    except Exception as e:
        print(f"❌ Report generator failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_pipeline_integration() 