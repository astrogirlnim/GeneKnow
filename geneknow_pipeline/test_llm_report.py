#!/usr/bin/env python3
"""
Test LLM-enhanced report generation.
"""

import sys
from pathlib import Path

# Add the parent directory to the path so we can import the module
sys.path.insert(0, str(Path(__file__).parent))

from nodes.report_generator import process as generate_report

def create_sample_state():
    """Create a sample state dict that mimics what the pipeline would pass."""
    return {
        "case_id": "BRCA_TEST_001",
        "structured_json": {
            "case_id": "BRCA_TEST_001",
            "processing_time_seconds": 45.2,
            "variant_count": 1250,
            "risk_scores": {
                "breast": 82.5,
                "ovarian": 45.2,
                "colon": 15.3,
                "lung": 8.7,
                "prostate": 3.2,
                "pancreatic": 2.1
            },
            "variants": [
                {
                    "gene": "BRCA1",
                    "position": "chr17:41246747",
                    "type": "frameshift",
                    "impact": "HIGH",
                    "clinvar_significance": "Pathogenic",
                    "cadd_score": 35.0,
                    "tcga_cancer_relevance": 0.95,
                    "pathway_damage_assessment": {
                        "pathways_affected": ["DNA repair", "Cell cycle"],
                        "damage_score": 0.89
                    }
                },
                {
                    "gene": "TP53",
                    "position": "chr17:7577121",
                    "type": "missense",
                    "impact": "MODERATE",
                    "clinvar_significance": "Likely pathogenic",
                    "cadd_score": 28.5,
                    "tcga_cancer_relevance": 0.75,
                    "pathway_damage_assessment": {
                        "pathways_affected": ["p53 signaling"],
                        "damage_score": 0.72
                    }
                }
            ],
            "metrics": {
                "variant_metrics": {
                    "total_variants": 1250,
                    "high_impact": 45,
                    "moderate_impact": 230,
                    "low_impact": 975
                },
                "confidence_metrics": {
                    "overall_confidence": 0.92,
                    "data_quality": 0.98
                }
            },
            "report_sections": {
                "overview": {
                    "title": "Analysis Overview",
                    "content": "Genomic analysis identified significant cancer risk factors.",
                    "severity": "high"
                },
                "key_findings": {
                    "title": "Key Findings",
                    "content": "BRCA1 pathogenic variant detected with high cancer relevance.",
                    "severity": "high"
                }
            }
        }
    }

def test_llm_report():
    print("=== Testing LLM-Enhanced Report Generation ===")
    
    # Create sample state
    state = create_sample_state()
    
    # Generate report (should now use LLM)
    result = generate_report(state)
    
    # Check results
    print(f"Original sections preserved: {'report_sections' in result}")
    print(f"Enhanced reports generated: {'enhanced_report_paths' in result}")
    print(f"Report info available: {'report_generator_info' in result}")
    
    if "enhanced_report_paths" in result:
        for fmt, path in result["enhanced_report_paths"].items():
            print(f"  {fmt}: {path}")
            
            # Read and display first 1000 chars of markdown
            if fmt == "markdown" and Path(path).exists():
                with open(path, 'r') as f:
                    content = f.read()
                    print(f"\n=== First 1000 chars of LLM-enhanced report ===")
                    print(content[:1000] + "..." if len(content) > 1000 else content)
                    print("=" * 50)
                    
                    # Check for key indicators
                    if "NON-LLM MODE" in content:
                        print("⚠️  Still in fallback mode")
                    else:
                        print("✅ LLM-enhanced report generated!")
                    
                    if "BRCA1" in content:
                        print("✅ High-risk variant mentioned")
                    if "Medical Glossary" in content:
                        print("✅ Glossary section found")
    
    if "report_generator_info" in result:
        info = result["report_generator_info"]
        print(f"\nReport Generator Info:")
        print(f"  Backend: {info.get('backend_used', 'unknown')}")
        print(f"  Model: {info.get('model_used', 'none')}")
        print(f"  LLM Enhanced: {info.get('llm_enhanced', False)}")
        print(f"  High-risk findings: {info.get('high_risk_findings_count', 0)}")

if __name__ == "__main__":
    test_llm_report() 