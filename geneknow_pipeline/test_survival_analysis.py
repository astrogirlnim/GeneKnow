#!/usr/bin/env python3
"""
Test script for survival analysis functionality.
Tests survival curve generation, data formatting, and API integration.
"""

import sys
import os
import json
from datetime import datetime
import logging

# Add the geneknow_pipeline directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import required modules
from nodes.survival_analyzer import process as survival_analyzer_process
from nodes.formatter import process as formatter_process
from nodes.pathway_burden import process as pathway_burden_process
from nodes.risk_model import process as risk_model_process

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def create_mock_state():
    """Create a mock state with realistic genomic data for testing."""
    return {
        "file_path": "test_sample.vcf",
        "file_type": "vcf",
        "file_metadata": {"qc_stats": {"variants_total": 150}},
        "variant_count": 150,
        "filtered_variants": [
            {
                "gene": "TP53",
                "variant_id": "chr17:7577121:C:T",
                "consequence": "missense_variant",
                "quality": 95,
                "depth": 45,
                "allele_freq": 0.52,
                "clinical_significance": "pathogenic",
                "cadd_phred": 32.5,
                "protein_change": "p.R175H"
            },
            {
                "gene": "KRAS",
                "variant_id": "chr12:25398284:C:T",
                "consequence": "missense_variant",
                "quality": 88,
                "depth": 38,
                "allele_freq": 0.48,
                "clinical_significance": "pathogenic",
                "cadd_phred": 28.7,
                "protein_change": "p.G12D"
            },
            {
                "gene": "BRCA1",
                "variant_id": "chr17:43094077:G:A",
                "consequence": "nonsense_variant",
                "quality": 92,
                "depth": 42,
                "allele_freq": 0.49,
                "clinical_significance": "pathogenic",
                "cadd_phred": 35.2,
                "protein_change": "p.R1751X"
            }
        ],
        "risk_scores": {
            "breast": 72.5,
            "lung": 68.3,
            "colon": 45.2,
            "prostate": 28.7,
            "blood": 15.3
        },
        "risk_genes": {
            "breast": ["BRCA1", "TP53"],
            "lung": ["TP53", "KRAS"],
            "colon": ["KRAS"],
            "prostate": [],
            "blood": []
        },
        "pathway_burden_results": {
            "dna_repair": {
                "burden_score": 0.85,
                "damaging_genes": ["BRCA1", "TP53"],
                "pathway_name": "DNA Repair"
            },
            "oncogenes": {
                "burden_score": 0.72,
                "damaging_genes": ["KRAS"],
                "pathway_name": "Oncogenes"
            },
            "tumor_suppressors": {
                "burden_score": 0.68,
                "damaging_genes": ["TP53"],
                "pathway_name": "Tumor Suppressors"
            }
        },
        "pathway_burden_summary": {
            "overall_burden_score": 0.75,
            "high_burden_pathways": ["dna_repair", "oncogenes", "tumor_suppressors"]
        },
        "variant_details": [],
        "tcga_matches": {},
        "cadd_stats": {},
        "clinvar_annotations": {},
        "prs_results": {},
        "completed_nodes": ["risk_model", "pathway_burden"],
        "pipeline_start_time": datetime.now(),
        "language": "en"
    }


def test_survival_analysis_generation():
    """Test survival analysis curve generation."""
    print("\n🧪 Testing Survival Analysis Generation")
    print("=" * 50)
    
    # Create mock state
    state = create_mock_state()
    
    # Process with survival analyzer
    result = survival_analyzer_process(state)
    
    # Check if survival analysis was generated
    survival_analysis = result.get("survival_analysis")
    
    if not survival_analysis:
        print("❌ No survival analysis generated")
        return False
    
    print("✅ Survival analysis generated successfully")
    
    # Check structure
    expected_keys = ["survival_curves", "prognostic_factors", "clinical_interpretation"]
    
    for key in expected_keys:
        if key not in survival_analysis:
            print(f"❌ Missing key: {key}")
            return False
        print(f"✅ Found key: {key}")
    
    # Check survival curves
    survival_curves = survival_analysis.get("survival_curves", {})
    if not survival_curves:
        print("❌ No survival curves generated")
        return False
    
    print(f"✅ Generated survival curves for {len(survival_curves)} cancer types")
    
    # Validate curve structure
    for cancer_type, curve_data in survival_curves.items():
        required_fields = ["time_points", "population_survival", "patient_survival", 
                         "median_survival", "hazard_ratio", "five_year_survival"]
        
        for field in required_fields:
            if field not in curve_data:
                print(f"❌ Missing field {field} in {cancer_type} curve")
                return False
        
        # Check data integrity
        time_points = curve_data["time_points"]
        pop_survival = curve_data["population_survival"]
        pat_survival = curve_data["patient_survival"]
        
        if len(time_points) != len(pop_survival) or len(time_points) != len(pat_survival):
            print(f"❌ Curve data length mismatch for {cancer_type}")
            return False
        
        print(f"✅ {cancer_type}: {len(time_points)} time points, HR={curve_data['hazard_ratio']}")
    
    return True


def test_formatter_integration():
    """Test that formatter properly includes survival analysis data."""
    print("\n🧪 Testing Formatter Integration")
    print("=" * 50)
    
    # Create mock state with survival analysis
    state = create_mock_state()
    state = survival_analyzer_process(state)
    
    # Process with formatter
    result = formatter_process(state)
    
    # Check if structured_json includes survival analysis
    structured_json = result.get("structured_json")
    
    if not structured_json:
        print("❌ No structured_json generated")
        return False
    
    survival_analysis = structured_json.get("survival_analysis")
    
    if not survival_analysis:
        print("❌ Survival analysis not included in structured_json")
        return False
    
    print("✅ Survival analysis included in structured_json")
    
    # Check if survival curves are present
    survival_curves = survival_analysis.get("survival_curves", {})
    if len(survival_curves) == 0:
        print("❌ No survival curves in formatted output")
        return False
    
    print(f"✅ Formatted output includes {len(survival_curves)} survival curves")
    
    # Validate data structure for frontend
    for cancer_type, curve_data in survival_curves.items():
        # Check if all required fields are present and properly formatted
        if not isinstance(curve_data.get("time_points"), list):
            print(f"❌ time_points not a list for {cancer_type}")
            return False
        
        if not isinstance(curve_data.get("population_survival"), list):
            print(f"❌ population_survival not a list for {cancer_type}")
            return False
        
        if not isinstance(curve_data.get("patient_survival"), list):
            print(f"❌ patient_survival not a list for {cancer_type}")
            return False
        
        print(f"✅ {cancer_type} curve data properly formatted for frontend")
    
    return True


def test_api_data_structure():
    """Test that the API data structure matches frontend expectations."""
    print("\n🧪 Testing API Data Structure")
    print("=" * 50)
    
    # Create complete pipeline state
    state = create_mock_state()
    state = survival_analyzer_process(state)
    state = formatter_process(state)
    
    # Extract the data that would be sent to frontend
    structured_json = state.get("structured_json", {})
    survival_analysis = structured_json.get("survival_analysis")
    
    if not survival_analysis:
        print("❌ No survival analysis in API response")
        return False
    
    # Check frontend-expected structure
    expected_structure = {
        "survival_curves": dict,
        "prognostic_factors": dict,
        "clinical_interpretation": dict
    }
    
    for field, expected_type in expected_structure.items():
        if field not in survival_analysis:
            print(f"❌ Missing field: {field}")
            return False
        
        if not isinstance(survival_analysis[field], expected_type):
            print(f"❌ Wrong type for {field}: expected {expected_type}, got {type(survival_analysis[field])}")
            return False
        
        print(f"✅ {field}: {expected_type.__name__}")
    
    # Check curve data structure matches frontend TypeScript interface
    survival_curves = survival_analysis["survival_curves"]
    for cancer_type, curve_data in survival_curves.items():
        required_fields = {
            "time_points": list,
            "population_survival": list,
            "patient_survival": list,
            "median_survival": dict,
            "hazard_ratio": (int, float),
            "five_year_survival": dict
        }
        
        for field, expected_type in required_fields.items():
            if field not in curve_data:
                print(f"❌ Missing field {field} in {cancer_type}")
                return False
            
            if not isinstance(curve_data[field], expected_type):
                print(f"❌ Wrong type for {cancer_type}.{field}: expected {expected_type}, got {type(curve_data[field])}")
                return False
        
        print(f"✅ {cancer_type} curve structure matches frontend expectations")
    
    return True


def test_clinical_recommendations():
    """Test clinical recommendations generation."""
    print("\n🧪 Testing Clinical Recommendations")
    print("=" * 50)
    
    # Create mock state
    state = create_mock_state()
    state = survival_analyzer_process(state)
    
    survival_analysis = state.get("survival_analysis")
    if not survival_analysis:
        print("❌ No survival analysis to test recommendations")
        return False
    
    clinical_interpretation = survival_analysis.get("clinical_interpretation")
    if not clinical_interpretation:
        print("❌ No clinical interpretation generated")
        return False
    
    recommendation = clinical_interpretation.get("recommendation")
    if not recommendation:
        print("❌ No clinical recommendation generated")
        return False
    
    print(f"✅ Clinical recommendation: {recommendation}")
    
    # Check prognostic factors
    prognostic_factors = survival_analysis.get("prognostic_factors")
    if not prognostic_factors:
        print("❌ No prognostic factors generated")
        return False
    
    positive_factors = prognostic_factors.get("positive_prognostic_factors", [])
    negative_factors = prognostic_factors.get("negative_prognostic_factors", [])
    
    print(f"✅ Positive prognostic factors: {len(positive_factors)}")
    print(f"✅ Negative prognostic factors: {len(negative_factors)}")
    
    return True


def test_end_to_end_pipeline():
    """Test the complete pipeline flow including survival analysis."""
    print("\n🧪 Testing End-to-End Pipeline")
    print("=" * 50)
    
    try:
        # Test importing the main pipeline
        from graph import run_pipeline
        
        print("✅ Pipeline import successful")
        
        # Note: We won't run the full pipeline as it requires actual files
        # but we can verify the survival analysis is properly integrated
        
        # Create a minimal test state
        state = create_mock_state()
        
        # Test the survival analysis flow
        state = survival_analyzer_process(state)
        state = formatter_process(state)
        
        # Verify final output structure
        structured_json = state.get("structured_json", {})
        
        if "survival_analysis" not in structured_json:
            print("❌ Survival analysis missing from final output")
            return False
        
        print("✅ Survival analysis present in final pipeline output")
        
        # Save test output for inspection
        with open("test_survival_output.json", "w") as f:
            json.dump(structured_json.get("survival_analysis"), f, indent=2)
        
        print("✅ Test output saved to test_survival_output.json")
        
        return True
        
    except Exception as e:
        print(f"❌ Pipeline test failed: {str(e)}")
        return False


def run_all_tests():
    """Run all survival analysis tests."""
    print("🧬 GeneKnow Survival Analysis Test Suite")
    print("=" * 60)
    
    tests = [
        ("Survival Analysis Generation", test_survival_analysis_generation),
        ("Formatter Integration", test_formatter_integration),
        ("API Data Structure", test_api_data_structure),
        ("Clinical Recommendations", test_clinical_recommendations),
        ("End-to-End Pipeline", test_end_to_end_pipeline)
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
                print(f"\n✅ {test_name}: PASSED")
            else:
                print(f"\n❌ {test_name}: FAILED")
        except Exception as e:
            print(f"\n❌ {test_name}: ERROR - {str(e)}")
    
    print("\n" + "=" * 60)
    print(f"📊 Test Results: {passed}/{total} passed ({passed/total*100:.1f}%)")
    
    if passed == total:
        print("🎉 All tests passed! Survival analysis is working correctly.")
    else:
        print("⚠️  Some tests failed. Review the output above.")
    
    return passed == total


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1) 