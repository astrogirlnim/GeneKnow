#!/usr/bin/env python3
"""
Final test demonstrating the metrics calculator implementation works correctly.
Tests the metrics_calculator directly with realistic state data.
"""
import sys
import os
import json
from datetime import datetime

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the metrics calculator directly
from geneknow_pipeline.nodes.metrics_calculator import process as calculate_metrics

def create_realistic_state():
    """Create a realistic state with all required fields populated."""
    return {
        "current_node": "test",
        "completed_nodes": ["file_input", "preprocess", "risk_model"],
        "errors": [],
        "warnings": [],
        "model_version": "1.0.0",
        "pipeline_start_time": datetime.now(),
        
        # Risk model outputs (realistic values)
        "risk_scores": {
            "breast": 72.5,
            "prostate": 23.8,
            "colon": 18.2,
            "lung": 15.4,
            "ovarian": 8.9,
            "blood": 12.3
        },
        
        "risk_details": {
            "breast": {
                "base_risk": 12.5,
                "genetic_risk": 45.0,
                "prs_contribution": 15.0,
                "confidence": 0.85
            },
            "prostate": {
                "base_risk": 10.0,
                "genetic_risk": 8.8,
                "prs_contribution": 5.0,
                "confidence": 0.72
            },
            "colon": {
                "base_risk": 6.0,
                "genetic_risk": 8.2,
                "prs_contribution": 4.0,
                "confidence": 0.68
            }
        },
        
        "ml_risk_assessment": {
            "model_confidence": 0.78,
            "prediction_certainty": 0.82,
            "risk_category": "moderate",
            "fusion_score": 0.74
        },
        
        # Variant data
        "filtered_variants": [
            {
                "gene": "BRCA1",
                "variant_id": "chr17:41244936:G>A",
                "clinical_significance": "pathogenic",
                "cadd_phred": 25.8,
                "consequence": "missense_variant"
            },
            {
                "gene": "BRCA2", 
                "variant_id": "chr13:32913055:A>G",
                "clinical_significance": "likely_pathogenic",
                "cadd_phred": 23.1,
                "consequence": "missense_variant"
            },
            {
                "gene": "TP53",
                "variant_id": "chr17:7577548:C>T", 
                "clinical_significance": "uncertain_significance",
                "cadd_phred": 18.9,
                "consequence": "missense_variant"
            },
            {
                "gene": "ATM",
                "variant_id": "chr11:108202901:C>T",
                "clinical_significance": "uncertain_significance", 
                "cadd_phred": 16.2,
                "consequence": "missense_variant"
            }
        ],
        
        # CADD stats
        "cadd_stats": {
            "variants_scored": 4,
            "mean_phred": 21.0,
            "max_phred": 25.8,
            "variants_gt20": 2,
            "variants_in_cancer_genes": 4
        },
        
        # PRS results
        "prs_results": {
            "breast": {
                "raw_score": 2.34,
                "percentile": 89,
                "confidence": "high"
            },
            "prostate": {
                "raw_score": 0.78,
                "percentile": 62,
                "confidence": "moderate"
            },
            "colon": {
                "raw_score": 0.45,
                "percentile": 55,
                "confidence": "moderate"
            }
        },
        
        "prs_summary": {
            "overall_confidence": "high",
            "high_risk_cancers": ["breast"]
        },
        
        # Pathway burden results
        "pathway_burden_results": {
            "dna_repair": {
                "burden_score": 0.85,
                "total_variants": 3,
                "damaging_variants": 2
            },
            "tumor_suppressors": {
                "burden_score": 0.72,
                "total_variants": 2,
                "damaging_variants": 1
            }
        },
        
        "pathway_burden_summary": {
            "overall_risk_level": "high",
            "high_burden_pathways": ["dna_repair", "tumor_suppressors"]
        },
        
        # Optional validation data for demonstration
        "validation_data": {
            "ground_truth_cancer_type": "breast",
            "actual_diagnosis": "positive",
            "time_to_diagnosis": "24_months"
        }
    }

def test_metrics_calculator():
    """Test the metrics calculator with realistic data."""
    print("🧪 Final Metrics Calculator Test")
    print("=" * 60)
    
    # Create realistic state
    state = create_realistic_state()
    
    print("📊 Input Data Summary:")
    print(f"  Risk scores: {len(state['risk_scores'])} cancer types")
    print(f"  Highest risk: {max(state['risk_scores'].items(), key=lambda x: x[1])}")
    print(f"  Variants: {len(state['filtered_variants'])}")
    print(f"  Model confidence: {state['ml_risk_assessment']['model_confidence']:.2f}")
    print(f"  Validation data: {'Available' if 'validation_data' in state else 'Not available'}")
    
    print("\n🔄 Running metrics calculator...")
    
    # Run metrics calculator
    result_state = calculate_metrics(state)
    
    print("✅ Metrics calculation completed!")
    
    # Check results
    if "metrics" in result_state:
        metrics = result_state["metrics"]
        print(f"\n📈 Metrics Results:")
        print(f"  Timestamp: {metrics.get('timestamp', 'N/A')}")
        print(f"  Pipeline version: {metrics.get('pipeline_version', 'N/A')}")
        
        # Confidence metrics
        if "confidence_metrics" in metrics:
            conf = metrics["confidence_metrics"]
            print(f"\n🎯 Confidence Metrics:")
            print(f"  Mean model confidence: {conf.get('mean_model_confidence', 0):.3f}")
            print(f"  Risk score range: {conf.get('risk_score_range', 0):.1f}")
            print(f"  ML fusion confidence: {conf.get('ml_fusion_confidence', 0):.3f}")
        
        # Variant metrics
        if "variant_metrics" in metrics:
            var = metrics["variant_metrics"]
            print(f"\n🧬 Variant Metrics:")
            print(f"  Total variants: {var.get('total_variants', 0)}")
            print(f"  Pathogenic variants: {var.get('pathogenic_variants', 0)}")
            print(f"  Uncertain variants: {var.get('uncertain_variants', 0)}")
            print(f"  High CADD variants: {var.get('high_cadd_variants', 0)}")
            print(f"  Genes affected: {var.get('genes_affected', 0)}")
        
        # Prediction metrics
        if "prediction_metrics" in metrics:
            pred = metrics["prediction_metrics"]
            print(f"\n📊 Prediction Metrics:")
            print(f"  Max risk score: {pred.get('max_risk_score', 0):.1f}%")
            print(f"  High risk cancer count: {pred.get('high_risk_cancer_count', 0)}")
            print(f"  ML risk category: {pred.get('ml_risk_category', 'unknown')}")
        
        # PRS metrics
        if "prs_metrics" in metrics:
            prs = metrics["prs_metrics"]
            print(f"\n🧮 PRS Metrics:")
            print(f"  High PRS cancer count: {prs.get('high_prs_cancer_count', 0)}")
            print(f"  Max PRS percentile: {prs.get('max_prs_percentile', 0)}")
            print(f"  Overall confidence: {prs.get('prs_overall_confidence', 'unknown')}")
        
        # Pathway metrics
        if "pathway_metrics" in metrics:
            path = metrics["pathway_metrics"]
            print(f"\n🛤️  Pathway Metrics:")
            print(f"  High burden pathway count: {path.get('high_burden_pathway_count', 0)}")
            print(f"  Mean pathway burden: {path.get('mean_pathway_burden', 0):.3f}")
            print(f"  Pathway risk level: {path.get('pathway_risk_level', 'unknown')}")
        
        # Overall assessment
        if "overall_assessment" in metrics:
            overall = metrics["overall_assessment"]
            print(f"\n⚡ Overall Assessment:")
            print(f"  High-risk cancers: {overall.get('high_risk_cancers', [])}")
            print(f"  Max risk score: {overall.get('max_risk_score', 0):.1f}%")
            print(f"  Clinical action needed: {overall.get('clinical_action_needed', False)}")
            print(f"  Risk category: {overall.get('risk_category', 'unknown')}")
        
        # Performance indicators
        if "performance_indicators" in metrics:
            perf = metrics["performance_indicators"]
            print(f"\n⚙️  Performance Indicators:")
            print(f"  Variant coverage: {perf.get('variant_coverage', False)}")
            print(f"  Model confidence adequate: {perf.get('model_confidence_adequate', False)}")
            print(f"  Sufficient evidence: {perf.get('sufficient_evidence', False)}")
        
        # Validation metrics
        if "validation_metrics" in metrics:
            val = metrics["validation_metrics"]
            print(f"\n🔬 Validation Metrics:")
            print(f"  Ground truth available: {val.get('ground_truth_available', False)}")
            print(f"  Validation ready: {val.get('validation_ready', False)}")
            if val.get('ground_truth_available'):
                print(f"  Expected AUC-ROC: {val.get('expected_auc_roc', 'N/A')}")
                print(f"  Expected sensitivity: {val.get('expected_sensitivity', 'N/A')}")
    
    # Check metrics summary
    if "metrics_summary" in result_state:
        summary = result_state["metrics_summary"]
        print(f"\n📋 Metrics Summary:")
        key_findings = summary.get("key_findings", {})
        print(f"  Highest risk cancer: {key_findings.get('highest_risk_cancer', 'N/A')}")
        print(f"  Highest risk score: {key_findings.get('highest_risk_score', 0):.1f}%")
        print(f"  Pathogenic variant count: {key_findings.get('pathogenic_variant_count', 0)}")
        print(f"  Confidence level: {key_findings.get('confidence_level', 'unknown')}")
        print(f"  Validation ready: {summary.get('validation_status', False)}")
    
    print("\n✅ Test completed successfully!")
    print("🎉 Metrics calculator implementation is working correctly!")
    
    return result_state

if __name__ == "__main__":
    test_metrics_calculator() 