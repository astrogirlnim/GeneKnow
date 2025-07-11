#!/usr/bin/env python3
"""
Test script to verify metrics calculator integration in the GeneKnow pipeline.
"""
import os
import sys
import json
import logging
from datetime import datetime
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent))

from graph import create_genomic_pipeline, run_pipeline

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def test_metrics_calculator():
    """Test the metrics calculator node with sample data."""
    
    # Use a test VCF file if available
    test_file = None
    test_paths = [
        "test_data/test_output.vcf",
        "data/test_output.vcf",
        "../test_data/test_output.vcf"
    ]
    
    for path in test_paths:
        if os.path.exists(path):
            test_file = path
            break
    
    if not test_file:
        logger.warning("No test VCF file found. Creating minimal test case...")
        # Create a minimal test VCF for demonstration
        test_file = "test_minimal.vcf"
        with open(test_file, 'w') as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("17\t43044295\trs80357506\tG\tA\t100\tPASS\tGENE=BRCA1\n")
            f.write("13\t32914437\trs80359550\tG\tA\t100\tPASS\tGENE=BRCA2\n")
    
    # Run the pipeline
    logger.info(f"Running pipeline with test file: {test_file}")
    
    user_preferences = {
        "patient_data": {
            "age": 45,
            "sex": "F"
        },
        "include_technical": True,
        "language": "en"
    }
    
    try:
        result = run_pipeline(test_file, user_preferences)
        
        # Check pipeline status
        logger.info(f"Pipeline status: {result['pipeline_status']}")
        logger.info(f"Completed nodes: {result['completed_nodes']}")
        
        # Check if metrics calculator ran
        if "metrics_calculator" in result['completed_nodes']:
            logger.info("✅ Metrics calculator successfully integrated!")
            
            # Display metrics results
            metrics = result.get("metrics", {})
            if metrics:
                logger.info("\n📊 Metrics Results:")
                
                # Confidence metrics
                conf_metrics = metrics.get("confidence_metrics", {})
                if conf_metrics:
                    logger.info("\nConfidence Metrics:")
                    logger.info(f"  Mean model confidence: {conf_metrics.get('mean_model_confidence', 0):.3f}")
                    logger.info(f"  ML fusion confidence: {conf_metrics.get('ml_fusion_confidence', 0):.3f}")
                    logger.info(f"  Risk category: {conf_metrics.get('ml_fusion_risk_category', 'unknown')}")
                
                # Variant metrics
                var_metrics = metrics.get("variant_metrics", {})
                if var_metrics:
                    logger.info("\nVariant Metrics:")
                    logger.info(f"  Total variants: {var_metrics.get('total_variants', 0)}")
                    logger.info(f"  Pathogenic variants: {var_metrics.get('pathogenic_variants', 0)}")
                    logger.info(f"  Uncertain variants: {var_metrics.get('uncertain_variants', 0)}")
                    logger.info(f"  Genes affected: {var_metrics.get('genes_affected', 0)}")
                    logger.info(f"  Cancer genes affected: {var_metrics.get('cancer_genes_affected', 0)}")
                
                # Overall assessment
                overall = metrics.get("overall_assessment", {})
                if overall:
                    logger.info("\nOverall Assessment:")
                    logger.info(f"  High-risk cancers: {', '.join(overall.get('high_risk_cancers', [])) or 'None'}")
                    logger.info(f"  Max risk score: {overall.get('max_risk_score', 0):.1f}%")
                    logger.info(f"  Risk category: {overall.get('risk_category', 'unknown')}")
                    logger.info(f"  Clinical action needed: {overall.get('clinical_action_needed', False)}")
                
                # Performance indicators
                perf_indicators = metrics.get("performance_indicators", {})
                if perf_indicators:
                    logger.info("\nModel Performance Indicators:")
                    for indicator, value in perf_indicators.items():
                        logger.info(f"  {indicator}: {value}")
                
                # Validation structure
                validation = metrics.get("validation_metrics", {})
                if validation:
                    logger.info("\nValidation Metrics Structure:")
                    logger.info(f"  Validation ready: {validation.get('validation_ready', False)}")
                    logger.info(f"  Ground truth available: {validation.get('ground_truth_available', False)}")
                    logger.info(f"  Note: {validation.get('validation_note', '')}")
                
            # Check metrics summary
            metrics_summary = result.get("metrics_summary", {})
            if metrics_summary:
                logger.info("\n📋 Metrics Summary:")
                key_findings = metrics_summary.get("key_findings", {})
                if key_findings:
                    logger.info(f"  Highest risk cancer: {key_findings.get('highest_risk_cancer', 'None')}")
                    logger.info(f"  Highest risk score: {key_findings.get('highest_risk_score', 0):.1f}%")
                    logger.info(f"  Pathogenic variant count: {key_findings.get('pathogenic_variant_count', 0)}")
                    logger.info(f"  Confidence level: {key_findings.get('confidence_level', 'unknown')}")
                
        else:
            logger.error("❌ Metrics calculator did not run!")
            
        # Check for errors
        if result.get("errors"):
            logger.error("\n⚠️ Errors encountered:")
            for error in result["errors"]:
                logger.error(f"  Node: {error.get('node', 'unknown')}")
                logger.error(f"  Error: {error.get('error', 'unknown error')}")
                
        # Save results for inspection
        output_file = "metrics_test_results.json"
        with open(output_file, 'w') as f:
            # Convert datetime objects to strings for JSON serialization
            json_safe_result = {}
            for key, value in result.items():
                if isinstance(value, datetime):
                    json_safe_result[key] = value.isoformat()
                elif key == "errors" and isinstance(value, list):
                    json_safe_result[key] = [
                        {k: v.isoformat() if isinstance(v, datetime) else v for k, v in err.items()}
                        for err in value
                    ]
                else:
                    json_safe_result[key] = value
                    
            json.dump(json_safe_result, f, indent=2)
        logger.info(f"\n💾 Full results saved to: {output_file}")
        
        # Clean up test file if we created it
        if test_file == "test_minimal.vcf":
            os.remove(test_file)
            
    except Exception as e:
        logger.error(f"Pipeline test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        
        # Clean up test file if we created it
        if test_file == "test_minimal.vcf" and os.path.exists(test_file):
            os.remove(test_file)


def test_metrics_calculation_logic():
    """Test the metrics calculation functions directly."""
    from nodes.metrics_calculator import (
        calculate_confidence_metrics,
        calculate_variant_metrics,
        calculate_prediction_metrics
    )
    
    logger.info("\n🧪 Testing metrics calculation logic...")
    
    # Create sample data
    risk_scores = {
        "breast": 65.0,
        "colon": 25.0,
        "lung": 15.0
    }
    
    risk_details = {
        "breast": {
            "model_confidence": 0.85,
            "gene_count": 3,
            "pathogenic_count": 2
        },
        "colon": {
            "model_confidence": 0.75,
            "gene_count": 1,
            "pathogenic_count": 0
        },
        "lung": {
            "model_confidence": 0.80,
            "gene_count": 0,
            "pathogenic_count": 0
        }
    }
    
    ml_assessment = {
        "aggregate_risk_score": 0.7,
        "risk_category": "high",
        "confidence": 0.82,
        "high_risk_variants": 2,
        "contributing_factors": {
            "cadd_score": 0.8,
            "clinvar_classification": 0.9
        }
    }
    
    # Test confidence metrics
    conf_metrics = calculate_confidence_metrics(risk_scores, risk_details, ml_assessment)
    logger.info("\nConfidence Metrics:")
    for key, value in conf_metrics.items():
        logger.info(f"  {key}: {value}")
    
    # Test prediction metrics
    pred_metrics = calculate_prediction_metrics(risk_scores, risk_details, ml_assessment)
    logger.info("\nPrediction Metrics:")
    for key, value in pred_metrics.items():
        if not key.startswith("factor_"):
            logger.info(f"  {key}: {value}")
    
    logger.info("\n✅ Metrics calculation logic test complete!")


if __name__ == "__main__":
    logger.info("🚀 Starting metrics calculator integration test...")
    
    # Test the full pipeline
    test_metrics_calculator()
    
    # Test individual functions
    test_metrics_calculation_logic()
    
    logger.info("\n🎉 All tests complete!") 