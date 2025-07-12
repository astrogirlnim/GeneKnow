#!/usr/bin/env python3
"""
Test script for SHAP validator node.
Tests the integration with ML fusion and validates the sanity rules.
"""
import numpy as np
from pathlib import Path
import sys
import logging

# Add parent directory to path
sys.path.append(str(Path(__file__).parent))

from nodes.shap_validator import process as shap_validator_process
from ml_models.fusion_layer import FusionLayer, StaticModelInputs

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def create_test_state_high_risk_valid():
    """Create test state with high risk driven by pathogenic variant (should PASS)."""
    # Create a trained fusion layer
    fusion_layer = FusionLayer(model_type="gradient_boosting")

    # Create training data
    training_data = []
    np.random.seed(42)

    # Add samples with pathogenic variants having high risk
    # Make pathogenic the dominant signal
    for _ in range(300):
        # Strong signal: pathogenic = high risk
        inputs = StaticModelInputs(
            prs_score=np.random.uniform(0.1, 0.4),  # Low to moderate PRS
            clinvar_classification="pathogenic",
            cadd_score=np.random.uniform(15, 25),  # Moderate CADD
            tcga_enrichment=np.random.uniform(1, 2),  # Low enrichment
            gene_burden_score=np.random.uniform(0, 1),  # Low burden
        )
        risk_score = np.random.uniform(0.75, 0.95)  # Very high risk
        training_data.append((inputs, risk_score))

    # Add benign samples with low risk
    for _ in range(300):
        inputs = StaticModelInputs(
            prs_score=np.random.uniform(0.1, 0.9),  # Any PRS
            clinvar_classification="benign",
            cadd_score=np.random.uniform(0, 35),  # Any CADD
            tcga_enrichment=np.random.uniform(0.5, 5),  # Any enrichment
            gene_burden_score=np.random.uniform(0, 3),  # Any burden
        )
        risk_score = np.random.uniform(0.05, 0.25)  # Low risk
        training_data.append((inputs, risk_score))

    # Add uncertain and not_found variants with medium risk
    for _ in range(200):
        inputs = StaticModelInputs(
            prs_score=np.random.uniform(0.3, 0.6),
            clinvar_classification=np.random.choice(["uncertain", "not_found"]),
            cadd_score=np.random.uniform(10, 20),
            tcga_enrichment=np.random.uniform(1, 2),
            gene_burden_score=np.random.uniform(0, 1),
        )
        risk_score = np.random.uniform(0.3, 0.5)  # Medium risk
        training_data.append((inputs, risk_score))

    # Train the model
    fusion_layer.train(training_data, validation_split=0)

    # Create test input with pathogenic variant
    test_inputs = [
        StaticModelInputs(
            prs_score=0.25,  # Low PRS (matching training range for pathogenic)
            clinvar_classification="pathogenic",
            cadd_score=20.0,  # Moderate CADD
            tcga_enrichment=1.5,  # Low enrichment
            gene_burden_score=0.5,  # Low burden
        )
    ]

    # Get feature matrix
    X_encoded = fusion_layer._encode_features(test_inputs)
    X_scaled = fusion_layer.scaler.transform(X_encoded)

    # Create state
    state = {
        "ml_fusion_model_instance": fusion_layer,
        "ml_fusion_feature_matrix": X_scaled,
        "risk_scores": {"aggregate_risk_score": 0.75},  # High risk
        "filtered_variants": [
            {"variant_id": "chr1:12345", "gene": "BRCA1", "clinvar_clinical_significance": "Pathogenic"}
        ],
    }

    return state


def create_test_state_high_risk_invalid():
    """Create test state with high risk NOT driven by pathogenic variant (should FLAG)."""
    # Similar setup but train model to give high risk for benign variants
    fusion_layer = FusionLayer(model_type="gradient_boosting")

    training_data = []
    np.random.seed(43)

    # Train model incorrectly - high risk for benign
    for _ in range(100):
        inputs = StaticModelInputs(
            prs_score=np.random.uniform(0.7, 0.9),
            clinvar_classification="benign",
            cadd_score=np.random.uniform(5, 15),
            tcga_enrichment=np.random.uniform(3, 6),
            gene_burden_score=np.random.uniform(2, 4),
        )
        risk_score = 0.8  # High risk for benign (incorrect)
        training_data.append((inputs, risk_score))

    fusion_layer.train(training_data, validation_split=0)

    # Test with benign variant
    test_inputs = [
        StaticModelInputs(
            prs_score=0.8, clinvar_classification="benign", cadd_score=10.0, tcga_enrichment=4.0, gene_burden_score=3.0
        )
    ]

    X_encoded = fusion_layer._encode_features(test_inputs)
    X_scaled = fusion_layer.scaler.transform(X_encoded)

    state = {
        "ml_fusion_model_instance": fusion_layer,
        "ml_fusion_feature_matrix": X_scaled,
        "risk_scores": {"aggregate_risk_score": 0.75},  # High risk
        "filtered_variants": [
            {"variant_id": "chr2:67890", "gene": "TP53", "clinvar_clinical_significance": "Benign"},
            {
                "variant_id": "chr17:41276045",
                "gene": "BRCA1",
                "clinvar_clinical_significance": "Pathogenic",  # Has pathogenic but model ignores it
            },
        ],
    }

    return state


def create_test_state_low_risk_with_pathogenic():
    """Create test state with low risk despite pathogenic variant (should FLAG)."""
    fusion_layer = FusionLayer(model_type="gradient_boosting")

    training_data = []
    np.random.seed(44)

    # Train model to give low risk even for pathogenic
    for _ in range(100):
        inputs = StaticModelInputs(
            prs_score=np.random.uniform(0.1, 0.3),
            clinvar_classification="pathogenic",
            cadd_score=np.random.uniform(20, 35),
            tcga_enrichment=np.random.uniform(0.5, 1.5),
            gene_burden_score=np.random.uniform(0, 1),
        )
        risk_score = 0.05  # Low risk despite pathogenic
        training_data.append((inputs, risk_score))

    fusion_layer.train(training_data, validation_split=0)

    test_inputs = [
        StaticModelInputs(
            prs_score=0.2,
            clinvar_classification="pathogenic",
            cadd_score=25.0,
            tcga_enrichment=1.0,
            gene_burden_score=0.5,
        )
    ]

    X_encoded = fusion_layer._encode_features(test_inputs)
    X_scaled = fusion_layer.scaler.transform(X_encoded)

    state = {
        "ml_fusion_model_instance": fusion_layer,
        "ml_fusion_feature_matrix": X_scaled,
        "risk_scores": {"aggregate_risk_score": 0.05},  # Low risk
        "filtered_variants": [
            {"variant_id": "chr3:11111", "gene": "MLH1", "clinvar_clinical_significance": "Pathogenic"}
        ],
    }

    return state


def test_shap_validator():
    """Run comprehensive tests on SHAP validator."""
    print("=" * 60)
    print("Testing SHAP Validator Node")
    print("=" * 60)

    # Test 1: High risk with pathogenic variant (should PASS)
    print("\n1. Testing high risk with pathogenic variant...")
    state1 = create_test_state_high_risk_valid()
    result1 = shap_validator_process(state1)

    print(f"   Status: {result1['shap_validation_status']}")
    print(f"   Reasons: {result1['shap_validation_reasons']}")
    print("   Top contributors:")
    for contrib in result1["shap_top_contributors"][:3]:
        print(f"     - {contrib['display_name']}: {contrib['direction']} risk (SHAP: {contrib['shap_value']:.3f})")

    # Debug: Print all feature importances
    print("\n   All feature importances:")
    for feature, importance in result1["shap_feature_importance"].items():
        print(f"     - {feature}: {importance:.3f}")

    assert result1["shap_validation_status"] == "PASS", "Should PASS for valid high risk"
    print("   ✅ Test passed!")

    # Test 2: High risk without pathogenic variant (should FLAG)
    print("\n2. Testing high risk without pathogenic variant...")
    state2 = create_test_state_high_risk_invalid()
    result2 = shap_validator_process(state2)

    print(f"   Status: {result2['shap_validation_status']}")
    print(f"   Reasons: {result2['shap_validation_reasons']}")
    print("   Top contributors:")
    for contrib in result2["shap_top_contributors"][:3]:
        print(f"     - {contrib['display_name']}: {contrib['direction']} risk")

    assert result2["shap_validation_status"] == "FLAG_FOR_REVIEW", "Should FLAG invalid high risk"
    print("   ✅ Test passed!")

    # Test 3: Low risk with pathogenic variant (should FLAG)
    print("\n3. Testing low risk with pathogenic variant...")
    state3 = create_test_state_low_risk_with_pathogenic()
    result3 = shap_validator_process(state3)

    print(f"   Status: {result3['shap_validation_status']}")
    print(f"   Reasons: {result3['shap_validation_reasons']}")
    print("   Top contributors:")
    for contrib in result3["shap_top_contributors"][:3]:
        print(f"     - {contrib['display_name']}: {contrib['direction']} risk")

    assert result3["shap_validation_status"] == "FLAG_FOR_REVIEW", "Should FLAG low risk with pathogenic"
    print("   ✅ Test passed!")

    # Test 4: Missing model (should handle gracefully)
    print("\n4. Testing missing model scenario...")
    state4 = {"risk_scores": {"aggregate_risk_score": 0.5}}
    result4 = shap_validator_process(state4)

    print(f"   Status: {result4['shap_validation_status']}")
    print(f"   Reasons: {result4['shap_validation_reasons']}")

    assert result4["shap_validation_status"] == "SKIPPED", "Should SKIP when model missing"
    print("   ✅ Test passed!")

    print("\n" + "=" * 60)
    print("All SHAP validator tests passed! ✅")
    print("=" * 60)


if __name__ == "__main__":
    test_shap_validator()
