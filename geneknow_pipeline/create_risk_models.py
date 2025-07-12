"""
Create simple risk prediction models for cancer types.
Uses logistic regression based on variant presence in key genes.
"""

import os
import pickle
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import json


def create_risk_models():
    """Create and save simple risk prediction models."""

    models_dir = "models"
    os.makedirs(models_dir, exist_ok=True)

    # Define feature genes for each cancer type
    cancer_genes = {
        "breast": ["BRCA1", "BRCA2", "TP53", "PIK3CA", "PALB2", "ATM", "CHEK2"],
        "colon": ["APC", "KRAS", "TP53", "PIK3CA", "SMAD4", "BRAF", "MSH2", "MLH1"],
        "lung": ["TP53", "KRAS", "EGFR", "ALK", "ROS1", "BRAF", "MET"],
        "prostate": ["AR", "PTEN", "TP53", "BRCA2", "ATM", "FOXA1", "SPOP"],
        "blood": ["JAK2", "FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2", "RUNX1"],
    }

    # Create synthetic training data
    # In production, this would come from real patient data
    np.random.seed(42)

    models = {}
    scalers = {}

    for cancer_type, genes in cancer_genes.items():
        print(f"\n🔨 Creating {cancer_type} cancer risk model...")

        # Generate synthetic training data
        n_samples = 1000
        n_features = len(genes) + 2  # genes + age + sex

        # Features: [gene1_variant, gene2_variant, ..., age, sex]
        X = np.zeros((n_samples, n_features))

        # Simulate variant presence (binary)
        for i in range(len(genes)):
            # Higher-risk genes have higher variant probability in cases
            if i < 3:  # First 3 genes are high-risk
                X[:, i] = np.random.binomial(1, 0.15, n_samples)
            else:
                X[:, i] = np.random.binomial(1, 0.05, n_samples)

        # Age (normalized 20-80)
        X[:, -2] = np.random.uniform(20, 80, n_samples) / 100

        # Sex (0=M, 1=F)
        X[:, -1] = np.random.binomial(1, 0.5, n_samples)

        # Generate labels based on risk factors
        # Simple model: more variants = higher risk
        risk_score = np.sum(X[:, :3] * [0.4, 0.3, 0.3], axis=1)  # High-risk genes
        risk_score += np.sum(X[:, 3:-2] * 0.1, axis=1)  # Other genes
        risk_score += X[:, -2] * 0.2  # Age factor

        # Add noise and create binary labels
        risk_score += np.random.normal(0, 0.1, n_samples)
        y = (risk_score > 0.3).astype(int)

        # Balance the dataset
        positive_ratio = np.mean(y)
        print(f"  Dataset balance: {positive_ratio:.1%} positive cases")

        # Train model
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        model = LogisticRegression(random_state=42, max_iter=1000)
        model.fit(X_scaled, y)

        # Evaluate on training data (in production, use cross-validation)
        accuracy = model.score(X_scaled, y)
        print(f"  Training accuracy: {accuracy:.1%}")

        # Save model and scaler
        model_path = os.path.join(models_dir, f"{cancer_type}_model.pkl")
        scaler_path = os.path.join(models_dir, f"{cancer_type}_scaler.pkl")

        with open(model_path, "wb") as f:
            pickle.dump(model, f)
        with open(scaler_path, "wb") as f:
            pickle.dump(scaler, f)

        models[cancer_type] = model
        scalers[cancer_type] = scaler

        print(f"  ✅ Saved model to {model_path}")

    # Save feature configuration
    config_path = os.path.join(models_dir, "model_config.json")
    config = {
        "cancer_genes": cancer_genes,
        "feature_order": {
            cancer: genes + ["age", "sex"] for cancer, genes in cancer_genes.items()
        },
        "model_type": "logistic_regression",
        "version": "1.0.0",
    }

    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)

    print(f"\n✅ Created risk models for {len(models)} cancer types")
    print(f"✅ Saved configuration to {config_path}")

    # Test prediction
    print("\n🧪 Testing breast cancer model:")
    test_features = np.zeros(len(cancer_genes["breast"]) + 2)
    test_features[0] = 1  # BRCA1 variant present
    test_features[-2] = 0.45  # Age 45
    test_features[-1] = 1  # Female

    test_scaled = scalers["breast"].transform([test_features])
    risk_prob = models["breast"].predict_proba(test_scaled)[0, 1]
    print(f"  BRCA1 variant carrier (45F): {risk_prob*100:.1f}% risk")

    return models_dir


if __name__ == "__main__":
    create_risk_models()
