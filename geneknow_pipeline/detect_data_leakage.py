#!/usr/bin/env python3
"""
Data Leakage Detection Script

Analyzes the suspicious perfect scores in our fusion layer training
to identify data leakage issues.
"""

import pandas as pd
import numpy as np
import sqlite3


def analyze_data_leakage():
    """Detect and explain the data leakage in our training."""

    print("🚨 GeneKnow Data Leakage Analysis")
    print("=" * 50)

    print("📊 Our Training Process:")
    print("1. Load variants from database")
    print("2. Generate static model features")
    print("3. Calculate risk targets using SAME features")
    print("4. Train model to predict targets from features")
    print("")

    print("🔍 Risk Target Calculation:")
    print("We calculate risk_score as:")
    print("  risk = 0.5 * clinvar_score")
    print("       + 0.2 * cadd_score")
    print("       + 0.15 * tcga_score")
    print("       + 0.1 * gene_burden_score")
    print("       + 0.05 * prs_score")
    print("")

    print("🎯 Model Input Features:")
    print("  - clinvar_classification (same as above)")
    print("  - cadd_score (same as above)")
    print("  - tcga_enrichment (same as above)")
    print("  - gene_burden_score (same as above)")
    print("  - prs_score (same as above)")
    print("")

    print("🚨 THE PROBLEM:")
    print("We are asking the model to predict a LINEAR COMBINATION")
    print("of its own input features! This is PERFECT DATA LEAKAGE.")
    print("")
    print("It's like asking someone to predict the result of:")
    print("  Y = 2*A + 3*B + 4*C")
    print("And then giving them A, B, C as inputs.")
    print("Of course they'll get it perfectly right!")
    print("")

    # Load some real data to demonstrate
    conn = sqlite3.connect("population_variants.db")

    # Get a small sample
    query = """
    SELECT clinical_significance, gnomad_af, consequence
    FROM population_variants
    WHERE clinical_significance IS NOT NULL
    LIMIT 5
    """

    sample_df = pd.read_sql_query(query, conn)
    print("📋 Sample Data:")
    print(sample_df)
    print("")

    # Demonstrate the leakage with a simple example
    print("💡 Simple Example of the Leakage:")
    print("If ClinVar = 'Pathogenic':")
    print("  -> Risk target gets +0.5 (50% weight)")
    print("  -> Model sees clinvar_pathogenic = 1 as input")
    print("  -> Model learns: if clinvar_pathogenic=1, predict ~0.5")
    print("")
    print("This is why ClinVar has 75.7% feature importance!")
    print("The model is just reconstructing our formula.")
    print("")

    print("✅ SOLUTIONS:")
    print("1. Use independent ground truth labels (not derived from features)")
    print("2. Use external validation dataset with real clinical outcomes")
    print("3. Predict actual clinical endpoints (cancer diagnosis, survival)")
    print("4. Use held-out features not included in target calculation")
    print("5. Train on one dataset, validate on completely different cohort")

    conn.close()


def demonstrate_perfect_prediction():
    """Show how we can perfectly predict our target."""

    print("\n" + "=" * 50)
    print("🔬 DEMONSTRATION: Perfect Prediction")
    print("=" * 50)

    # Create simple example
    np.random.seed(42)

    # Simulate our exact process
    print("Creating simulated data using our EXACT process:")

    n_samples = 1000

    # Generate features
    prs = np.random.beta(2, 5, n_samples)
    clinvar = np.random.choice(["pathogenic", "benign", "uncertain"], n_samples, p=[0.2, 0.6, 0.2])
    cadd = np.random.exponential(10, n_samples)
    tcga = np.random.lognormal(0, 1, n_samples)
    burden = np.random.poisson(2, n_samples)

    # Calculate target using our EXACT formula
    risk_targets = np.zeros(n_samples)

    for i in range(n_samples):
        risk = 0.0

        # ClinVar contribution (50% weight) - SAME LOGIC AS TRAINING
        if clinvar[i] == "pathogenic":
            risk += 0.5 * 1.0
        elif clinvar[i] == "benign":
            risk += 0.5 * 0.0
        elif clinvar[i] == "uncertain":
            risk += 0.5 * 0.3

        # CADD contribution (20% weight)
        if cadd[i] > 25:
            risk += 0.2 * 1.0
        elif cadd[i] > 15:
            risk += 0.2 * 0.6
        elif cadd[i] > 10:
            risk += 0.2 * 0.3

        # Other contributions...
        risk += 0.15 * min(tcga[i] / 10, 1.0)  # TCGA
        risk += 0.1 * min(burden[i] / 10, 1.0)  # Burden
        risk += 0.05 * prs[i]  # PRS

        risk_targets[i] = min(risk, 1.0)

    print(f"Generated {n_samples} samples")
    print("Risk target stats:")
    print(f"  Mean: {np.mean(risk_targets):.3f}")
    print(f"  Std: {np.std(risk_targets):.3f}")
    print(f"  Range: {np.min(risk_targets):.3f} - {np.max(risk_targets):.3f}")
    print("")

    print("Now if we train a model using prs, clinvar, cadd, tcga, burden")
    print("to predict these risk_targets...")
    print("")
    print("💥 OF COURSE IT WILL BE PERFECT!")
    print("The model just needs to learn our deterministic formula!")
    print("")
    print("This explains:")
    print("  ✅ Perfect AUC scores (1.0000)")
    print("  ✅ Perfect accuracy (1.0000)")
    print("  ✅ Near-perfect R² (0.9999)")
    print("  ✅ Tiny MSE (0.000007)")
    print("  ✅ ClinVar dominance (75.7% importance)")


if __name__ == "__main__":
    analyze_data_leakage()
    demonstrate_perfect_prediction()

    print("\n" + "=" * 50)
    print("🚨 CONCLUSION: MASSIVE DATA LEAKAGE DETECTED")
    print("=" * 50)
    print("Our 'perfect' results are completely invalid.")
    print("We need to redesign our training approach.")
    print("The model learned to reconstruct our target formula,")
    print("not to make real genomic risk predictions.")
