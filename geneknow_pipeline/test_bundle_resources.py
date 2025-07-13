#!/usr/bin/env python3
"""
Test script to verify all required resources are available for production bundle.
Run this to ensure the bundle has everything needed for the pipeline to work.
"""

import os
import sys
import sqlite3
import pickle
from pathlib import Path

def check_resource(path, description, required=True):
    """Check if a resource exists and report status."""
    exists = os.path.exists(path)
    status = "[OK]" if exists else ("[FAIL]" if required else "[WARN]")
    
    if exists:
        size = os.path.getsize(path)
        if size > 1024 * 1024:  # > 1MB
            size_str = f"{size / (1024 * 1024):.1f}MB"
        elif size > 1024:  # > 1KB
            size_str = f"{size / 1024:.1f}KB"
        else:
            size_str = f"{size}B"
        print(f"{status} {description}: {path} ({size_str})")
    else:
        print(f"{status} {description}: {path} (NOT FOUND)")
    
    return exists

def check_database(path, description):
    """Check database and report table count."""
    if not os.path.exists(path):
        print(f"[FAIL] {description}: {path} (NOT FOUND)")
        return False
    
    try:
        conn = sqlite3.connect(path)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = cursor.fetchall()
        
        # Get row count for main table
        main_table = tables[0][0] if tables else None
        row_count = 0
        if main_table:
            cursor.execute(f"SELECT COUNT(*) FROM {main_table}")
            row_count = cursor.fetchone()[0]
        
        conn.close()
        
        size = os.path.getsize(path) / (1024 * 1024)  # MB
        print(f"✅ {description}: {path} ({size:.1f}MB, {len(tables)} tables, {row_count:,} rows)")
        return True
    except Exception as e:
        print(f"❌ {description}: {path} (ERROR: {e})")
        return False

def main():
    """Check all required resources."""
    print("[CHECK] GeneKnow Bundle Resource Verification\n")
    
    # Get base directory (either bundle or development)
    base_dir = Path(__file__).parent
    print(f"Base directory: {base_dir}\n")
    
    all_good = True
    
    # 1. Check ML Fusion Models
    print("📊 ML Fusion Models:")
    fusion_models = [
        ("ml_models/best_fusion_model.pkl", "Main fusion model", True),
        ("ml_models/fusion_gradient_boosting.pkl", "Gradient boosting variant", False),
        ("ml_models/fusion_random_forest.pkl", "Random forest variant", False),
        ("ml_models/fusion_linear.pkl", "Linear variant", False),
    ]
    
    for path, desc, required in fusion_models:
        if not check_resource(base_dir / path, desc, required) and required:
            all_good = False
    
    # Check alternative fusion models
    print("\n📊 Alternative Fusion Models:")
    alt_models = [
        "best_fusion_model_FIXED.pkl",
        "best_fusion_model_real_data.pkl",
        "fusion_gradient_boosting_FIXED.pkl",
        "fusion_gradient_boosting_real_data.pkl",
        "fusion_linear_FIXED.pkl",
        "fusion_linear_real_data.pkl",
        "fusion_random_forest_FIXED.pkl",
        "fusion_random_forest_real_data.pkl"
    ]
    
    for model in alt_models:
        check_resource(base_dir / model, f"Alternative: {model}", False)
    
    # 2. Check No-Leakage Models
    print("\n🛡️ No-Leakage Models:")
    no_leak_models = [
        ("ml_models_no_leakage/best_model.pkl", "Best model", True),
        ("ml_models_no_leakage/gradient_boosting_class_weight.pkl", "Gradient boosting", True),
        ("ml_models_no_leakage/logistic_regression_class_weight.pkl", "Logistic regression", True),
        ("ml_models_no_leakage/random_forest_class_weight.pkl", "Random forest", True),
        ("ml_models_no_leakage/svm_class_weight.pkl", "SVM", True),
        ("ml_models_no_leakage/naive_bayes_class_weight.pkl", "Naive Bayes", True),
        ("ml_models_no_leakage/feature_columns.pkl", "Feature columns", True),
        ("ml_models_no_leakage/label_encoders.pkl", "Label encoders", True),
        ("ml_models_no_leakage/scaler.pkl", "Scaler", True),
        ("ml_models_no_leakage/model_metadata.json", "Model metadata", True),
    ]
    
    for path, desc, required in no_leak_models:
        if not check_resource(base_dir / path, desc, required) and required:
            all_good = False
    
    # 3. Check Databases
    print("\n🗄️ Databases:")
    databases = [
        ("population_variants.db", "Population variants database"),
        ("clinvar_annotations.db", "ClinVar annotations database"),
        ("prs_snps.db", "PRS SNPs database"),
    ]
    
    for db_path, desc in databases:
        if not check_database(base_dir / db_path, desc):
            all_good = False
    
    # 4. Check Other Resources
    print("\n📁 Other Resources:")
    other_resources = [
        ("test_reference/test_genome.fa", "Test reference genome", False),
        ("create_simple_fusion_model.py", "Fusion model creator", True),
        ("create_all_fusion_models.py", "All models creator", False),
        ("create_population_database.py", "Database creator", True),
    ]
    
    for path, desc, required in other_resources:
        if not check_resource(base_dir / path, desc, required) and required:
            all_good = False
    
    # Summary
    print("\n" + "="*50)
    if all_good:
        print("✅ All required resources are present!")
        print("The bundle is ready for production deployment.")
    else:
        print("❌ Some required resources are missing!")
        print("Run the bundle script again or create missing resources.")
        sys.exit(1)

if __name__ == "__main__":
    main() 