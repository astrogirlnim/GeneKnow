#!/usr/bin/env python3
"""
Train ML Models WITHOUT Data Leakage
Quick test to see real performance on genomic features only.
"""

import logging
import sys
import os
<<<<<<< HEAD
=======
from datetime import datetime
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from ml_trainer import GenomicMLTrainer
from ml_feature_extractor_no_leakage import GenomicFeatureExtractorNoLeakage

<<<<<<< HEAD

def main():
    """Train ML models WITHOUT clinical significance features (no data leakage)."""

    # Setup logging
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )
    logger = logging.getLogger(__name__)

=======
def main():
    """Train ML models WITHOUT clinical significance features (no data leakage)."""
    
    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    print("🧬 GeneKnow ML Training - NO DATA LEAKAGE VERSION")
    print("=" * 60)
    print("Training on genomic features ONLY (no clinical significance)")
    print("=" * 60)
<<<<<<< HEAD

=======
    
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    # Check if database exists
    if not os.path.exists("population_variants.db"):
        print("❌ Database not found!")
        print("Please run: python create_population_database.py --cancer-genes-only")
        return
<<<<<<< HEAD

=======
    
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    try:
        # Initialize trainer with no-leakage feature extractor
        trainer = GenomicMLTrainer(random_state=42)
        trainer.feature_extractor = GenomicFeatureExtractorNoLeakage()
<<<<<<< HEAD

        # Load and prepare data (NO clinical significance features)
        logger.info("Loading training data WITHOUT clinical significance features...")
        X, y = trainer.feature_extractor.load_training_data()

        print("\n📊 Training Dataset Overview (NO LEAKAGE):")
        print(f"   Total variants: {len(X):,}")
        print(
            f"   Features: {X.shape[1]} (reduced from ~25 to exclude clinical features)"
        )
        print(f"   Pathogenic variants: {sum(y):,} ({sum(y)/len(y)*100:.1f}%)")
        print(
            f"   Benign variants: {len(y)-sum(y):,} ({(len(y)-sum(y))/len(y)*100:.1f}%)"
        )

        # Show features being used
        print("\n🔬 Features Used for Training:")
        for i, feature in enumerate(X.columns, 1):
            print(f"   {i:2d}. {feature}")

        # Check that no clinical features are present
        clinical_features = [col for col in X.columns if "clinvar" in col.lower()]
        if clinical_features:
            print(f"⚠️  WARNING: Clinical features still present: {clinical_features}")
        else:
            print("✅ Confirmed: NO clinical significance features (no data leakage)")

        # Prepare for training
        logger.info("Preparing data for ML training...")
        X_train, X_test, y_train, y_test = (
            trainer.feature_extractor.prepare_for_training(X, y)
        )

        print("\n🎯 Training/Test Split:")
        print(f"   Training samples: {len(X_train):,}")
        print(f"   Test samples: {len(X_test):,}")

        # Train models (faster - just test key algorithms)
        print("\n🚀 Training ML Models (genomic features only)...")
        print("   Testing: Random Forest, Gradient Boosting, Logistic Regression")

        # Quick training with fewer models for faster testing
        results = trainer.train_models(
            X_train,
            X_test,
            y_train,
            y_test,
            sampling_methods=["class_weight"],  # Just class weighting for speed
            cv_folds=3,  # Reduced CV folds for faster training
        )

        # Display realistic results
        print("\n📈 REALISTIC ML Performance (NO DATA LEAKAGE):")
        print(f"   Best model: {trainer.best_model_name}")

=======
        
        # Load and prepare data (NO clinical significance features)
        logger.info("Loading training data WITHOUT clinical significance features...")
        X, y = trainer.feature_extractor.load_training_data()
        
        print(f"\n📊 Training Dataset Overview (NO LEAKAGE):")
        print(f"   Total variants: {len(X):,}")
        print(f"   Features: {X.shape[1]} (reduced from ~25 to exclude clinical features)")
        print(f"   Pathogenic variants: {sum(y):,} ({sum(y)/len(y)*100:.1f}%)")
        print(f"   Benign variants: {len(y)-sum(y):,} ({(len(y)-sum(y))/len(y)*100:.1f}%)")
        
        # Show features being used
        print(f"\n🔬 Features Used for Training:")
        for i, feature in enumerate(X.columns, 1):
            print(f"   {i:2d}. {feature}")
        
        # Check that no clinical features are present
        clinical_features = [col for col in X.columns if 'clinvar' in col.lower()]
        if clinical_features:
            print(f"⚠️  WARNING: Clinical features still present: {clinical_features}")
        else:
            print(f"✅ Confirmed: NO clinical significance features (no data leakage)")
        
        # Prepare for training
        logger.info("Preparing data for ML training...")
        X_train, X_test, y_train, y_test = trainer.feature_extractor.prepare_for_training(X, y)
        
        print(f"\n🎯 Training/Test Split:")
        print(f"   Training samples: {len(X_train):,}")
        print(f"   Test samples: {len(X_test):,}")
        
        # Train models (faster - just test key algorithms)
        print(f"\n🚀 Training ML Models (genomic features only)...")
        print("   Testing: Random Forest, Gradient Boosting, Logistic Regression")
        
        # Quick training with fewer models for faster testing
        results = trainer.train_models(
            X_train, X_test, y_train, y_test,
            sampling_methods=['class_weight'],  # Just class weighting for speed
            cv_folds=3  # Reduced CV folds for faster training
        )
        
        # Display realistic results
        print(f"\n📈 REALISTIC ML Performance (NO DATA LEAKAGE):")
        print(f"   Best model: {trainer.best_model_name}")
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        # Find best performance
        best_performance = None
        for sampling_method, method_results in results.items():
            for model_name, metrics in method_results.items():
<<<<<<< HEAD
                if (
                    best_performance is None
                    or metrics["roc_auc"] > best_performance["roc_auc"]
                ):
                    best_performance = metrics
                    best_performance["full_name"] = f"{model_name}_{sampling_method}"

        if best_performance:
            print(
                f"   ROC AUC: {best_performance['roc_auc']:.3f} (realistic genomic performance)"
            )
            print(f"   Balanced Accuracy: {best_performance['balanced_accuracy']:.3f}")
            print(f"   F1 Score: {best_performance['f1_score']:.3f}")
            print(
                f"   Matthews Correlation: {best_performance['matthews_corrcoef']:.3f}"
            )

            # Interpret performance
            auc = best_performance["roc_auc"]
=======
                if best_performance is None or metrics['roc_auc'] > best_performance['roc_auc']:
                    best_performance = metrics
                    best_performance['full_name'] = f"{model_name}_{sampling_method}"
        
        if best_performance:
            print(f"   ROC AUC: {best_performance['roc_auc']:.3f} (realistic genomic performance)")
            print(f"   Balanced Accuracy: {best_performance['balanced_accuracy']:.3f}")
            print(f"   F1 Score: {best_performance['f1_score']:.3f}")
            print(f"   Matthews Correlation: {best_performance['matthews_corrcoef']:.3f}")
            
            # Interpret performance
            auc = best_performance['roc_auc']
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
            if auc > 0.85:
                interpretation = "🎉 Excellent performance!"
            elif auc > 0.75:
                interpretation = "✅ Good performance for genomics!"
            elif auc > 0.65:
                interpretation = "👍 Decent performance, could be improved"
            else:
                interpretation = "📈 Needs more feature engineering"
<<<<<<< HEAD

            print(f"   Interpretation: {interpretation}")

=======
            
            print(f"   Interpretation: {interpretation}")
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
        # Save models
        logger.info("Saving trained models...")
        output_dir = "ml_models_no_leakage"
        trainer.save_models(output_dir)
<<<<<<< HEAD

        print("\n✅ Training Complete (NO DATA LEAKAGE)!")
        print(f"   Models saved to: {output_dir}/")
        print(f"   Best model: {output_dir}/best_model.pkl")

        # Show feature importance if available
        if hasattr(trainer.best_model, "feature_importances_"):
            print("\n🔍 Top 10 Most Important Features:")
            importances = trainer.best_model.feature_importances_
            feature_names = trainer.feature_extractor.feature_columns

            # Get top 10 features
            indices = sorted(
                range(len(importances)), key=lambda i: importances[i], reverse=True
            )[:10]

            for i, idx in enumerate(indices, 1):
                feat_name = (
                    feature_names[idx] if idx < len(feature_names) else f"feature_{idx}"
                )
                importance = importances[idx]
                print(f"   {i:2d}. {feat_name:<30} ({importance:.3f})")

        # Compare to expected performance
        print("\n📚 Performance Context:")
        print(
            "   - ROC AUC 0.60-0.70: Typical for genomic prediction without clinical data"
        )
        print("   - ROC AUC 0.70-0.80: Good genomic prediction")
        print("   - ROC AUC 0.80+: Excellent genomic prediction")
        print(f"   - Your result: {best_performance['roc_auc']:.3f}")

        # Next steps
        print("\n🔬 Next Steps:")
        print("   1. This model can now be used for REAL pathogenicity prediction")
        print("   2. Clinical significance can be used for validation/confidence")
        print("   3. Integration into pipeline maintains both ML + clinical insights")

=======
        
        print(f"\n✅ Training Complete (NO DATA LEAKAGE)!")
        print(f"   Models saved to: {output_dir}/")
        print(f"   Best model: {output_dir}/best_model.pkl")
        
        # Show feature importance if available
        if hasattr(trainer.best_model, 'feature_importances_'):
            print(f"\n🔍 Top 10 Most Important Features:")
            importances = trainer.best_model.feature_importances_
            feature_names = trainer.feature_extractor.feature_columns
            
            # Get top 10 features
            indices = sorted(range(len(importances)), key=lambda i: importances[i], reverse=True)[:10]
            
            for i, idx in enumerate(indices, 1):
                feat_name = feature_names[idx] if idx < len(feature_names) else f"feature_{idx}"
                importance = importances[idx]
                print(f"   {i:2d}. {feat_name:<30} ({importance:.3f})")
        
        # Compare to expected performance
        print(f"\n📚 Performance Context:")
        print(f"   - ROC AUC 0.60-0.70: Typical for genomic prediction without clinical data")
        print(f"   - ROC AUC 0.70-0.80: Good genomic prediction")
        print(f"   - ROC AUC 0.80+: Excellent genomic prediction")
        print(f"   - Your result: {best_performance['roc_auc']:.3f}")
        
        # Next steps
        print(f"\n🔬 Next Steps:")
        print(f"   1. This model can now be used for REAL pathogenicity prediction")
        print(f"   2. Clinical significance can be used for validation/confidence")
        print(f"   3. Integration into pipeline maintains both ML + clinical insights")
        
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    except Exception as e:
        logger.error(f"Training failed: {e}")
        print(f"❌ Training failed: {e}")
        return 1
<<<<<<< HEAD

=======
    
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
    return 0


if __name__ == "__main__":
<<<<<<< HEAD
    exit(main())
=======
    exit(main()) 
>>>>>>> 2c09325 (ML Fusion Model Integration & Pipeline Enhancements (#32))
