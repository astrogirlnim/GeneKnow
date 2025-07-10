#!/usr/bin/env python3
"""
Train ML Models for GeneKnow Risk Fusion
Quick start script to train ML models using existing database.
"""

import logging
import os
import sys
from datetime import datetime

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from ml_trainer import GenomicMLTrainer
from ml_feature_extractor import GenomicFeatureExtractor

def main():
    """Train ML models for genomic risk prediction."""
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('ml_training.log'),
            logging.StreamHandler()
        ]
    )
    
    logger = logging.getLogger(__name__)
    
    print("🧬 GeneKnow ML Model Training")
    print("=" * 50)
    
    # Check if database exists
    if not os.path.exists("population_variants.db"):
        print("❌ Database not found!")
        print("Please run: python create_population_database.py --cancer-genes-only")
        print("This will create the training database from ClinVar data.")
        return
    
    try:
        # Initialize components
        logger.info("Initializing ML training components...")
        trainer = GenomicMLTrainer(random_state=42)
        
        # Load and prepare data
        logger.info("Loading training data from database...")
        X, y = trainer.feature_extractor.load_training_data()
        
        print(f"📊 Training Dataset Overview:")
        print(f"   Total variants: {len(X):,}")
        print(f"   Features: {X.shape[1]}")
        print(f"   Pathogenic variants: {sum(y):,} ({sum(y)/len(y)*100:.1f}%)")
        print(f"   Benign variants: {len(y)-sum(y):,} ({(len(y)-sum(y))/len(y)*100:.1f}%)")
        print(f"   Class imbalance ratio: {(len(y)-sum(y))/sum(y):.1f}:1")
        
        # Prepare for training
        logger.info("Preparing data for ML training...")
        X_train, X_test, y_train, y_test = trainer.feature_extractor.prepare_for_training(X, y)
        
        print(f"\n🎯 Training/Test Split:")
        print(f"   Training samples: {len(X_train):,}")
        print(f"   Test samples: {len(X_test):,}")
        
        # Train models with different approaches
        print(f"\n🚀 Training ML Models...")
        print("   This may take several minutes...")
        
        results = trainer.train_models(
            X_train, X_test, y_train, y_test,
            sampling_methods=['none', 'smote', 'class_weight'],
            cv_folds=5
        )
        
        # Display results
        print(f"\n📈 Training Results:")
        print(f"   Best model: {trainer.best_model_name}")
        
        # Find best performance
        best_performance = None
        for sampling_method, method_results in results.items():
            for model_name, metrics in method_results.items():
                if best_performance is None or metrics['roc_auc'] > best_performance['roc_auc']:
                    best_performance = metrics
                    best_performance['full_name'] = f"{model_name}_{sampling_method}"
        
        if best_performance:
            print(f"   Best ROC AUC: {best_performance['roc_auc']:.3f}")
            print(f"   Best Balanced Accuracy: {best_performance['balanced_accuracy']:.3f}")
            print(f"   Best F1 Score: {best_performance['f1_score']:.3f}")
            print(f"   Matthews Correlation: {best_performance['matthews_corrcoef']:.3f}")
        
        # Save models
        logger.info("Saving trained models...")
        trainer.save_models()
        
        # Generate plots
        logger.info("Generating performance plots...")
        try:
            trainer.plot_results()
        except Exception as e:
            logger.warning(f"Plot generation failed: {e}")
        
        print(f"\n✅ Training Complete!")
        print(f"   Models saved to: ml_models/")
        print(f"   Best model: ml_models/best_model.pkl")
        print(f"   Training log: ml_training.log")
        
        # Show next steps
        print(f"\n🔬 Next Steps:")
        print(f"   1. Review training results in ml_models/model_metadata.json")
        print(f"   2. Test the model with: python test_ml_predictions.py")
        print(f"   3. Integrate into pipeline by updating nodes/feature_vector_builder.py")
        
    except Exception as e:
        logger.error(f"Training failed: {e}")
        print(f"❌ Training failed: {e}")
        print(f"Check ml_training.log for details")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main()) 