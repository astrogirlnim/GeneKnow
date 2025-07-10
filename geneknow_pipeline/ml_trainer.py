#!/usr/bin/env python3
"""
ML Trainer for GeneKnow Risk Fusion
Handles class imbalance and trains multiple models for genomic risk prediction.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
import logging
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import (
    classification_report, confusion_matrix, roc_auc_score, 
    precision_recall_curve, average_precision_score, roc_curve,
    balanced_accuracy_score, matthews_corrcoef, f1_score
)
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.utils.class_weight import compute_class_weight
from imblearn.over_sampling import SMOTE, ADASYN
from imblearn.under_sampling import RandomUnderSampler
from imblearn.combine import SMOTEENN
import pickle
import os
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import json

from ml_feature_extractor import GenomicFeatureExtractor

logger = logging.getLogger(__name__)

class GenomicMLTrainer:
    """
    Train ML models for genomic risk prediction with proper handling of class imbalance.
    """
    
    def __init__(self, random_state: int = 42):
        self.random_state = random_state
        self.models = {}
        self.best_model = None
        self.best_model_name = None
        self.feature_extractor = GenomicFeatureExtractor()
        self.training_history = {}
        
    def _get_models(self, class_weight: Optional[Dict] = None) -> Dict:
        """Get a dictionary of models to train."""
        models = {
            'logistic_regression': LogisticRegression(
                random_state=self.random_state,
                class_weight=class_weight or 'balanced',
                max_iter=1000,
                solver='liblinear'
            ),
            'random_forest': RandomForestClassifier(
                n_estimators=100,
                random_state=self.random_state,
                class_weight=class_weight or 'balanced',
                max_depth=10,
                min_samples_split=5,
                min_samples_leaf=2
            ),
            'gradient_boosting': GradientBoostingClassifier(
                n_estimators=100,
                random_state=self.random_state,
                learning_rate=0.1,
                max_depth=5,
                min_samples_split=10,
                min_samples_leaf=4
            ),
            'svm': SVC(
                random_state=self.random_state,
                class_weight=class_weight or 'balanced',
                probability=True,
                kernel='rbf',
                C=1.0
            ),
            'naive_bayes': GaussianNB()
        }
        return models
    
    def _handle_class_imbalance(self, X: np.ndarray, y: np.ndarray, 
                               method: str = 'smote') -> Tuple[np.ndarray, np.ndarray]:
        """Handle class imbalance using various sampling techniques."""
        
        logger.info(f"Original class distribution: {np.bincount(y)}")
        
        if method == 'smote':
            sampler = SMOTE(random_state=self.random_state, k_neighbors=5)
        elif method == 'adasyn':
            sampler = ADASYN(random_state=self.random_state)
        elif method == 'smoteenn':
            sampler = SMOTEENN(random_state=self.random_state)
        elif method == 'undersample':
            sampler = RandomUnderSampler(random_state=self.random_state)
        else:
            logger.info("No sampling applied")
            return X, y
        
        try:
            X_resampled, y_resampled = sampler.fit_resample(X, y)
            logger.info(f"Resampled class distribution: {np.bincount(y_resampled)}")
            return X_resampled, y_resampled
        except Exception as e:
            logger.warning(f"Sampling failed: {e}. Using original data.")
            return X, y
    
    def _calculate_class_weights(self, y: np.ndarray) -> Dict:
        """Calculate class weights for imbalanced learning."""
        classes = np.unique(y)
        class_weights = compute_class_weight('balanced', classes=classes, y=y)
        return dict(zip(classes, class_weights))
    
    def _evaluate_model(self, model, X_test: np.ndarray, y_test: np.ndarray, 
                       model_name: str) -> Dict:
        """Comprehensive model evaluation."""
        # Predictions
        y_pred = model.predict(X_test)
        y_pred_proba = model.predict_proba(X_test)[:, 1]
        
        # Metrics
        metrics = {
            'accuracy': model.score(X_test, y_test),
            'balanced_accuracy': balanced_accuracy_score(y_test, y_pred),
            'f1_score': f1_score(y_test, y_pred),
            'roc_auc': roc_auc_score(y_test, y_pred_proba),
            'average_precision': average_precision_score(y_test, y_pred_proba),
            'matthews_corrcoef': matthews_corrcoef(y_test, y_pred)
        }
        
        # Classification report
        report = classification_report(y_test, y_pred, output_dict=True)
        metrics['classification_report'] = report
        
        # Confusion matrix
        cm = confusion_matrix(y_test, y_pred)
        metrics['confusion_matrix'] = cm.tolist()
        
        # Feature importance (if available)
        if hasattr(model, 'feature_importances_'):
            feature_importance = model.feature_importances_
            metrics['feature_importance'] = feature_importance.tolist()
        elif hasattr(model, 'coef_'):
            feature_importance = np.abs(model.coef_[0])
            metrics['feature_importance'] = feature_importance.tolist()
        
        logger.info(f"{model_name} - ROC AUC: {metrics['roc_auc']:.3f}, "
                   f"Balanced Accuracy: {metrics['balanced_accuracy']:.3f}, "
                   f"F1: {metrics['f1_score']:.3f}")
        
        return metrics
    
    def train_models(self, X_train: np.ndarray, X_test: np.ndarray, 
                    y_train: np.ndarray, y_test: np.ndarray,
                    sampling_methods: List[str] = ['none', 'smote', 'class_weight'],
                    cv_folds: int = 5) -> Dict:
        """Train multiple models with different class imbalance strategies."""
        
        results = {}
        
        for sampling_method in sampling_methods:
            logger.info(f"\n{'='*50}")
            logger.info(f"Training with sampling method: {sampling_method}")
            logger.info(f"{'='*50}")
            
            # Prepare training data
            if sampling_method == 'class_weight':
                X_train_balanced = X_train.copy()
                y_train_balanced = y_train.copy()
                class_weights = self._calculate_class_weights(y_train)
                models = self._get_models(class_weight=class_weights)
            elif sampling_method == 'none':
                X_train_balanced = X_train.copy()
                y_train_balanced = y_train.copy()
                models = self._get_models()
            else:
                X_train_balanced, y_train_balanced = self._handle_class_imbalance(
                    X_train, y_train, method=sampling_method
                )
                models = self._get_models()
            
            # Train each model
            method_results = {}
            for model_name, model in models.items():
                logger.info(f"Training {model_name}...")
                
                try:
                    # Cross-validation
                    cv_scores = cross_val_score(
                        model, X_train_balanced, y_train_balanced,
                        cv=StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=self.random_state),
                        scoring='roc_auc'
                    )
                    
                    # Train final model
                    model.fit(X_train_balanced, y_train_balanced)
                    
                    # Evaluate
                    metrics = self._evaluate_model(model, X_test, y_test, model_name)
                    metrics['cv_scores'] = cv_scores.tolist()
                    metrics['cv_mean'] = cv_scores.mean()
                    metrics['cv_std'] = cv_scores.std()
                    
                    # Store model and results
                    full_model_name = f"{model_name}_{sampling_method}"
                    self.models[full_model_name] = model
                    method_results[model_name] = metrics
                    
                    logger.info(f"  CV ROC AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
                    
                except Exception as e:
                    logger.error(f"Training failed for {model_name}: {e}")
                    continue
            
            results[sampling_method] = method_results
        
        # Find best model
        best_score = 0
        for sampling_method, method_results in results.items():
            for model_name, metrics in method_results.items():
                score = metrics['roc_auc']
                if score > best_score:
                    best_score = score
                    self.best_model_name = f"{model_name}_{sampling_method}"
                    self.best_model = self.models[self.best_model_name]
        
        logger.info(f"\nBest model: {self.best_model_name} (ROC AUC: {best_score:.3f})")
        
        self.training_history = results
        return results
    
    def plot_results(self, save_dir: str = "ml_models"):
        """Plot training results and model comparisons."""
        if not self.training_history:
            logger.warning("No training history found. Run train_models first.")
            return
        
        os.makedirs(save_dir, exist_ok=True)
        
        # ROC AUC comparison
        plt.figure(figsize=(12, 8))
        
        methods = []
        models = []
        scores = []
        
        for sampling_method, method_results in self.training_history.items():
            for model_name, metrics in method_results.items():
                methods.append(sampling_method)
                models.append(model_name)
                scores.append(metrics['roc_auc'])
        
        # Create comparison plot
        comparison_df = pd.DataFrame({
            'Sampling Method': methods,
            'Model': models,
            'ROC AUC': scores
        })
        
        plt.subplot(2, 2, 1)
        sns.barplot(data=comparison_df, x='Model', y='ROC AUC', hue='Sampling Method')
        plt.title('ROC AUC Comparison')
        plt.xticks(rotation=45)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Feature importance for best model
        if self.best_model and hasattr(self.best_model, 'feature_importances_'):
            plt.subplot(2, 2, 2)
            importances = self.best_model.feature_importances_
            feature_names = self.feature_extractor.feature_columns[:len(importances)]
            
            # Top 20 features
            top_indices = np.argsort(importances)[-20:]
            top_importances = importances[top_indices]
            top_names = [feature_names[i] for i in top_indices]
            
            plt.barh(range(len(top_importances)), top_importances)
            plt.yticks(range(len(top_importances)), top_names)
            plt.title(f'Top 20 Features - {self.best_model_name}')
            plt.xlabel('Importance')
        
        # Class distribution
        plt.subplot(2, 2, 3)
        # This would show the class distribution - simplified for now
        plt.bar(['Benign', 'Pathogenic'], [147704, 52296])
        plt.title('Class Distribution')
        plt.ylabel('Count')
        
        # Cross-validation scores
        plt.subplot(2, 2, 4)
        if self.best_model_name in self.training_history:
            # Get CV scores from best model
            sampling_method = self.best_model_name.split('_')[-1]
            model_name = '_'.join(self.best_model_name.split('_')[:-1])
            cv_scores = self.training_history[sampling_method][model_name]['cv_scores']
            
            plt.boxplot(cv_scores)
            plt.title(f'CV Scores - {self.best_model_name}')
            plt.ylabel('ROC AUC')
        
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, 'training_results.png'), dpi=300, bbox_inches='tight')
        plt.show()
    
    def save_models(self, save_dir: str = "ml_models"):
        """Save trained models and metadata."""
        os.makedirs(save_dir, exist_ok=True)
        
        # Save best model
        if self.best_model:
            with open(os.path.join(save_dir, 'best_model.pkl'), 'wb') as f:
                pickle.dump(self.best_model, f)
            
            # Save best model metadata
            metadata = {
                'best_model_name': self.best_model_name,
                'training_date': datetime.now().isoformat(),
                'training_history': self.training_history
            }
            
            with open(os.path.join(save_dir, 'model_metadata.json'), 'w') as f:
                json.dump(metadata, f, indent=2)
        
        # Save all models
        for model_name, model in self.models.items():
            with open(os.path.join(save_dir, f'{model_name}.pkl'), 'wb') as f:
                pickle.dump(model, f)
        
        # Save feature extractor
        self.feature_extractor.save_preprocessors(save_dir)
        
        logger.info(f"Models saved to {save_dir}")
    
    def load_model(self, model_path: str = "ml_models/best_model.pkl"):
        """Load a saved model."""
        with open(model_path, 'rb') as f:
            self.best_model = pickle.load(f)
        
        # Load metadata if available
        metadata_path = os.path.join(os.path.dirname(model_path), 'model_metadata.json')
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
                self.best_model_name = metadata.get('best_model_name')
                self.training_history = metadata.get('training_history', {})
        
        # Load feature extractor
        self.feature_extractor.load_preprocessors(os.path.dirname(model_path))
        
        logger.info(f"Model loaded: {self.best_model_name}")
    
    def predict_risk(self, variants: List[Dict]) -> Dict:
        """Predict pathogenicity risk for new variants."""
        if not self.best_model:
            raise ValueError("No trained model available. Train a model first.")
        
        # Extract features
        feature_df = self.feature_extractor.extract_features(variants)
        
        # Handle categorical encoding
        for col in ['gene']:
            if col in feature_df.columns and col in self.feature_extractor.label_encoders:
                le = self.feature_extractor.label_encoders[col]
                # Handle unseen categories
                feature_df[col] = feature_df[col].astype(str)
                feature_df[col] = feature_df[col].apply(
                    lambda x: le.transform([x])[0] if x in le.classes_ else -1
                )
        
        # Scale features
        if self.feature_extractor.scaler:
            X_scaled = self.feature_extractor.scaler.transform(feature_df)
        else:
            X_scaled = feature_df.values
        
        # Predict
        predictions = self.best_model.predict(X_scaled)
        probabilities = self.best_model.predict_proba(X_scaled)
        
        # Format results
        results = {
            'predictions': predictions.tolist(),
            'probabilities': probabilities.tolist(),
            'risk_scores': probabilities[:, 1].tolist(),  # Probability of pathogenic
            'model_used': self.best_model_name
        }
        
        return results


def main():
    """Train and evaluate genomic risk prediction models."""
    logging.basicConfig(level=logging.INFO)
    
    # Initialize trainer
    trainer = GenomicMLTrainer()
    
    # Load and prepare data
    logger.info("Loading training data...")
    X, y = trainer.feature_extractor.load_training_data()
    
    # Prepare for training
    X_train, X_test, y_train, y_test = trainer.feature_extractor.prepare_for_training(X, y)
    
    # Train models
    logger.info("Training models...")
    results = trainer.train_models(
        X_train, X_test, y_train, y_test,
        sampling_methods=['none', 'smote', 'class_weight'],
        cv_folds=5
    )
    
    # Plot results
    trainer.plot_results()
    
    # Save models
    trainer.save_models()
    
    logger.info("Training completed successfully!")


if __name__ == "__main__":
    main() 