#!/usr/bin/env python3
"""
ML Fusion Layer for GeneKnow Risk Assessment

This module implements the fusion layer that combines outputs from the 5 static models:
1. PRS (Polygenic Risk Score)
2. ClinVar (Pathogenic/Benign classifications)
3. CADD (Deleteriousness scores)
4. TCGA (Tumor enrichment matching)
5. Gene/Pathway Burden (Damaging variants in key genes)

The fusion layer learns optimal weights to combine these pre-computed scores
into a final cancer risk assessment.
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import accuracy_score, roc_auc_score, mean_squared_error
import pickle
import json
from typing import Dict, List, Tuple, Any, Optional
from dataclasses import dataclass
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class StaticModelInputs:
    """Data class for static model inputs to the fusion layer."""
    prs_score: float
    clinvar_classification: str  # 'pathogenic', 'benign', 'uncertain', 'not_found'
    cadd_score: float
    tcga_enrichment: float
    gene_burden_score: float
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for processing."""
        return {
            'prs_score': self.prs_score,
            'clinvar_classification': self.clinvar_classification,
            'cadd_score': self.cadd_score,
            'tcga_enrichment': self.tcga_enrichment,
            'gene_burden_score': self.gene_burden_score
        }

@dataclass
class FusionOutput:
    """Data class for fusion layer output."""
    risk_score: float
    confidence: float
    contributing_factors: Dict[str, float]
    risk_category: str  # 'low', 'moderate', 'high', 'very_high'
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'risk_score': self.risk_score,
            'confidence': self.confidence,
            'contributing_factors': self.contributing_factors,
            'risk_category': self.risk_category
        }

class FusionLayer:
    """
    ML Fusion Layer that combines static model outputs into final risk assessment.
    
    This is the meta-learning component that learns optimal weights for combining:
    - PRS (inherited background risk)
    - ClinVar (known pathogenic variants)
    - CADD (functional impact predictions)
    - TCGA (tumor enrichment evidence)
    - Gene burden (pathway-level risk)
    """
    
    def __init__(self, model_type: str = 'gradient_boosting'):
        """
        Initialize the fusion layer.
        
        Args:
            model_type: Type of ML model ('gradient_boosting', 'random_forest', 'linear')
        """
        self.model_type = model_type
        self.model = None
        self.scaler = StandardScaler()
        self.label_encoder = LabelEncoder()
        self.feature_names = ['prs_score', 'cadd_score', 'tcga_enrichment', 'gene_burden_score']
        self.clinvar_categories = ['pathogenic', 'benign', 'uncertain', 'not_found']
        self.is_trained = False
        
        # Risk thresholds for categorization
        self.risk_thresholds = {
            'low': 0.25,
            'moderate': 0.5,
            'high': 0.75,
            'very_high': 1.0
        }
        
        # Initialize model based on type
        if model_type == 'gradient_boosting':
            self.model = GradientBoostingRegressor(
                n_estimators=100,
                learning_rate=0.1,
                max_depth=3,
                random_state=42
            )
        elif model_type == 'random_forest':
            self.model = RandomForestRegressor(
                n_estimators=100,
                max_depth=5,
                random_state=42
            )
        elif model_type == 'linear':
            # Use LinearRegression for continuous targets
            self.model = LinearRegression()
        else:
            raise ValueError(f"Unknown model type: {model_type}")
    
    def _encode_features(self, inputs: List[StaticModelInputs]) -> np.ndarray:
        """
        Encode static model inputs into feature vectors.
        
        Args:
            inputs: List of StaticModelInputs
            
        Returns:
            Encoded feature matrix
        """
        # Convert to DataFrame for easier processing
        df = pd.DataFrame([inp.to_dict() for inp in inputs])
        
        # One-hot encode ClinVar classifications
        clinvar_encoded = pd.get_dummies(df['clinvar_classification'], prefix='clinvar')
        
        # Ensure all expected columns are present
        for cat in self.clinvar_categories:
            col_name = f'clinvar_{cat}'
            if col_name not in clinvar_encoded.columns:
                clinvar_encoded[col_name] = 0
        
        # Reorder columns consistently
        clinvar_cols = [f'clinvar_{cat}' for cat in self.clinvar_categories]
        clinvar_encoded = clinvar_encoded[clinvar_cols]
        
        # Combine numerical features with encoded categorical
        numerical_features = df[self.feature_names]
        features = pd.concat([numerical_features, clinvar_encoded], axis=1)
        
        return features.values
    
    def train(self, training_data: List[Tuple[StaticModelInputs, float]], 
              validation_split: float = 0.2) -> Dict[str, Any]:
        """
        Train the fusion layer on static model outputs and risk labels.
        
        Args:
            training_data: List of (StaticModelInputs, risk_score) pairs
            validation_split: Fraction of data to use for validation
            
        Returns:
            Training metrics and results
        """
        logger.info(f"Training fusion layer with {len(training_data)} samples")
        
        # Extract inputs and targets
        inputs = [item[0] for item in training_data]
        targets = np.array([item[1] for item in training_data])
        
        # Encode features
        X = self._encode_features(inputs)
        
        # Scale features
        X_scaled = self.scaler.fit_transform(X)
        
        # Split data only if validation_split > 0
        if validation_split > 0:
            X_train, X_val, y_train, y_val = train_test_split(
                X_scaled, targets, test_size=validation_split, random_state=42
            )
        else:
            # Use all data for training when validation_split is 0
            X_train, X_val, y_train, y_val = X_scaled, X_scaled, targets, targets
        
        # Train model
        self.model.fit(X_train, y_train)
        self.is_trained = True
        
        # Evaluate
        train_pred = self.model.predict(X_train)
        val_pred = self.model.predict(X_val)
        
        train_mse = mean_squared_error(y_train, train_pred)
        val_mse = mean_squared_error(y_val, val_pred)
        
        # Cross-validation
        cv_scores = cross_val_score(self.model, X_scaled, targets, cv=5, scoring='neg_mean_squared_error')
        
        # Feature importance (if available)
        feature_importance = None
        if hasattr(self.model, 'feature_importances_'):
            feature_names = self.feature_names + [f'clinvar_{cat}' for cat in self.clinvar_categories]
            feature_importance = dict(zip(feature_names, self.model.feature_importances_))
        
        results = {
            'train_mse': train_mse,
            'val_mse': val_mse,
            'cv_mse_mean': -cv_scores.mean(),
            'cv_mse_std': cv_scores.std(),
            'feature_importance': feature_importance,
            'n_samples': len(training_data),
            'model_type': self.model_type
        }
        
        logger.info(f"Training completed. Validation MSE: {val_mse:.4f}")
        return results
    
    def predict(self, inputs: StaticModelInputs) -> FusionOutput:
        """
        Make risk prediction using the fusion layer.
        
        Args:
            inputs: Static model inputs
            
        Returns:
            Fusion output with risk score and metadata
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before making predictions")
        
        # Encode features
        X = self._encode_features([inputs])
        X_scaled = self.scaler.transform(X)
        
        # Make prediction
        risk_score = self.model.predict(X_scaled)[0]
        
        # Clip to valid range
        risk_score = np.clip(risk_score, 0.0, 1.0)
        
        # Calculate confidence (simplified - could be improved)
        confidence = 0.8  # Placeholder - could use prediction interval or ensemble variance
        
        # Determine risk category
        risk_category = 'low'
        for category, threshold in self.risk_thresholds.items():
            if risk_score <= threshold:
                risk_category = category
                break
        
        # Calculate contributing factors (feature importance weighted by input values)
        contributing_factors = {}
        if hasattr(self.model, 'feature_importances_'):
            feature_names = self.feature_names + [f'clinvar_{cat}' for cat in self.clinvar_categories]
            feature_values = X_scaled[0]
            
            for i, (name, importance) in enumerate(zip(feature_names, self.model.feature_importances_)):
                contributing_factors[name] = float(importance * abs(feature_values[i]))
        
        return FusionOutput(
            risk_score=float(risk_score),
            confidence=float(confidence),
            contributing_factors=contributing_factors,
            risk_category=risk_category
        )
    
    def save_model(self, filepath: str) -> None:
        """Save the trained fusion layer model."""
        if not self.is_trained:
            raise ValueError("Cannot save untrained model")
        
        model_data = {
            'model': self.model,
            'scaler': self.scaler,
            'model_type': self.model_type,
            'feature_names': self.feature_names,
            'clinvar_categories': self.clinvar_categories,
            'risk_thresholds': self.risk_thresholds,
            'is_trained': self.is_trained
        }
        
        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)
        
        logger.info(f"Model saved to {filepath}")
    
    def load_model(self, filepath: str) -> None:
        """Load a trained fusion layer model."""
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)
        
        self.model = model_data['model']
        self.scaler = model_data['scaler']
        self.model_type = model_data['model_type']
        self.feature_names = model_data['feature_names']
        self.clinvar_categories = model_data['clinvar_categories']
        self.risk_thresholds = model_data['risk_thresholds']
        self.is_trained = model_data['is_trained']
        
        logger.info(f"Model loaded from {filepath}")

def create_synthetic_training_data(n_samples: int = 1000) -> List[Tuple[StaticModelInputs, float]]:
    """
    Create synthetic training data for testing the fusion layer.
    
    This simulates the outputs that would come from the 5 static models.
    In production, this would be replaced with real data from your pipeline.
    """
    np.random.seed(42)
    training_data = []
    
    for _ in range(n_samples):
        # Generate synthetic static model outputs
        prs_score = np.random.beta(2, 5)  # Most people have low PRS
        
        # ClinVar classification (weighted toward benign/uncertain)
        clinvar_weights = [0.05, 0.7, 0.2, 0.05]  # pathogenic, benign, uncertain, not_found
        clinvar_classification = np.random.choice(
            ['pathogenic', 'benign', 'uncertain', 'not_found'],
            p=clinvar_weights
        )
        
        # CADD score (0-50, higher = more deleterious)
        cadd_score = np.random.exponential(5)
        cadd_score = np.clip(cadd_score, 0, 50)
        
        # TCGA enrichment (fold-change in tumor frequency)
        tcga_enrichment = np.random.lognormal(0, 1)
        tcga_enrichment = np.clip(tcga_enrichment, 0.1, 20)
        
        # Gene burden score (number of damaging variants in key genes)
        gene_burden_score = np.random.poisson(1)
        gene_burden_score = np.clip(gene_burden_score, 0, 10)
        
        # Create synthetic risk score based on logical combination
        risk_score = 0.1  # Base risk
        
        # PRS contribution
        risk_score += prs_score * 0.3
        
        # ClinVar contribution
        if clinvar_classification == 'pathogenic':
            risk_score += 0.5
        elif clinvar_classification == 'benign':
            risk_score -= 0.2
        
        # CADD contribution
        if cadd_score > 20:
            risk_score += 0.3
        elif cadd_score > 15:
            risk_score += 0.1
        
        # TCGA contribution
        if tcga_enrichment > 5:
            risk_score += 0.2
        elif tcga_enrichment > 2:
            risk_score += 0.1
        
        # Gene burden contribution
        risk_score += gene_burden_score * 0.05
        
        # Add some noise and clip to valid range
        risk_score += np.random.normal(0, 0.1)
        risk_score = np.clip(risk_score, 0.0, 1.0)
        
        # Create input object
        inputs = StaticModelInputs(
            prs_score=float(prs_score),
            clinvar_classification=clinvar_classification,
            cadd_score=float(cadd_score),
            tcga_enrichment=float(tcga_enrichment),
            gene_burden_score=float(gene_burden_score)
        )
        
        training_data.append((inputs, risk_score))
    
    return training_data

if __name__ == "__main__":
    # Example usage and testing
    print("🔬 Testing Fusion Layer Implementation")
    
    # Create synthetic training data
    print("📊 Generating synthetic training data...")
    training_data = create_synthetic_training_data(n_samples=2000)
    
    # Initialize and train fusion layer
    print("🎯 Training fusion layer...")
    fusion = FusionLayer(model_type='gradient_boosting')
    results = fusion.train(training_data)
    
    # Print results
    print(f"✅ Training completed:")
    print(f"  Validation MSE: {results['val_mse']:.4f}")
    print(f"  CV MSE: {results['cv_mse_mean']:.4f} ± {results['cv_mse_std']:.4f}")
    
    if results['feature_importance']:
        print("📈 Feature Importance:")
        for feature, importance in sorted(results['feature_importance'].items(), 
                                        key=lambda x: x[1], reverse=True):
            print(f"  {feature}: {importance:.3f}")
    
    # Test prediction
    print("\n🧪 Testing prediction...")
    test_input = StaticModelInputs(
        prs_score=0.8,
        clinvar_classification='pathogenic',
        cadd_score=25.0,
        tcga_enrichment=3.0,
        gene_burden_score=2.0
    )
    
    prediction = fusion.predict(test_input)
    print(f"Risk Score: {prediction.risk_score:.3f}")
    print(f"Risk Category: {prediction.risk_category}")
    print(f"Confidence: {prediction.confidence:.3f}")
    
    # Save model
    fusion.save_model('fusion_model.pkl')
    print("💾 Model saved successfully!") 