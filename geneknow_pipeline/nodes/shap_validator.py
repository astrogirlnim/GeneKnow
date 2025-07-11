"""
SHAP Validator Node for LangGraph genomic pipeline.
Performs automated validation of ML predictions using SHAP explanations.
"""
import logging
import numpy as np
from typing import Dict, Any, List, Tuple
import shap
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LinearRegression

# Configure logging
logger = logging.getLogger(__name__)

# Feature display names for user-friendly output
FEATURE_DISPLAY_NAMES = {
    'prs_score': 'Polygenic Risk Score',
    'cadd_score': 'CADD Deleteriousness Score',
    'tcga_enrichment': 'TCGA Tumor Enrichment',
    'gene_burden_score': 'Gene/Pathway Burden Score',
    'clinvar_pathogenic': 'ClinVar Pathogenic Variant',
    'clinvar_benign': 'ClinVar Benign Variant',
    'clinvar_uncertain': 'ClinVar Uncertain Significance',
    'clinvar_not_found': 'No ClinVar Evidence'
}

# Validation thresholds
HIGH_RISK_THRESHOLD = 0.6
LOW_RISK_THRESHOLD = 0.1
MIN_PATHOGENIC_CONTRIBUTION = 0.15  # Minimum SHAP contribution expected from pathogenic variants


class SHAPValidator:
    """
    Validates ML predictions using SHAP explanations and predefined sanity rules.
    """
    
    def __init__(self):
        self.explainer = None
        self.validation_rules = {
            'high_risk_sanity': self._validate_high_risk,
            'low_risk_sanity': self._validate_low_risk,
            'consistency_check': self._validate_consistency
        }
    
    def _create_explainer(self, model: Any, feature_matrix: np.ndarray) -> shap.Explainer:
        """
        Create appropriate SHAP explainer based on model type.
        
        Args:
            model: Trained ML model
            feature_matrix: Feature matrix used for training
            
        Returns:
            SHAP explainer instance
        """
        if isinstance(model, (GradientBoostingRegressor, RandomForestRegressor)):
            # Use TreeExplainer for tree-based models (fastest)
            return shap.TreeExplainer(model)
        elif isinstance(model, LinearRegression):
            # Use LinearExplainer for linear models
            return shap.LinearExplainer(model, feature_matrix)
        else:
            # Fallback to KernelExplainer (slower but universal)
            logger.warning(f"Using KernelExplainer for model type: {type(model).__name__}")
            return shap.KernelExplainer(model.predict, feature_matrix)
    
    def _get_top_contributors(self, 
                            shap_values: np.ndarray, 
                            feature_names: List[str], 
                            n_top: int = 3) -> List[Dict[str, Any]]:
        """
        Get top contributing features based on SHAP values.
        
        Args:
            shap_values: SHAP values for all features
            feature_names: List of feature names
            n_top: Number of top contributors to return
            
        Returns:
            List of top contributors with their contributions
        """
        # Get absolute SHAP values for ranking
        abs_shap_values = np.abs(shap_values)
        
        # Get indices of top contributors
        top_indices = np.argsort(abs_shap_values)[-n_top:][::-1]
        
        contributors = []
        for idx in top_indices:
            feature_name = feature_names[idx]
            display_name = FEATURE_DISPLAY_NAMES.get(feature_name, feature_name)
            
            contributors.append({
                'feature': feature_name,
                'display_name': display_name,
                'shap_value': float(shap_values[idx]),
                'abs_contribution': float(abs_shap_values[idx]),
                'direction': 'increases' if shap_values[idx] > 0 else 'decreases'
            })
        
        return contributors
    
    def _validate_high_risk(self, 
                          risk_score: float, 
                          shap_values: np.ndarray,
                          feature_names: List[str],
                          variant_data: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """
        Validate high risk predictions.
        
        Rule: High risk (>0.6) should have pathogenic variant in top 3 contributors.
        
        Returns:
            (passed, reasons) tuple
        """
        if risk_score <= HIGH_RISK_THRESHOLD:
            return True, []
        
        reasons = []
        top_contributors = self._get_top_contributors(shap_values, feature_names, n_top=3)
        
        # Check if any pathogenic feature is in top contributors
        pathogenic_features = ['clinvar_pathogenic']
        has_pathogenic_contributor = any(
            contrib['feature'] in pathogenic_features and contrib['shap_value'] > MIN_PATHOGENIC_CONTRIBUTION
            for contrib in top_contributors
        )
        
        if not has_pathogenic_contributor:
            # Check if pathogenic variant exists in data
            has_pathogenic_variant = any(
                (v.get('clinvar_clinical_significance') or '').lower() in ['pathogenic', 'likely pathogenic']
                for v in variant_data.get('filtered_variants', [])
            )
            
            if has_pathogenic_variant:
                reasons.append(
                    f"The AI predicted HIGH RISK ({risk_score*100:.0f}%) but this appears to be based on "
                    f"indirect factors ({', '.join(c['display_name'] for c in top_contributors[:2])}) "
                    f"rather than known disease-causing mutations. The prediction may be less reliable."
                )
            else:
                # High risk without pathogenic variants might be valid (polygenic risk)
                logger.info("High risk score based on polygenic/burden factors - this may be valid")
        
        return len(reasons) == 0, reasons
    
    def _validate_low_risk(self,
                         risk_score: float,
                         shap_values: np.ndarray,
                         feature_names: List[str],
                         variant_data: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """
        Validate low risk predictions.
        
        Rule: Low risk (<0.1) should flag if pathogenic variants present.
        
        Returns:
            (passed, reasons) tuple
        """
        if risk_score >= LOW_RISK_THRESHOLD:
            return True, []
        
        reasons = []
        
        # Check for pathogenic variants in the data
        pathogenic_variants = [
            v for v in variant_data.get('filtered_variants', [])
            if (v.get('clinvar_clinical_significance') or '').lower() in ['pathogenic', 'likely pathogenic']
        ]
        
        if pathogenic_variants:
            # Get SHAP value for pathogenic feature
            pathogenic_idx = feature_names.index('clinvar_pathogenic') if 'clinvar_pathogenic' in feature_names else None
            
            if pathogenic_idx is not None:
                pathogenic_shap = shap_values[pathogenic_idx]
                
                if pathogenic_shap < MIN_PATHOGENIC_CONTRIBUTION:
                    variant_genes = [v.get('gene', 'Unknown') for v in pathogenic_variants[:3]]
                    reasons.append(
                        f"The AI predicted LOW RISK ({risk_score*100:.0f}%) even though "
                        f"{len(pathogenic_variants)} known disease-causing mutation(s) were found in: {', '.join(variant_genes)}. "
                        f"The risk score may be incorrectly low."
                    )
        
        return len(reasons) == 0, reasons
    
    def _validate_consistency(self,
                            risk_score: float,
                            shap_values: np.ndarray,
                            feature_names: List[str],
                            variant_data: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """
        Check for consistency in SHAP explanations.
        
        Returns:
            (passed, reasons) tuple
        """
        reasons = []
        
        # Check for negative SHAP values for pathogenic features
        pathogenic_idx = feature_names.index('clinvar_pathogenic') if 'clinvar_pathogenic' in feature_names else None
        if pathogenic_idx is not None and shap_values[pathogenic_idx] < -0.05:
            reasons.append(
                f"Model logic error: Disease-causing mutations are being counted as protective factors. "
                f"This violates basic genetic principles and suggests the model's reasoning is flawed."
            )
        
        # Check for very high positive SHAP values for benign features
        benign_idx = feature_names.index('clinvar_benign') if 'clinvar_benign' in feature_names else None
        if benign_idx is not None and shap_values[benign_idx] > 0.2:
            reasons.append(
                f"Model logic error: Harmless genetic variants are being counted as major risk factors. "
                f"This suggests the model's internal logic may be compromised."
            )
        
        return len(reasons) == 0, reasons
    
    def validate_prediction(self,
                          model: Any,
                          feature_matrix: np.ndarray,
                          feature_names: List[str],
                          risk_scores: Dict[str, float],
                          variant_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate ML predictions using SHAP explanations.
        
        Args:
            model: Trained ML model
            feature_matrix: Feature matrix for predictions
            feature_names: List of feature names
            risk_scores: Risk scores from the model
            variant_data: Full variant data from pipeline
            
        Returns:
            Validation results with status and explanations
        """
        try:
            # Create SHAP explainer
            explainer = self._create_explainer(model, feature_matrix)
            
            # Calculate SHAP values
            logger.info("Calculating SHAP values...")
            shap_values = explainer(feature_matrix).values
            
            # Handle multi-output models
            if len(shap_values.shape) > 2:
                # Average across outputs for overall importance
                shap_values = np.mean(shap_values, axis=0)
            
            # If we have multiple samples, focus on the first one
            # (In production, we'd validate each sample)
            if len(shap_values.shape) > 1:
                shap_values = shap_values[0]
            
            # Get aggregate risk score
            aggregate_risk = risk_scores.get('aggregate_risk_score', 0.0)
            
            # Run validation rules
            all_passed = True
            all_reasons = []
            rule_results = {}
            
            for rule_name, rule_func in self.validation_rules.items():
                passed, reasons = rule_func(aggregate_risk, shap_values, feature_names, variant_data)
                rule_results[rule_name] = {
                    'passed': passed,
                    'reasons': reasons
                }
                if not passed:
                    all_passed = False
                    all_reasons.extend(reasons)
            
            # Get top contributors
            top_contributors = self._get_top_contributors(shap_values, feature_names, n_top=3)
            
            # Determine validation status
            validation_status = "PASS" if all_passed else "FLAG_FOR_REVIEW"
            
            # Create validation summary
            validation_summary = {
                'status': validation_status,
                'risk_score': aggregate_risk,
                'top_contributors': top_contributors,
                'validation_reasons': all_reasons,
                'rule_results': rule_results,
                'shap_values': shap_values.tolist(),
                'feature_names': feature_names,
                'model_type': type(model).__name__
            }
            
            logger.info(f"SHAP validation completed: {validation_status}")
            if not all_passed:
                logger.warning(f"Validation issues: {'; '.join(all_reasons)}")
            
            return validation_summary
            
        except Exception as e:
            logger.error(f"Error in SHAP validation: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            
            # Return a safe default if validation fails
            return {
                'status': 'ERROR',
                'error': str(e),
                'risk_score': risk_scores.get('aggregate_risk_score', 0.0),
                'top_contributors': [],
                'validation_reasons': [f"SHAP validation error: {str(e)}"]
            }


def process(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    LangGraph node function for SHAP validation.
    
    Args:
        state: Current pipeline state
        
    Returns:
        Updated state with SHAP validation results
    """
    logger.info("🔍 Starting SHAP validation")
    
    # Check if we have required inputs
    ml_fusion_model = state.get('ml_fusion_model_instance')
    ml_fusion_features = state.get('ml_fusion_feature_matrix')
    ml_fusion_results = state.get('ml_fusion_results', {})
    risk_scores = state.get('risk_scores', {})
    
    if not ml_fusion_model or ml_fusion_features is None:
        logger.warning("Missing ML fusion model or features, skipping SHAP validation")
        return {
            'shap_validation_status': 'SKIPPED',
            'shap_validation_reasons': ['ML fusion model or features not available'],
            'shap_top_contributors': [],
            'shap_feature_importance': {},
            'shap_validation_details': {}
        }
    
    # Get the actual model from fusion layer
    if hasattr(ml_fusion_model, 'model'):
        model = ml_fusion_model.model
        feature_names = (ml_fusion_model.feature_names + 
                        [f'clinvar_{cat}' for cat in ml_fusion_model.clinvar_categories])
    else:
        logger.error("ML fusion model instance doesn't have expected structure")
        return {
            'shap_validation_status': 'ERROR',
            'shap_validation_reasons': ['Invalid model structure'],
            'shap_top_contributors': [],
            'shap_feature_importance': {},
            'shap_validation_details': {}
        }
    
    # Run SHAP validation
    validator = SHAPValidator()
    validation_results = validator.validate_prediction(
        model=model,
        feature_matrix=ml_fusion_features,
        feature_names=feature_names,
        risk_scores=risk_scores,
        variant_data=state
    )
    
    # Extract key results for state
    return {
        'shap_validation_status': validation_results['status'],
        'shap_validation_reasons': validation_results.get('validation_reasons', []),
        'shap_top_contributors': validation_results.get('top_contributors', []),
        'shap_feature_importance': {
            name: float(value) 
            for name, value in zip(
                validation_results.get('feature_names', []),
                validation_results.get('shap_values', [])
            )
        },
        'shap_validation_details': validation_results
    } 