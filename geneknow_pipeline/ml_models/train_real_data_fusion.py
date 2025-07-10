#!/usr/bin/env python3
"""
Real Data Fusion Layer Training

This script trains the ML fusion layer using actual variant data from our database
instead of synthetic data. It implements proper data splitting, feature generation
from all 5 static models, and comprehensive evaluation.
"""

import sqlite3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import mean_squared_error, r2_score, roc_auc_score, accuracy_score
from sklearn.metrics import classification_report, confusion_matrix
from fusion_layer import FusionLayer, StaticModelInputs, FusionOutput
import pickle
import json
from datetime import datetime
from typing import Dict, List, Tuple, Any
import logging
import warnings
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RealDataFusionTrainer:
    """
    Fusion layer trainer using real genomic variant data.
    
    FIXED: No data leakage! Uses ClinVar as independent target.
    
    Implements 4 static model features:
    1. PRS (Polygenic Risk Score)
    2. CADD (Deleteriousness Scores)
    3. TCGA (Tumor Enrichment)
    4. Gene Burden (Pathway Analysis)
    
    Target: ClinVar pathogenicity (independent ground truth)
    """
    
    def __init__(self, db_path: str = "population_variants.db"):
        self.db_path = db_path
        self.cancer_gene_groups = self._define_cancer_gene_groups()
        self.feature_columns = []
        
    def _define_cancer_gene_groups(self) -> Dict[str, List[str]]:
        """Define cancer gene groups for gene burden analysis."""
        return {
            'tumor_suppressor': ['TP53', 'RB1', 'APC', 'BRCA1', 'BRCA2', 'PTEN', 'VHL', 'STK11'],
            'oncogenes': ['KRAS', 'BRAF', 'PIK3CA', 'EGFR', 'MYC', 'RAS'],
            'dna_repair': ['BRCA1', 'BRCA2', 'ATM', 'PALB2', 'CHEK2', 'MLH1', 'MSH2', 'MSH6', 'PMS2'],
            'cell_cycle': ['TP53', 'RB1', 'CDKN2A', 'ATM', 'CHEK2'],
            'mismatch_repair': ['MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM'],
            'breast_cancer': ['BRCA1', 'BRCA2', 'PALB2', 'ATM', 'CHEK2', 'TP53', 'PTEN'],
            'colon_cancer': ['APC', 'KRAS', 'TP53', 'PIK3CA', 'SMAD4', 'BRAF', 'MLH1', 'MSH2'],
            'lung_cancer': ['TP53', 'KRAS', 'EGFR', 'STK11', 'KEAP1', 'NF1', 'RB1'],
            'hematologic': ['JAK2', 'FLT3', 'NPM1', 'DNMT3A', 'TET2', 'IDH1', 'IDH2']
        }
    
    def load_database_variants(self) -> pd.DataFrame:
        """Load variants from the database with all available annotations."""
        logger.info("Loading variants from database...")
        
        conn = sqlite3.connect(self.db_path)
        
        # Main query to get population variants with TCGA enrichment
        query = """
        SELECT 
            pv.chrom, pv.pos, pv.ref, pv.alt, pv.gene,
            pv.gnomad_af, pv.clinical_significance, pv.is_pathogenic,
            pv.consequence, pv.review_status,
            COALESCE(MAX(tv.enrichment_score), 1.0) as tcga_enrichment_score,
            COALESCE(MAX(tv.tumor_frequency), 0.0) as tcga_tumor_frequency,
            COALESCE(MAX(tv.sample_count), 0) as tcga_sample_count
        FROM population_variants pv
        LEFT JOIN tcga_variants tv ON 
            pv.chrom = tv.chrom AND pv.pos = tv.pos 
            AND pv.ref = tv.ref AND pv.alt = tv.alt
        WHERE pv.gene IS NOT NULL AND pv.clinical_significance IS NOT NULL
        GROUP BY pv.chrom, pv.pos, pv.ref, pv.alt
        """
        
        df = pd.read_sql_query(query, conn)
        conn.close()
        
        logger.info(f"Loaded {len(df):,} variants from database")
        logger.info(f"Clinical significance distribution:")
        clin_counts = df['clinical_significance'].value_counts()
        for clin, count in clin_counts.head(10).items():
            logger.info(f"  {clin}: {count:,}")
        
        return df
    
    def generate_static_model_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate features from all 5 static models.
        
        Returns DataFrame with static model features for fusion layer.
        """
        logger.info("Generating static model features...")
        
        features_df = df.copy()
        
        # 1. PRS (Polygenic Risk Score) - Simplified calculation
        logger.info("  Generating PRS scores...")
        # Use allele frequency as proxy for population-level risk contribution
        # In reality, this would use GWAS effect sizes
        features_df['prs_score'] = self._calculate_prs_scores(features_df)
        
        # 2. ClinVar Clinical Significance - EXCLUDED FROM FEATURES!
        # We use this as our target, not as an input feature
        logger.info("  Excluding ClinVar from features (using as target)...")
        
        # 3. CADD Scores (simulate based on consequence and gene importance)
        logger.info("  Generating CADD scores...")
        features_df['cadd_score'] = self._generate_cadd_scores(features_df)
        
        # 4. TCGA Enrichment (already joined from database)
        logger.info("  Processing TCGA enrichment...")
        features_df['tcga_enrichment'] = features_df['tcga_enrichment_score']
        
        # 5. Gene Burden Scores
        logger.info("  Calculating gene burden scores...")
        features_df['gene_burden_score'] = self._calculate_gene_burden_scores(features_df)
        
        # Create target variable (risk score)
        features_df['risk_score'] = self._calculate_risk_targets(features_df)
        
        logger.info(f"Generated features for {len(features_df):,} variants")
        
        return features_df
    
    def _calculate_prs_scores(self, df: pd.DataFrame) -> np.ndarray:
        """Calculate PRS scores based on allele frequency and gene importance."""
        prs_scores = np.zeros(len(df))
        
        # Base score from allele frequency (rare variants get higher base risk)
        af = df['gnomad_af'].fillna(0.0)
        prs_scores = 0.1 + (1 - af) * 0.4  # Rare variants = higher base risk
        
        # Gene-specific multipliers
        for gene_group, genes in self.cancer_gene_groups.items():
            mask = df['gene'].isin(genes)
            if gene_group in ['tumor_suppressor', 'dna_repair']:
                prs_scores[mask] *= 1.5  # Higher risk for key cancer genes
            elif gene_group in ['oncogenes']:
                prs_scores[mask] *= 1.3
        
        # Add some population-level variation
        np.random.seed(42)
        prs_scores += np.random.beta(2, 5, len(df)) * 0.2
        
        return np.clip(prs_scores, 0.0, 1.0)
    
    # _process_clinvar_classifications method removed - we use ClinVar as target, not feature
    
    def _generate_cadd_scores(self, df: pd.DataFrame) -> np.ndarray:
        """Generate CADD-like scores based on consequence and gene context."""
        cadd_scores = np.zeros(len(df))
        
        # Base scores by consequence type
        consequence_scores = {
            'stop_gained': 35, 'stop_lost': 35, 'start_lost': 30,
            'frameshift': 30, 'missense': 20, 'splice': 25,
            'synonymous': 5, 'utr': 8, 'intronic': 3
        }
        
        for i, consequence in enumerate(df['consequence']):
            base_score = 10  # Default
            
            if pd.isna(consequence):
                consequence = ''
            consequence_lower = str(consequence).lower()
            
            for cons_type, score in consequence_scores.items():
                if cons_type in consequence_lower:
                    base_score = score
                    break
            
            # Gene-specific adjustments
            gene = df.iloc[i]['gene']
            if gene in ['TP53', 'BRCA1', 'BRCA2']:
                base_score *= 1.3  # Higher scores for critical genes
            elif gene in self.cancer_gene_groups['tumor_suppressor']:
                base_score *= 1.2
            
            # Allele frequency adjustment
            af = df.iloc[i]['gnomad_af']
            if pd.notna(af) and af < 0.001:
                base_score += 3  # Rare variants get bonus
                
            cadd_scores[i] = base_score
        
        # Add noise and clip
        np.random.seed(42)
        cadd_scores += np.random.normal(0, 2, len(df))
        
        return np.clip(cadd_scores, 0, 50)
    
    def _calculate_gene_burden_scores(self, df: pd.DataFrame) -> np.ndarray:
        """Calculate gene burden scores based on damaging variants in pathways."""
        burden_scores = np.zeros(len(df))
        
        # Count variants per gene (higher burden = more damaging variants)
        gene_counts = df['gene'].value_counts()
        
        for i, gene in enumerate(df['gene']):
            score = 0
            
            # Base score from gene variant count (normalized)
            if gene in gene_counts:
                score += min(gene_counts[gene] / 1000, 5)  # Max 5 points
            
            # Pathway membership bonuses
            for group_name, genes in self.cancer_gene_groups.items():
                if gene in genes:
                    if group_name == 'dna_repair':
                        score += 3
                    elif group_name == 'tumor_suppressor':
                        score += 2
                    else:
                        score += 1
            
            burden_scores[i] = score
        
        return np.clip(burden_scores, 0, 10)
    
    def _calculate_risk_targets(self, df: pd.DataFrame) -> np.ndarray:
        """
        Calculate training targets using INDEPENDENT ground truth from ClinVar.
        
        FIXED: No more data leakage! We use ClinVar clinical significance as 
        the independent target and EXCLUDE it from model features.
        
        Target: Binary pathogenicity from ClinVar expert curation
        - Pathogenic/Likely_pathogenic: 1 (high risk)
        - Benign/Likely_benign: 0 (low risk)  
        - Uncertain: 0.5 (moderate risk)
        - Others: 0.1 (low risk)
        """
        risk_scores = np.zeros(len(df))
        
        for i in range(len(df)):
            clin_sig = str(df.iloc[i]['clinical_significance']).lower()
            
            # Use ClinVar as INDEPENDENT ground truth target
            if 'pathogenic' in clin_sig and 'likely' not in clin_sig:
                risk_scores[i] = 1.0  # Definitely pathogenic
            elif 'likely_pathogenic' in clin_sig:
                risk_scores[i] = 0.8  # Likely pathogenic
            elif 'benign' in clin_sig:
                risk_scores[i] = 0.0  # Definitely benign
            elif 'uncertain' in clin_sig:
                risk_scores[i] = 0.5  # Unknown - let model decide
            else:
                risk_scores[i] = 0.1  # Default low risk
        
        return risk_scores
    
    def create_stratified_splits(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Create stratified train/validation/test splits.
        
        Stratifies primarily by clinical significance with backup random sampling.
        """
        logger.info("Creating stratified data splits...")
        
        # Create stratification variables
        df = df.copy()
        
        # Simplified stratification - use original clinical significance
        # Group rare categories together
        def simplify_clinical_significance(clin_sig):
            clin_lower = str(clin_sig).lower()
            if 'pathogenic' in clin_lower:
                return 'pathogenic'
            elif 'benign' in clin_lower:
                return 'benign'
            elif 'uncertain' in clin_lower:
                return 'uncertain'
            else:
                return 'other'
        
        df['stratify_simple'] = df['clinical_significance'].apply(simplify_clinical_significance)
        
        # Check stratification distribution
        strat_counts = df['stratify_simple'].value_counts()
        logger.info("Stratification distribution:")
        for strat, count in strat_counts.items():
            logger.info(f"  {strat}: {count:,}")
        
        # Only stratify if all groups have sufficient samples
        min_samples = strat_counts.min()
        if min_samples >= 10:
            # Use stratified split
            try:
                # First split: 80% train+val, 20% test
                train_val, test = train_test_split(
                    df, 
                    test_size=0.2, 
                    random_state=42,
                    stratify=df['stratify_simple']
                )
                
                # Second split: 60% train, 20% val from the remaining 80%
                train, val = train_test_split(
                    train_val,
                    test_size=0.25,  # 0.25 of 80% = 20% of total
                    random_state=42,
                    stratify=train_val['stratify_simple']
                )
                logger.info("Using stratified splits")
            except ValueError as e:
                logger.warning(f"Stratified split failed: {e}")
                logger.info("Falling back to random splits")
                # Fall back to random splits
                train_val, test = train_test_split(df, test_size=0.2, random_state=42)
                train, val = train_test_split(train_val, test_size=0.25, random_state=42)
        else:
            logger.warning(f"Insufficient samples for stratification (min: {min_samples})")
            logger.info("Using random splits")
            # Use random splits
            train_val, test = train_test_split(df, test_size=0.2, random_state=42)
            train, val = train_test_split(train_val, test_size=0.25, random_state=42)
        
        logger.info(f"Data splits created:")
        logger.info(f"  Training: {len(train):,} variants ({len(train)/len(df)*100:.1f}%)")
        logger.info(f"  Validation: {len(val):,} variants ({len(val)/len(df)*100:.1f}%)")
        logger.info(f"  Test: {len(test):,} variants ({len(test)/len(df)*100:.1f}%)")
        
        return train, val, test
    
    def prepare_fusion_inputs(self, df: pd.DataFrame) -> List[Tuple[StaticModelInputs, float]]:
        """Convert DataFrame to fusion layer input format (NO CLINVAR FEATURES)."""
        fusion_data = []
        
        for _, row in df.iterrows():
            inputs = StaticModelInputs(
                prs_score=float(row['prs_score']),
                clinvar_classification='not_found',  # Exclude ClinVar from features!
                cadd_score=float(row['cadd_score']),
                tcga_enrichment=float(row['tcga_enrichment']),
                gene_burden_score=float(row['gene_burden_score'])
            )
            
            risk_score = float(row['risk_score'])
            fusion_data.append((inputs, risk_score))
        
        return fusion_data
    
    def train_and_evaluate_models(self, train_df: pd.DataFrame, val_df: pd.DataFrame, test_df: pd.DataFrame):
        """Train multiple fusion layer models and evaluate performance."""
        logger.info("Training fusion layer models...")
        
        # Prepare data
        train_data = self.prepare_fusion_inputs(train_df)
        val_data = self.prepare_fusion_inputs(val_df)
        test_data = self.prepare_fusion_inputs(test_df)
        
        # Model types to train
        model_types = ['gradient_boosting', 'random_forest', 'linear']
        results = {}
        
        # Train each model type
        for model_type in model_types:
            logger.info(f"\n🎯 Training {model_type} model...")
            
            # Initialize and train
            fusion = FusionLayer(model_type=model_type)
            training_results = fusion.train(train_data, validation_split=0.0)  # We have separate val set
            
            # Evaluate on validation set
            val_predictions = []
            val_targets = []
            
            for inputs, target in val_data:
                pred = fusion.predict(inputs)
                val_predictions.append(pred.risk_score)
                val_targets.append(target)
            
            val_predictions = np.array(val_predictions)
            val_targets = np.array(val_targets)
            
            # Calculate metrics
            val_mse = mean_squared_error(val_targets, val_predictions)
            val_r2 = r2_score(val_targets, val_predictions)
            
            # Binary classification metrics
            val_binary_targets = (val_targets > 0.5).astype(int)
            val_binary_preds = (val_predictions > 0.5).astype(int)
            val_accuracy = accuracy_score(val_binary_targets, val_binary_preds)
            val_auc = roc_auc_score(val_binary_targets, val_predictions)
            
            # Store results
            results[model_type] = {
                'model': fusion,
                'training_results': training_results,
                'val_mse': val_mse,
                'val_r2': val_r2,
                'val_accuracy': val_accuracy,
                'val_auc': val_auc,
                'val_predictions': val_predictions,
                'val_targets': val_targets
            }
            
            logger.info(f"  Validation MSE: {val_mse:.4f}")
            logger.info(f"  Validation R²: {val_r2:.4f}")
            logger.info(f"  Validation AUC: {val_auc:.4f}")
            
            # Save model
            model_path = f'fusion_{model_type}_FIXED.pkl'
            fusion.save_model(model_path)
            logger.info(f"  Model saved to {model_path}")
        
        # Select best model based on validation AUC
        best_model_type = max(results.keys(), key=lambda k: results[k]['val_auc'])
        best_model = results[best_model_type]['model']
        
        logger.info(f"\n🏆 Best model: {best_model_type} (AUC: {results[best_model_type]['val_auc']:.4f})")
        
        # Final evaluation on test set (ONLY ONCE!)
        logger.info("\n🧪 Final evaluation on held-out test set...")
        test_predictions = []
        test_targets = []
        
        for inputs, target in test_data:
            pred = best_model.predict(inputs)
            test_predictions.append(pred.risk_score)
            test_targets.append(target)
        
        test_predictions = np.array(test_predictions)
        test_targets = np.array(test_targets)
        
        # Calculate final metrics
        test_mse = mean_squared_error(test_targets, test_predictions)
        test_r2 = r2_score(test_targets, test_predictions)
        
        test_binary_targets = (test_targets > 0.5).astype(int)
        test_binary_preds = (test_predictions > 0.5).astype(int)
        test_accuracy = accuracy_score(test_binary_targets, test_binary_preds)
        test_auc = roc_auc_score(test_binary_targets, test_predictions)
        
        # Save best model as default
        best_model.save_model('best_fusion_model_FIXED.pkl')
        
        # Create comprehensive results summary
        final_results = {
            'best_model_type': best_model_type,
            'training_date': datetime.now().isoformat(),
            'data_stats': {
                'total_variants': len(train_df) + len(val_df) + len(test_df),
                'train_size': len(train_df),
                'val_size': len(val_df),
                'test_size': len(test_df)
            },
            'all_models': {k: {
                'val_mse': v['val_mse'],
                'val_r2': v['val_r2'],
                'val_accuracy': v['val_accuracy'],
                'val_auc': v['val_auc']
            } for k, v in results.items()},
            'final_test_results': {
                'mse': test_mse,
                'r2': test_r2,
                'accuracy': test_accuracy,
                'auc': test_auc
            }
        }
        
        # Save results
        with open('real_data_training_results_FIXED.json', 'w') as f:
            json.dump(final_results, f, indent=2)
        
        # Print final summary
        print("\n" + "="*60)
        print("🎉 REAL DATA TRAINING COMPLETE!")
        print("="*60)
        print(f"Best Model: {best_model_type}")
        print(f"Training Data: {len(train_df):,} variants")
        print(f"Validation Data: {len(val_df):,} variants")
        print(f"Test Data: {len(test_df):,} variants")
        print("\nFinal Test Performance:")
        print(f"  MSE: {test_mse:.4f}")
        print(f"  R²: {test_r2:.4f}")
        print(f"  Accuracy: {test_accuracy:.4f}")
        print(f"  AUC: {test_auc:.4f}")
        print(f"\nModel saved as: best_fusion_model_FIXED.pkl")
        print(f"Results saved as: real_data_training_results_FIXED.json")
        
        return results, final_results

def main():
    """Main training function - FIXED VERSION (no data leakage)."""
    print("🧬 GeneKnow Real Data Fusion Layer Training - FIXED!")
    print("🔧 Uses ClinVar as independent target, not as feature")
    print("=" * 60)
    
    # Initialize trainer
    trainer = RealDataFusionTrainer()
    
    # Load real data from database
    df = trainer.load_database_variants()
    
    # Generate static model features
    df_with_features = trainer.generate_static_model_features(df)
    
    # Create stratified splits
    train_df, val_df, test_df = trainer.create_stratified_splits(df_with_features)
    
    # Train and evaluate models
    model_results, final_results = trainer.train_and_evaluate_models(train_df, val_df, test_df)
    
    print(f"\n✅ Training completed successfully!")
    
    return model_results, final_results

if __name__ == "__main__":
    results, final = main() 