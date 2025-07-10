#!/usr/bin/env python3
"""
ML Feature Extractor (No Data Leakage Version)
Excludes clinical significance features to avoid data leakage in training.
"""

import sqlite3
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
import logging
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split
import pickle
import os

logger = logging.getLogger(__name__)

class GenomicFeatureExtractorNoLeakage:
    """
    Extract ML-ready features from the GeneKnow database WITHOUT clinical significance.
    
    Features extracted (NO DATA LEAKAGE):
    1. Population frequency (gnomAD) - objective measurement
    2. Gene-level features (cancer gene groups) - biological knowledge
    3. Variant consequence severity - functional prediction
    4. TCGA cancer enrichment - cancer relevance
    5. CADD scores - deleteriousness prediction
    
    EXCLUDED to avoid leakage:
    - clinical_significance (this is what we're trying to predict!)
    - is_pathogenic (this is our target label)
    """
    
    def __init__(self, db_path: str = "population_variants.db"):
        self.db_path = db_path
        self.feature_columns = []
        self.label_encoders = {}
        self.scaler = None
        self.cancer_gene_groups = self._define_cancer_gene_groups()
        
    def _define_cancer_gene_groups(self) -> Dict[str, List[str]]:
        """Define cancer gene groups for feature engineering."""
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
    
    def extract_features(self, variant_batch: List[Dict]) -> pd.DataFrame:
        """Extract features from a batch of variants (NO clinical significance)."""
        features = []
        
        for variant in variant_batch:
            feature_vector = self._extract_single_variant_features(variant)
            features.append(feature_vector)
        
        return pd.DataFrame(features)
    
    def _extract_single_variant_features(self, variant: Dict) -> Dict:
        """Extract features for a single variant (NO DATA LEAKAGE)."""
        features = {}
        
        # 1. Population frequency features (KEEP - objective measurement)
        gnomad_af = variant.get('gnomad_af', 0.0)
        features['population_frequency'] = gnomad_af
        features['is_common_variant'] = 1 if gnomad_af > 0.01 else 0
        features['is_rare_variant'] = 1 if gnomad_af < 0.001 else 0
        features['log_population_frequency'] = np.log10(gnomad_af + 1e-8)
        
        # 2. Gene-level features (KEEP - biological knowledge)
        gene = variant.get('gene', 'Unknown')
        features['gene'] = gene
        
        # Cancer gene group memberships
        for group_name, genes in self.cancer_gene_groups.items():
            features[f'in_{group_name}'] = 1 if gene in genes else 0
        
        # 3. Variant consequence severity (KEEP - functional prediction)
        consequence = variant.get('consequence', '').lower()
        consequence_severity = self._get_consequence_severity(consequence)
        features['consequence_severity'] = consequence_severity
        features['is_truncating'] = 1 if any(term in consequence for term in ['nonsense', 'frameshift', 'stop_gained', 'stop_lost']) else 0
        features['is_missense'] = 1 if 'missense' in consequence else 0
        features['is_synonymous'] = 1 if 'synonymous' in consequence else 0
        
        # 4. REMOVED - Clinical significance features (CAUSES DATA LEAKAGE)
        # These are derived from the same source as our target labels!
        # ❌ features['clinvar_pathogenic'] = ...
        # ❌ features['clinvar_likely_pathogenic'] = ...  
        # ❌ features['clinvar_benign'] = ...
        # ❌ features['clinvar_likely_benign'] = ...
        # ❌ features['clinvar_uncertain'] = ...
        
        # 5. Review status quality (KEEP - data quality measure)
        review_status = variant.get('review_status', '').lower()
        features['review_status_quality'] = self._get_review_status_score(review_status)
        
        # 6. TCGA enrichment features (KEEP - cancer relevance)
        tcga_enrichment = variant.get('tcga_enrichment_score', 0.0)
        features['tcga_enrichment'] = tcga_enrichment
        features['in_tcga_database'] = 1 if tcga_enrichment > 0 else 0
        
        # 7. CADD scores (KEEP - functional prediction)
        cadd_phred = variant.get('cadd_phred', 0.0)
        features['cadd_phred'] = cadd_phred
        features['cadd_high_impact'] = 1 if cadd_phred > 20 else 0
        features['cadd_moderate_impact'] = 1 if 10 <= cadd_phred <= 20 else 0
        
        return features
    
    def _get_consequence_severity(self, consequence: str) -> float:
        """Convert consequence to numeric severity score."""
        severity_map = {
            'transcript_ablation': 1.0,
            'splice_acceptor_variant': 1.0,
            'splice_donor_variant': 1.0,
            'stop_gained': 1.0,
            'frameshift_variant': 1.0,
            'stop_lost': 0.9,
            'start_lost': 0.9,
            'transcript_amplification': 0.8,
            'inframe_insertion': 0.7,
            'inframe_deletion': 0.7,
            'missense_variant': 0.6,
            'protein_altering_variant': 0.6,
            'splice_region_variant': 0.5,
            'incomplete_terminal_codon_variant': 0.5,
            'stop_retained_variant': 0.4,
            'synonymous_variant': 0.2,
            'coding_sequence_variant': 0.2,
            'mature_mirna_variant': 0.2,
            '5_prime_utr_variant': 0.1,
            '3_prime_utr_variant': 0.1,
            'non_coding_transcript_exon_variant': 0.1,
            'intron_variant': 0.05,
            'nmd_transcript_variant': 0.05,
            'non_coding_transcript_variant': 0.05,
            'upstream_gene_variant': 0.02,
            'downstream_gene_variant': 0.02,
            'tfbs_ablation': 0.8,
            'tfbs_amplification': 0.8,
            'tf_binding_site_variant': 0.3,
            'regulatory_region_ablation': 0.8,
            'regulatory_region_amplification': 0.8,
            'feature_elongation': 0.5,
            'regulatory_region_variant': 0.3,
            'feature_truncation': 0.5,
            'intergenic_variant': 0.01
        }
        
        # Handle multiple consequences
        if '&' in consequence:
            consequences = consequence.split('&')
            return max(severity_map.get(c.strip(), 0.1) for c in consequences)
        
        return severity_map.get(consequence, 0.1)
    
    def _get_review_status_score(self, review_status: str) -> float:
        """Convert review status to quality score."""
        if 'practice_guideline' in review_status:
            return 1.0
        elif 'reviewed_by_expert_panel' in review_status:
            return 0.9
        elif 'criteria_provided' in review_status and 'multiple_submitters' in review_status:
            return 0.8
        elif 'criteria_provided' in review_status and 'single_submitter' in review_status:
            return 0.6
        elif 'no_assertion' in review_status:
            return 0.2
        else:
            return 0.3
    
    def load_training_data(self, cancer_specific: Optional[str] = None) -> Tuple[pd.DataFrame, pd.Series]:
        """Load training data from database (NO CLINICAL SIGNIFICANCE FEATURES)."""
        conn = sqlite3.connect(self.db_path)
        
        # Base query - NO clinical significance in SELECT
        query = """
        SELECT 
            pv.chrom, pv.pos, pv.ref, pv.alt, pv.gene,
            pv.gnomad_af, pv.is_pathogenic,
            pv.consequence, pv.review_status
        FROM population_variants pv
        WHERE pv.gene IS NOT NULL
        """
        
        # Add TCGA enrichment if available
        tcga_query = """
        SELECT 
            pv.chrom, pv.pos, pv.ref, pv.alt, pv.gene,
            pv.gnomad_af, pv.is_pathogenic,
            pv.consequence, pv.review_status,
            COALESCE(MAX(tv.enrichment_score), 0) as tcga_enrichment_score
        FROM population_variants pv
        LEFT JOIN tcga_variants tv ON 
            pv.chrom = tv.chrom AND pv.pos = tv.pos 
            AND pv.ref = tv.ref AND pv.alt = tv.alt
        WHERE pv.gene IS NOT NULL
        GROUP BY pv.chrom, pv.pos, pv.ref, pv.alt
        """
        
        # Check if TCGA table exists
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='tcga_variants'")
        has_tcga = cursor.fetchone() is not None
        
        if has_tcga:
            query = tcga_query
        
        # Add cancer-specific filtering if requested
        if cancer_specific:
            cancer_genes = []
            for group_name, genes in self.cancer_gene_groups.items():
                if cancer_specific.lower() in group_name:
                    cancer_genes.extend(genes)
            
            if cancer_genes:
                gene_list = "', '".join(cancer_genes)
                query += f" AND pv.gene IN ('{gene_list}')"
        
        # Load data
        df = pd.read_sql_query(query, conn)
        conn.close()
        
        logger.info(f"Loaded {len(df)} variants for training (NO clinical significance features)")
        
        # Convert to feature vectors (without clinical significance)
        variants = df.to_dict('records')
        feature_df = self.extract_features(variants)
        
        # Prepare labels
        labels = df['is_pathogenic']
        
        # Store feature columns
        self.feature_columns = feature_df.columns.tolist()
        
        logger.info(f"Feature columns (no leakage): {len(self.feature_columns)} features")
        logger.info(f"Features: {self.feature_columns}")
        
        return feature_df, labels
    
    def prepare_for_training(self, X: pd.DataFrame, y: pd.Series, 
                           test_size: float = 0.2, 
                           random_state: int = 42) -> Tuple:
        """Prepare data for ML training with proper encoding and scaling."""
        
        # Handle categorical variables
        categorical_cols = ['gene']
        for col in categorical_cols:
            if col in X.columns:
                le = LabelEncoder()
                X[col] = le.fit_transform(X[col].astype(str))
                self.label_encoders[col] = le
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_state, 
            stratify=y  # Maintain class balance
        )
        
        # Scale features
        self.scaler = StandardScaler()
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        logger.info(f"Training set: {len(X_train)} samples")
        logger.info(f"Test set: {len(X_test)} samples")
        logger.info(f"Feature dimensions: {X_train_scaled.shape[1]}")
        
        return X_train_scaled, X_test_scaled, y_train, y_test
    
    def save_preprocessors(self, output_dir: str = "ml_models_no_leakage"):
        """Save preprocessing components."""
        os.makedirs(output_dir, exist_ok=True)
        
        # Save scaler
        if self.scaler:
            with open(os.path.join(output_dir, 'scaler.pkl'), 'wb') as f:
                pickle.dump(self.scaler, f)
        
        # Save label encoders
        if self.label_encoders:
            with open(os.path.join(output_dir, 'label_encoders.pkl'), 'wb') as f:
                pickle.dump(self.label_encoders, f)
        
        # Save feature columns
        with open(os.path.join(output_dir, 'feature_columns.pkl'), 'wb') as f:
            pickle.dump(self.feature_columns, f)
    
    def load_preprocessors(self, model_dir: str = "ml_models_no_leakage"):
        """Load preprocessing components."""
        # Load scaler
        scaler_path = os.path.join(model_dir, 'scaler.pkl')
        if os.path.exists(scaler_path):
            with open(scaler_path, 'rb') as f:
                self.scaler = pickle.load(f)
        
        # Load label encoders
        encoders_path = os.path.join(model_dir, 'label_encoders.pkl')
        if os.path.exists(encoders_path):
            with open(encoders_path, 'rb') as f:
                self.label_encoders = pickle.load(f)
        
        # Load feature columns
        features_path = os.path.join(model_dir, 'feature_columns.pkl')
        if os.path.exists(features_path):
            with open(features_path, 'rb') as f:
                self.feature_columns = pickle.load(f)


def main():
    """Test the feature extractor without data leakage."""
    extractor = GenomicFeatureExtractorNoLeakage()
    
    # Load training data
    X, y = extractor.load_training_data()
    
    print("=== NO DATA LEAKAGE VERSION ===")
    print("Feature columns:", X.columns.tolist())
    print("Feature shape:", X.shape)
    print("Label distribution:", y.value_counts())
    
    # Show that clinical significance features are excluded
    clinical_features = [col for col in X.columns if 'clinvar' in col.lower()]
    print(f"Clinical significance features: {clinical_features}")  # Should be empty!
    
    # Prepare for training
    X_train, X_test, y_train, y_test = extractor.prepare_for_training(X, y)
    
    print("Training data prepared successfully!")
    print(f"Training set shape: {X_train.shape}")
    print(f"Test set shape: {X_test.shape}")


if __name__ == "__main__":
    main() 