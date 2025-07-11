import React from 'react';
import ConfidenceCheck from './ConfidenceCheck';
import type { SHAPValidation } from '../api/geneknowPipeline';

// Mock data for testing all states
const mockValidations: Record<string, SHAPValidation> = {
  flagForReview: {
    status: 'FLAG_FOR_REVIEW',
    reasons: [
      'The AI predicted HIGH RISK (75%) but this appears to be based on indirect factors (Gene/Pathway Burden Score, TCGA Tumor Enrichment) rather than known disease-causing mutations. The prediction may be less reliable.'
    ],
    top_contributors: [
      {
        feature: 'gene_burden_score',
        display_name: 'Gene/Pathway Burden Score',
        shap_value: 0.25,
        abs_contribution: 0.25,
        direction: 'increases'
      },
      {
        feature: 'tcga_enrichment',
        display_name: 'TCGA Tumor Enrichment',
        shap_value: 0.18,
        abs_contribution: 0.18,
        direction: 'increases'
      },
      {
        feature: 'prs_score',
        display_name: 'Polygenic Risk Score',
        shap_value: 0.12,
        abs_contribution: 0.12,
        direction: 'increases'
      }
    ],
    feature_importance: {
      'gene_burden_score': 0.25,
      'tcga_enrichment': 0.18,
      'prs_score': 0.12,
      'clinvar_pathogenic': 0.05,
      'clinvar_benign': -0.08,
      'cadd_score': 0.10
    },
    details: {
      status: 'FLAG_FOR_REVIEW',
      risk_score: 0.75,
      top_contributors: [],
      validation_reasons: [],
      rule_results: {},
      shap_values: [],
      feature_names: [],
      model_type: 'GradientBoostingRegressor'
    }
  },

  pass: {
    status: 'PASS',
    reasons: [],
    top_contributors: [
      {
        feature: 'clinvar_pathogenic',
        display_name: 'ClinVar Pathogenic Variant',
        shap_value: 0.45,
        abs_contribution: 0.45,
        direction: 'increases'
      },
      {
        feature: 'cadd_score',
        display_name: 'CADD Deleteriousness Score',
        shap_value: 0.22,
        abs_contribution: 0.22,
        direction: 'increases'
      },
      {
        feature: 'gene_burden_score',
        display_name: 'Gene/Pathway Burden Score',
        shap_value: 0.15,
        abs_contribution: 0.15,
        direction: 'increases'
      }
    ],
    feature_importance: {
      'clinvar_pathogenic': 0.45,
      'cadd_score': 0.22,
      'gene_burden_score': 0.15,
      'tcga_enrichment': 0.08,
      'prs_score': 0.10
    },
    details: {
      status: 'PASS',
      risk_score: 0.82,
      top_contributors: [],
      validation_reasons: [],
      rule_results: {},
      shap_values: [],
      feature_names: [],
      model_type: 'GradientBoostingRegressor'
    }
  },

  error: {
    status: 'ERROR',
    reasons: ['SHAP validation error: Model structure incompatible with explainer'],
    top_contributors: [],
    feature_importance: {},
    details: {
      status: 'ERROR',
      risk_score: 0.0,
      top_contributors: [],
      validation_reasons: ['SHAP validation error: Model structure incompatible with explainer'],
      rule_results: {},
      shap_values: [],
      feature_names: [],
      model_type: 'Unknown'
    }
  },

  skipped: {
    status: 'SKIPPED',
    reasons: ['ML fusion model or features not available'],
    top_contributors: [],
    feature_importance: {},
    details: {
      status: 'SKIPPED',
      risk_score: 0.0,
      top_contributors: [],
      validation_reasons: ['ML fusion model or features not available'],
      rule_results: {},
      shap_values: [],
      feature_names: [],
      model_type: 'None'
    }
  }
};

const ConfidenceCheckTest: React.FC = () => {
  const handleNavigate = () => {
    alert('Would navigate to In-Depth Analysis');
  };

  return (
    <div style={{ padding: '2rem', maxWidth: '800px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2rem', fontWeight: 'bold', marginBottom: '2rem', color: '#111827' }}>
        Confidence Check Component Test
      </h1>
      
      <div style={{ display: 'flex', flexDirection: 'column', gap: '3rem' }}>
        
        {/* Dashboard Summary Views */}
        <section>
          <h2 style={{ fontSize: '1.5rem', fontWeight: '600', marginBottom: '1rem', color: '#111827' }}>
            Dashboard Summary Views
          </h2>
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', gap: '1rem' }}>
            {Object.entries(mockValidations).map(([key, validation]) => (
              <div key={key} style={{
                padding: '1rem',
                background: '#FFFFFF',
                border: '1px solid #E5E7EB',
                borderRadius: '0.5rem',
                boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
              }}>
                <h3 style={{ fontSize: '1rem', fontWeight: '600', marginBottom: '0.5rem', textTransform: 'capitalize' }}>
                  {key.replace(/([A-Z])/g, ' $1')} State
                </h3>
                <ConfidenceCheck 
                  validation={validation}
                  onNavigateToDetail={handleNavigate}
                  isDetailed={false}
                />
              </div>
            ))}
          </div>
        </section>

        {/* Detailed Clinical Alerts Views */}
        <section>
          <h2 style={{ fontSize: '1.5rem', fontWeight: '600', marginBottom: '1rem', color: '#111827' }}>
            Clinical Alerts Detailed Views
          </h2>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
            {Object.entries(mockValidations).map(([key, validation]) => (
              <div key={key} style={{
                padding: '1rem',
                background: '#F9FAFB',
                border: '1px solid #E5E7EB',
                borderRadius: '0.5rem'
              }}>
                <h3 style={{ fontSize: '1rem', fontWeight: '600', marginBottom: '1rem', textTransform: 'capitalize' }}>
                  {key.replace(/([A-Z])/g, ' $1')} State - Detailed View
                </h3>
                <ConfidenceCheck 
                  validation={validation}
                  isDetailed={true}
                />
              </div>
            ))}
          </div>
        </section>

        {/* Null/Empty State */}
        <section>
          <h2 style={{ fontSize: '1.5rem', fontWeight: '600', marginBottom: '1rem', color: '#111827' }}>
            Null State (Component Hidden)
          </h2>
          <div style={{
            padding: '1rem',
            background: '#F3F4F6',
            border: '1px solid #E5E7EB',
            borderRadius: '0.5rem',
            textAlign: 'center'
          }}>
            <ConfidenceCheck validation={null} isDetailed={false} />
            <p style={{ color: '#6B7280', fontStyle: 'italic' }}>
              Component returns null when validation data is not available
            </p>
          </div>
        </section>

      </div>
    </div>
  );
};

export default ConfidenceCheckTest; 