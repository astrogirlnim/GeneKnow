import React from 'react';
import type { PipelineResult } from '../api/geneknowPipeline';

// Static model validation metrics from training
const MODEL_VALIDATION_METRICS = {
  auc_roc: 0.7628,
  f1_score: 0.5712,
  matthews_corrcoef: 0.4270,
  accuracy: 0.5712,
  r2_score: 0.4270,
  mse: 0.07947,
  training_date: "2025-07-10",
  total_variants: 200000,
  best_model: "gradient_boosting"
};

interface ModelPerformanceTabProps {
  pipelineResults?: PipelineResult;
}

// Icon component for tooltips
const InformationCircleIcon = ({ className = "w-5 h-5" }) => (
  <svg className={className} xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth="2" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
  </svg>
);

// Tooltip component
const Tooltip = ({ content, link }: { content: string; link?: string }) => (
  <div style={{
    position: 'absolute',
    bottom: '100%',
    marginBottom: '0.5rem',
    width: '16rem',
    padding: '0.75rem',
    background: '#1F2937',
    color: '#FFFFFF',
    fontSize: '0.875rem',
    borderRadius: '0.5rem',
    boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
    opacity: 0,
    transition: 'opacity 300ms ease',
    zIndex: 10,
    pointerEvents: 'none'
  }}
  className="group-hover:opacity-100">
    <p>{content}</p>
    {link && (
      <a href={link} target="_blank" rel="noopener noreferrer" style={{ color: '#60A5FA', textDecoration: 'underline', marginTop: '0.25rem', display: 'inline-block' }}>
        Learn more
      </a>
    )}
  </div>
);

// Metric card component
const MetricCard = ({ 
  title, 
  value, 
  unit, 
  description, 
  isReference = false 
}: { 
  title: string; 
  value: string | number; 
  unit: string; 
  description: string;
  isReference?: boolean;
}) => (
  <div style={{
    background: isReference ? '#F0F9FF' : '#FFFFFF',
    padding: '1rem',
    borderRadius: '0.75rem',
    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
    border: isReference ? '1px solid #BAE6FD' : '1px solid #E5E7EB',
    minWidth: 0,
  }}>
    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
      <h4 style={{ 
        fontWeight: '600', 
        color: isReference ? '#0369A1' : '#4B5563',
        fontSize: '0.875rem'
      }}>
        {title}
      </h4>
      <div style={{ position: 'relative' }} className="group">
        <InformationCircleIcon className="w-4 h-4 cursor-pointer text-gray-400" />
        <Tooltip content={description} />
      </div>
    </div>
    <p style={{ 
      fontSize: '1.5rem', 
      fontWeight: 'bold', 
      color: isReference ? '#0369A1' : '#111827',
      marginTop: '0.5rem',
      overflow: 'hidden',
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap'
    }}>
      {value} <span style={{ fontSize: '1rem', fontWeight: '500', color: '#6B7280' }}>{unit}</span>
    </p>
  </div>
);

// Risk category distribution component
const RiskCategoryDistribution = ({ 
  categories 
}: { 
  categories: Record<string, number> 
}) => {
  const total = Object.values(categories).reduce((sum, count) => sum + count, 0);
  
  return (
    <div style={{
      background: '#FFFFFF',
      padding: '1.5rem',
      borderRadius: '0.75rem',
      boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
      border: '1px solid #E5E7EB'
    }}>
      <h4 style={{ 
        fontWeight: '600', 
        color: '#111827', 
        marginBottom: '1rem',
        fontSize: '1.125rem'
      }}>
        Risk Category Distribution
      </h4>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(100px, 1fr))', gap: '1rem' }}>
        {Object.entries(categories).map(([category, count]) => {
          const percentage = total > 0 ? (count / total * 100).toFixed(1) : '0';
          const color = category === 'high' ? '#DC2626' : 
                       category === 'moderate' ? '#F59E0B' : 
                       category === 'low' ? '#059669' : '#6B7280';
          
          return (
            <div key={category} style={{ textAlign: 'center' }}>
              <div style={{
                width: '60px',
                height: '60px',
                borderRadius: '50%',
                background: color,
                color: '#FFFFFF',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                fontSize: '1.25rem',
                fontWeight: 'bold',
                margin: '0 auto 0.5rem'
              }}>
                {count}
              </div>
              <p style={{ 
                fontSize: '0.875rem', 
                fontWeight: '600', 
                color: '#374151',
                textTransform: 'capitalize'
              }}>
                {category}
              </p>
              <p style={{ 
                fontSize: '0.75rem', 
                color: '#6B7280'
              }}>
                {percentage}%
              </p>
            </div>
          );
        })}
      </div>
    </div>
  );
};

const ModelPerformanceTab: React.FC<ModelPerformanceTabProps> = ({ pipelineResults }) => {
  console.log('🏆 ModelPerformanceTab: Rendering with pipelineResults:', !!pipelineResults);
  console.log('🏆 ModelPerformanceTab: Metrics available:', !!pipelineResults?.structured_json?.metrics);
  
  // Extract metrics from pipeline results
  const metrics = pipelineResults?.structured_json?.metrics;
  
  console.log('🏆 ModelPerformanceTab: Full metrics object:', metrics);
  console.log('🏆 ModelPerformanceTab: Confidence metrics:', metrics?.confidence_metrics);
  console.log('🏆 ModelPerformanceTab: Variant metrics:', metrics?.variant_metrics);
  console.log('🏆 ModelPerformanceTab: Integration metrics:', metrics?.integration_metrics);
  console.log('🏆 ModelPerformanceTab: Performance indicators:', metrics?.performance_indicators);

  // Extract user analysis metrics
  const userMetrics = React.useMemo(() => {
    if (!metrics) return null;
    
    const confidence = metrics.confidence_metrics || {};
    const variant = metrics.variant_metrics || {};
    const integration = metrics.integration_metrics || {};
    const performance = metrics.performance_indicators || {};
    
    console.log('🏆 ModelPerformanceTab: Processing confidence metrics:', confidence);
    console.log('🏆 ModelPerformanceTab: Processing variant metrics:', variant);
    console.log('🏆 ModelPerformanceTab: Processing integration metrics:', integration);
    console.log('🏆 ModelPerformanceTab: Processing performance metrics:', performance);
    
    return {
      confidence: confidence,
      variant: variant,
      integration: integration,
      performance: performance
    };
  }, [metrics]);

  // Calculate risk category distribution
  const riskCategories = React.useMemo(() => {
    if (!pipelineResults?.risk_scores) return { low: 0, moderate: 0, high: 0 };
    
    const scores = Object.values(pipelineResults.risk_scores);
    const categories = { low: 0, moderate: 0, high: 0 };
    
    scores.forEach(score => {
      if (score >= 50) categories.high++;
      else if (score >= 20) categories.moderate++;
      else categories.low++;
    });
    
    console.log('🏆 ModelPerformanceTab: Risk categories calculated:', categories);
    return categories;
  }, [pipelineResults?.risk_scores]);

  return (
    <div style={{
      padding: '1.5rem',
      background: '#FFFFFF',
      borderRadius: '0.75rem',
      boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
      border: '1px solid #E5E7EB'
    }}>
      <div style={{ marginBottom: '2rem' }}>
        <h3 style={{ 
          fontWeight: '600', 
          fontSize: '1.5rem', 
          color: '#111827',
          marginBottom: '0.5rem'
        }}>
          Model Performance
        </h3>
        <p style={{ 
          color: '#6B7280', 
          fontSize: '0.875rem' 
        }}>
          Analysis metrics for your sample and model validation reference data
        </p>
      </div>

      {/* Your Analysis Metrics */}
      <div style={{ marginBottom: '2rem' }}>
        <h4 style={{ 
          fontWeight: '600', 
          fontSize: '1.25rem', 
          color: '#111827',
          marginBottom: '1rem'
        }}>
          Your Analysis Metrics
        </h4>
        
        {userMetrics ? (
          <>
            {/* Model Confidence Section */}
            <div style={{ marginBottom: '1.5rem' }}>
              <h5 style={{ 
                fontWeight: '600', 
                fontSize: '1rem', 
                color: '#374151',
                marginBottom: '0.75rem'
              }}>
                Model Confidence
              </h5>
              <div style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
                gap: '1rem',
                marginBottom: '1rem'
              }}>
                <MetricCard
                  title="Mean Confidence"
                  value={userMetrics.confidence.mean_model_confidence?.toFixed(3) || 'N/A'}
                  unit=""
                  description="Average confidence of model predictions across all variants"
                />
                <MetricCard
                  title="Min Confidence"
                  value={userMetrics.confidence.min_model_confidence?.toFixed(3) || 'N/A'}
                  unit=""
                  description="Lowest confidence score for any prediction"
                />
                <MetricCard
                  title="Max Confidence"
                  value={userMetrics.confidence.max_model_confidence?.toFixed(3) || 'N/A'}
                  unit=""
                  description="Highest confidence score for any prediction"
                />
                <MetricCard
                  title="ML Fusion Confidence"
                  value={userMetrics.confidence.ml_fusion_confidence?.toFixed(3) || 'N/A'}
                  unit=""
                  description="Confidence score from the ML fusion layer"
                />
              </div>
            </div>

            {/* Risk Score Distribution */}
            <div style={{ marginBottom: '1.5rem' }}>
              <h5 style={{ 
                fontWeight: '600', 
                fontSize: '1rem', 
                color: '#374151',
                marginBottom: '0.75rem'
              }}>
                Risk Score Distribution
              </h5>
              <div style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
                gap: '1rem',
                marginBottom: '1rem'
              }}>
                <MetricCard
                  title="Mean Risk Score"
                  value={userMetrics.confidence.risk_score_mean?.toFixed(3) || 'N/A'}
                  unit=""
                  description="Average risk score across all cancer types"
                />
                <MetricCard
                  title="Risk Score StdDev"
                  value={userMetrics.confidence.risk_score_std?.toFixed(3) || 'N/A'}
                  unit=""
                  description="Standard deviation of risk scores"
                />
                <MetricCard
                  title="Max Risk Score"
                  value={userMetrics.confidence.max_risk_score?.toFixed(3) || 'N/A'}
                  unit=""
                  description="Highest risk score detected"
                />
                <MetricCard
                  title="Coefficient of Variation"
                  value={userMetrics.confidence.risk_score_cv?.toFixed(3) || 'N/A'}
                  unit=""
                  description="Risk score variability (std/mean)"
                />
              </div>
            </div>

            {/* Variant Metrics */}
            <div style={{ marginBottom: '1.5rem' }}>
              <h5 style={{ 
                fontWeight: '600', 
                fontSize: '1rem', 
                color: '#374151',
                marginBottom: '0.75rem'
              }}>
                Variant Metrics
              </h5>
              <div style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
                gap: '1rem',
                marginBottom: '1rem'
              }}>
                <MetricCard
                  title="Total Variants"
                  value={userMetrics.variant.total_variants || 'N/A'}
                  unit="variants"
                  description="Total number of genetic variants identified"
                />
                <MetricCard
                  title="Pathogenic Variants"
                  value={userMetrics.variant.pathogenic_variants || 'N/A'}
                  unit="variants"
                  description="Number of variants classified as pathogenic"
                />
                <MetricCard
                  title="Benign Variants"
                  value={userMetrics.variant.benign_variants || 'N/A'}
                  unit="variants"
                  description="Number of variants classified as benign"
                />
                <MetricCard
                  title="Uncertain Variants"
                  value={userMetrics.variant.uncertain_variants || 'N/A'}
                  unit="variants"
                  description="Number of variants with uncertain significance"
                />
                <MetricCard
                  title="Mean CADD Score"
                  value={userMetrics.variant.mean_cadd_score?.toFixed(1) || 'N/A'}
                  unit=""
                  description="Average CADD deleteriousness score"
                />
                <MetricCard
                  title="Max CADD Score"
                  value={userMetrics.variant.max_cadd_score?.toFixed(1) || 'N/A'}
                  unit=""
                  description="Highest CADD score detected"
                />
              </div>
            </div>

            {/* PRS and Pathway Metrics */}
            <div style={{ marginBottom: '1.5rem' }}>
              <h5 style={{ 
                fontWeight: '600', 
                fontSize: '1rem', 
                color: '#374151',
                marginBottom: '0.75rem'
              }}>
                PRS and Pathway Burden
              </h5>
              <div style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
                gap: '1rem',
                marginBottom: '1rem'
              }}>
                                 <MetricCard
                   title="PRS Confidence"
                   value={(userMetrics.integration as any).prs_confidence || 'N/A'}
                   unit=""
                   description="Polygenic Risk Score confidence level"
                 />
                 <MetricCard
                   title="High-Risk Cancers"
                   value={(userMetrics.integration as any).prs_high_risk_cancers?.length || 0}
                   unit="types"
                   description="Number of cancer types with high PRS risk"
                 />
                 <MetricCard
                   title="Pathway Burden Score"
                   value={(userMetrics.integration as any).pathway_burden_score?.toFixed(3) || 'N/A'}
                   unit=""
                   description="Overall pathway disruption burden"
                 />
                 <MetricCard
                   title="High-Burden Pathways"
                   value={(userMetrics.integration as any).high_burden_pathways?.length || 0}
                   unit="pathways"
                   description="Number of pathways with high burden"
                 />
              </div>
            </div>

            {/* Risk Category Distribution */}
            <RiskCategoryDistribution categories={riskCategories} />
          </>
        ) : (
          <div style={{
            padding: '2rem',
            textAlign: 'center',
            background: '#F9FAFB',
            borderRadius: '0.5rem',
            border: '1px solid #E5E7EB'
          }}>
            <p style={{ color: '#374151', fontWeight: '600' }}>
              No Analysis Metrics Available
            </p>
            <p style={{ color: '#6B7280', fontSize: '0.875rem', marginTop: '0.5rem' }}>
              Metrics are only available for completed pipeline analyses with structured results.
            </p>
          </div>
        )}
      </div>

      {/* Model Validation (Reference Only) */}
      <div style={{ marginBottom: '1rem' }}>
        <h4 style={{ 
          fontWeight: '600', 
          fontSize: '1.25rem', 
          color: '#111827',
          marginBottom: '0.5rem'
        }}>
          Model Validation (Reference Only)
        </h4>
        <p style={{ 
          color: '#6B7280', 
          fontSize: '0.875rem',
          marginBottom: '1rem',
          fontStyle: 'italic'
        }}>
          These metrics are from model validation on public datasets and do not reflect your individual results.
        </p>
        
        <div style={{
          display: 'grid',
          gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
          gap: '1rem'
        }}>
          <MetricCard
            title="AUC-ROC"
            value={MODEL_VALIDATION_METRICS.auc_roc.toFixed(3)}
            unit=""
            description="Area Under the Curve - Receiver Operating Characteristic"
            isReference={true}
          />
          <MetricCard
            title="F1-Score"
            value={MODEL_VALIDATION_METRICS.f1_score.toFixed(3)}
            unit=""
            description="Harmonic mean of precision and recall"
            isReference={true}
          />
          <MetricCard
            title="Matthews Correlation"
            value={MODEL_VALIDATION_METRICS.matthews_corrcoef.toFixed(3)}
            unit=""
            description="Correlation coefficient between predicted and actual classifications"
            isReference={true}
          />
          <MetricCard
            title="Accuracy"
            value={MODEL_VALIDATION_METRICS.accuracy.toFixed(3)}
            unit=""
            description="Proportion of correct predictions"
            isReference={true}
          />
          <MetricCard
            title="R² Score"
            value={MODEL_VALIDATION_METRICS.r2_score.toFixed(3)}
            unit=""
            description="Coefficient of determination"
            isReference={true}
          />
          <MetricCard
            title="Mean Squared Error"
            value={MODEL_VALIDATION_METRICS.mse.toFixed(5)}
            unit=""
            description="Average squared difference between predicted and actual values"
            isReference={true}
          />
        </div>
        
        <div style={{
          marginTop: '1rem',
          padding: '1rem',
          background: '#EFF6FF',
          borderRadius: '0.5rem',
          border: '1px solid #BFDBFE'
        }}>
          <p style={{ 
            color: '#1E40AF', 
            fontSize: '0.875rem',
            fontWeight: '500'
          }}>
            <strong>Training Details:</strong> Model: {MODEL_VALIDATION_METRICS.best_model}, 
            Variants: {MODEL_VALIDATION_METRICS.total_variants.toLocaleString()}, 
            Date: {MODEL_VALIDATION_METRICS.training_date}
          </p>
        </div>
      </div>
    </div>
  );
};

export default ModelPerformanceTab; 