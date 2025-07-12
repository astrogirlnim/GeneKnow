import React, { useState, useRef, useEffect } from 'react';
import type { PipelineResult } from '../api/geneknowPipeline';

// Static model validation metrics from actual training results
// These metrics are from real model training/validation on public datasets
const MODEL_VALIDATION_METRICS = {
  auc_roc: 0.7628,
  accuracy: 0.5712,
  r2_score: 0.4270,
  mse: 0.07947,
  training_date: "2025-07-10",
  total_variants: 200000,
  best_model: "gradient_boosting",
  // Note: F1-score and Matthews correlation coefficient not available in training results
  // These would require additional calculation during training
  f1_score: undefined as number | undefined,
  matthews_corrcoef: undefined as number | undefined
};

interface ModelPerformanceTabProps {
  pipelineResults?: PipelineResult;
}

// Icon component for tooltips with hover state management
const InformationCircleIcon = React.forwardRef<HTMLDivElement, { 
  style?: React.CSSProperties; 
  onMouseEnter?: (e: React.MouseEvent<HTMLDivElement>) => void;
  onMouseLeave?: (e: React.MouseEvent<HTMLDivElement>) => void;
}>(({ style, onMouseEnter, onMouseLeave }, ref) => (
  <div 
    ref={ref}
    style={{
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      width: '1rem',
      height: '1rem',
      borderRadius: '50%',
      backgroundColor: '#9CA3AF',
      color: '#FFFFFF',
      fontSize: '0.75rem',
      fontWeight: 'bold',
      cursor: 'pointer',
      transition: 'all 0.2s ease',
      ...style
    }}
    onMouseEnter={onMouseEnter}
    onMouseLeave={onMouseLeave}
  >
    i
  </div>
));

// Enhanced tooltip component with smart positioning (adapted from Clinical View)
const SmartTooltip = ({ content, isVisible, triggerRef }: { 
  content: string; 
  isVisible: boolean; 
  triggerRef?: React.RefObject<HTMLDivElement | null> | null;
}) => {
  const [position, setPosition] = useState<'right' | 'above' | 'below'>('right');
  const tooltipRef = useRef<HTMLDivElement>(null);
  
  useEffect(() => {
    if (isVisible && triggerRef?.current && tooltipRef.current) {
      const triggerRect = triggerRef.current.getBoundingClientRect();
      const viewportWidth = window.innerWidth;
      const viewportHeight = window.innerHeight;
      
      // Check if there's enough space to the right
      const spaceToRight = viewportWidth - triggerRect.right;
      const tooltipWidth = 288; // 18rem = 288px
      
      // Check if there's enough space above/below
      const spaceAbove = triggerRect.top;
      const spaceBelow = viewportHeight - triggerRect.bottom;
      const tooltipHeight = 100; // Approximate height
      
      // Determine best position
      if (spaceToRight >= tooltipWidth + 16) {
        setPosition('right');
      } else if (spaceAbove >= tooltipHeight + 16) {
        setPosition('above');
      } else if (spaceBelow >= tooltipHeight + 16) {
        setPosition('below');
      } else {
        setPosition('right');
      }
    }
  }, [isVisible, triggerRef]);
  
  const getTooltipStyles = () => {
    const baseStyles = {
      position: 'absolute' as const,
      width: '18rem',
      padding: '0.75rem',
      background: '#1F2937',
      color: '#FFFFFF',
      fontSize: '0.75rem',
      borderRadius: '0.5rem',
      boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
      opacity: isVisible ? 1 : 0,
      visibility: isVisible ? 'visible' as const : 'hidden' as const,
      transition: 'opacity 300ms ease, visibility 300ms ease',
      zIndex: 1000,
      pointerEvents: 'none' as const,
      lineHeight: '1.4',
      border: '1px solid #374151'
    };
    
    switch (position) {
      case 'right':
        return {
          ...baseStyles,
          left: 'calc(100% + 0.5rem)',
          top: '50%',
          transform: 'translateY(-50%)',
        };
      case 'above':
        return {
          ...baseStyles,
          left: '50%',
          bottom: 'calc(100% + 0.5rem)',
          transform: 'translateX(-50%)',
        };
      case 'below':
        return {
          ...baseStyles,
          left: '50%',
          top: 'calc(100% + 0.5rem)',
          transform: 'translateX(-50%)',
        };
      default:
        return baseStyles;
    }
  };
  
  const getArrowStyles = () => {
    const baseArrowStyles = {
      position: 'absolute' as const,
      width: '0',
      height: '0',
    };
    
    switch (position) {
      case 'right':
        return {
          ...baseArrowStyles,
          left: '-0.5rem',
          top: '50%',
          transform: 'translateY(-50%)',
          borderTop: '0.5rem solid transparent',
          borderBottom: '0.5rem solid transparent',
          borderRight: '0.5rem solid #1F2937'
        };
      case 'above':
        return {
          ...baseArrowStyles,
          left: '50%',
          bottom: '-0.5rem',
          transform: 'translateX(-50%)',
          borderLeft: '0.5rem solid transparent',
          borderRight: '0.5rem solid transparent',
          borderTop: '0.5rem solid #1F2937'
        };
      case 'below':
        return {
          ...baseArrowStyles,
          left: '50%',
          top: '-0.5rem',
          transform: 'translateX(-50%)',
          borderLeft: '0.5rem solid transparent',
          borderRight: '0.5rem solid transparent',
          borderBottom: '0.5rem solid #1F2937'
        };
      default:
        return baseArrowStyles;
    }
  };
  
  return (
    <div ref={tooltipRef} style={getTooltipStyles()}>
      <p style={{ color: '#D1D5DB', marginBottom: '0' }}>{content}</p>
      <div style={getArrowStyles()}></div>
    </div>
  );
};

// Metric card component with enhanced tooltip system
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
}) => {
  const [showTooltip, setShowTooltip] = useState(false);
  const triggerRef = useRef<HTMLDivElement>(null);

  return (
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
        <div style={{ position: 'relative' }}>
          <InformationCircleIcon 
            ref={triggerRef}
            style={{
              backgroundColor: showTooltip ? '#6B7280' : '#9CA3AF'
            }}
            onMouseEnter={() => setShowTooltip(true)}
            onMouseLeave={() => setShowTooltip(false)}
          />
          <SmartTooltip 
            content={description} 
            isVisible={showTooltip} 
            triggerRef={triggerRef}
          />
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
};

// Risk category distribution component
const RiskCategoryDistribution = ({ 
  categories 
}: { 
  categories: Record<string, number> 
}) => {
  const total = Object.values(categories).reduce((sum, count) => sum + count, 0);
  const [showTooltip, setShowTooltip] = useState(false);
  const triggerRef = useRef<HTMLDivElement>(null);
  
  return (
    <div style={{
      background: '#FFFFFF',
      padding: '1.5rem',
      borderRadius: '0.75rem',
      boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
      border: '1px solid #E5E7EB'
    }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
        <h4 style={{ 
          fontWeight: '600', 
          color: '#111827', 
          fontSize: '1.125rem'
        }}>
          Risk Category Distribution
        </h4>
        <div style={{ position: 'relative' }}>
          <InformationCircleIcon 
            ref={triggerRef}
            style={{
              backgroundColor: showTooltip ? '#6B7280' : '#9CA3AF'
            }}
            onMouseEnter={() => setShowTooltip(true)}
            onMouseLeave={() => setShowTooltip(false)}
          />
          <SmartTooltip 
            content="Distribution of cancer risk predictions across risk categories. Low risk represents baseline or below-average risk, moderate risk indicates elevated risk requiring monitoring, and high risk suggests significant predisposition warranting enhanced screening and preventive measures. This visualization helps understand your overall risk profile." 
            isVisible={showTooltip} 
            triggerRef={triggerRef}
          />
        </div>
      </div>
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
                  value={userMetrics.confidence.mean_model_confidence?.toFixed(3) ?? 'N/A'}
                  unit=""
                  description="Average confidence level of the AI model's predictions across all cancer types analyzed. Values range from 0-1, where higher values indicate the model is more certain about its risk assessments. This metric helps assess the reliability of your overall risk profile. Scores above 0.7 are considered high confidence."
                />
                <MetricCard
                  title="Min Confidence"
                  value={userMetrics.confidence.min_model_confidence?.toFixed(3) ?? 'N/A'}
                  unit=""
                  description="Lowest confidence score among all cancer risk predictions. This represents the cancer type where the model has the least certainty in its assessment. Lower values may indicate insufficient genetic information or complex variant patterns requiring additional clinical evaluation. Values below 0.5 suggest caution in interpretation."
                />
                <MetricCard
                  title="Max Confidence"
                  value={userMetrics.confidence.max_model_confidence?.toFixed(3) ?? 'N/A'}
                  unit=""
                  description="Highest confidence score among all cancer risk predictions. This represents the cancer type where the model has the strongest certainty in its assessment, typically based on clear genetic signals and well-characterized pathogenic variants. High values (>0.8) indicate strong genetic evidence."
                />
                <MetricCard
                  title="ML Fusion Confidence"
                  value={userMetrics.confidence.ml_fusion_confidence?.toFixed(3) ?? 'N/A'}
                  unit=""
                  description="Overall confidence score from the machine learning fusion model that combines multiple data sources including genetic variants, pathway analysis, and polygenic risk scores. This represents the model's confidence in the integrated risk assessment across all cancer types. The fusion model provides more robust predictions than individual models."
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
                  value={userMetrics.confidence.risk_score_mean?.toFixed(3) ?? 'N/A'}
                  unit=""
                  description="Average cancer risk score across all cancer types analyzed. Risk scores represent the estimated percentage likelihood of developing cancer based on your genetic variants. Higher scores indicate elevated risk compared to the general population. Scores are calculated using machine learning models trained on cancer genomics data."
                />
                <MetricCard
                  title="Risk Score StdDev"
                  value={userMetrics.confidence.risk_score_std?.toFixed(3) ?? 'N/A'}
                  unit=""
                  description="Standard deviation of risk scores across cancer types, measuring how much risk scores vary from the average. Higher values indicate greater variability in risk between different cancer types. Low standard deviation suggests consistent risk levels across cancers, while high values indicate specific high-risk cancers."
                />
                <MetricCard
                  title="Max Risk Score"
                  value={userMetrics.confidence.max_risk_score?.toFixed(3) ?? 'N/A'}
                  unit=""
                  description="Highest cancer risk score detected among all cancer types analyzed. This represents your greatest genetic predisposition to any single cancer type. Values above 50 are considered high risk and may warrant enhanced screening or preventive measures. This score identifies your primary cancer risk concern."
                />
                <MetricCard
                  title="Coefficient of Variation"
                  value={userMetrics.confidence.risk_score_cv?.toFixed(3) ?? 'N/A'}
                  unit=""
                  description="Coefficient of variation (standard deviation divided by mean) for risk scores, providing a normalized measure of risk variability. Values above 0.5 indicate high variability between cancer types, suggesting specific genetic predispositions rather than general cancer risk. Lower values indicate more uniform risk distribution."
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
                  value={userMetrics.variant.total_variants ?? 'N/A'}
                  unit="variants"
                  description="Total number of genetic variants identified in your genomic data after quality filtering. This includes single nucleotide variants (SNVs), insertions/deletions (indels), and other genetic changes. Higher numbers indicate more comprehensive genomic coverage but don't necessarily correlate with increased cancer risk."
                />
                <MetricCard
                  title="Pathogenic Variants"
                  value={userMetrics.variant.pathogenic_variants ?? 'N/A'}
                  unit="variants"
                  description="Number of variants classified as pathogenic or likely pathogenic according to clinical databases like ClinVar. These variants have established evidence for causing disease or significantly increasing cancer risk. Even one pathogenic variant in a cancer-associated gene can substantially affect risk and screening recommendations."
                />
                <MetricCard
                  title="Benign Variants"
                  value={userMetrics.variant.benign_variants ?? 'N/A'}
                  unit="variants"
                  description="Number of variants classified as benign or likely benign, meaning they have no known harmful effects on health. These variants represent normal genetic variation and don't contribute to cancer risk. A higher proportion of benign variants generally indicates a lower genetic risk profile."
                />
                <MetricCard
                  title="Uncertain Variants"
                  value={userMetrics.variant.uncertain_variants ?? 'N/A'}
                  unit="variants"
                  description="Number of variants with uncertain significance (VUS) where the clinical impact is unknown or conflicting. These variants require careful interpretation and may be reclassified as more data becomes available. They are typically not used for clinical decision-making but are monitored for future updates."
                />
                <MetricCard
                  title="Mean CADD Score"
                  value={userMetrics.variant.mean_cadd_score?.toFixed(1) ?? 'N/A'}
                  unit=""
                  description="Average CADD (Combined Annotation Dependent Depletion) score across all variants. CADD scores predict the deleteriousness of genetic variants on a scale from 0-40+. Scores above 20 are considered highly deleterious, 15-20 moderately deleterious, and below 15 less likely to be harmful. Higher average scores suggest more potentially damaging variants."
                />
                <MetricCard
                  title="Max CADD Score"
                  value={userMetrics.variant.max_cadd_score?.toFixed(1) ?? 'N/A'}
                  unit=""
                  description="Highest CADD score detected among all variants, representing the most potentially damaging genetic change identified. Scores above 30 are extremely rare and highly deleterious, 20-30 are considered pathogenic, and 15-20 are moderately damaging. This metric identifies your most concerning genetic variant from a functional impact perspective."
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
                   value={(userMetrics.integration as any).prs_confidence ?? 'N/A'}
                   unit=""
                   description="Polygenic Risk Score confidence level based on SNP coverage and data quality. PRS combines many small genetic effects to estimate cancer risk. High confidence indicates good coverage of known risk variants. Low confidence suggests limited genetic information and results should be interpreted cautiously."
                 />
                 <MetricCard
                   title="High-Risk Cancers"
                   value={(userMetrics.integration as any).prs_high_risk_cancers?.length ?? 0}
                   unit="types"
                   description="Number of cancer types where your polygenic risk score falls in the high-risk category (typically >95th percentile). These cancers show elevated genetic predisposition based on cumulative effects of multiple common variants. Enhanced screening may be recommended for these cancer types."
                 />
                 <MetricCard
                   title="Pathway Burden Score"
                   value={(userMetrics.integration as any).pathway_burden_score?.toFixed(3) ?? 'N/A'}
                   unit=""
                   description="Overall pathway disruption burden score representing the cumulative impact of genetic variants on biological pathways. Scores range from 0-1, where higher values indicate more pathway disruption. Values above 0.5 suggest significant pathway dysfunction that may contribute to cancer risk."
                 />
                 <MetricCard
                   title="High-Burden Pathways"
                   value={(userMetrics.integration as any).high_burden_pathways?.length ?? 0}
                   unit="pathways"
                   description="Number of biological pathways with high disruption burden from genetic variants. These pathways show significant dysfunction that may contribute to cancer development. Common high-burden pathways include DNA repair, cell cycle control, and tumor suppression. Multiple disrupted pathways indicate increased cancer risk."
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
            description="Area Under the Receiver Operating Characteristic curve measuring the model's ability to distinguish between cancer and non-cancer cases. Values range from 0-1, where 0.5 is random chance and 1.0 is perfect classification. Our model achieves 0.763, indicating good discriminative ability. Values above 0.7 are considered clinically useful."
            isReference={true}
          />
          {MODEL_VALIDATION_METRICS.f1_score !== undefined && (
            <MetricCard
              title="F1-Score"
              value={MODEL_VALIDATION_METRICS.f1_score!.toFixed(3)}
              unit=""
              description="Harmonic mean of precision and recall, providing a balanced measure of model performance. F1-score considers both false positives and false negatives, making it ideal for medical applications where both types of errors are costly. Values closer to 1.0 indicate better overall performance."
              isReference={true}
            />
          )}
          {MODEL_VALIDATION_METRICS.matthews_corrcoef !== undefined && (
            <MetricCard
              title="Matthews Correlation"
              value={MODEL_VALIDATION_METRICS.matthews_corrcoef!.toFixed(3)}
              unit=""
              description="Matthews Correlation Coefficient providing a balanced measure of classification quality even with imbalanced datasets. Values range from -1 to +1, where +1 indicates perfect prediction, 0 indicates random prediction, and -1 indicates completely wrong prediction. This metric is particularly reliable for medical classification tasks."
              isReference={true}
            />
          )}
          <MetricCard
            title="Accuracy"
            value={MODEL_VALIDATION_METRICS.accuracy.toFixed(3)}
            unit=""
            description="Proportion of correct predictions out of total predictions made during validation. While accuracy provides an overall performance measure, it can be misleading with imbalanced datasets. Our model achieves 57.1% accuracy, which is meaningful when combined with other metrics like AUC-ROC and F1-score."
            isReference={true}
          />
          <MetricCard
            title="R² Score"
            value={MODEL_VALIDATION_METRICS.r2_score.toFixed(3)}
            unit=""
            description="Coefficient of determination measuring how well the model explains the variance in cancer risk predictions. Values range from 0-1, where 1.0 indicates perfect prediction of risk variance. Our R² of 0.427 means the model explains 42.7% of the variance in cancer risk, which is substantial for complex genetic risk prediction."
            isReference={true}
          />
          <MetricCard
            title="Mean Squared Error"
            value={MODEL_VALIDATION_METRICS.mse.toFixed(5)}
            unit=""
            description="Average squared difference between predicted and actual risk values, measuring prediction accuracy. Lower values indicate better performance. MSE is sensitive to outliers and provides a comprehensive view of prediction errors. This metric helps assess the precision of continuous risk score predictions."
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