import React from 'react';
import { useLocation, useNavigate, useSearchParams } from 'react-router-dom';
import Layout from '../components/Layout';
import ConfidenceCheck from '../components/ConfidenceCheck';
import MarkdownRenderer from '../components/MarkdownRenderer';
import AnalysisDisclaimer from '../components/AnalysisDisclaimer';
import ModelPerformanceTab from '../components/ModelPerformanceTab';
import type { PipelineResult, SHAPValidation } from '../api/geneknowPipeline';
import { invoke } from '@tauri-apps/api/core';

// Type definitions
interface MetricData {
  id: number;
  title: string;
  value: string | number;
  unit: string;
  tooltipContent: {
    content: string;
    link?: string;
  };
}

interface DisplayData {
  probability: number;
  hazardScore: string;
  patient?: {
    name: string;
    condition: string;
  };
  riskLevel?: string;
  condition?: string;
  otherMetrics: MetricData[];
  topCancerRisks: Array<[string, number]>;
}

// Icon components
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

// Risk Probability Card
const ProbabilityCard = ({ value, tooltipContent }: { value: number; tooltipContent: { content: string; link?: string } }) => {
  const getProbabilityColor = (prob: number) => {
    if (prob >= 75) return { bg: '#FEF2F2', border: '#FECACA', text: '#991B1B' };
    if (prob >= 20) return { bg: '#FFFBEB', border: '#FDE68A', text: '#92400E' };
    return { bg: '#F0FDF4', border: '#BBF7D0', text: '#166534' };
  };

  const colors = getProbabilityColor(value);

  return (
    <div style={{
      padding: '1.5rem',
      borderRadius: '0.75rem',
      border: `1px solid ${colors.border}`,
      backgroundColor: colors.bg,
      color: colors.text
    }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <h3 style={{ fontWeight: '600', fontSize: '1.125rem' }}>Highest Cancer Risk</h3>
        <div style={{ position: 'relative' }} className="group">
          <InformationCircleIcon className="w-5 h-5 cursor-pointer text-gray-500" />
          <Tooltip content={tooltipContent.content} link={tooltipContent.link} />
        </div>
      </div>
      <p style={{ fontSize: '3rem', fontWeight: 'bold', marginTop: '0.5rem' }}>{value.toFixed(1)}%</p>
      <p style={{ fontSize: '0.875rem', marginTop: '0.25rem' }}>Highest individual cancer type risk detected.</p>
    </div>
  );
};

// Hazard Score Card - Enhanced visual presentation
const HazardScoreCard = ({ 
  value, 
  tooltipContent, 
  shapValidation
}: { 
  value: string; 
  tooltipContent: { content: string; link?: string };
  shapValidation?: SHAPValidation | null;
}) => {
  const numValue = parseFloat(value);
  const getHazardColor = (score: number) => {
    if (score >= 2.0) return { bg: '#FEF2F2', border: '#FECACA', text: '#991B1B', level: 'High Risk' };
    if (score >= 1.0) return { bg: '#FFFBEB', border: '#FDE68A', text: '#92400E', level: 'Moderate Risk' };
    return { bg: '#F0FDF4', border: '#BBF7D0', text: '#166534', level: 'Low Risk' };
  };

  const colors = getHazardColor(numValue);

  return (
    <div style={{
      padding: '1.5rem',
      borderRadius: '0.75rem',
      border: `1px solid ${colors.border}`,
      backgroundColor: colors.bg,
      color: colors.text
    }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
        <h3 style={{ fontWeight: '600', fontSize: '1.125rem' }}>Hazard Score</h3>
        <div style={{ position: 'relative' }} className="group">
          <InformationCircleIcon className="w-5 h-5 cursor-pointer text-gray-500" />
          <Tooltip content={tooltipContent.content} link={tooltipContent.link} />
        </div>
      </div>
      <div style={{ display: 'flex', alignItems: 'baseline', gap: '0.5rem', marginTop: '0.5rem' }}>
        <p style={{ fontSize: '3rem', fontWeight: 'bold' }}>{value}</p>
        <p style={{ fontSize: '1.5rem', fontWeight: '500' }}>×</p>
      </div>
      <p style={{ fontSize: '0.875rem', marginTop: '0.25rem' }}>{colors.level} compared to baseline</p>
      <div style={{
        marginTop: '1rem',
        height: '8px',
        backgroundColor: colors.bg,
        borderRadius: '4px',
        overflow: 'hidden',
        border: `1px solid ${colors.border}`
      }}>
        <div style={{
          height: '100%',
          width: `${Math.min(numValue * 25, 100)}%`,
          backgroundColor: colors.text,
          transition: 'width 500ms ease'
        }} />
      </div>
      
      {/* Confidence Check Summary */}
      <ConfidenceCheck 
        validation={shapValidation ?? null} 
        isDetailed={false}
      />
    </div>
  );
};

// Metric Card
const MetricCard = ({ title, value, unit, tooltipContent }: { 
  title: string; 
  value: string | number; 
  unit: string; 
  tooltipContent: { content: string; link?: string } 
}) => (
  <div style={{
    background: '#FFFFFF',
    padding: '1rem',
    borderRadius: '0.75rem',
    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
    border: '1px solid #E5E7EB',
    minWidth: 0, // Enable text truncation
  }}>
    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
      <h4 style={{ fontWeight: '600', color: '#4B5563' }}>{title}</h4>
      <div style={{ position: 'relative' }} className="group">
        <InformationCircleIcon className="w-5 h-5 cursor-pointer text-gray-400" />
        <Tooltip content={tooltipContent.content} link={tooltipContent.link} />
      </div>
    </div>
    <p style={{ 
      fontSize: '1.875rem', 
      fontWeight: 'bold', 
      color: '#111827', 
      marginTop: '0.5rem',
      overflow: 'hidden',
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap'
    }}>
      {value} <span style={{ fontSize: '1.25rem', fontWeight: '500', color: '#6B7280' }}>{unit}</span>
    </p>
  </div>
);

// Mock data functions removed - only real pipeline data is used for confidence checks

// Mock data sets for different risk levels - completely anonymous
const mockDataSets = {
  high: {
    probability: 82,
    hazardScore: "2.4",
    riskLevel: 'High Risk',
    condition: 'Hereditary Breast and Ovarian Cancer Syndrome',
    topCancerRisks: [['breast', 82], ['ovarian', 45]] as Array<[string, number]>,
    otherMetrics: [
      { 
        id: 1, 
        title: "Key Gene Variant", 
        value: "BRCA1", 
        unit: "(c.5266dupC)", 
        tooltipContent: { 
          content: "The most significant gene variant identified in the analysis.", 
          link: "#" 
        } 
      },
      { 
        id: 2, 
        title: "Somatic Mutations", 
        value: 28, 
        unit: "found", 
        tooltipContent: { 
          content: "Number of cancer-related somatic mutations detected.", 
          link: "#" 
        } 
      },
      { 
        id: 3, 
        title: "Tumor Mutational Burden", 
        value: 15.8, 
        unit: "muts/Mb", 
        tooltipContent: { 
          content: "A measure of the total number of mutations per megabase of DNA.", 
          link: "#" 
        } 
      },
      { 
        id: 4, 
        title: "Variant Allele Freq.", 
        value: "67.3", 
        unit: "%", 
        tooltipContent: { 
          content: "The percentage of sequence reads that match a specific variant.", 
          link: "#" 
        } 
      },
      { 
        id: 5, 
        title: "Ploidy", 
        value: 2.8, 
        unit: "", 
        tooltipContent: { 
          content: "The average number of chromosome sets in a cell.", 
          link: "#" 
        } 
      },
      { 
        id: 6, 
        title: "Loss of Heterozygosity", 
        value: "32%", 
        unit: "of genome", 
        tooltipContent: { 
          content: "The percentage of the genome where one parental copy of a gene is lost.", 
          link: "#" 
        } 
      },
      { 
        id: 7, 
        title: "Clonal Hematopoiesis", 
        value: "Detected", 
        unit: "", 
        tooltipContent: { 
          content: "Presence of somatic mutations in blood cells.", 
          link: "#" 
        } 
      },
      { 
        id: 8, 
        title: "Data Quality Score", 
        value: 97.2, 
        unit: "Q-Score", 
        tooltipContent: { 
          content: "An overall score representing the quality of the input sequencing data.", 
          link: "#" 
        } 
      },
    ]
  },
  medium: {
    probability: 45,
    hazardScore: "1.8",
    riskLevel: 'Medium Risk',
    condition: 'Lynch Syndrome',
    topCancerRisks: [['colon', 45], ['endometrial', 28]] as Array<[string, number]>,
    otherMetrics: [
      { 
        id: 1, 
        title: "Key Gene Variant", 
        value: "MLH1", 
        unit: "(c.1558G>A)", 
        tooltipContent: { 
          content: "The most significant gene variant identified in the analysis.", 
          link: "#" 
        } 
      },
      { 
        id: 2, 
        title: "Somatic Mutations", 
        value: 16, 
        unit: "found", 
        tooltipContent: { 
          content: "Number of cancer-related somatic mutations detected.", 
          link: "#" 
        } 
      },
      { 
        id: 3, 
        title: "Tumor Mutational Burden", 
        value: 11.2, 
        unit: "muts/Mb", 
        tooltipContent: { 
          content: "A measure of the total number of mutations per megabase of DNA.", 
          link: "#" 
        } 
      },
      { 
        id: 4, 
        title: "Variant Allele Freq.", 
        value: "52.7", 
        unit: "%", 
        tooltipContent: { 
          content: "The percentage of sequence reads that match a specific variant.", 
          link: "#" 
        } 
      },
      { 
        id: 5, 
        title: "Ploidy", 
        value: 2.3, 
        unit: "", 
        tooltipContent: { 
          content: "The average number of chromosome sets in a cell.", 
          link: "#" 
        } 
      },
      { 
        id: 6, 
        title: "Loss of Heterozygosity", 
        value: "21%", 
        unit: "of genome", 
        tooltipContent: { 
          content: "The percentage of the genome where one parental copy of a gene is lost.", 
          link: "#" 
        } 
      },
      { 
        id: 7, 
        title: "Clonal Hematopoiesis", 
        value: "Not Detected", 
        unit: "", 
        tooltipContent: { 
          content: "Presence of somatic mutations in blood cells.", 
          link: "#" 
        } 
      },
      { 
        id: 8, 
        title: "Data Quality Score", 
        value: 98.9, 
        unit: "Q-Score", 
        tooltipContent: { 
          content: "An overall score representing the quality of the input sequencing data.", 
          link: "#" 
        } 
      },
    ]
  },
  low: {
    probability: 15,
    hazardScore: "0.9",
    riskLevel: 'Low Risk',
    condition: 'Li-Fraumeni Syndrome',
    topCancerRisks: [] as Array<[string, number]>,
    otherMetrics: [
      { 
        id: 1, 
        title: "Key Gene Variant", 
        value: "TP53", 
        unit: "(c.743G>A)", 
        tooltipContent: { 
          content: "The most significant gene variant identified in the analysis.", 
          link: "#" 
        } 
      },
      { 
        id: 2, 
        title: "Somatic Mutations", 
        value: 8, 
        unit: "found", 
        tooltipContent: { 
          content: "Number of cancer-related somatic mutations detected.", 
          link: "#" 
        } 
      },
      { 
        id: 3, 
        title: "Tumor Mutational Burden", 
        value: 5.3, 
        unit: "muts/Mb", 
        tooltipContent: { 
          content: "A measure of the total number of mutations per megabase of DNA.", 
          link: "#" 
        } 
      },
      { 
        id: 4, 
        title: "Variant Allele Freq.", 
        value: "31.4", 
        unit: "%", 
        tooltipContent: { 
          content: "The percentage of sequence reads that match a specific variant.", 
          link: "#" 
        } 
      },
      { 
        id: 5, 
        title: "Ploidy", 
        value: 2.0, 
        unit: "", 
        tooltipContent: { 
          content: "The average number of chromosome sets in a cell.", 
          link: "#" 
        } 
      },
      { 
        id: 6, 
        title: "Loss of Heterozygosity", 
        value: "8%", 
        unit: "of genome", 
        tooltipContent: { 
          content: "The percentage of the genome where one parental copy of a gene is lost.", 
          link: "#" 
        } 
      },
      { 
        id: 7, 
        title: "Clonal Hematopoiesis", 
        value: "Not Detected", 
        unit: "", 
        tooltipContent: { 
          content: "Presence of somatic mutations in blood cells.", 
          link: "#" 
        } 
      },
      { 
        id: 8, 
        title: "Data Quality Score", 
        value: 99.1, 
        unit: "Q-Score", 
        tooltipContent: { 
          content: "An overall score representing the quality of the input sequencing data.", 
          link: "#" 
        } 
      },
    ]
  },
  skipped: {
    probability: 8,
    hazardScore: "0.6",
    riskLevel: 'Low Risk',
    condition: 'Rare Variant Profile',
    topCancerRisks: [['pancreatic', 8], ['sarcoma', 5]] as Array<[string, number]>,
    otherMetrics: [
      { 
        id: 1, 
        title: "Key Gene Variant", 
        value: "CDKN2A", 
        unit: "(c.442G>T)", 
        tooltipContent: { 
          content: "The most significant gene variant identified in the analysis.", 
          link: "#" 
        } 
      },
      { 
        id: 2, 
        title: "Somatic Mutations", 
        value: 3, 
        unit: "found", 
        tooltipContent: { 
          content: "Number of cancer-related somatic mutations detected.", 
          link: "#" 
        } 
      },
      { 
        id: 3, 
        title: "Tumor Mutational Burden", 
        value: 2.8, 
        unit: "muts/Mb", 
        tooltipContent: { 
          content: "A measure of the total number of mutations per megabase of DNA.", 
          link: "#" 
        } 
      },
      { 
        id: 4, 
        title: "Variant Allele Freq.", 
        value: "28.9", 
        unit: "%", 
        tooltipContent: { 
          content: "The percentage of sequence reads that match a specific variant.", 
          link: "#" 
        } 
      },
      { 
        id: 5, 
        title: "Ploidy", 
        value: 2.1, 
        unit: "", 
        tooltipContent: { 
          content: "The average number of chromosome sets in a cell.", 
          link: "#" 
        } 
      },
      { 
        id: 6, 
        title: "Loss of Heterozygosity", 
        value: "3%", 
        unit: "of genome", 
        tooltipContent: { 
          content: "The percentage of the genome where one parental copy of a gene is lost.", 
          link: "#" 
        } 
      },
      { 
        id: 7, 
        title: "Clonal Hematopoiesis", 
        value: "Not Detected", 
        unit: "", 
        tooltipContent: { 
          content: "Presence of somatic mutations in blood cells.", 
          link: "#" 
        } 
      },
      { 
        id: 8, 
        title: "Data Quality Score", 
        value: 95.8, 
        unit: "Q-Score", 
        tooltipContent: { 
          content: "An overall score representing the quality of the input sequencing data.", 
          link: "#" 
        } 
      },
    ]
  }
};

const baseTooltips = {
  probability: { 
    content: "This shows the highest individual cancer type risk detected in your genetic analysis. Each cancer type has its own separate risk percentage.", 
    link: "#" 
  },
  hazardScore: { 
    content: "The hazard score compares your risk to a baseline population over time.", 
    link: "#" 
  },
};

const DashboardPage: React.FC = () => {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const location = useLocation();
  const [activeTab, setActiveTab] = React.useState<'dashboard' | 'performance' | 'reports'>('dashboard');
  const [markdownContent, setMarkdownContent] = React.useState<string>('');
  const [loadingMarkdown] = React.useState<boolean>(false);
  
  // Check if we have real results from the pipeline
  const pipelineResults = location.state?.results as PipelineResult | undefined;
  const fileName = location.state?.fileName as string | undefined;
  const mockRiskLevel = location.state?.mockRiskLevel as string | undefined;
  
  // Get the risk level from URL parameters or state, default to 'low' if not specified
  const riskLevel = searchParams.get('risk') || mockRiskLevel || 'low';
  const mockData = mockDataSets[riskLevel as keyof typeof mockDataSets] || mockDataSets.low;
  
  // Use real SHAP validation from pipeline results only - no mock data fallback
  const getConfidenceCheckValidation = () => {
    // Use real SHAP validation if available
    if (pipelineResults?.structured_json?.shap_validation) {
      return pipelineResults.structured_json.shap_validation;
    }
    
    // If no pipeline results or no SHAP validation, return SKIPPED status
    return {
      status: 'SKIPPED' as const,
      reasons: ['No pipeline results available'],
      top_contributors: [],
      feature_importance: {},
      details: {
        status: 'SKIPPED' as const,
        risk_score: 0.0,
        top_contributors: [],
        validation_reasons: ['No pipeline results available'],
        rule_results: {},
        shap_values: [],
        feature_names: [],
        model_type: 'None'
      }
    };
  };
  
  const confidenceCheckValidation = getConfidenceCheckValidation();

  // Use real data if available, otherwise use mock data
  const displayData: DisplayData = React.useMemo(() => {
    if (pipelineResults) {
      // DEBUG: Log the actual pipeline results structure
      console.log('🔍 DASHBOARD DEBUG: pipelineResults received:', pipelineResults);
      console.log('🔍 DASHBOARD DEBUG: pipelineResults.risk_scores:', pipelineResults.risk_scores);
      console.log('🔍 DASHBOARD DEBUG: pipelineResults.variant_count:', pipelineResults.variant_count);
      console.log('🔍 DASHBOARD DEBUG: pipelineResults keys:', Object.keys(pipelineResults));
      
      // Extract the highest risk score for display
      const riskScores = Object.entries(pipelineResults.risk_scores || {});
      console.log('🔍 DASHBOARD DEBUG: riskScores entries:', riskScores);
      
      const highestRisk = riskScores.reduce((prev, [cancer, score]) => 
        score > prev.score ? { cancer, score } : prev,
        { cancer: '', score: 0 }
      );
      
      console.log('🔍 DASHBOARD DEBUG: highestRisk calculated:', highestRisk);
      
      // Use pipeline results directly (already in percentage format)
      const probability = parseFloat(highestRisk.score.toFixed(1));
      const hazardScore = highestRisk.score ? (highestRisk.score / 100 * 3).toFixed(1) : "0.0"; // Convert from percentage to hazard score scale
      
      // Extract key metrics from the pipeline results
      const variantCount = pipelineResults.variant_count || 0;
      
      // Sort cancer risks by score and filter out baseline risks
      const BASELINE_THRESHOLD = 5.0; // Only show risks above 5%
      const topCancerRisks = riskScores
        .filter(([, score]) => score > BASELINE_THRESHOLD)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 10); // Show top 10 risks above threshold
      
      const otherMetrics: MetricData[] = [
        {
          id: 1,
          title: "Total Variants",
          value: variantCount,
          unit: "found",
          tooltipContent: {
            content: "Total number of genetic variants identified in the sample.",
            link: "#"
          }
        },
        {
          id: 2,
          title: "Processing Time",
          value: pipelineResults.processing_time_seconds ? pipelineResults.processing_time_seconds.toFixed(1) : "N/A",
          unit: pipelineResults.processing_time_seconds ? "seconds" : "",
          tooltipContent: {
            content: "Time taken to process and analyze the genomic data.",
            link: "#"
          }
        }
      ];
      
      return {
        probability,
        hazardScore,
        patient: { 
          name: fileName || 'Analyzed Sample', 
          condition: `${highestRisk.cancer} - ${probability.toFixed(1)}% risk` 
        },
        otherMetrics,
        topCancerRisks // New field for cancer risks
      };
    } else {
      // Use mock data
      return mockData;
    }
  }, [pipelineResults, fileName, mockData]);

  const handleNewAnalysis = () => {
    navigate('/upload');
  };

  const handleClinicalView = () => {
    // Navigate to clinical view with the results or risk level
    if (pipelineResults) {
      navigate('/clinical', { state: { results: pipelineResults, fileName } });
    } else {
      navigate(`/clinical?risk=${riskLevel}`);
    }
  };

  // Generate markdown content from pipeline results
  const generateMarkdownFromResults = React.useCallback((results: PipelineResult): string => {
    const currentDate = new Date().toLocaleDateString();
    const currentTime = new Date().toLocaleTimeString();
    
    let markdown = `**Generated:** ${currentDate} at ${currentTime}  
**Analysis Type:** In-Depth Genomic Analysis  
**File Analyzed:** ${fileName || 'Unknown'}

---

## Executive Summary

This genomic analysis has been completed using advanced bioinformatics pipelines. The analysis includes comprehensive variant assessment, population frequency analysis, and risk stratification based on current scientific literature.

## Risk Assessment

### Key Findings

Based on the genomic analysis, the following risk assessments have been determined:

`;

    // Add report sections if available
    if (results.report_sections) {
      Object.entries(results.report_sections).forEach(([, section]) => {
        const severityEmoji = section.severity === 'high' ? '🔴' : 
                             section.severity === 'medium' ? '🟡' : '🟢';
        
        markdown += `### ${severityEmoji} ${section.title}

${section.content}

`;
        
        if (section.technical_details) {
          markdown += `**Technical Details:** ${section.technical_details}

`;
        }
      });
    }

    // Add risk information if available
    if (displayData.topCancerRisks && displayData.topCancerRisks.length > 0) {
      markdown += `## Cancer Risk Stratification

The analysis evaluated risk across multiple cancer types using polygenic risk scores and variant pathogenicity assessment.

| Cancer Type | Risk Level | Percentage |
|-------------|------------|------------|
`;
      
      displayData.topCancerRisks.forEach(([cancer, score]) => {
        const riskLevel = score >= 30 ? 'High' : score >= 15 ? 'Moderate' : 'Slightly Elevated';
        markdown += `| ${cancer.charAt(0).toUpperCase() + cancer.slice(1)} | ${riskLevel} | ${score.toFixed(1)}% |
`;
      });
      
      markdown += `
`;
    }

    // Add variant information
    markdown += `## Variant Analysis

### Significant Variants

All variants have been evaluated for:
- **Research Significance**: Using ClinVar database annotations
- **Population Frequency**: Compared against 1000 Genomes and gnomAD databases  
- **Pathogenicity Prediction**: CADD scores and ensemble predictions
- **Literature Evidence**: Cross-referenced with TCGA and published studies

`;

    // Add metrics if available
    if (displayData.otherMetrics && displayData.otherMetrics.length > 0) {
      markdown += `### Analysis Metrics

`;
      displayData.otherMetrics.forEach(metric => {
        markdown += `- **${metric.title}**: ${metric.value} ${metric.unit}
`;
      });
      markdown += `
`;
    }

    markdown += `## Recommendations

### Recommended Follow-up

1. **Routine Screening**: Continue standard screening protocols for your age group
2. **Family History**: Consider family history in risk assessment
3. **Lifestyle Factors**: Maintain healthy lifestyle choices

### Genetic Counseling

Consider genetic counseling if:
- Family history of cancer is present
- Questions about risk interpretation arise
- Additional testing is desired

## Technical Details

### Analysis Pipeline

- **Variant Calling**: Advanced genomic analysis pipeline
- **Quality Control**: Comprehensive filtering applied
- **Annotation**: Multi-database cross-reference
- **Risk Modeling**: Statistical analysis approach

### Data Sources

- **ClinVar**: Research variant interpretation
- **TCGA**: Cancer genomics reference
- **1000 Genomes**: Population frequency data
- **PRS Catalog**: Polygenic risk scores

## Glossary

**Variant**: A genetic difference from the reference genome  
**ClinVar**: Database of genetic variants and their research significance  
**CADD Score**: Combined Annotation Dependent Depletion score for variant pathogenicity  
**Polygenic Risk Score (PRS)**: Combined effect of multiple genetic variants on disease risk  

---

*This report is for informational purposes only and should not replace professional medical advice.*`;

    return markdown;
  }, [fileName, displayData]);

  // Auto-generate markdown content when we have pipeline results but no enhanced report
  React.useEffect(() => {
    if (activeTab === 'reports' && pipelineResults) {
      // Reset markdown content when switching to reports tab to ensure fresh loading
      if (!markdownContent) {
        // First priority: Use the in-memory report content if available (HIPAA compliant)
        if (pipelineResults.enhanced_report_content?.markdown) {
          console.log('Using in-memory enhanced report content');
          setMarkdownContent(pipelineResults.enhanced_report_content.markdown);
        } else {
          // Fallback: Use the old generateMarkdownFromResults if no enhanced report exists
          console.log('No enhanced report content found, using generated markdown');
          const generatedMarkdown = generateMarkdownFromResults(pipelineResults);
          setMarkdownContent(generatedMarkdown);
        }
      }
    }
  }, [activeTab, pipelineResults, markdownContent, generateMarkdownFromResults]);

  // Generate PDF from markdown content using pandoc
  const downloadPDF = React.useCallback(async () => {
    const content = markdownContent || pipelineResults?.enhanced_report_content?.markdown || '';
    if (!content) {
      console.error('No content available for PDF generation');
      return;
    }

    try {
      // Create complete markdown content with header and footer
      const currentDate = new Date().toLocaleDateString();
      const filename = `genomic-report-${fileName || 'analysis'}-${new Date().toISOString().split('T')[0]}.pdf`;
      
      const completeMarkdown = `---
title: "Genomic Analysis Report"
author: "GeneKnow Platform"
date: "${currentDate}"
geometry: margin=2cm
fontsize: 11pt
papersize: a4
documentclass: article
header-includes:
  - \\usepackage{fancyhdr}
  - \\pagestyle{fancy}
  - \\fancyfoot[C]{This report is for informational purposes only and should not replace professional medical advice.}
---

${content}`;

      // Use Tauri to convert markdown to PDF with pandoc
      try {
        const savedPath = await invoke<string>('convert_markdown_to_pdf', {
          markdownContent: completeMarkdown,
          filename: filename
        });
        
        console.log('PDF saved successfully to:', savedPath);
        console.log('PDF generated using pandoc');
      } catch (error) {
        // Fallback: If pandoc command is not available, show helpful error
        if (error && typeof error === 'string' && error.includes('pandoc')) {
          alert('Pandoc is not installed. Please install pandoc to generate PDF reports.\n\nOn macOS: brew install pandoc\nOn Windows: Download from https://pandoc.org/installing.html\nOn Linux: sudo apt-get install pandoc');
        } else if (error !== 'Save cancelled by user') {
          console.error('Failed to generate PDF with pandoc:', error);
          alert('Failed to generate PDF: ' + error);
        }
      }
    } catch (error) {
      console.error('Error preparing PDF generation:', error);
      alert('Failed to generate PDF. Please try again.');
    }
  }, [markdownContent, pipelineResults, fileName]);
  return (
    <Layout>
      <section style={{ 
        background: '#F9FAFB',
        minHeight: 'calc(100vh - 4rem)',
        padding: '3rem 0'
      }}>
        <div style={{ 
          maxWidth: '1200px',
          margin: '0 auto',
          padding: '0 1.5rem'
        }}>
          <div style={{ 
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            marginBottom: '2rem'
          }}>
            <div>
              <h2 style={{
                fontSize: '1.875rem',
                fontWeight: 'bold',
                color: '#111827'
              }}>
                Analysis Dashboard
              </h2>
              {pipelineResults && (
                <p style={{
                  fontSize: '0.875rem',
                  color: '#374151',
                  marginTop: '0.25rem'
                }}>
                  Analysis completed.
                </p>
              )}
            </div>
            <div style={{ 
              display: 'flex',
              alignItems: 'center',
              gap: '1rem'
            }}>
              <button
                onClick={handleClinicalView}
                style={{
                  padding: '0.5rem 1.25rem',
                  color: '#2563EB',
                  background: '#DBEAFE',
                  border: 'none',
                  borderRadius: '0.5rem',
                  cursor: 'pointer',
                  transition: 'all 200ms ease',
                  fontSize: '0.875rem',
                  fontWeight: '500'
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.background = '#BFDBFE';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.background = '#DBEAFE';
                }}
              >
                In-Depth Analysis
              </button>
              <button
                onClick={handleNewAnalysis}
                style={{
                  padding: '0.5rem 1.25rem',
                  color: '#FFFFFF',
                  background: '#2563EB',
                  border: 'none',
                  borderRadius: '0.5rem',
                  cursor: 'pointer',
                  transition: 'all 200ms ease',
                  fontSize: '0.875rem',
                  fontWeight: '500'
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.background = '#1D4ED8';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.background = '#2563EB';
                }}
              >
                Run New Analysis
              </button>
            </div>
          </div>



          {/* Tab Navigation */}
          <div style={{
            display: 'flex',
            borderBottom: '1px solid #E5E7EB',
            marginBottom: '2rem'
          }}>
            <button
              onClick={() => setActiveTab('dashboard')}
              style={{
                padding: '1rem 1.5rem',
                background: 'none',
                border: 'none',
                borderBottom: activeTab === 'dashboard' ? '2px solid #2563EB' : '2px solid transparent',
                color: activeTab === 'dashboard' ? '#2563EB' : '#6B7280',
                fontWeight: activeTab === 'dashboard' ? '600' : '500',
                cursor: 'pointer',
                fontSize: '1rem',
                transition: 'all 200ms ease'
              }}
              onMouseEnter={(e) => {
                if (activeTab !== 'dashboard') {
                  e.currentTarget.style.color = '#374151';
                }
              }}
              onMouseLeave={(e) => {
                if (activeTab !== 'dashboard') {
                  e.currentTarget.style.color = '#6B7280';
                }
              }}
            >
              Analysis Dashboard
            </button>
            <button
              onClick={() => setActiveTab('performance')}
              style={{
                padding: '1rem 1.5rem',
                background: 'none',
                border: 'none',
                borderBottom: activeTab === 'performance' ? '2px solid #2563EB' : '2px solid transparent',
                color: activeTab === 'performance' ? '#2563EB' : '#6B7280',
                fontWeight: activeTab === 'performance' ? '600' : '500',
                cursor: 'pointer',
                fontSize: '1rem',
                transition: 'all 200ms ease'
              }}
              onMouseEnter={(e) => {
                if (activeTab !== 'performance') {
                  e.currentTarget.style.color = '#374151';
                }
              }}
              onMouseLeave={(e) => {
                if (activeTab !== 'performance') {
                  e.currentTarget.style.color = '#6B7280';
                }
              }}
            >
              Model Performance
            </button>
            <button
              onClick={() => setActiveTab('reports')}
              style={{
                padding: '1rem 1.5rem',
                background: 'none',
                border: 'none',
                borderBottom: activeTab === 'reports' ? '2px solid #2563EB' : '2px solid transparent',
                color: activeTab === 'reports' ? '#2563EB' : '#6B7280',
                fontWeight: activeTab === 'reports' ? '600' : '500',
                cursor: 'pointer',
                fontSize: '1rem',
                transition: 'all 200ms ease'
              }}
              onMouseEnter={(e) => {
                if (activeTab !== 'reports') {
                  e.currentTarget.style.color = '#374151';
                }
              }}
              onMouseLeave={(e) => {
                if (activeTab !== 'reports') {
                  e.currentTarget.style.color = '#6B7280';
                }
              }}
            >
              Report
              {pipelineResults?.report_generator_info?.llm_enhanced && (
                <span style={{
                  marginLeft: '0.5rem',
                  padding: '0.125rem 0.375rem',
                  background: '#10B981',
                  color: '#FFFFFF',
                  fontSize: '0.75rem',
                  borderRadius: '9999px',
                  fontWeight: '500'
                }}>
                  AI
                </span>
              )}
            </button>
          </div>

          {/* Tab Content */}
          {activeTab === 'dashboard' && (
            <>
              {/* Headline Metrics */}
              <div style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(350px, 1fr))',
                gap: '1.5rem',
                marginBottom: '2rem'
              }}>
                <ProbabilityCard 
                  value={displayData.probability} 
                  tooltipContent={baseTooltips.probability} 
                />
                <HazardScoreCard 
                  value={displayData.hazardScore} 
                  tooltipContent={baseTooltips.hazardScore}
                  shapValidation={confidenceCheckValidation as unknown as SHAPValidation}
                />
              </div>



              {/* Cancer Risk Assessment - Only show if we have real pipeline results */}
              {pipelineResults && (
                <div style={{
                  marginBottom: '2rem',
                  padding: '1.5rem',
                  background: '#FFFFFF',
                  borderRadius: '0.75rem',
                  boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
                  border: '1px solid #E5E7EB'
                }}>
                  <h3 style={{ 
                    fontWeight: '600', 
                    fontSize: '1.25rem', 
                    marginBottom: '1rem',
                    color: '#111827'
                  }}>
                    Cancer Risk Assessment
                  </h3>
                  
                  {displayData.topCancerRisks && displayData.topCancerRisks.length > 0 ? (
                    <>
                      <p style={{ 
                        color: '#6B7280', 
                        marginBottom: '1rem',
                        fontSize: '0.875rem'
                      }}>
                        Showing cancer types with elevated risk (above 5% baseline)
                      </p>
                      <div style={{
                        display: 'grid',
                        gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
                        gap: '1rem'
                      }}>
                        {displayData.topCancerRisks.map(([cancer, score]: [string, number]) => {
                          const getRiskColor = (risk: number) => {
                            if (risk >= 30) return '#DC2626'; // Red for high risk
                            if (risk >= 15) return '#F59E0B'; // Amber for medium risk
                            return '#059669'; // Green for lower risk
                          };
                          
                          const getRiskLevel = (risk: number) => {
                            if (risk >= 30) return 'High';
                            if (risk >= 15) return 'Moderate';
                            return 'Slightly Elevated';
                          };
                          
                          return (
                            <div key={cancer} style={{
                              padding: '1rem',
                              borderRadius: '0.5rem',
                              background: '#F9FAFB',
                              border: '1px solid #E5E7EB',
                              display: 'flex',
                              justifyContent: 'space-between',
                              alignItems: 'center'
                            }}>
                              <div>
                                <p style={{ 
                                  fontWeight: '600',
                                  color: '#374151',
                                  textTransform: 'capitalize',
                                  marginBottom: '0.25rem'
                                }}>
                                  {cancer}
                                </p>
                                <p style={{ 
                                  fontSize: '0.75rem',
                                  color: getRiskColor(score)
                                }}>
                                  {getRiskLevel(score)} Risk
                                </p>
                              </div>
                              <p style={{ 
                                fontSize: '1.5rem',
                                fontWeight: 'bold',
                                color: getRiskColor(score)
                              }}>
                                {score.toFixed(1)}%
                              </p>
                            </div>
                          );
                        })}
                      </div>
                      {displayData.topCancerRisks.length <= 2 && (
                        <p style={{ 
                          color: '#6B7280', 
                          marginTop: '1rem',
                          fontSize: '0.875rem',
                          fontStyle: 'italic'
                        }}>
                          Other cancer types are within normal baseline risk levels.
                        </p>
                      )}
                    </>
                  ) : (
                    <div style={{
                      textAlign: 'center',
                      padding: '2rem',
                      background: '#F0FDF4',
                      borderRadius: '0.5rem',
                      border: '1px solid #BBF7D0'
                    }}>
                      <div style={{ 
                        fontSize: '3rem',
                        marginBottom: '0.5rem'
                      }}>
                        ✅
                      </div>
                      <p style={{ 
                        color: '#166534',
                        fontWeight: '600',
                        fontSize: '1.125rem'
                      }}>
                        No Elevated Cancer Risks Detected
                      </p>
                      <p style={{ 
                        color: '#15803D',
                        marginTop: '0.5rem',
                        fontSize: '0.875rem'
                      }}>
                        All cancer risk scores are within normal baseline levels (below 5%)
                      </p>
                    </div>
                  )}
                </div>
              )}

              {/* Other Metrics Grid */}
              <div style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))',
                gap: '1.5rem'
              }}>
                {displayData.otherMetrics.map((metric: MetricData) => (
                  <MetricCard 
                    key={metric.id}
                    title={metric.title} 
                    value={metric.value} 
                    unit={metric.unit} 
                    tooltipContent={metric.tooltipContent} 
                  />
                ))}
              </div>
            </>
          )}

          {/* Model Performance Tab Content */}
          {activeTab === 'performance' && (
            <ModelPerformanceTab pipelineResults={pipelineResults} />
          )}

          {/* Reports Tab Content */}
          {activeTab === 'reports' && (
            <div style={{
              padding: '1.5rem',
              background: '#FFFFFF',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB'
            }}>
              <div style={{ 
                display: 'flex', 
                justifyContent: 'space-between', 
                alignItems: 'center', 
                marginBottom: '1.5rem' 
              }}>
                <h3 style={{ 
                  fontWeight: '600', 
                  fontSize: '1.25rem', 
                  color: '#111827',
                  margin: 0
                }}>
                  In-Depth Reports
                </h3>
                
{/* Download PDF Button */}
                {(markdownContent || pipelineResults?.enhanced_report_content?.markdown) && (
                  <button
                    onClick={downloadPDF}
                    style={{
                      padding: '0.5rem 1rem',
                      background: '#2563EB',
                      color: '#FFFFFF',
                      border: 'none',
                      borderRadius: '0.375rem',
                      fontSize: '0.875rem',
                      fontWeight: '500',
                      cursor: 'pointer',
                      transition: 'all 200ms ease',
                      display: 'flex',
                      alignItems: 'center',
                      gap: '0.5rem'
                    }}
                    onMouseEnter={(e) => {
                      e.currentTarget.style.background = '#1D4ED8';
                    }}
                    onMouseLeave={(e) => {
                      e.currentTarget.style.background = '#2563EB';
                    }}
                  >
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                      <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/>
                      <polyline points="7,10 12,15 17,10"/>
                      <line x1="12" y1="15" x2="12" y2="3"/>
                    </svg>
                    Download PDF
                  </button>
                )}

              </div>

              {pipelineResults ? (
                <>
                  {/* Main Report Content */}
                  <div style={{
                    padding: '1rem',
                    background: '#FFFFFF',
                    borderRadius: '0.375rem',
                    border: '1px solid #D1D5DB',
                    maxHeight: '80vh',
                    overflowY: 'auto',
                    marginBottom: '1.5rem'
                  }}>
                    {loadingMarkdown ? (
                      <div style={{ textAlign: 'center', padding: '2rem' }}>
                        <div style={{
                          display: 'inline-block',
                          width: '2rem',
                          height: '2rem',
                          border: '3px solid #E5E7EB',
                          borderTop: '3px solid #2563EB',
                          borderRadius: '50%',
                          animation: 'spin 1s linear infinite'
                        }} />
                        <p style={{ marginTop: '1rem', color: '#6B7280' }}>Loading report...</p>
                      </div>
                    ) : markdownContent ? (
                      <MarkdownRenderer content={markdownContent} />
                    ) : pipelineResults.enhanced_report_content?.markdown ? (
                      <MarkdownRenderer content={pipelineResults.enhanced_report_content.markdown} />
                    ) : (
                      <div style={{ textAlign: 'center', padding: '2rem' }}>
                        <p style={{ color: '#6B7280', fontStyle: 'italic' }}>
                          Generating report from analysis results...
                        </p>
                      </div>
                    )}
                  </div>

                  {/* Report Generation Info - Only show if available */}
                  {pipelineResults.report_generator_info && (
                    <div style={{
                      padding: '1rem',
                      background: '#F9FAFB',
                      borderRadius: '0.5rem',
                      marginBottom: '1.5rem',
                      border: '1px solid #E5E7EB'
                    }}>
                      <h4 style={{ fontWeight: '600', marginBottom: '0.5rem', color: '#374151' }}>
                        Report Generation Details
                      </h4>
                      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '1rem' }}>
                        <div>
                          <p style={{ fontSize: '0.875rem', color: '#6B7280' }}>Backend</p>
                          <p style={{ fontWeight: '600', color: '#374151' }}>
                            {pipelineResults.report_generator_info.backend_used}
                            {pipelineResults.report_generator_info.llm_enhanced && (
                              <span style={{
                                marginLeft: '0.5rem',
                                padding: '0.125rem 0.375rem',
                                background: '#10B981',
                                color: '#FFFFFF',
                                fontSize: '0.75rem',
                                borderRadius: '9999px'
                              }}>
                                LLM Enhanced
                              </span>
                            )}
                          </p>
                        </div>
                        {pipelineResults.report_generator_info.model_used && (
                          <div>
                            <p style={{ fontSize: '0.875rem', color: '#6B7280' }}>Model</p>
                            <p style={{ fontWeight: '600', color: '#374151' }}>
                              {pipelineResults.report_generator_info.model_used}
                            </p>
                          </div>
                        )}
                        <div>
                          <p style={{ fontSize: '0.875rem', color: '#6B7280' }}>High-Risk Findings</p>
                          <p style={{ fontWeight: '600', color: '#374151' }}>
                            {pipelineResults.report_generator_info.high_risk_findings_count}
                          </p>
                        </div>
                      </div>
                    </div>
                  )}
                </>
              ) : (
                <div style={{
                  padding: '2rem',
                  textAlign: 'center',
                  background: '#F3F4F6',
                  borderRadius: '0.5rem',
                  border: '1px solid #D1D5DB'
                }}>
                  <p style={{ color: '#374151', fontWeight: '600' }}>
                    No report data available
                  </p>
                  <p style={{ color: '#6B7280', fontSize: '0.875rem', marginTop: '0.5rem' }}>
                    Reports are only available for completed pipeline analyses.
                  </p>
                </div>
              )}
            </div>
          )}
          
          {/* Analysis Disclaimer */}
          <AnalysisDisclaimer />
        </div>
      </section>
    </Layout>
  );
};

export default DashboardPage; 