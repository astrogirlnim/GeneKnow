import React from 'react';
import { useNavigate, useSearchParams, useLocation } from 'react-router-dom';
import Layout from '../components/Layout';
import type { PipelineResult } from '../api/geneknowPipeline';

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
const HazardScoreCard = ({ value, tooltipContent }: { value: string; tooltipContent: { content: string; link?: string } }) => {
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
    </div>
  );
};

// Dynamic font sizing for text that needs to fit
const getDynamicFontSize = (text: string, maxWidth: number = 200) => {
  const baseSize = 30; // Base font size in pixels
  const minSize = 12; // Minimum readable size
  const charWidth = 0.6; // Approximate character width ratio
  
  const textWidth = text.length * baseSize * charWidth;
  if (textWidth <= maxWidth) return baseSize;
  
  const scaledSize = Math.floor(maxWidth / (text.length * charWidth));
  return Math.max(scaledSize, minSize);
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

// Filename Card with dynamic font sizing
const FilenameCard = ({ filename, tooltipContent }: { 
  filename: string; 
  tooltipContent: { content: string; link?: string } 
}) => {
  const fontSize = React.useMemo(() => {
    // Calculate font size based on filename length
    if (filename.length <= 15) return '1.875rem';
    if (filename.length <= 25) return '1.5rem';
    if (filename.length <= 35) return '1.25rem';
    if (filename.length <= 45) return '1rem';
    return '0.875rem'; // Minimum size for very long filenames
  }, [filename]);

  return (
    <div style={{
      background: '#FFFFFF',
      padding: '1rem',
      borderRadius: '0.75rem',
      boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
      border: '1px solid #E5E7EB',
      minWidth: 0,
    }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
        <h4 style={{ fontWeight: '600', color: '#4B5563' }}>File Analyzed</h4>
        <div style={{ position: 'relative' }} className="group">
          <InformationCircleIcon className="w-5 h-5 cursor-pointer text-gray-400" />
          <Tooltip content={tooltipContent.content} link={tooltipContent.link} />
        </div>
      </div>
      <p style={{ 
        fontSize,
        fontWeight: 'bold', 
        color: '#111827', 
        marginTop: '0.5rem',
        wordBreak: 'break-all',
        lineHeight: '1.2'
      }}>
        {filename}
      </p>
    </div>
  );
};

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
  
  // Check if we have real results from the pipeline
  const pipelineResults = location.state?.results as PipelineResult | undefined;
  const fileName = location.state?.fileName as string | undefined;
  const mockRiskLevel = location.state?.mockRiskLevel as string | undefined;
  
  // Get the risk level from URL parameters or state, default to 'low' if not specified
  const riskLevel = searchParams.get('risk') || mockRiskLevel || 'low';
  const mockData = mockDataSets[riskLevel as keyof typeof mockDataSets] || mockDataSets.low;
  
  // Use real data if available, otherwise use mock data
  const displayData: DisplayData = React.useMemo(() => {
    if (pipelineResults) {
      // Extract the highest risk score for display
      const riskScores = Object.entries(pipelineResults.risk_scores || {});
      const highestRisk = riskScores.reduce((prev, [cancer, score]) => 
        score > prev.score ? { cancer, score } : prev,
        { cancer: '', score: 0 }
      );
      
      // Use pipeline results directly (already in percentage format)
      const probability = parseFloat(highestRisk.score.toFixed(1));
      const hazardScore = highestRisk.score ? (highestRisk.score / 100 * 3).toFixed(1) : "0.0"; // Convert from percentage to hazard score scale
      
      // Extract key metrics from the pipeline results
      const keyVariants = pipelineResults.variants?.slice(0, 3) || [];
      const variantCount = pipelineResults.variant_count || 0;
      
      // Sort cancer risks by score and filter out baseline risks
      const BASELINE_THRESHOLD = 5.0; // Only show risks above 5%
      const topCancerRisks = riskScores
        .filter(([_, score]) => score > BASELINE_THRESHOLD)
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
      
      // Add variant details if available
      if (keyVariants.length > 0) {
        keyVariants.forEach((variant, index) => {
          otherMetrics.push({
            id: 4 + index,
            title: `Variant ${index + 1}`,
            value: variant.gene,
            unit: variant.type,
            tooltipContent: {
              content: `${variant.impact} impact variant at position ${variant.position}`,
              link: "#"
            }
          });
        });
      }
      
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
                  Analysis completed for: {fileName}
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
                Clinical View
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

          {/* Report Sections from Pipeline Results */}
          {pipelineResults && pipelineResults.report_sections && (
            <div style={{
              marginBottom: '2rem',
              padding: '1.5rem',
              background: '#FFFFFF',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB'
            }}>
              <h3 style={{ fontWeight: '600', fontSize: '1.25rem', marginBottom: '1rem', color: '#111827' }}>
                Analysis Report
              </h3>
              {Object.entries(pipelineResults.report_sections || {}).map(([key, section]) => (
                <div key={key} style={{ marginBottom: '1rem' }}>
                  <h4 style={{ 
                    fontWeight: '600', 
                    fontSize: '1rem',
                    color: section.severity === 'high' ? '#DC2626' : 
                           section.severity === 'medium' ? '#F59E0B' : '#059669',
                    marginBottom: '0.5rem'
                  }}>
                    {section.title}
                  </h4>
                  <p style={{ color: '#4B5563', lineHeight: '1.5' }}>{section.content}</p>
                  {section.technical_details && (
                    <details style={{ marginTop: '0.5rem' }}>
                      <summary style={{ cursor: 'pointer', color: '#6B7280', fontSize: '0.875rem' }}>
                        Technical Details
                      </summary>
                      <p style={{ marginTop: '0.5rem', fontSize: '0.875rem', color: '#6B7280' }}>
                        {section.technical_details}
                      </p>
                    </details>
                  )}
                </div>
              ))}
            </div>
          )}

          {/* Comprehensive Metrics Section */}
          {pipelineResults && pipelineResults.structured_json && pipelineResults.structured_json.metrics && (
            <div style={{
              marginBottom: '2rem',
              padding: '1.5rem',
              background: '#FFFFFF',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB'
            }}>
              <h3 style={{ fontWeight: '600', fontSize: '1.25rem', marginBottom: '1rem' }}>
                📊 Comprehensive Metrics
              </h3>
              
              {/* Metrics Summary */}
              {pipelineResults.structured_json.metrics_summary && (
                <div style={{ 
                  marginBottom: '1.5rem', 
                  padding: '1rem',
                  background: '#F9FAFB',
                  borderRadius: '0.5rem',
                  border: '1px solid #E5E7EB'
                }}>
                  <h4 style={{ fontWeight: '600', fontSize: '1rem', marginBottom: '0.5rem' }}>
                    Summary
                  </h4>
                  <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '1rem' }}>
                    <div>
                      <span style={{ fontWeight: '500' }}>Highest Risk Cancer: </span>
                      <span style={{ color: '#DC2626', fontWeight: '600' }}>
                        {pipelineResults.structured_json.metrics_summary.highest_risk_cancer}
                      </span>
                    </div>
                    <div>
                      <span style={{ fontWeight: '500' }}>Risk Score: </span>
                      <span style={{ color: '#DC2626', fontWeight: '600' }}>
                        {pipelineResults.structured_json.metrics_summary.highest_risk_score}%
                      </span>
                    </div>
                    <div>
                      <span style={{ fontWeight: '500' }}>Pathogenic Variants: </span>
                      <span style={{ color: '#EF4444', fontWeight: '600' }}>
                        {pipelineResults.structured_json.metrics_summary.pathogenic_variant_count}
                      </span>
                    </div>
                    <div>
                      <span style={{ fontWeight: '500' }}>Confidence Level: </span>
                      <span style={{ 
                        color: pipelineResults.structured_json.metrics_summary.confidence_level === 'high' ? '#22C55E' : 
                               pipelineResults.structured_json.metrics_summary.confidence_level === 'medium' ? '#F59E0B' : '#EF4444',
                        fontWeight: '600'
                      }}>
                        {pipelineResults.structured_json.metrics_summary.confidence_level}
                      </span>
                    </div>
                  </div>
                </div>
              )}

              {/* Detailed Metrics */}
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', gap: '1rem' }}>
                
                {/* Variant Metrics */}
                <div style={{ 
                  padding: '1rem',
                  background: '#F0F9FF',
                  borderRadius: '0.5rem',
                  border: '1px solid #BAE6FD'
                }}>
                  <h4 style={{ fontWeight: '600', fontSize: '1rem', marginBottom: '0.5rem', color: '#0369A1' }}>
                    🧬 Variant Analysis
                  </h4>
                  <div style={{ fontSize: '0.875rem', lineHeight: '1.4' }}>
                    <div>Total Variants: <strong>{pipelineResults.structured_json.metrics.variant_metrics.total_variants}</strong></div>
                    <div>Pathogenic: <strong style={{ color: '#DC2626' }}>{pipelineResults.structured_json.metrics.variant_metrics.pathogenic_variants}</strong></div>
                    <div>Uncertain: <strong style={{ color: '#F59E0B' }}>{pipelineResults.structured_json.metrics.variant_metrics.uncertain_variants}</strong></div>
                    <div>High CADD Score: <strong>{pipelineResults.structured_json.metrics.variant_metrics.high_cadd_variants}</strong></div>
                    <div>Genes Affected: <strong>{pipelineResults.structured_json.metrics.variant_metrics.genes_affected}</strong></div>
                    {pipelineResults.structured_json.metrics.variant_metrics.mean_cadd_score > 0 && (
                      <div>Mean CADD Score: <strong>{pipelineResults.structured_json.metrics.variant_metrics.mean_cadd_score.toFixed(2)}</strong></div>
                    )}
                  </div>
                </div>

                {/* Confidence Metrics */}
                <div style={{ 
                  padding: '1rem',
                  background: '#F0FDF4',
                  borderRadius: '0.5rem',
                  border: '1px solid #BBF7D0'
                }}>
                  <h4 style={{ fontWeight: '600', fontSize: '1rem', marginBottom: '0.5rem', color: '#166534' }}>
                    🎯 Model Confidence
                  </h4>
                  <div style={{ fontSize: '0.875rem', lineHeight: '1.4' }}>
                    <div>Model Confidence: <strong>{(pipelineResults.structured_json.metrics.confidence_metrics.mean_model_confidence * 100).toFixed(1)}%</strong></div>
                    <div>ML Fusion Confidence: <strong>{(pipelineResults.structured_json.metrics.confidence_metrics.ml_fusion_confidence * 100).toFixed(1)}%</strong></div>
                    <div>Risk Score Range: <strong>{pipelineResults.structured_json.metrics.confidence_metrics.risk_score_range.toFixed(1)}%</strong></div>
                  </div>
                </div>

                {/* PRS Metrics */}
                <div style={{ 
                  padding: '1rem',
                  background: '#FEF3C7',
                  borderRadius: '0.5rem',
                  border: '1px solid #FCD34D'
                }}>
                  <h4 style={{ fontWeight: '600', fontSize: '1rem', marginBottom: '0.5rem', color: '#92400E' }}>
                    🧮 Polygenic Risk Score
                  </h4>
                  <div style={{ fontSize: '0.875rem', lineHeight: '1.4' }}>
                    <div>Max PRS Percentile: <strong>{pipelineResults.structured_json.metrics.prs_metrics.max_prs_percentile}%</strong></div>
                    <div>Mean PRS Score: <strong>{pipelineResults.structured_json.metrics.prs_metrics.mean_prs_score.toFixed(3)}</strong></div>
                    <div>Overall Confidence: <strong>{pipelineResults.structured_json.metrics.prs_metrics.overall_confidence}</strong></div>
                  </div>
                </div>

                {/* Pathway Metrics */}
                <div style={{ 
                  padding: '1rem',
                  background: '#FDF4FF',
                  borderRadius: '0.5rem',
                  border: '1px solid #E9D5FF'
                }}>
                  <h4 style={{ fontWeight: '600', fontSize: '1rem', marginBottom: '0.5rem', color: '#7C3AED' }}>
                    🛤️ Pathway Analysis
                  </h4>
                  <div style={{ fontSize: '0.875rem', lineHeight: '1.4' }}>
                    <div>High Burden Pathways: <strong>{pipelineResults.structured_json.metrics.pathway_metrics.high_burden_pathway_count}</strong></div>
                    <div>Mean Pathway Burden: <strong>{pipelineResults.structured_json.metrics.pathway_metrics.mean_pathway_burden.toFixed(3)}</strong></div>
                    <div>Risk Level: <strong>{pipelineResults.structured_json.metrics.pathway_metrics.pathway_risk_level}</strong></div>
                  </div>
                </div>

                {/* Overall Assessment */}
                <div style={{ 
                  padding: '1rem',
                  background: '#FEF2F2',
                  borderRadius: '0.5rem',
                  border: '1px solid #FECACA'
                }}>
                  <h4 style={{ fontWeight: '600', fontSize: '1rem', marginBottom: '0.5rem', color: '#DC2626' }}>
                    ⚡ Overall Assessment
                  </h4>
                  <div style={{ fontSize: '0.875rem', lineHeight: '1.4' }}>
                    <div>High-Risk Cancers: <strong>{pipelineResults.structured_json.metrics.overall_assessment.high_risk_cancers.join(', ')}</strong></div>
                    <div>Max Risk Score: <strong>{pipelineResults.structured_json.metrics.overall_assessment.max_risk_score.toFixed(1)}%</strong></div>
                    <div>Clinical Action Needed: <strong style={{ color: pipelineResults.structured_json.metrics.overall_assessment.clinical_action_needed ? '#DC2626' : '#22C55E' }}>
                      {pipelineResults.structured_json.metrics.overall_assessment.clinical_action_needed ? 'Yes' : 'No'}
                    </strong></div>
                    <div>Risk Category: <strong>{pipelineResults.structured_json.metrics.overall_assessment.risk_category}</strong></div>
                  </div>
                </div>

                {/* Performance Indicators */}
                <div style={{ 
                  padding: '1rem',
                  background: '#F8FAFC',
                  borderRadius: '0.5rem',
                  border: '1px solid #E2E8F0'
                }}>
                  <h4 style={{ fontWeight: '600', fontSize: '1rem', marginBottom: '0.5rem', color: '#475569' }}>
                    ⚙️ Performance Indicators
                  </h4>
                  <div style={{ fontSize: '0.875rem', lineHeight: '1.4' }}>
                    <div>Variant Coverage: <strong style={{ color: pipelineResults.structured_json.metrics.performance_indicators.variant_coverage ? '#22C55E' : '#EF4444' }}>
                      {pipelineResults.structured_json.metrics.performance_indicators.variant_coverage ? 'Adequate' : 'Limited'}
                    </strong></div>
                    <div>Model Confidence: <strong style={{ color: pipelineResults.structured_json.metrics.performance_indicators.model_confidence_adequate ? '#22C55E' : '#EF4444' }}>
                      {pipelineResults.structured_json.metrics.performance_indicators.model_confidence_adequate ? 'Adequate' : 'Limited'}
                    </strong></div>
                    <div>Evidence Sufficiency: <strong style={{ color: pipelineResults.structured_json.metrics.performance_indicators.sufficient_evidence ? '#22C55E' : '#EF4444' }}>
                      {pipelineResults.structured_json.metrics.performance_indicators.sufficient_evidence ? 'Sufficient' : 'Insufficient'}
                    </strong></div>
                  </div>
                </div>
              </div>

              {/* Validation Metrics (if available) */}
              {pipelineResults.structured_json.metrics.validation_metrics.ground_truth_available && (
                <div style={{ 
                  marginTop: '1rem',
                  padding: '1rem',
                  background: '#F3F4F6',
                  borderRadius: '0.5rem',
                  border: '1px solid #D1D5DB'
                }}>
                  <h4 style={{ fontWeight: '600', fontSize: '1rem', marginBottom: '0.5rem', color: '#374151' }}>
                    🔬 Validation Metrics
                  </h4>
                  <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '1rem', fontSize: '0.875rem' }}>
                    {pipelineResults.structured_json.metrics.validation_metrics.auc_roc && (
                      <div>AUC-ROC: <strong>{pipelineResults.structured_json.metrics.validation_metrics.auc_roc.toFixed(3)}</strong></div>
                    )}
                    {pipelineResults.structured_json.metrics.validation_metrics.sensitivity && (
                      <div>Sensitivity: <strong>{pipelineResults.structured_json.metrics.validation_metrics.sensitivity.toFixed(3)}</strong></div>
                    )}
                    {pipelineResults.structured_json.metrics.validation_metrics.specificity && (
                      <div>Specificity: <strong>{pipelineResults.structured_json.metrics.validation_metrics.specificity.toFixed(3)}</strong></div>
                    )}
                    {pipelineResults.structured_json.metrics.validation_metrics.f1_score && (
                      <div>F1 Score: <strong>{pipelineResults.structured_json.metrics.validation_metrics.f1_score.toFixed(3)}</strong></div>
                    )}
                  </div>
                </div>
              )}
            </div>
          )}

          {/* Headline Metrics */}
          <div style={{
            display: 'grid',
            gridTemplateColumns: fileName ? 'repeat(auto-fit, minmax(320px, 1fr))' : 'repeat(auto-fit, minmax(350px, 1fr))',
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
            />
            {fileName && (
              <FilenameCard 
                filename={fileName} 
                tooltipContent={{ content: `Full filename: ${fileName}`, link: "#" }} 
              />
            )}
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
                    {displayData.topCancerRisks.map(([cancer, score]: [string, number], index: number) => {
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
        </div>
      </section>
    </Layout>
  );
};

export default DashboardPage; 