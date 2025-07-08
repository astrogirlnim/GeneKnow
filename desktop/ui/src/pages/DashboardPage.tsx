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
        <h3 style={{ fontWeight: '600', fontSize: '1.125rem' }}>Risk Probability</h3>
        <div style={{ position: 'relative' }} className="group">
          <InformationCircleIcon className="w-5 h-5 cursor-pointer text-gray-500" />
          <Tooltip content={tooltipContent.content} link={tooltipContent.link} />
        </div>
      </div>
      <p style={{ fontSize: '3rem', fontWeight: 'bold', marginTop: '0.5rem' }}>{value}%</p>
      <p style={{ fontSize: '0.875rem', marginTop: '0.25rem' }}>Probability of having this cancer type.</p>
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
    border: '1px solid #E5E7EB'
  }}>
    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
      <h4 style={{ fontWeight: '600', color: '#4B5563' }}>{title}</h4>
      <div style={{ position: 'relative' }} className="group">
        <InformationCircleIcon className="w-5 h-5 cursor-pointer text-gray-400" />
        <Tooltip content={tooltipContent.content} link={tooltipContent.link} />
      </div>
    </div>
    <p style={{ fontSize: '1.875rem', fontWeight: 'bold', color: '#111827', marginTop: '0.5rem' }}>
      {value} <span style={{ fontSize: '1.25rem', fontWeight: '500', color: '#6B7280' }}>{unit}</span>
    </p>
  </div>
);

// Mock data sets for different risk levels
const mockDataSets = {
  high: {
    probability: 82,
    hazardScore: 2.4,
    patient: { name: 'Emma Rodriguez', condition: 'Hereditary Breast and Ovarian Cancer Syndrome' },
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
    hazardScore: 1.8,
    patient: { name: 'David Kim', condition: 'Lynch Syndrome' },
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
    hazardScore: 0.9,
    patient: { name: 'Sarah Johnson', condition: 'Li-Fraumeni Syndrome' },
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
    content: "This score represents the likelihood of a specific cancer type based on the genetic markers found.", 
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
  const displayData = React.useMemo(() => {
    if (pipelineResults) {
      // Extract the highest risk score for display
      const riskScores = Object.entries(pipelineResults.risk_scores);
      const highestRisk = riskScores.reduce((prev, [cancer, score]) => 
        score > prev.score ? { cancer, score } : prev,
        { cancer: '', score: 0 }
      );
      
      // Convert pipeline results to display format
      const probability = Math.round(highestRisk.score * 100);
      const hazardScore = (highestRisk.score * 3).toFixed(1); // Convert to hazard score scale
      
      // Extract key metrics from the pipeline results
      const keyVariants = pipelineResults.variants?.slice(0, 3) || [];
      const variantCount = pipelineResults.variant_count || 0;
      
      const otherMetrics: MetricData[] = [
        {
          id: 1,
          title: "Highest Risk Cancer",
          value: highestRisk.cancer || "None detected",
          unit: "",
          tooltipContent: {
            content: "The cancer type with the highest predicted risk based on genetic analysis.",
            link: "#"
          }
        },
        {
          id: 2,
          title: "Total Variants",
          value: variantCount,
          unit: "found",
          tooltipContent: {
            content: "Total number of genetic variants identified in the sample.",
            link: "#"
          }
        },
        {
          id: 3,
          title: "Processing Time",
          value: pipelineResults.processing_time_seconds.toFixed(1),
          unit: "seconds",
          tooltipContent: {
            content: "Time taken to process and analyze the genomic data.",
            link: "#"
          }
        },
        {
          id: 4,
          title: "File Analyzed",
          value: fileName || "Unknown",
          unit: "",
          tooltipContent: {
            content: "The genomic file that was analyzed.",
            link: "#"
          }
        }
      ];
      
      // Add variant details if available
      if (keyVariants.length > 0) {
        keyVariants.forEach((variant, index) => {
          otherMetrics.push({
            id: 5 + index,
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
      
      // Add risk scores for other cancer types
      let metricId = 8;
      riskScores.forEach(([cancer, score]) => {
        if (cancer !== highestRisk.cancer) {
          otherMetrics.push({
            id: metricId++,
            title: `${cancer} Risk`,
            value: (score * 100).toFixed(1),
            unit: "%",
            tooltipContent: {
              content: `Predicted risk score for ${cancer}.`,
              link: "#"
            }
          });
        }
      });
      
      return {
        probability,
        hazardScore,
        patient: { 
          name: fileName || 'Analyzed Sample', 
          condition: `${highestRisk.cancer} - ${probability}% risk` 
        },
        otherMetrics
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
                  color: '#6B7280',
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
              <h3 style={{ fontWeight: '600', fontSize: '1.25rem', marginBottom: '1rem' }}>
                Analysis Report
              </h3>
              {Object.entries(pipelineResults.report_sections).map(([key, section]) => (
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

          {/* Headline Metrics */}
          <div style={{
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))',
            gap: '1.5rem',
            marginBottom: '2rem'
          }}>
            <ProbabilityCard 
              value={displayData.probability} 
              tooltipContent={baseTooltips.probability} 
            />
            <MetricCard 
              title="Hazard Score" 
              value={displayData.hazardScore} 
              unit="" 
              tooltipContent={baseTooltips.hazardScore} 
            />
          </div>

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