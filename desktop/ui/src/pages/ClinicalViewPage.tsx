import React, { useState } from 'react';
import { useNavigate, useSearchParams } from 'react-router-dom';
import Layout from '../components/Layout';

// Type definitions
interface Alert {
  type: 'critical' | 'warning' | 'info' | 'success';
  title: string;
  desc: string;
}

// Mock data sets for different risk levels - completely anonymous
const mockDataSets = {
  high: {
    riskLevel: 'High Risk',
    riskScore: '82/100',
    condition: 'Hereditary Breast and Ovarian Cancer Syndrome',
    details: 'Family History: Breast Cancer<br/>Referral: Oncology<br/>Previous Tests: BRCA1/2 Panel',
    alerts: [
      {
        type: 'critical' as const,
        title: 'BRCA1 Pathogenic Variant',
        desc: 'c.5266dupC - Frameshift mutation detected'
      },
      {
        type: 'warning' as const,
        title: 'High Risk Classification',
        desc: 'Immediate genetic counseling recommended'
      }
    ]
  },
  medium: {
    riskLevel: 'Medium Risk',
    riskScore: '45/100',
    condition: 'Lynch Syndrome',
    details: 'Family History: Colorectal Cancer<br/>Referral: Oncology<br/>Previous Tests: MSI-H positive',
    alerts: [
      {
        type: 'warning' as const,
        title: 'MLH1 Pathogenic Variant',
        desc: 'c.1989-1G>A - Splice site mutation'
      },
      {
        type: 'info' as const,
        title: 'Screening Recommendations',
        desc: 'Enhanced colonoscopy surveillance required'
      }
    ]
  },
  low: {
    riskLevel: 'Low Risk',
    riskScore: '15/100',
    condition: 'Li-Fraumeni Syndrome',
    details: 'Family History: Multiple Sarcomas<br/>Referral: Genetics<br/>Previous Tests: TP53 Sequencing',
    alerts: [
      {
        type: 'info' as const,
        title: 'TP53 Variant of Uncertain Significance',
        desc: 'c.743G>A - Missense mutation'
      },
      {
        type: 'success' as const,
        title: 'Low Risk Assessment',
        desc: 'Routine follow-up recommended'
      }
    ]
  }
};

const ClinicalViewPage: React.FC = () => {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const [activeTab, setActiveTab] = useState('analysis');
  
  // Get the risk level from URL parameters, default to 'low' if not specified
  const riskLevel = searchParams.get('risk') || 'low';
  const currentData = mockDataSets[riskLevel as keyof typeof mockDataSets] || mockDataSets.low;

  const renderTabContent = () => {
    switch (activeTab) {
      case 'analysis':
        return (
          <div style={{ 
            display: 'flex',
            flexDirection: 'column',
            alignItems: 'center',
            padding: '3rem 2rem',
            textAlign: 'center'
          }}>
            <h2 style={{ 
              color: '#111827',
              fontSize: '1.5rem',
              fontWeight: '600',
              marginBottom: '1rem'
            }}>
              Genomic Analysis Complete
            </h2>
            <p style={{ 
              color: '#4B5563',
              fontSize: '1rem',
              marginBottom: '2rem',
              maxWidth: '600px'
            }}>
              Comprehensive genomic analysis has been processed successfully.
            </p>
            
            <div style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB',
              maxWidth: '600px',
              width: '100%',
              textAlign: 'left'
            }}>
              <h3 style={{ 
                color: '#111827',
                fontSize: '1.125rem',
                fontWeight: '600',
                marginBottom: '1rem'
              }}>
                Analysis Summary
              </h3>
              <div style={{ color: '#4B5563', fontSize: '0.875rem', lineHeight: '1.5' }}>
                <p style={{ marginBottom: '0.5rem' }}>
                  <strong style={{ color: '#111827' }}>Risk Level:</strong> {currentData?.riskLevel}
                </p>
                <p style={{ marginBottom: '0.5rem' }}>
                  <strong style={{ color: '#111827' }}>Condition:</strong> {currentData?.condition}
                </p>
                <p style={{ marginBottom: '0.5rem' }}>
                  <strong style={{ color: '#111827' }}>Risk Score:</strong> {currentData?.riskScore}
                </p>
                <p style={{ marginBottom: '0.5rem' }}>
                  <strong style={{ color: '#111827' }}>Status:</strong> Analysis Complete
                </p>
              </div>
              
              <button 
                style={{
                  marginTop: '1.5rem',
                  padding: '0.5rem 1rem',
                  background: '#DBEAFE',
                  color: '#2563EB',
                  border: '1px solid #BFDBFE',
                  borderRadius: '0.375rem',
                  fontSize: '0.875rem',
                  fontWeight: '500',
                  cursor: 'pointer',
                  transition: 'all 200ms ease'
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.background = '#BFDBFE';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.background = '#DBEAFE';
                }}
                onClick={() => alert('Detailed Analysis\n\nThis would display:\n• Complete variant analysis\n• Statistical significance\n• Population frequencies\n• Clinical interpretations\n• Pathway involvement\n• Treatment recommendations')}
              >
                View Detailed Analysis
              </button>
            </div>
          </div>
        );
      default:
        return (
          <div style={{ 
            display: 'flex',
            flexDirection: 'column',
            alignItems: 'center',
            justifyContent: 'center',
            padding: '3rem 2rem',
            textAlign: 'center',
            minHeight: '400px'
          }}>
            <h3 style={{ 
              marginBottom: '1rem',
              color: '#111827',
              fontSize: '1.25rem',
              fontWeight: '600'
            }}>
              Analysis Module Available
            </h3>
            <p style={{ 
              textAlign: 'center',
              maxWidth: '400px',
              lineHeight: '1.6',
              marginBottom: '2rem',
              color: '#4B5563'
            }}>
              {activeTab === 'variants' && 'Comprehensive variant heatmap and detailed genomic analysis are now available for review.'}
              {activeTab === 'pathways' && 'Biological pathway enrichment analysis has been completed and is ready for clinical interpretation.'}
              {activeTab === 'clinical' && 'Comprehensive clinical report has been generated with actionable insights and recommendations.'}
              {activeTab === 'family' && 'Family pedigree analysis and genetic testing recommendations are available for review.'}
            </p>
            <button 
              style={{
                padding: '0.5rem 1rem',
                background: '#DBEAFE',
                color: '#2563EB',
                border: '1px solid #BFDBFE',
                borderRadius: '0.375rem',
                fontSize: '0.875rem',
                fontWeight: '500',
                cursor: 'pointer',
                transition: 'all 200ms ease'
              }}
              onMouseEnter={(e) => {
                e.currentTarget.style.background = '#BFDBFE';
              }}
              onMouseLeave={(e) => {
                e.currentTarget.style.background = '#DBEAFE';
              }}
              onClick={() => alert(`${activeTab} analysis data is now available for detailed clinical review and interpretation.`)}
            >
              View {activeTab} Analysis
            </button>
          </div>
        );
    }
  };

  return (
    <Layout>
      <section style={{ 
        background: '#F9FAFB',
        minHeight: 'calc(100vh - 4rem)',
        padding: '2rem 0'
      }}>
        <div style={{ 
          maxWidth: '1400px',
          margin: '0 auto',
          padding: '0 1.5rem'
        }}>
          {/* Header */}
          <div style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            marginBottom: '2rem'
          }}>
            <div>
              <h1 style={{
                fontSize: '1.875rem',
                fontWeight: 'bold',
                color: '#111827',
                marginBottom: '0.5rem'
              }}>
                Clinical Genomics Dashboard
              </h1>
              <p style={{
                color: '#4B5563',
                fontSize: '1rem'
              }}>
                Comprehensive genomic analysis • {currentData.riskLevel}
              </p>
            </div>
            <div style={{ display: 'flex', gap: '1rem', alignItems: 'center' }}>
              <button 
                onClick={() => navigate(`/dashboard?risk=${riskLevel}`)}
                style={{
                  padding: '0.5rem 1rem',
                  background: '#E5E7EB',
                  color: '#374151',
                  border: 'none',
                  borderRadius: '0.375rem',
                  fontSize: '0.875rem',
                  fontWeight: '500',
                  cursor: 'pointer',
                  transition: 'all 200ms ease'
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.background = '#D1D5DB';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.background = '#E5E7EB';
                }}
              >
                Back to Dashboard
              </button>
              <div style={{
                background: '#22C55E',
                color: '#FFFFFF',
                padding: '0.5rem 1rem',
                borderRadius: '0.375rem',
                fontSize: '0.875rem',
                fontWeight: '500'
              }}>
                Analysis Complete
              </div>
            </div>
          </div>

          <div style={{ display: 'grid', gridTemplateColumns: '350px 1fr', gap: '2rem' }}>
            {/* Sidebar */}
            <div style={{
              background: '#FFFFFF',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB',
              overflow: 'hidden'
            }}>
              {/* Risk Info */}
              <div style={{
                background: '#F9FAFB',
                padding: '1.5rem',
                borderBottom: '1px solid #E5E7EB'
              }}>
                <div style={{ 
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  marginBottom: '0.5rem',
                  color: '#111827'
                }}>
                  {currentData.riskLevel}
                </div>
                <div style={{ 
                  fontSize: '0.875rem',
                  lineHeight: '1.4',
                  color: '#4B5563'
                }}>
                  <span dangerouslySetInnerHTML={{ __html: currentData.details }} />
                </div>
                <div style={{
                  background: '#DBEAFE',
                  color: '#2563EB',
                  padding: '0.5rem',
                  borderRadius: '0.375rem',
                  marginTop: '1rem',
                  fontSize: '0.875rem',
                  fontWeight: '500',
                  textAlign: 'center'
                }}>
                  Risk Score: {currentData.riskScore}
                </div>
              </div>

              {/* Navigation */}
              <div style={{ padding: '1rem' }}>
                <div style={{
                  fontSize: '0.75rem',
                  fontWeight: '600',
                  textTransform: 'uppercase',
                  letterSpacing: '0.05em',
                  color: '#4B5563',
                  marginBottom: '0.75rem'
                }}>
                  Clinical Workflow
                </div>
                
                {[
                  { id: 'analysis', label: 'Genomic Analysis' },
                  { id: 'variants', label: 'Variant Heatmap' },
                  { id: 'pathways', label: 'Pathway Analysis' },
                  { id: 'clinical', label: 'Clinical Report' },
                  { id: 'family', label: 'Family Analysis' }
                ].map(item => (
                  <div 
                    key={item.id}
                    style={{
                      display: 'flex',
                      alignItems: 'center',
                      gap: '0.75rem',
                      padding: '0.75rem',
                      marginBottom: '0.25rem',
                      borderRadius: '0.375rem',
                      cursor: 'pointer',
                      transition: 'all 200ms ease',
                      fontWeight: '500',
                      fontSize: '0.875rem',
                      color: activeTab === item.id ? '#FFFFFF' : '#111827',
                      background: activeTab === item.id ? '#2563EB' : 'transparent'
                    }}
                    onMouseEnter={(e) => {
                      if (activeTab !== item.id) {
                        e.currentTarget.style.background = '#F3F4F6';
                      }
                    }}
                    onMouseLeave={(e) => {
                      if (activeTab !== item.id) {
                        e.currentTarget.style.background = 'transparent';
                      }
                    }}
                    onClick={() => setActiveTab(item.id)}
                  >
                    <span>{item.label}</span>
                  </div>
                ))}
              </div>

              {/* Alerts */}
              <div style={{ padding: '1rem', borderTop: '1px solid #E5E7EB' }}>
                <div style={{
                  fontSize: '0.75rem',
                  fontWeight: '600',
                  textTransform: 'uppercase',
                  letterSpacing: '0.05em',
                  color: '#4B5563',
                  marginBottom: '0.75rem'
                }}>
                  Clinical Alerts
                </div>
                
                {currentData.alerts.map((alert: Alert, index: number) => (
                  <div 
                    key={index}
                    style={{
                      background: alert.type === 'critical' ? '#FEF2F2' : 
                                 alert.type === 'warning' ? '#FFFBEB' : 
                                 alert.type === 'info' ? '#EFF6FF' : '#F0FDF4',
                      border: `1px solid ${alert.type === 'critical' ? '#FECACA' : 
                                           alert.type === 'warning' ? '#FDE68A' : 
                                           alert.type === 'info' ? '#BFDBFE' : '#BBF7D0'}`,
                      borderLeft: `3px solid ${alert.type === 'critical' ? '#EF4444' : 
                                               alert.type === 'warning' ? '#F59E0B' : 
                                               alert.type === 'info' ? '#2563EB' : '#22C55E'}`,
                      padding: '0.75rem',
                      marginBottom: '0.5rem',
                      borderRadius: '0.375rem',
                      fontSize: '0.75rem',
                      cursor: 'pointer'
                    }}
                  >
                    <div style={{
                      display: 'flex',
                      alignItems: 'center',
                      gap: '0.5rem',
                      marginBottom: '0.25rem'
                    }}>
                      <div style={{
                        width: '0.5rem',
                        height: '0.5rem',
                        borderRadius: '50%',
                        background: alert.type === 'critical' ? '#EF4444' : 
                                   alert.type === 'warning' ? '#F59E0B' : 
                                   alert.type === 'info' ? '#2563EB' : '#22C55E'
                      }}></div>
                      <strong style={{ fontSize: '0.8rem', color: '#111827' }}>{alert.title}</strong>
                    </div>
                    <p style={{ fontSize: '0.75rem', opacity: 0.8, color: '#4B5563' }}>{alert.desc}</p>
                  </div>
                ))}
              </div>
            </div>

            {/* Main Content */}
            <div style={{
              background: '#FFFFFF',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB',
              minHeight: '600px'
            }}>
              {renderTabContent()}
            </div>
          </div>
        </div>
      </section>
    </Layout>
  );
};

export default ClinicalViewPage; 