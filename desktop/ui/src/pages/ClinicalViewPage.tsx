import React, { useState } from 'react';
import { useNavigate, useSearchParams } from 'react-router-dom';

// Type definitions
interface Alert {
  type: 'critical' | 'warning' | 'info' | 'success';
  title: string;
  desc: string;
}



// Mock data sets for different risk levels
const mockDataSets = {
  high: {
    name: 'Emma Rodriguez',
    age: 38,
    sex: 'Female',
    condition: 'Hereditary Breast and Ovarian Cancer Syndrome',
    riskScore: '82/100 (High)',
    details: 'Age: 38 • Female<br/>Family History: Breast Cancer<br/>Referral: Oncology<br/>Previous Tests: BRCA1/2 Panel',
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
    name: 'David Kim',
    age: 42,
    sex: 'Male',
    condition: 'Lynch Syndrome',
    riskScore: '45/100 (Medium)',
    details: 'Age: 42 • Male<br/>Family History: Colorectal Cancer<br/>Referral: Oncology<br/>Previous Tests: MSI-H positive',
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
    name: 'Sarah Johnson',
    age: 29,
    sex: 'Female',
    condition: 'Li-Fraumeni Syndrome',
    riskScore: '15/100 (Low)',
    details: 'Age: 29 • Female<br/>Family History: Multiple Sarcomas<br/>Referral: Genetics<br/>Previous Tests: TP53 Sequencing',
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
  const currentPatient = mockDataSets[riskLevel as keyof typeof mockDataSets] || mockDataSets.low;

  const renderTabContent = () => {
    switch (activeTab) {
      case 'analysis':
        return (
          <div style={{ textAlign: 'center', padding: '40px' }}>
            <div style={{ fontSize: '48px', marginBottom: '20px' }}>🧬</div>
            <h2 style={{ color: '#1e293b', marginBottom: '16px' }}>Genomic Analysis Complete</h2>
            <p style={{ color: '#64748b', marginBottom: '24px' }}>Analysis for {currentPatient?.name} has been processed successfully.</p>
            <div style={{ background: 'white', padding: '24px', borderRadius: '12px', boxShadow: '0 4px 12px rgba(0, 0, 0, 0.05)', maxWidth: '600px', margin: '0 auto' }}>
              <h3 style={{ color: '#1e293b', marginBottom: '16px' }}>Analysis Summary</h3>
              <div style={{ textAlign: 'left', color: '#64748b' }}>
                <p><strong>Patient:</strong> {currentPatient?.name}</p>
                <p><strong>Condition:</strong> {currentPatient?.condition}</p>
                <p><strong>Risk Score:</strong> {currentPatient?.riskScore}</p>
                <p><strong>Status:</strong> Analysis Complete</p>
              </div>
              <button 
                style={{ 
                  background: '#f0f9ff', 
                  color: '#2563eb', 
                  border: '1px solid #bfdbfe', 
                  padding: '8px 16px', 
                  borderRadius: '8px', 
                  fontSize: '14px', 
                  fontWeight: '600', 
                  cursor: 'pointer', 
                  marginTop: '20px' 
                }}
                onClick={() => alert('🔍 Detailed Analysis\n\nThis would show:\n• Complete variant analysis\n• Statistical significance\n• Population frequencies\n• Clinical interpretations\n• Pathway involvement\n• Treatment recommendations')}
              >
                🔍 View Detailed Analysis
              </button>
            </div>
          </div>
        );
      default:
        return (
          <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', flexDirection: 'column', color: '#64748b' }}>
            <div style={{ fontSize: '64px', marginBottom: '24px', color: '#2563eb' }}>🔥</div>
            <h3 style={{ marginBottom: '12px', color: '#1e293b' }}>Analysis Available</h3>
            <p style={{ textAlign: 'center', maxWidth: '400px', lineHeight: 1.6, marginBottom: '24px' }}>
              {activeTab === 'variants' && 'Variant heatmap and detailed analysis now available'}
              {activeTab === 'pathways' && 'Biological pathway enrichment analysis completed'}
              {activeTab === 'clinical' && 'Comprehensive clinical report generated'}
              {activeTab === 'family' && 'Family pedigree and testing recommendations available'}
            </p>
            <button 
              style={{ 
                background: '#f0f9ff', 
                color: '#2563eb', 
                border: '1px solid #bfdbfe', 
                padding: '8px 16px', 
                borderRadius: '8px', 
                fontSize: '14px', 
                fontWeight: '600', 
                cursor: 'pointer' 
              }}
              onClick={() => alert(`📊 ${activeTab} data is now available for detailed analysis!`)}
            >
              📊 View {activeTab}
            </button>
          </div>
        );
    }
  };

  return (
    <div style={{ 
      fontFamily: "'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif",
      background: 'linear-gradient(135deg, #1e3c72 0%, #2a5298 100%)',
      color: '#2c3e50',
      height: '100vh',
      overflow: 'hidden'
    }}>
      <div style={{ 
        height: '100vh',
        display: 'flex',
        flexDirection: 'column',
        background: '#f8fafc',
        borderRadius: '16px',
        margin: '8px',
        boxShadow: '0 20px 40px rgba(0, 0, 0, 0.15)',
        overflow: 'hidden',
        position: 'relative'
      }}>
        {/* Medical Header */}
        <div style={{ 
          background: 'linear-gradient(135deg, #2563eb 0%, #1e40af 100%)',
          color: 'white',
          padding: '16px 24px',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          boxShadow: '0 8px 25px rgba(0, 0, 0, 0.1)',
          position: 'relative',
          zIndex: 200
        }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '16px' }}>
            <div style={{ 
              fontSize: '24px', 
              color: '#2563eb', 
              display: 'flex', 
              alignItems: 'center', 
              justifyContent: 'center', 
              width: '32px', 
              height: '32px' 
            }}>
              <svg viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg" style={{ width: '28px', height: '28px', stroke: 'currentColor', strokeWidth: '2', fill: 'none', strokeLinecap: 'round', strokeLinejoin: 'round' }}>
                <path d="M4 3c4 0 8 2 12 6" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" fill="none"/>
                <path d="M4 7c4 0 8 2 12 6" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" fill="none"/>
                <path d="M4 11c4 0 8 2 12 6" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" fill="none"/>
                <path d="M4 15c4 0 8 2 12 6" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" fill="none"/>
                <path d="M20 3c-4 0-8 2-12 6" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" fill="none"/>
                <path d="M20 7c-4 0-8 2-12 6" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" fill="none"/>
                <path d="M20 11c-4 0-8 2-12 6" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" fill="none"/>
                <path d="M20 15c-4 0-8 2-12 6" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" fill="none"/>
                <line x1="6" y1="5" x2="18" y2="5" stroke="currentColor" strokeWidth="2" strokeLinecap="round"/>
                <line x1="8" y1="9" x2="16" y2="9" stroke="currentColor" strokeWidth="2" strokeLinecap="round"/>
                <line x1="10" y1="13" x2="14" y2="13" stroke="currentColor" strokeWidth="2" strokeLinecap="round"/>
                <line x1="6" y1="17" x2="18" y2="17" stroke="currentColor" strokeWidth="2" strokeLinecap="round"/>
                <circle cx="4" cy="3" r="2" fill="currentColor"/>
                <circle cx="20" cy="3" r="2" fill="currentColor"/>
                <circle cx="4" cy="21" r="2" fill="currentColor"/>
                <circle cx="20" cy="21" r="2" fill="currentColor"/>
              </svg>
            </div>
            <div>
              <div style={{ fontSize: '26px', fontWeight: '700', letterSpacing: '-0.5px' }}>GeneKnow</div>
              <div style={{ fontSize: '14px', opacity: 0.9, fontWeight: '400' }}>Advanced Clinical Genomics Platform</div>
            </div>
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: '16px' }}>
            <div style={{ 
              background: 'rgba(16, 185, 129, 0.2)',
              padding: '12px 20px',
              borderRadius: '25px',
              fontSize: '14px',
              backdropFilter: 'blur(15px)',
              border: '1px solid rgba(16, 185, 129, 0.3)',
              color: '#10b981'
            }}>
              ✅ Analysis Complete • Clinical View Active
            </div>
          </div>
        </div>

        <div style={{ flex: 1, display: 'flex', overflow: 'hidden' }}>
          {/* Medical Sidebar */}
          <div style={{ 
            width: '350px',
            background: '#ffffff',
            borderRight: '1px solid #e2e8f0',
            display: 'flex',
            flexDirection: 'column',
            boxShadow: '0 8px 25px rgba(0, 0, 0, 0.1)',
            zIndex: 50,
            position: 'relative',
            height: '100%'
          }}>
            <div style={{ 
              padding: '18px',
              background: 'linear-gradient(135deg, #2563eb 0%, #1e40af 100%)',
              color: 'white',
              margin: '12px 16px 8px 16px',
              borderRadius: '12px',
              boxShadow: '0 8px 25px rgba(0, 0, 0, 0.1)',
              position: 'relative',
              overflow: 'hidden',
              zIndex: 2
            }}>
              <div style={{ fontSize: '17px', fontWeight: '700', marginBottom: '8px', position: 'relative', zIndex: 2 }}>
                {currentPatient.name}
              </div>
              <div style={{ fontSize: '12px', opacity: 0.95, lineHeight: 1.4, position: 'relative', zIndex: 2 }}>
                <span dangerouslySetInnerHTML={{ __html: currentPatient.details }} />
              </div>
              <div style={{ 
                background: 'rgba(255, 255, 255, 0.2)',
                padding: '5px 10px',
                borderRadius: '8px',
                marginTop: '8px',
                fontSize: '11px',
                fontWeight: '600',
                position: 'relative',
                zIndex: 2
              }}>
                <strong>Risk Score: {currentPatient.riskScore}</strong>
              </div>
            </div>

            <div style={{ 
              padding: '12px 16px 8px 16px',
              borderBottom: '1px solid #e2e8f0',
              position: 'relative',
              zIndex: 2
            }}>
              <div style={{ 
                fontSize: '10px',
                fontWeight: '700',
                textTransform: 'uppercase',
                letterSpacing: '1px',
                color: '#64748b',
                marginBottom: '10px'
              }}>Clinical Workflow</div>
              
              {[
                { id: 'analysis', icon: '🧬', label: 'Genomic Analysis' },
                { id: 'variants', icon: '🔥', label: 'Variant Heatmap' },
                { id: 'pathways', icon: '🔬', label: 'Pathway Analysis' },
                { id: 'clinical', icon: '📋', label: 'Clinical Report' },
                { id: 'family', icon: '👥', label: 'Family Analysis' }
              ].map(item => (
                <div 
                  key={item.id}
                  style={{ 
                    display: 'flex',
                    alignItems: 'center',
                    gap: '8px',
                    padding: '10px 12px',
                    marginBottom: '4px',
                    borderRadius: '8px',
                    cursor: 'pointer',
                    transition: 'all 0.3s ease',
                    fontWeight: '500',
                    position: 'relative',
                    overflow: 'hidden',
                    fontSize: '13px',
                    ...(activeTab === item.id && { background: 'linear-gradient(135deg, #2563eb 0%, #1e40af 100%)', color: 'white', boxShadow: '0 8px 25px rgba(0, 0, 0, 0.1)', transform: 'translateX(6px)' })
                  }}
                  onClick={() => setActiveTab(item.id)}
                >
                  <span style={{ fontSize: '16px', width: '18px', textAlign: 'center' }}>{item.icon}</span>
                  <span>{item.label}</span>
                </div>
              ))}
            </div>

            <div style={{ 
              flex: 1,
              padding: '12px 16px 16px 16px',
              overflowY: 'auto',
              position: 'relative',
              zIndex: 2,
              display: 'flex',
              flexDirection: 'column'
            }}>
              <div style={{ 
                fontSize: '10px',
                fontWeight: '700',
                textTransform: 'uppercase',
                letterSpacing: '1px',
                color: '#64748b',
                marginBottom: '8px'
              }}>Clinical Alerts</div>
              
              <div style={{ flex: 1, overflowY: 'auto' }}>
                {currentPatient.alerts.map((alert: Alert, index: number) => (
                  <div 
                    key={index}
                    style={{ 
                      background: alert.type === 'critical' ? '#fff5f5' : alert.type === 'warning' ? '#fffbeb' : alert.type === 'info' ? '#eff6ff' : '#f0fdf4',
                      border: `1px solid ${alert.type === 'critical' ? '#fecaca' : alert.type === 'warning' ? '#fcd34d' : alert.type === 'info' ? '#93c5fd' : '#86efac'}`,
                      borderLeft: `3px solid ${alert.type === 'critical' ? '#ef4444' : alert.type === 'warning' ? '#f59e0b' : alert.type === 'info' ? '#2563eb' : '#10b981'}`,
                      padding: '10px',
                      marginBottom: '6px',
                      borderRadius: '6px',
                      fontSize: '12px',
                      cursor: 'pointer',
                      lineHeight: 1.3
                    }}
                  >
                    <div style={{ 
                      width: '10px', 
                      height: '10px', 
                      borderRadius: '50%', 
                      display: 'inline-block', 
                      marginRight: '6px',
                      background: alert.type === 'critical' ? '#ef4444' : alert.type === 'warning' ? '#f59e0b' : '#10b981'
                    }}></div>
                    <strong style={{ display: 'block', marginBottom: '4px', fontSize: '13px' }}>{alert.title}</strong>
                    <small style={{ fontSize: '11px', opacity: 0.9 }}>{alert.desc}</small>
                  </div>
                ))}
              </div>
            </div>
          </div>

          {/* Content Area */}
          <div style={{ flex: 1, display: 'flex', flexDirection: 'column', overflow: 'hidden', position: 'relative' }}>
            <div style={{ 
              background: 'white',
              padding: '32px 32px 24px',
              borderBottom: '1px solid #e2e8f0',
              boxShadow: '0 4px 12px rgba(0, 0, 0, 0.05)',
              position: 'relative',
              zIndex: 2
            }}>
              <div style={{ 
                fontSize: '32px',
                fontWeight: '700',
                color: '#1e293b',
                marginBottom: '8px',
                background: 'linear-gradient(135deg, #1e293b 0%, #2563eb 100%)',
                WebkitBackgroundClip: 'text',
                WebkitTextFillColor: 'transparent',
                backgroundClip: 'text'
              }}>🧬 Genomic Analysis Dashboard</div>
              <div style={{ color: '#64748b', fontSize: '16px', marginBottom: '20px' }}>
                Comprehensive genomic analysis for {currentPatient.name}
              </div>
              <div style={{ display: 'flex', gap: '16px', alignItems: 'center', flexWrap: 'wrap' }}>
                <button 
                  onClick={() => navigate(`/dashboard?risk=${riskLevel}`)}
                  style={{ 
                    background: '#f3f4f6',
                    color: '#374151',
                    border: '1px solid #d1d5db',
                    padding: '8px 16px',
                    borderRadius: '8px',
                    fontSize: '14px',
                    fontWeight: '600',
                    cursor: 'pointer'
                  }}
                >
                  ← Back to Dashboard
                </button>
                {[
                  '📄 Export Report',
                  '📤 Share Findings', 
                  '📅 Schedule Consult',
                  '📊 View History'
                ].map(label => (
                  <button 
                    key={label}
                    style={{ 
                      background: '#f0f9ff',
                      color: '#2563eb',
                      border: '1px solid #bfdbfe',
                      padding: '8px 16px',
                      borderRadius: '8px',
                      fontSize: '14px',
                      fontWeight: '600',
                      cursor: 'pointer'
                    }}
                  >
                    {label}
                  </button>
                ))}
              </div>
            </div>

            <div style={{ 
              flex: 1,
              padding: '32px',
              overflowY: 'auto',
              background: '#fafbfc',
              position: 'relative',
              zIndex: 2
            }}>
              {renderTabContent()}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ClinicalViewPage; 