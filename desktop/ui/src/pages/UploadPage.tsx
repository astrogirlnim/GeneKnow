import React, { useState, useRef } from 'react';

// Type definitions
interface Alert {
  type: 'critical' | 'warning' | 'info' | 'success';
  title: string;
  desc: string;
}

interface PatientData {
  name: string;
  age: number;
  sex: string;
  condition: string;
  riskScore: string;
  details: string;
  alerts: Alert[];
}

interface FileInfo {
  name: string;
  size: string;
  type: string;
}

const UploadPage: React.FC = () => {
  const [uploadStatus, setUploadStatus] = useState('ready');
  const [isDashboardEnabled, setIsDashboardEnabled] = useState(false);
  const [currentPatient, setCurrentPatient] = useState<PatientData | null>(null);
  const [activeTab, setActiveTab] = useState('analysis');
  const [uploadProgress, setUploadProgress] = useState(0);
  const [isUploading, setIsUploading] = useState(false);
  const [showFileInfo, setShowFileInfo] = useState(false);
  const [fileInfo, setFileInfo] = useState<FileInfo | null>(null);
  const [isDragOver, setIsDragOver] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Mock genome datasets
  const mockGenomeData = {
    patient1: {
      name: "Emma Rodriguez",
      age: 34,
      sex: "Female",
      condition: "BRCA1/BRCA2 Positive",
      riskScore: "95/100 (Extremely High)",
      details: "Age: 34 • Female<br>Family History: Breast Cancer (Maternal)<br>Referral: Genetic Counseling<br>Previous Tests: None",
      alerts: [
        { type: "critical", title: "🚨 Pathogenic BRCA1 Variant", desc: "c.5266dupC - Immediate action required" },
        { type: "critical", title: "⚠️ Pathogenic BRCA2 Variant", desc: "c.9976A>T - High penetrance" },
        { type: "info", title: "📊 Family History Match", desc: "Consistent with maternal lineage" },
        { type: "success", title: "✅ Quality Control", desc: "All metrics passed (>99% coverage)" }
      ]
    },
    patient2: {
      name: "David Kim",
      age: 42,
      sex: "Male",
      condition: "Lynch Syndrome",
      riskScore: "78/100 (High)",
      details: "Age: 42 • Male<br>Family History: Colorectal Cancer<br>Referral: Oncology<br>Previous Tests: MSI-H positive",
      alerts: [
        { type: "critical", title: "🚨 MLH1 Pathogenic Variant", desc: "c.1989-1G>A - Splice site mutation" },
        { type: "warning", title: "⚠️ MSH2 VUS", desc: "c.2634+1G>T - Needs monitoring" },
        { type: "info", title: "📊 Population Frequency", desc: "Variant rare in Asian populations" },
        { type: "success", title: "✅ Microsatellite Analysis", desc: "MSI-H confirmed" }
      ]
    },
    patient3: {
      name: "Sarah Johnson",
      age: 28,
      sex: "Female",
      condition: "Li-Fraumeni Syndrome",
      riskScore: "88/100 (Very High)",
      details: "Age: 28 • Female<br>Family History: Multiple Cancers<br>Referral: Pediatric Oncology<br>Previous Tests: TP53 screening",
      alerts: [
        { type: "critical", title: "🚨 TP53 Pathogenic Variant", desc: "c.742C>T - Guardian of genome" },
        { type: "warning", title: "⚠️ Early Onset Risk", desc: "Childhood cancer surveillance needed" },
        { type: "info", title: "📊 Penetrance Data", desc: "90% lifetime cancer risk" },
        { type: "success", title: "✅ Functional Analysis", desc: "DNA binding domain affected" }
      ]
    }
  } as const;

  const handleFileSelect = (files: FileList) => {
    if (files.length > 0) {
      const file = files[0];
      setFileInfo({
        name: file.name,
        size: formatFileSize(file.size),
        type: file.type || 'Unknown'
      });
      setShowFileInfo(true);
      uploadFile(file);
    }
  };

  const uploadFile = (file: File) => {
    setIsUploading(true);
    setUploadStatus('uploading');
    setUploadProgress(0);

    // Simulate upload progress
    const uploadInterval = setInterval(() => {
      setUploadProgress(prev => {
        const newProgress = prev + Math.random() * 15;
        if (newProgress >= 100) {
          clearInterval(uploadInterval);
          completeUpload(file);
          return 100;
        }
        return newProgress;
      });
    }, 200);
  };

  const completeUpload = (file: File) => {
    setUploadStatus('complete');
    
    setTimeout(() => {
      const mockData = {
        name: file.name.replace(/\.[^/.]+$/, ""),
        age: Math.floor(Math.random() * 50) + 20,
        sex: Math.random() > 0.5 ? "Female" : "Male",
        condition: "File-based Analysis",
        riskScore: Math.floor(Math.random() * 50) + 30 + "/100",
        details: `File: ${file.name}<br>Size: ${formatFileSize(file.size)}<br>Type: Genomic Analysis<br>Status: Processing Complete`,
        alerts: [
          { type: "success", title: "✅ File Processed", desc: "Genome data successfully analyzed" },
          { type: "info", title: "📊 Quality Control", desc: "All metrics within normal range" },
          { type: "warning", title: "⚠️ Analysis Ready", desc: "Results available for review" }
        ]
      };
      
      activateDashboard(mockData);
    }, 2000);
  };

  const loadMockData = (mockType: keyof typeof mockGenomeData) => {
    setUploadStatus('uploading');
    setIsUploading(true);
    
    setTimeout(() => {
      setUploadStatus('complete');
      setIsUploading(false);
      activateDashboard(mockGenomeData[mockType]);
    }, 1500);
  };

  const activateDashboard = (data: PatientData) => {
    setCurrentPatient(data);
    setIsDashboardEnabled(true);
    setUploadStatus('complete');
  };

  const formatFileSize = (bytes: number): string => {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
  };

  const handleDragEnter = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragOver(true);
  };

  const handleDragLeave = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragOver(false);
  };

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragOver(false);
    const files = e.dataTransfer.files;
    handleFileSelect(files);
  };

  const renderTabContent = () => {
    if (!isDashboardEnabled) {
      return (
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', flexDirection: 'column', color: '#64748b' }}>
          <div style={{ fontSize: '64px', marginBottom: '24px', opacity: 0.5 }}>🧬</div>
          <h3 style={{ marginBottom: '12px', color: '#1e293b' }}>Genomic Analysis Ready</h3>
          <p style={{ textAlign: 'center', maxWidth: '400px', lineHeight: 1.6 }}>
            Upload your genome data file or select from our mock datasets to begin comprehensive genomic analysis and clinical interpretation.
          </p>
        </div>
      );
    }

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
        position: 'relative',
        filter: isDashboardEnabled ? 'none' : 'blur(8px)',
        opacity: isDashboardEnabled ? 1 : 0.3,
        pointerEvents: isDashboardEnabled ? 'auto' : 'none'
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
              background: 'rgba(255, 255, 255, 0.15)',
              padding: '12px 20px',
              borderRadius: '25px',
              fontSize: '14px',
              backdropFilter: 'blur(15px)',
              border: '1px solid rgba(255, 255, 255, 0.2)',
              ...(uploadStatus === 'uploading' && { background: 'rgba(245, 158, 11, 0.2)', borderColor: '#f59e0b', color: '#f59e0b' }),
              ...(uploadStatus === 'complete' && { background: 'rgba(16, 185, 129, 0.2)', borderColor: '#10b981', color: '#10b981' })
            }}>
              {uploadStatus === 'ready' && '📁 Ready for genome upload'}
              {uploadStatus === 'uploading' && '📤 Uploading genome data...'}
              {uploadStatus === 'complete' && '✅ Upload complete • Processing...'}
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
                {currentPatient ? currentPatient.name : 'Upload Genome Data'}
              </div>
              <div style={{ fontSize: '12px', opacity: 0.95, lineHeight: 1.4, position: 'relative', zIndex: 2 }}>
                {currentPatient ? 
                  <span dangerouslySetInnerHTML={{ __html: currentPatient.details }} /> :
                  'Select a genome file or use mock data to begin analysis'
                }
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
                <strong>
                  {currentPatient ? `Risk Score: ${currentPatient.riskScore}` : 'Status: Awaiting Data'}
                </strong>
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
                    cursor: isDashboardEnabled ? 'pointer' : 'default',
                    transition: 'all 0.3s ease',
                    fontWeight: '500',
                    position: 'relative',
                    overflow: 'hidden',
                    fontSize: '13px',
                    ...(activeTab === item.id && { background: 'linear-gradient(135deg, #2563eb 0%, #1e40af 100%)', color: 'white', boxShadow: '0 8px 25px rgba(0, 0, 0, 0.1)', transform: 'translateX(6px)' })
                  }}
                  onClick={() => isDashboardEnabled && setActiveTab(item.id)}
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
                {currentPatient ? currentPatient.alerts.map((alert: Alert, index: number) => (
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
                )) : (
                  <>
                    <div style={{ background: '#eff6ff', border: '1px solid #93c5fd', borderLeft: '3px solid #2563eb', padding: '10px', marginBottom: '6px', borderRadius: '6px', fontSize: '12px', lineHeight: 1.3 }}>
                      <div style={{ width: '10px', height: '10px', borderRadius: '50%', display: 'inline-block', marginRight: '6px', background: '#10b981' }}></div>
                      <strong style={{ display: 'block', marginBottom: '4px', fontSize: '13px' }}>📁 System Status</strong>
                      <small style={{ fontSize: '11px', opacity: 0.9 }}>Ready for genome data upload</small>
                    </div>
                    <div style={{ background: '#fffbeb', border: '1px solid #fcd34d', borderLeft: '3px solid #f59e0b', padding: '10px', marginBottom: '6px', borderRadius: '6px', fontSize: '12px', lineHeight: 1.3 }}>
                      <div style={{ width: '10px', height: '10px', borderRadius: '50%', display: 'inline-block', marginRight: '6px', background: '#f59e0b' }}></div>
                      <strong style={{ display: 'block', marginBottom: '4px', fontSize: '13px' }}>🔬 Analysis Pending</strong>
                      <small style={{ fontSize: '11px', opacity: 0.9 }}>Awaiting genomic data for processing</small>
                    </div>
                  </>
                )}
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
                Upload genome data to begin comprehensive analysis
              </div>
              <div style={{ display: 'flex', gap: '16px', alignItems: 'center', flexWrap: 'wrap' }}>
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
                      cursor: isDashboardEnabled ? 'pointer' : 'not-allowed',
                      opacity: isDashboardEnabled ? 1 : 0.6
                    }}
                    disabled={!isDashboardEnabled}
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

      {/* File Upload Overlay */}
      {!isDashboardEnabled && (
        <div style={{ 
          position: 'fixed',
          top: 0,
          left: 0,
          right: 0,
          bottom: 0,
          background: 'rgba(0, 0, 0, 0.4)',
          backdropFilter: 'blur(8px)',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          zIndex: 1000
        }}>
          <div 
            style={{ 
              background: 'rgba(255, 255, 255, 0.95)',
              backdropFilter: 'blur(20px)',
              borderRadius: '24px',
              padding: '40px',
              maxWidth: '550px',
              width: '85%',
              textAlign: 'center',
              boxShadow: '0 25px 50px rgba(0, 0, 0, 0.25)',
              border: '2px solid rgba(255, 255, 255, 0.3)',
              position: 'relative',
              overflow: 'hidden',
              ...(isDragOver && { borderColor: '#2563eb', background: 'rgba(240, 249, 255, 0.98)', transform: 'scale(1.02)', boxShadow: '0 30px 60px rgba(37, 99, 235, 0.3)' })
            }}
            onDragEnter={handleDragEnter}
            onDragLeave={handleDragLeave}
            onDragOver={handleDragOver}
            onDrop={handleDrop}
          >
            <div style={{ position: 'relative', zIndex: 10 }}>
              <div style={{ fontSize: '64px', color: '#2563eb', marginBottom: '20px' }}>📁</div>
              <h2 style={{ fontSize: '24px', fontWeight: '700', color: '#1e293b', marginBottom: '10px' }}>Upload Genome Data</h2>
              <p style={{ fontSize: '15px', color: '#64748b', marginBottom: '28px', lineHeight: 1.5 }}>
                Select a VCF, FASTQ, or genome analysis file to begin processing.<br/>
                Supported formats: .vcf, .vcf.gz, .fastq, .fq, .bam, .sam
              </p>
              
              <input 
                ref={fileInputRef}
                type="file" 
                accept=".vcf,.vcf.gz,.fastq,.fq,.bam,.sam" 
                onChange={(e) => e.target.files && handleFileSelect(e.target.files)}
                style={{ position: 'absolute', opacity: 0, pointerEvents: 'none' }}
              />
              
              <button 
                style={{ 
                  background: 'linear-gradient(135deg, #2563eb 0%, #1e40af 100%)',
                  color: 'white',
                  border: 'none',
                  padding: '14px 28px',
                  borderRadius: '12px',
                  fontSize: '15px',
                  fontWeight: '600',
                  cursor: isUploading ? 'not-allowed' : 'pointer',
                  marginBottom: '20px',
                  position: 'relative',
                  overflow: 'hidden',
                  opacity: isUploading ? 0.6 : 1
                }}
                onClick={() => fileInputRef.current?.click()}
                disabled={isUploading}
              >
                {isUploading ? '📤 Uploading...' : '📤 Choose Genome File'}
              </button>
              
              {isUploading && (
                <div style={{ width: '100%', height: '8px', background: '#e2e8f0', borderRadius: '4px', overflow: 'hidden', margin: '24px 0' }}>
                  <div style={{ height: '100%', background: 'linear-gradient(135deg, #2563eb 0%, #1e40af 100%)', width: `${uploadProgress}%`, transition: 'width 0.3s ease', borderRadius: '4px' }}></div>
                </div>
              )}
              
              {isUploading && (
                <div style={{ width: '40px', height: '40px', border: '4px solid #e2e8f0', borderTop: '4px solid #2563eb', borderRadius: '50%', animation: 'spin 1s linear infinite', margin: '20px auto' }}>
                  <style>
                    {`
                      @keyframes spin {
                        0% { transform: rotate(0deg); }
                        100% { transform: rotate(360deg); }
                      }
                    `}
                  </style>
                </div>
              )}
              
              {uploadStatus === 'complete' && !isUploading && (
                <div style={{ width: '60px', height: '60px', borderRadius: '50%', background: '#10b981', position: 'relative', margin: '20px auto', display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                  <span style={{ color: 'white', fontSize: '24px', fontWeight: 'bold' }}>✓</span>
                </div>
              )}
              
              {showFileInfo && fileInfo && (
                <div style={{ background: '#f8fafc', padding: '16px', borderRadius: '8px', marginTop: '16px', borderLeft: '4px solid #2563eb' }}>
                  <h4 style={{ color: '#1e293b', marginBottom: '8px', fontSize: '16px' }}>File Information</h4>
                  <p style={{ color: '#64748b', fontSize: '14px', marginBottom: '4px' }}>Name: {fileInfo.name}</p>
                  <p style={{ color: '#64748b', fontSize: '14px', marginBottom: '4px' }}>Size: {fileInfo.size}</p>
                  <p style={{ color: '#64748b', fontSize: '14px', marginBottom: '4px' }}>Type: {fileInfo.type}</p>
                </div>
              )}
            </div>
            
            <div style={{ borderTop: '1px solid rgba(226, 232, 240, 0.6)', paddingTop: '24px', marginTop: '24px' }}>
              <h3 style={{ fontSize: '16px', fontWeight: '600', color: '#1e293b', marginBottom: '16px' }}>📊 Or Use Mock Genome Data</h3>
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(140px, 1fr))', gap: '12px' }}>
                {[
                  { id: 'patient1', icon: '👩‍⚕️', name: 'Emma Rodriguez', desc: 'BRCA1/2 Positive\nHigh-risk profile' },
                  { id: 'patient2', icon: '👨‍⚕️', name: 'David Kim', desc: 'Lynch Syndrome\nColorectal cancer risk' },
                  { id: 'patient3', icon: '👵', name: 'Sarah Johnson', desc: 'TP53 Mutation\nLi-Fraumeni syndrome' }
                ].map(patient => (
                  <div 
                    key={patient.id}
                    style={{ 
                      background: 'rgba(248, 250, 252, 0.8)',
                      backdropFilter: 'blur(10px)',
                      border: '1.5px solid rgba(226, 232, 240, 0.8)',
                      borderRadius: '10px',
                      padding: '16px',
                      textAlign: 'center',
                      cursor: 'pointer',
                      transition: 'all 0.3s ease',
                      position: 'relative',
                      overflow: 'hidden'
                    }}
                    onMouseEnter={(e) => {
                      e.currentTarget.style.transform = 'translateY(-3px)';
                      e.currentTarget.style.boxShadow = '0 8px 25px rgba(0, 0, 0, 0.15)';
                      e.currentTarget.style.borderColor = '#2563eb';
                      e.currentTarget.style.background = 'rgba(255, 255, 255, 0.95)';
                    }}
                    onMouseLeave={(e) => {
                      e.currentTarget.style.transform = 'translateY(0)';
                      e.currentTarget.style.boxShadow = 'none';
                      e.currentTarget.style.borderColor = 'rgba(226, 232, 240, 0.8)';
                      e.currentTarget.style.background = 'rgba(248, 250, 252, 0.8)';
                    }}
                    onClick={() => loadMockData(patient.id as keyof typeof mockGenomeData)}
                  >
                    <div style={{ fontSize: '28px', marginBottom: '10px', color: '#2563eb' }}>{patient.icon}</div>
                    <div style={{ fontSize: '13px', fontWeight: '600', color: '#1e293b', marginBottom: '6px' }}>{patient.name}</div>
                    <div style={{ fontSize: '11px', color: '#64748b', lineHeight: 1.3, whiteSpace: 'pre-line' }}>{patient.desc}</div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default UploadPage; 