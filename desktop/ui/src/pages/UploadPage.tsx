import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import Layout from '../components/Layout';
import { useGeneKnowTauri } from '../api/geneknowTauri';
import type { JobProgress } from '../api/geneknowPipeline';
import { invoke } from '@tauri-apps/api/core';
import { apiConfig } from '../api/apiConfig';

interface ReportConfig {
  backend: 'ollama' | 'huggingface' | 'none'
  model_name: string | null
  style: 'clinician' | 'technical' | 'patient'
}

// Icon components
const DocumentIcon = () => (
  <svg style={{ width: '2.5rem', height: '2.5rem', color: '#9CA3AF' }} xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth="1" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
  </svg>
);

// Loading Spinner component
const LoadingSpinner = () => (
  <div style={{
    width: '3rem',
    height: '3rem',
    border: '4px solid #E5E7EB',
    borderTop: '4px solid #2563EB',
    borderRadius: '50%',
    animation: 'spin 1s linear infinite'
  }}>
    <style>
      {`
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
      `}
    </style>
  </div>
);

// Loading Modal component
const LoadingModal: React.FC<{ isVisible: boolean; currentStep?: string }> = ({ isVisible, currentStep }) => {
  if (!isVisible) return null;

  return (
    <div style={{
      position: 'fixed',
      top: 0,
      left: 0,
      right: 0,
      bottom: 0,
      backgroundColor: 'rgba(0, 0, 0, 0.5)',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      zIndex: 1000
    }}>
      <div style={{
        backgroundColor: 'rgba(255, 255, 255, 0.92)',
        borderRadius: '0.75rem',
        padding: '4rem 3rem',
        maxWidth: '800px',
        width: '90%',
        minHeight: '500px',
        textAlign: 'center',
        boxShadow: '0 10px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)',
        display: 'flex',
        flexDirection: 'column',
        justifyContent: 'center'
      }}>
        <div style={{ 
          marginBottom: '2rem',
          display: 'flex',
          justifyContent: 'center',
          alignItems: 'center'
        }}>
          <LoadingSpinner />
        </div>
        <h3 style={{
          fontSize: '1.25rem',
          fontWeight: '600',
          color: '#111827',
          marginBottom: '1rem'
        }}>
          Analyzing Your Genomic Data
        </h3>
        <p style={{
          fontSize: '0.875rem',
          color: '#6B7280',
          marginBottom: '0.5rem'
        }}>
          {currentStep || 'Processing your file...'}
        </p>
        <p style={{
          fontSize: '0.75rem',
          color: '#9CA3AF'
        }}>
          This may take a few moments. Please don't close the application.
        </p>
      </div>
    </div>
  );
};

// Mock Test Case Card component
interface MockTestCaseProps {
  riskLevel: string;
  condition: string;
  description: string;
  onClick: () => void;
}

const MockTestCase: React.FC<MockTestCaseProps> = ({ riskLevel, condition, description, onClick }) => (
  <div
    onClick={onClick}
    style={{
      background: '#FFFFFF',
      padding: '1.25rem',
      borderRadius: '0.75rem',
      boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
      border: '1px solid #E5E7EB',
      cursor: 'pointer',
      transition: 'all 200ms ease',
      textAlign: 'center',
      minWidth: '220px',
      maxWidth: '250px',
      flex: '1'
    }}
    onMouseEnter={(e) => {
      e.currentTarget.style.transform = 'translateY(-2px)';
      e.currentTarget.style.boxShadow = '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)';
    }}
    onMouseLeave={(e) => {
      e.currentTarget.style.transform = 'translateY(0)';
      e.currentTarget.style.boxShadow = '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)';
    }}
  >
    <h4 style={{ fontWeight: '600', color: '#111827', marginBottom: '0.5rem', fontSize: '0.95rem' }}>{riskLevel}</h4>
    <p style={{ fontSize: '0.85rem', color: '#4B5563', marginBottom: '0.5rem', fontWeight: '500' }}>{condition}</p>
    <p style={{ fontSize: '0.75rem', color: '#6B7280', lineHeight: '1.25' }}>{description}</p>
  </div>
);

const UploadPage: React.FC = () => {
  const navigate = useNavigate();
  const [file, setFile] = useState<File | null>(null);
  const [filePath, setFilePath] = useState<string | null>(null);
  const [isDragging, setIsDragging] = useState(false);
  const [isProcessing, setIsProcessing] = useState(false);
  const [progress, setProgress] = useState<JobProgress | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [reportConfig, setReportConfig] = useState<ReportConfig>({
    backend: 'ollama',
    model_name: null,
    style: 'clinician'
  });
  const [showGuidance, setShowGuidance] = useState(false);
  
  const { processAndWait, ensureApiRunning } = useGeneKnowTauri();

  useEffect(() => {
    const fetchReportConfig = async () => {
      try {
        const response = await fetch(`${apiConfig.getBaseUrl()}/api/report-generator/config`);
        if (response.ok) {
          const config = await response.json();
          setReportConfig({
            backend: config.backend || 'ollama',
            model_name: config.model_name,
            style: config.style || 'clinician'
          });
          setShowGuidance(config.backend === 'none');
        } else {
          throw new Error('Failed to fetch config');
        }
      } catch (err) {
        console.error('Failed to fetch report config:', err);
        setShowGuidance(true); // Show guidance if config fetch fails
      }
    };
    fetchReportConfig();
  }, []);

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files[0]) {
      handleFileSelect(e.target.files[0]);
    }
  };

  const handleFileSelect = async (selectedFile: File) => {
    try {
      setFile(selectedFile);
      setError(null);
      
      // Read the file content
      const arrayBuffer = await selectedFile.arrayBuffer();
      const fileContent = new Uint8Array(arrayBuffer);
      
      // Save the file temporarily on the backend
      const tempPath = await invoke<string>('save_temp_file', {
        fileName: selectedFile.name,
        fileContent: Array.from(fileContent)
      });
      
      setFilePath(tempPath);
    } catch (err) {
      console.error('Error handling file:', err);
      setError('Failed to handle file');
    }
  };

  const handleDragEvents = (e: React.DragEvent, dragging: boolean) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(dragging);
  };

  const handleDrop = (e: React.DragEvent) => {
    handleDragEvents(e, false);
    if (e.dataTransfer.files && e.dataTransfer.files[0]) {
      handleFileSelect(e.dataTransfer.files[0]);
    }
  };

  const handleBackClick = () => {
    navigate('/');
  };

  const handleStartAnalysis = async () => {
    if (!file && !filePath) {
      setError('Please select a file first');
      return;
    }

    setIsProcessing(true);
    setError(null);
    setProgress(null);

    try {
      // Ensure API server is running
      await ensureApiRunning();

      // For desktop app, we need a file path
      if (!filePath) {
        setError('Please use the file browser to select a file');
        setIsProcessing(false);
        return;
      }

      // Process the file with progress tracking
      const result = await processAndWait(
        filePath,
        {
          language: 'en',
          include_technical: true,
          patient_data: {
            age: 45,
            sex: 'F',
            family_history: false
          }
        },
        (jobProgress) => {
          setProgress(jobProgress);
        }
      );
      
      if (!result) {
        throw new Error('No results returned from processing');
      }
      
      // Clean up the temporary file
      try {
        await invoke('delete_temp_file', { filePath });
      } catch (cleanupError) {
        console.warn('Failed to clean up temporary file:', cleanupError);
      }
      
      // Navigate to dashboard with results
      navigate('/dashboard', { 
        state: { 
          results: result,
          fileName: file?.name 
        } 
      });
      
    } catch (err) {
      console.error('Processing error:', err);
      
      const errorMessage = err instanceof Error ? err.message : 'Failed to process file';
      setError(`Processing failed: ${errorMessage}`);
      
      // Clean up the temporary file on error too
      if (filePath) {
        try {
          await invoke('delete_temp_file', { filePath });
        } catch (cleanupError) {
          console.warn('Failed to clean up temporary file:', cleanupError);
        }
      }
      
      setIsProcessing(false);
    }
  };

  const handleMockDataSelect = (riskLevel: string) => {
    // Navigate to dashboard with the risk level as a URL parameter
    navigate(`/dashboard?risk=${riskLevel}`);
  };



  return (
    <Layout>
      <section style={{ 
        background: '#F9FAFB',
        minHeight: 'calc(100vh - 4rem)',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        padding: '1rem 0'
      }}>
        <div style={{ 
          maxWidth: '50rem',
          width: '100%',
          padding: '0 1.5rem',
          textAlign: 'center'
        }}>
          <div>
            <h2 style={{
              fontSize: '1.875rem',
              fontWeight: 'bold',
              color: '#111827',
              marginBottom: '0.5rem'
            }}>
              Upload Your Genomic Data File
            </h2>
            <p style={{
              fontSize: '1.125rem',
              color: '#4B5563',
              marginBottom: '1.5rem'
            }}>
              Select your file to begin the analysis. Your data remains on your computer at all times.
            </p>

            {/* Error display */}
            {error && (
              <div style={{
                padding: '1rem',
                marginBottom: '1rem',
                backgroundColor: '#FEE2E2',
                border: '1px solid #FCA5A5',
                borderRadius: '0.5rem',
                color: '#B91C1C',
                fontSize: '0.875rem'
              }}>
                {error}
              </div>
            )}

            {/* Loading Modal */}
            <LoadingModal isVisible={isProcessing} currentStep={progress?.current_step} />

            <div
              style={{
                padding: '2rem 3rem',
                border: `2px dashed ${isDragging ? '#2563EB' : '#D1D5DB'}`,
                borderRadius: '0.75rem',
                backgroundColor: isDragging ? '#EFF6FF' : '#FFFFFF',
                transition: 'all 200ms ease',
                cursor: isProcessing ? 'not-allowed' : 'pointer',
                marginBottom: '1.5rem',
                maxWidth: '600px',
                margin: '0 auto 1.5rem',
                display: 'flex',
                flexDirection: 'column',
                alignItems: 'center',
                opacity: isProcessing ? 0.6 : 1
              }}
              onDragEnter={(e) => !isProcessing && handleDragEvents(e, true)}
              onDragLeave={(e) => !isProcessing && handleDragEvents(e, false)}
              onDragOver={(e) => !isProcessing && handleDragEvents(e, true)}
              onDrop={(e) => !isProcessing && handleDrop(e)}
              onClick={() => !isProcessing && document.getElementById('file-upload')?.click()}
            >
              <input
                type="file"
                id="file-upload"
                style={{ display: 'none' }}
                onChange={handleFileChange}
                accept=".fastq,.vcf,.bam,.maf,.gz"
                disabled={isProcessing}
              />
              <DocumentIcon />
              {file ? (
                <div style={{ marginTop: '1rem', textAlign: 'center' }}>
                  <p style={{ fontWeight: '600', color: '#111827' }}>Selected File:</p>
                  <p style={{ fontSize: '0.875rem', color: '#2563EB' }}>{file.name}</p>
                  {filePath && (
                    <p style={{ fontSize: '0.75rem', color: '#6B7280', marginTop: '0.25rem' }}>
                      Path: {filePath}
                    </p>
                  )}
                </div>
              ) : (
                <>
                  <p style={{ 
                    marginTop: '0.5rem',
                    fontWeight: '600',
                    color: '#2563EB'
                  }}>
                    Click to browse files
                  </p>
                  <p style={{ 
                    fontSize: '0.875rem',
                    color: '#6B7280'
                  }}>
                    Supported formats: FASTQ, VCF, BAM, MAF
                  </p>
                </>
              )}
            </div>

            <div style={{ 
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              gap: '1rem',
              marginBottom: '2rem'
            }}>
              <button
                onClick={handleBackClick}
                disabled={isProcessing}
                style={{
                  padding: '0.75rem 1.5rem',
                  color: '#374151',
                  background: '#E5E7EB',
                  border: 'none',
                  borderRadius: '0.5rem',
                  cursor: isProcessing ? 'not-allowed' : 'pointer',
                  transition: 'all 200ms ease',
                  fontSize: '1rem',
                  fontWeight: 'bold',
                  opacity: isProcessing ? 0.6 : 1
                }}
                onMouseEnter={(e) => {
                  if (!isProcessing) {
                    e.currentTarget.style.background = '#D1D5DB';
                  }
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.background = '#E5E7EB';
                }}
              >
                Back
              </button>

              <button
                onClick={handleStartAnalysis}
                disabled={(!file && !filePath) || isProcessing}
                style={{
                  padding: '0.75rem 1.5rem',
                  fontWeight: 'bold',
                  color: '#FFFFFF',
                  background: file && !isProcessing ? '#2563EB' : '#9CA3AF',
                  border: 'none',
                  borderRadius: '0.5rem',
                  cursor: file && !isProcessing ? 'pointer' : 'not-allowed',
                  transition: 'all 200ms ease',
                  boxShadow: file && !isProcessing ? '0 1px 2px 0 rgba(0, 0, 0, 0.05)' : 'none',
                  fontSize: '1rem'
                }}
                onMouseEnter={(e) => {
                  if (file && !isProcessing) {
                    e.currentTarget.style.background = '#1D4ED8';
                  }
                }}
                onMouseLeave={(e) => {
                  if (file && !isProcessing) {
                    e.currentTarget.style.background = '#2563EB';
                  }
                }}
              >
                {isProcessing ? 'Processing...' : 'Start Analysis'}
              </button>
            </div>

            {/* No LLM Guidance */}
            {showGuidance && reportConfig.backend === 'none' && (
              <div style={{ 
                marginTop: '1rem',
                padding: '1rem',
                background: '#FEF3C7',
                borderRadius: '0.5rem',
                border: '1px solid #F59E0B',
                textAlign: 'center',
                maxWidth: '600px',
                margin: '1rem auto 0'
              }}>
                <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', gap: '0.5rem', marginBottom: '0.5rem' }}>
                  <svg style={{ width: '1.25rem', height: '1.25rem', color: '#D97706' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" d="M9.594 3.94c.09-.542.56-.94 1.11-.94h2.593c.55 0 1.02.398 1.11.94l.213 1.281c.063.374.313.686.645.87.074.04.147.083.22.127.324.196.72.257 1.075.124l1.217-.456a1.125 1.125 0 011.37.49l1.296 2.247a1.125 1.125 0 01-.26 1.431l-1.003.827c-.293.24-.438.613-.431.992a6.759 6.759 0 010 .255c-.007.378.138.75.43.99l1.005.828c.424.35.534.954.26 1.43l-1.298 2.247a1.125 1.125 0 01-1.369.491l-1.217-.456c-.355-.133-.75-.072-1.076.124a6.57 6.57 0 01-.22.128c-.331.183-.581.495-.644.869l-.213 1.28c-.09.543-.56.941-1.11.941h-2.594c-.55 0-1.02-.398-1.11-.94l-.213-1.281c-.062-.374-.312-.686-.644-.87a6.52 6.52 0 01-.22-.127c-.325-.196-.72-.257-1.076-.124l-1.217.456a1.125 1.125 0 01-1.369-.49l-1.297-2.247a1.125 1.125 0 01.26-1.431l1.004-.827c.292-.24.437-.613.43-.992a6.932 6.932 0 010-.255c.007-.378-.138-.75-.43-.99l-1.004-.828a1.125 1.125 0 01-.26-1.43l1.297-2.247a1.125 1.125 0 011.37-.491l1.216.456c.356.133.751.072 1.076-.124.072-.044.146-.087.22-.128.332-.183.582-.495.644-.869l.214-1.281z" />
                    <path strokeLinecap="round" strokeLinejoin="round" d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
                  </svg>
                  <span style={{ fontWeight: '600', color: '#92400E', fontSize: '0.875rem' }}>
                    Template-Based Reports
                  </span>
                </div>
                <p style={{ fontSize: '0.75rem', color: '#92400E', margin: '0 0 0.75rem' }}>
                  Your reports will use standard templates without AI enhancement. For more comprehensive and personalized reports, please configure your AI models.
                </p>
                <button
                  onClick={() => navigate('/settings')}
                  style={{
                    padding: '0.5rem 1rem',
                    fontSize: '0.75rem',
                    fontWeight: '600',
                    color: '#92400E',
                    background: '#FFFFFF',
                    border: '1px solid #F59E0B',
                    borderRadius: '0.375rem',
                    cursor: 'pointer',
                    transition: 'all 200ms ease'
                  }}
                  onMouseEnter={(e) => {
                    e.currentTarget.style.background = '#FEF3C7';
                  }}
                  onMouseLeave={(e) => {
                    e.currentTarget.style.background = '#FFFFFF';
                  }}
                >
                  Configure AI Models
                </button>
              </div>
            )}
          </div>

          {/* Mock Test Cases Section */}
          <div style={{
            padding: '1.5rem 0',
            borderTop: '1px solid #E5E7EB'
          }}>
            <div style={{ 
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              gap: '0.5rem',
              marginBottom: '1.5rem'
            }}>
              <svg style={{ width: '1.5rem', height: '1.5rem', color: '#2563EB' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
                <path strokeLinecap="round" strokeLinejoin="round" d="M3 13.125C3 12.504 3.504 12 4.125 12h2.25c.621 0 1.125.504 1.125 1.125v6.75C7.5 20.496 6.996 21 6.375 21h-2.25A1.125 1.125 0 013 19.875v-6.75zM9.75 8.625c0-.621.504-1.125 1.125-1.125h2.25c.621 0 1.125.504 1.125 1.125v11.25c0 .621-.504 1.125-1.125 1.125h-2.25a1.125 1.125 0 01-1.125-1.125V8.625zM16.5 4.125c0-.621.504-1.125 1.125-1.125h2.25C20.496 3 21 3.504 21 4.125v15.75c0 .621-.504 1.125-1.125 1.125h-2.25a1.125 1.125 0 01-1.125-1.125V4.125z" />
              </svg>
              <h3 style={{
                fontSize: '1.25rem',
                fontWeight: 'bold',
                color: '#111827'
              }}>
                Or Use Test Cases
              </h3>
            </div>

            <div style={{
              display: 'flex',
              justifyContent: 'center',
              gap: '1.5rem',
              flexWrap: 'wrap',
              maxWidth: '900px',
              margin: '0 auto'
            }}>
              <MockTestCase
                riskLevel="High Risk"
                condition="Hereditary Breast and Ovarian Cancer"
                description="BRCA1/2 Positive - High-risk genetic profile with pathogenic variants"
                onClick={() => handleMockDataSelect('high')}
              />
              <MockTestCase
                riskLevel="Medium Risk"
                condition="Lynch Syndrome"
                description="Colorectal cancer predisposition with MSI-H positive markers"
                onClick={() => handleMockDataSelect('medium')}
              />
              <MockTestCase
                riskLevel="Low Risk"
                condition="Li-Fraumeni Syndrome"
                description="TP53 variant with low penetrance - protective factors identified"
                onClick={() => handleMockDataSelect('low')}
              />
            </div>
          </div>
        </div>
      </section>
    </Layout>
  );
};

export default UploadPage; 