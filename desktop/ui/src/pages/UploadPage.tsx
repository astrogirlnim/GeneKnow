import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import Layout from '../components/Layout';

// Icon components
const DocumentIcon = () => (
  <svg style={{ width: '2.5rem', height: '2.5rem', color: '#9CA3AF' }} xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth="1" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
  </svg>
);

const InformationCircleIcon = ({ className = "w-5 h-5" }) => (
  <svg className={className} xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth="2" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
  </svg>
);

// Mock Patient Card component
interface MockPatientCardProps {
  emoji: string;
  name: string;
  condition: string;
  onClick: () => void;
}

const MockPatientCard: React.FC<MockPatientCardProps> = ({ emoji, name, condition, onClick }) => (
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
    <div style={{ fontSize: '1.75rem', marginBottom: '0.5rem' }}>{emoji}</div>
    <h4 style={{ fontWeight: '600', color: '#111827', marginBottom: '0.5rem', fontSize: '0.95rem' }}>{name}</h4>
    <p style={{ fontSize: '0.8rem', color: '#6B7280', lineHeight: '1.25' }}>{condition}</p>
  </div>
);

const UploadPage: React.FC = () => {
  const navigate = useNavigate();
  const [file, setFile] = useState<File | null>(null);
  const [isDragging, setIsDragging] = useState(false);

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files[0]) {
      setFile(e.target.files[0]);
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
      setFile(e.dataTransfer.files[0]);
    }
  };

  const handleBackClick = () => {
    navigate('/');
  };

  const handleStartAnalysis = () => {
    // TODO: Implement analysis logic
    console.log('Starting analysis with file:', file?.name);
    // For now, just show an alert
    alert('Analysis started! This will be implemented to process the file.');
  };

  const handleMockDataSelect = (patientType: string) => {
    // TODO: Implement mock data selection logic
    console.log('Selected mock patient:', patientType);
    alert(`Mock analysis started for ${patientType}! This will be implemented to process mock data.`);
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

            <div
              style={{
                padding: '2rem 3rem',
                border: `2px dashed ${isDragging ? '#2563EB' : '#D1D5DB'}`,
                borderRadius: '0.75rem',
                backgroundColor: isDragging ? '#EFF6FF' : '#FFFFFF',
                transition: 'all 200ms ease',
                cursor: 'pointer',
                marginBottom: '1.5rem',
                maxWidth: '600px',
                margin: '0 auto 1.5rem',
                display: 'flex',
                flexDirection: 'column',
                alignItems: 'center'
              }}
              onDragEnter={(e) => handleDragEvents(e, true)}
              onDragLeave={(e) => handleDragEvents(e, false)}
              onDragOver={(e) => handleDragEvents(e, true)}
              onDrop={handleDrop}
              onClick={() => document.getElementById('file-upload')?.click()}
            >
              <input
                type="file"
                id="file-upload"
                style={{ display: 'none' }}
                onChange={handleFileChange}
                accept=".fastq,.vcf,.bam"
              />
              <DocumentIcon />
              {file ? (
                <div style={{ marginTop: '1rem', textAlign: 'left' }}>
                  <p style={{ fontWeight: '600', color: '#111827' }}>Selected File:</p>
                  <p style={{ fontSize: '0.875rem', color: '#2563EB' }}>{file.name}</p>
                  <p style={{ fontSize: '0.75rem', color: '#6B7280' }}>
                    Size: {(file.size / 1024 / 1024).toFixed(2)} MB
                  </p>
                </div>
              ) : (
                <>
                  <p style={{ 
                    marginTop: '0.5rem',
                    fontWeight: '600',
                    color: '#2563EB'
                  }}>
                    Click to browse or drag & drop
                  </p>
                  <p style={{ 
                    fontSize: '0.875rem',
                    color: '#6B7280'
                  }}>
                    Supported formats: FASTQ, VCF, BAM
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
                style={{
                  padding: '0.75rem 1.5rem',
                  color: '#374151',
                  background: '#E5E7EB',
                  border: 'none',
                  borderRadius: '0.5rem',
                  cursor: 'pointer',
                  transition: 'all 200ms ease',
                  fontSize: '1rem',
                  fontWeight: 'bold'
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.background = '#D1D5DB';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.background = '#E5E7EB';
                }}
              >
                Back
              </button>
              <button
                onClick={handleStartAnalysis}
                disabled={!file}
                style={{
                  padding: '0.75rem 1.5rem',
                  fontWeight: 'bold',
                  color: '#FFFFFF',
                  background: file ? '#2563EB' : '#9CA3AF',
                  border: 'none',
                  borderRadius: '0.5rem',
                  cursor: file ? 'pointer' : 'not-allowed',
                  transition: 'all 200ms ease',
                  boxShadow: file ? '0 1px 2px 0 rgba(0, 0, 0, 0.05)' : 'none',
                  fontSize: '1rem'
                }}
                onMouseEnter={(e) => {
                  if (file) {
                    e.currentTarget.style.background = '#1D4ED8';
                  }
                }}
                onMouseLeave={(e) => {
                  if (file) {
                    e.currentTarget.style.background = '#2563EB';
                  }
                }}
              >
                Start Analysis
              </button>
            </div>
          </div>

          {/* Mock Genome Data Section */}
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
                Or Use Mock Genome Data
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
              <MockPatientCard
                emoji="👩"
                name="Emma Rodriguez"
                condition="BRCA1/2 Positive High-risk profile"
                onClick={() => handleMockDataSelect('emma')}
              />
              <MockPatientCard
                emoji="👨"
                name="David Kim"
                condition="Lynch Syndrome Colorectal cancer risk"
                onClick={() => handleMockDataSelect('david')}
              />
              <MockPatientCard
                emoji="👩"
                name="Sarah Johnson"
                condition="TP53 Mutation Li-Fraumeni syndrome"
                onClick={() => handleMockDataSelect('sarah')}
              />
            </div>
          </div>
        </div>
      </section>
    </Layout>
  );
};

export default UploadPage; 