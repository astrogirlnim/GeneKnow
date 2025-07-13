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

interface UploadOverlayProps {
  onUploadComplete: (data: PatientData) => void;
}

const UploadOverlay: React.FC<UploadOverlayProps> = ({ onUploadComplete }) => {
  const [uploadStatus, setUploadStatus] = useState('ready');
  const [uploadProgress, setUploadProgress] = useState(0);
  const [isUploading, setIsUploading] = useState(false);
  const [showFileInfo, setShowFileInfo] = useState(false);
  const [fileInfo, setFileInfo] = useState<FileInfo | null>(null);
  const [isDragOver, setIsDragOver] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Mock genome datasets
  const mockGenomeData: Record<string, PatientData> = {
    patient1: {
      name: "Emma Rodriguez",
      age: 34,
      sex: "Female",
      condition: "BRCA1/BRCA2 Positive",
      riskScore: "95/100 (Extremely High)",
      details: "Analysis Type: Genomic Variant Analysis<br>Method: Machine Learning Risk Assessment<br>Data Source: Uploaded VCF File",
      alerts: [
        { type: "critical" as const, title: "🚨 Pathogenic BRCA1 Variant", desc: "c.5266dupC - Immediate action required" },
        { type: "critical" as const, title: "⚠️ Pathogenic BRCA2 Variant", desc: "c.9976A>T - High penetrance" },
        { type: "info" as const, title: "📊 Family History Match", desc: "Consistent with maternal lineage" },
        { type: "success" as const, title: "✅ Quality Control", desc: "All metrics passed (>99% coverage)" }
      ]
    },
    patient2: {
      name: "David Kim",
      age: 42,
      sex: "Male",
      condition: "Lynch Syndrome",
      riskScore: "78/100 (High)",
      details: "Analysis Type: Genomic Variant Analysis<br>Method: Machine Learning Risk Assessment<br>Data Source: Uploaded VCF File",
      alerts: [
        { type: "critical" as const, title: "🚨 MLH1 Pathogenic Variant", desc: "c.1989-1G>A - Splice site mutation" },
        { type: "warning" as const, title: "⚠️ MSH2 VUS", desc: "c.2634+1G>T - Needs monitoring" },
        { type: "info" as const, title: "📊 Population Frequency", desc: "Variant rare in Asian populations" },
        { type: "success" as const, title: "✅ Microsatellite Analysis", desc: "MSI-H confirmed" }
      ]
    },
    patient3: {
      name: "Sarah Johnson",
      age: 28,
      sex: "Female",
      condition: "Li-Fraumeni Syndrome",
      riskScore: "88/100 (Very High)",
      details: "Analysis Type: Genomic Variant Analysis<br>Method: Machine Learning Risk Assessment<br>Data Source: Uploaded VCF File",
      alerts: [
        { type: "critical" as const, title: "🚨 TP53 Pathogenic Variant", desc: "c.742C>T - Guardian of genome" },
        { type: "warning" as const, title: "⚠️ Early Onset Risk", desc: "Childhood cancer surveillance needed" },
        { type: "info" as const, title: "📊 Penetrance Data", desc: "90% lifetime cancer risk" },
        { type: "success" as const, title: "✅ Functional Analysis", desc: "DNA binding domain affected" }
      ]
    }
  };

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
      const mockData: PatientData = {
        name: file.name.replace(/\.[^/.]+$/, ""),
        age: Math.floor(Math.random() * 50) + 20,
        sex: Math.random() > 0.5 ? "Female" : "Male",
        condition: "File-based Analysis",
        riskScore: Math.floor(Math.random() * 50) + 30 + "/100",
        details: `File: ${file.name}<br>Size: ${formatFileSize(file.size)}<br>Type: Genomic Analysis<br>Status: Processing Complete`,
        alerts: [
          { type: "success" as const, title: "✅ File Processed", desc: "Genome data successfully analyzed" },
          { type: "info" as const, title: "📊 Quality Control", desc: "All metrics within normal range" },
          { type: "warning" as const, title: "⚠️ Analysis Ready", desc: "Results available for review" }
        ]
      };
      
      onUploadComplete(mockData);
    }, 2000);
  };

  const loadMockData = (mockType: keyof typeof mockGenomeData) => {
    setUploadStatus('uploading');
    setIsUploading(true);
    
    setTimeout(() => {
      setUploadStatus('complete');
      setIsUploading(false);
      onUploadComplete(mockGenomeData[mockType]);
    }, 1500);
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

  return (
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
  );
};

export default UploadOverlay; 