import React from 'react';

const AnalysisDisclaimer: React.FC = () => {
  return (
    <div style={{
      marginTop: '3rem',
      padding: '1.25rem',
      background: '#F9FAFB',
      borderRadius: '0.75rem',
      border: '1px solid #E5E7EB',
      boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)'
    }}>
      <div style={{
        display: 'flex',
        alignItems: 'flex-start',
        gap: '0.75rem'
      }}>
        <div style={{
          width: '20px',
          height: '20px',
          borderRadius: '50%',
          background: '#6B7280',
          color: '#FFFFFF',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          fontSize: '0.875rem',
          fontWeight: 'bold',
          flexShrink: 0,
          marginTop: '0.125rem'
        }}>
          i
        </div>
        
        <div style={{ flex: 1 }}>
          <h3 style={{
            fontSize: '1rem',
            fontWeight: '600',
            color: '#6B7280',
            margin: '0 0 0.5rem 0'
          }}>
            Analysis Limitations
          </h3>
          
          <p style={{
            fontSize: '0.875rem',
            lineHeight: '1.6',
            color: '#9CA3AF',
            margin: 0
          }}>
            This analysis is based solely on the provided genomic data and does not account for other critical factors like family history or lifestyle. The risk model is trained on currently known genetic associations and may not accurately assess risk from novel or uncharacterized variants. Always consult a healthcare professional.
          </p>
        </div>
      </div>
    </div>
  );
};

export default AnalysisDisclaimer; 