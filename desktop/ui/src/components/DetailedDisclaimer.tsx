import React, { useState } from 'react';

const ChevronDownIcon = () => (
  <svg width="20" height="20" viewBox="0 0 20 20" fill="currentColor">
    <path fillRule="evenodd" d="M5.293 7.293a1 1 0 011.414 0L10 10.586l3.293-3.293a1 1 0 111.414 1.414l-4 4a1 1 0 01-1.414 0l-4-4a1 1 0 010-1.414z" clipRule="evenodd" />
  </svg>
);

const DetailedDisclaimer: React.FC = () => {
  const [isExpanded, setIsExpanded] = useState(false);

  return (
    <div style={{
      marginTop: '3rem',
      background: '#F9FAFB',
      borderRadius: '0.75rem',
      border: '1px solid #E5E7EB',
      boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)'
    }}>
      <button
        onClick={() => setIsExpanded(!isExpanded)}
        style={{
          width: '100%',
          padding: '1.25rem',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          background: 'none',
          border: 'none',
          cursor: 'pointer',
          textAlign: 'left'
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '0.75rem' }}>
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
          }}>
            i
          </div>
          <h3 style={{
            fontSize: '1rem',
            fontWeight: '600',
            color: '#374151',
            margin: 0
          }}>
            Methodology & Limitations
          </h3>
        </div>
        <div style={{ transition: 'transform 200ms ease', transform: isExpanded ? 'rotate(180deg)' : 'rotate(0deg)' }}>
          <ChevronDownIcon />
        </div>
      </button>

      {isExpanded && (
        <div style={{
          padding: '0 1.25rem 1.25rem 1.25rem',
          borderTop: '1px solid #E5E7EB',
          color: '#4B5563',
          fontSize: '0.875rem',
          lineHeight: '1.6'
        }}>
          <p style={{ marginTop: '1.25rem', marginBottom: '1.5rem' }}>
            This genomic analysis was performed using a proprietary model trained on variant data from 1000 Genomes and gnomAD. The risk assessment is derived from currently known pathogenic and likely-pathogenic variants associated with specific cancer types.
          </p>
          <p style={{ marginBottom: '1rem', fontWeight: '600', color: '#374151' }}>
            Please consider the following limitations:
          </p>
          <ul style={{ paddingLeft: '1.25rem', margin: 0 }}>
            <li style={{ marginBottom: '1rem' }}>
              <strong>Scope of Analysis:</strong> This report is based exclusively on the variant calls from the provided data file. It does not incorporate essential clinical information such as patient demographics, family history, lifestyle, or environmental factors, all of which are critical for a comprehensive risk assessment.
            </li>
            <li style={{ marginBottom: '1rem' }}>
              <strong>Model Knowledge:</strong> The underlying risk model is limited to genetic associations that are well-documented in our training data. Its ability to accurately assess risk from novel, rare, or uncharacterized genetic variants is not validated.
            </li>
            <li style={{ marginBottom: '1rem' }}>
              <strong>Data Quality:</strong> The accuracy of this report is directly dependent on the quality of the input data. Sequencing artifacts, low read depth, or errors in the provided VCF/BAM/FASTQ file can impact the results.
            </li>
          </ul>
          <p style={{ 
            marginTop: '1.5rem', 
            paddingTop: '1rem',
            borderTop: '1px solid #E5E7EB',
            fontStyle: 'italic',
            fontSize: '0.8rem',
            color: '#6B7280'
          }}>
            For Investigational Use Only. All findings must be correlated with clinical findings and confirmed by an independent, certified laboratory.
          </p>
        </div>
      )}
    </div>
  );
};

export default DetailedDisclaimer; 