import React, { useState } from 'react';
import type { SHAPValidation, SHAPContributor } from '../api/geneknowPipeline';

interface ConfidenceCheckProps {
  validation: SHAPValidation | null;
  onNavigateToDetail?: () => void;
  isDetailed?: boolean;
}

const ConfidenceCheck: React.FC<ConfidenceCheckProps> = ({ 
  validation, 
  onNavigateToDetail, 
  isDetailed = false 
}) => {
  const [isAcknowledged, setIsAcknowledged] = useState(false);
  const [showTooltip, setShowTooltip] = useState(false);
  const [showDetails, setShowDetails] = useState(false);

  if (!validation) {
    return null;
  }

  const getStateConfig = (status: string) => {
    switch (status) {
      case 'FLAG_FOR_REVIEW':
        return {
          icon: (
            <div style={{
              width: '0.5rem',
              height: '0.5rem',
              borderRadius: '50%',
              background: '#F59E0B'
            }}></div>
          ),
          title: 'Confidence Check: Review Required',
          bgColor: '#FFFBEB',
          borderColor: '#FDE68A',
          textColor: '#92400E',
          buttonColor: '#F59E0B',
          buttonText: 'Review & Acknowledge',
          buttonStyle: 'primary' as const
        };
      case 'PASS':
        return {
          icon: (
            <div style={{
              width: '0.5rem',
              height: '0.5rem',
              borderRadius: '50%',
              background: '#22C55E'
            }}></div>
          ),
          title: 'Confidence Check: Passed',
          bgColor: '#F0FDF4',
          borderColor: '#BBF7D0',
          textColor: '#166534',
          buttonColor: '#22C55E',
          buttonText: 'Acknowledge',
          buttonStyle: 'secondary' as const
        };
      case 'ERROR':
        return {
          icon: (
            <div style={{
              width: '0.5rem',
              height: '0.5rem',
              borderRadius: '50%',
              background: '#EF4444'
            }}></div>
          ),
          title: 'Confidence Check: Error',
          bgColor: '#FEF2F2',
          borderColor: '#FECACA',
          textColor: '#991B1B',
          buttonColor: '#EF4444',
          buttonText: 'Acknowledge & Continue',
          buttonStyle: 'secondary' as const
        };
      case 'SKIPPED':
        return {
          icon: (
            <div style={{
              width: '0.5rem',
              height: '0.5rem',
              borderRadius: '50%',
              background: '#6B7280'
            }}></div>
          ),
          title: 'Confidence Check: Not Applicable',
          bgColor: '#F9FAFB',
          borderColor: '#E5E7EB',
          textColor: '#4B5563',
          buttonColor: '#6B7280',
          buttonText: 'Acknowledge',
          buttonStyle: 'secondary' as const
        };
      default:
        return {
          icon: (
            <div style={{
              width: '0.5rem',
              height: '0.5rem',
              borderRadius: '50%',
              background: '#6B7280'
            }}></div>
          ),
          title: 'Confidence Check: Unknown',
          bgColor: '#F9FAFB',
          borderColor: '#E5E7EB',
          textColor: '#4B5563',
          buttonColor: '#6B7280',
          buttonText: 'Acknowledge',
          buttonStyle: 'secondary' as const
        };
    }
  };

  const config = getStateConfig(validation.status);

  const handleAcknowledge = () => {
    console.log(`Confidence Check acknowledged: ${validation.status}`);
    setIsAcknowledged(true);
    setShowDetails(false); // Hide details when acknowledging
  };

  const handleShowDetails = () => {
    setShowDetails(true);
  };

  const handleNavigateToDetail = () => {
    if (onNavigateToDetail) {
      onNavigateToDetail();
    }
  };

  const renderDescription = () => {
    switch (validation.status) {
      case 'FLAG_FOR_REVIEW':
        return (
          <div>
            <p style={{ fontSize: '0.875rem', color: config.textColor, marginBottom: '0.75rem' }}>
              {validation.reasons.length > 0 ? validation.reasons[0] : 
                'The AI\'s prediction may need additional review based on the genetic evidence.'}
            </p>
            {validation.top_contributors.length > 0 && (
              <div style={{ marginTop: '0.75rem' }}>
                <p style={{ fontSize: '0.75rem', fontWeight: '600', color: config.textColor, marginBottom: '0.5rem' }}>
                  Top Contributing Factors:
                </p>
                <ul style={{ fontSize: '0.75rem', color: config.textColor, listStyle: 'none', padding: 0 }}>
                  {validation.top_contributors.slice(0, 3).map((contributor, index) => (
                    <li key={index} style={{ marginBottom: '0.25rem' }}>
                      • {contributor.display_name} ({contributor.direction} risk)
                    </li>
                  ))}
                </ul>
              </div>
            )}
          </div>
        );
      case 'PASS':
        return (
          <p style={{ fontSize: '0.875rem', color: config.textColor }}>
            The AI's reasoning aligns with the primary genomic evidence. The risk calculation is based on 
            clinically significant factors.
          </p>
        );
      case 'ERROR':
        return (
          <p style={{ fontSize: '0.875rem', color: config.textColor }}>
            A technical issue prevented the confidence check from completing. Your primary analysis is 
            unaffected, but please exercise extra caution when interpreting results.
          </p>
        );
      case 'SKIPPED':
        return (
          <p style={{ fontSize: '0.875rem', color: config.textColor }}>
            The confidence check was not designed for this type of analysis. No additional review 
            is required from a validation standpoint.
          </p>
        );
      default:
        return (
          <p style={{ fontSize: '0.875rem', color: config.textColor }}>
            Unknown validation status. Please contact support if this issue persists.
          </p>
        );
    }
  };

  const renderTooltip = () => (
    <div 
      style={{
        position: 'absolute',
        left: isDetailed ? 'calc(100% + 0.5rem)' : '50%',
        top: isDetailed ? '50%' : 'calc(100% + 0.5rem)',
        transform: isDetailed ? 'translateY(-50%)' : 'translateX(-50%)',
        width: '16rem',
        padding: '0.75rem',
        background: '#1F2937',
        color: '#FFFFFF',
        fontSize: '0.75rem',
        borderRadius: '0.5rem',
        boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
        opacity: showTooltip ? 1 : 0,
        visibility: showTooltip ? 'visible' : 'hidden',
        transition: 'opacity 300ms ease, visibility 300ms ease',
        zIndex: 1000,
        pointerEvents: 'none',
        lineHeight: '1.4',
        border: '1px solid #374151'
      }}
    >
      This automated check verifies that the model's risk calculation is based on the most 
      clinically significant genomic evidence.
      
      {/* Arrow */}
      <div style={{
        position: 'absolute',
        [isDetailed ? 'left' : 'top']: isDetailed ? '-0.5rem' : '-0.5rem',
        [isDetailed ? 'top' : 'left']: isDetailed ? '50%' : '50%',
        transform: isDetailed ? 'translateY(-50%)' : 'translateX(-50%)',
        width: '0',
        height: '0',
        [isDetailed ? 'borderTop' : 'borderLeft']: '0.5rem solid transparent',
        [isDetailed ? 'borderBottom' : 'borderRight']: '0.5rem solid transparent',
        [isDetailed ? 'borderRight' : 'borderBottom']: '0.5rem solid #1F2937'
      }}></div>
    </div>
  );

  // Summary view for dashboard
  if (!isDetailed) {
    return (
      <div 
        style={{
          cursor: 'default',
          padding: '0.75rem',
          marginTop: '0.75rem',
          background: config.bgColor,
          border: `1px solid ${config.borderColor}`,
          borderRadius: '0.5rem',
          fontSize: '0.875rem',
          color: config.textColor,
          transition: 'all 200ms ease',
          position: 'relative'
        }}
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
          <span style={{ fontSize: '1rem' }}>{config.icon}</span>
          <span style={{ fontWeight: '600' }}>{config.title}</span>
          <div 
            style={{ position: 'relative', marginLeft: 'auto' }}
            onMouseEnter={() => setShowTooltip(true)}
            onMouseLeave={() => setShowTooltip(false)}
          >
            <div
              style={{
                width: '1rem',
                height: '1rem',
                backgroundColor: '#E5E7EB',
                borderRadius: '50%',
                border: '1px solid #D1D5DB',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                fontSize: '0.625rem',
                fontWeight: '600',
                color: '#6B7280',
                cursor: 'help',
                opacity: 0.7,
                transition: 'all 200ms ease'
              }}
            >
              i
            </div>
            {renderTooltip()}
          </div>
        </div>
        
        <div style={{ marginTop: '0.5rem', fontSize: '0.75rem' }}>
          See In-Depth Analysis for details.
        </div>
      </div>
    );
  }

  // Detailed view for Clinical Alerts
  return (
    <div 
      style={{
        background: isAcknowledged && !showDetails ? '#F0FDF4' : config.bgColor, // Green background when reviewed and not showing details
        border: `1px solid ${isAcknowledged && !showDetails ? '#BBF7D0' : config.borderColor}`, // Green border when reviewed and not showing details
        borderLeft: `3px solid ${isAcknowledged && !showDetails ? '#22C55E' : config.buttonColor}`, // Green left border when reviewed and not showing details
        padding: '1rem',
        marginBottom: '1rem',
        borderRadius: '0.5rem',
        fontSize: '0.875rem',
        position: 'relative',
        transition: 'all 300ms ease'
      }}
    >
      <div style={{
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between',
        marginBottom: '0.75rem'
      }}>
        <div style={{
          display: 'flex',
          alignItems: 'center',
          gap: '0.5rem'
        }}>
          <span style={{ fontSize: '1.125rem' }}>
            {isAcknowledged && !showDetails ? (
              <div style={{
                width: '0.5rem',
                height: '0.5rem',
                borderRadius: '50%',
                background: '#22C55E'
              }}></div>
            ) : config.icon}
          </span>
          <strong style={{ fontSize: '0.875rem', color: isAcknowledged && !showDetails ? '#166534' : '#111827' }}>
            {isAcknowledged && !showDetails ? 'Confidence Check: Reviewed' : config.title}
          </strong>
        </div>
        
        <div 
          style={{ position: 'relative', zIndex: 10, overflow: 'visible' }}
          onMouseEnter={() => setShowTooltip(true)}
          onMouseLeave={() => setShowTooltip(false)}
        >
          <div
            style={{
              width: '1rem',
              height: '1rem',
              backgroundColor: '#E5E7EB',
              borderRadius: '50%',
              border: '1px solid #D1D5DB',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              fontSize: '0.625rem',
              fontWeight: '600',
              color: '#6B7280',
              cursor: 'help',
              transition: 'all 200ms ease'
            }}
            onMouseEnter={(e) => {
              e.currentTarget.style.backgroundColor = '#D1D5DB';
              e.currentTarget.style.color = '#374151';
              e.currentTarget.style.borderColor = '#9CA3AF';
            }}
            onMouseLeave={(e) => {
              e.currentTarget.style.backgroundColor = '#E5E7EB';
              e.currentTarget.style.color = '#6B7280';
              e.currentTarget.style.borderColor = '#D1D5DB';
            }}
          >
            i
          </div>
          {renderTooltip()}
        </div>
      </div>
      
      <div style={{ marginBottom: '1rem' }}>
        {isAcknowledged && !showDetails ? (
          <p style={{ fontSize: '0.875rem', color: '#4B5563', fontWeight: 'normal' }}>
            This analysis has been reviewed and confirmed by the user. The AI's risk assessment 
            and contributing factors have been acknowledged.
          </p>
        ) : (
          renderDescription()
        )}
      </div>
      
      <div style={{ display: 'flex', justifyContent: 'flex-end', gap: '0.5rem' }}>
        {isAcknowledged && !showDetails ? (
          <button
            onClick={handleShowDetails}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: '#FFFFFF',
              color: '#6B7280',
              border: '1px solid #D1D5DB',
              borderRadius: '0.375rem',
              fontSize: '0.875rem',
              fontWeight: '500',
              cursor: 'pointer',
              transition: 'all 200ms ease'
            }}
            onMouseEnter={(e) => {
              e.currentTarget.style.backgroundColor = '#F9FAFB';
              e.currentTarget.style.borderColor = '#9CA3AF';
            }}
            onMouseLeave={(e) => {
              e.currentTarget.style.backgroundColor = '#FFFFFF';
              e.currentTarget.style.borderColor = '#D1D5DB';
            }}
          >
            See Details
          </button>
        ) : (
          <button
            onClick={handleAcknowledge}
            style={{
              padding: '0.5rem 1rem',
              backgroundColor: config.buttonStyle === 'primary' ? config.buttonColor : '#FFFFFF',
              color: config.buttonStyle === 'primary' ? '#FFFFFF' : config.buttonColor,
              border: config.buttonStyle === 'primary' ? 'none' : `1px solid ${config.buttonColor}`,
              borderRadius: '0.375rem',
              fontSize: '0.875rem',
              fontWeight: '500',
              cursor: 'pointer',
              transition: 'all 200ms ease'
            }}
            onMouseEnter={(e) => {
              if (config.buttonStyle === 'primary') {
                e.currentTarget.style.backgroundColor = config.buttonColor + 'DD';
              } else {
                e.currentTarget.style.backgroundColor = config.buttonColor + '10';
              }
            }}
            onMouseLeave={(e) => {
              if (config.buttonStyle === 'primary') {
                e.currentTarget.style.backgroundColor = config.buttonColor;
              } else {
                e.currentTarget.style.backgroundColor = '#FFFFFF';
              }
            }}
          >
            {config.buttonText}
          </button>
        )}
      </div>
    </div>
  );
};

export default ConfidenceCheck; 