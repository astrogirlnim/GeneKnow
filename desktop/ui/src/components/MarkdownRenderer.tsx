import React from 'react';
import ReactMarkdown from 'react-markdown';

interface MarkdownRendererProps {
  content: string;
  className?: string;
}

const MarkdownRenderer: React.FC<MarkdownRendererProps> = ({ content, className = '' }) => {
  return (
    <div className={`markdown-content ${className}`} style={{
      fontFamily: 'ui-serif, Georgia, Cambria, "Times New Roman", Times, serif',
      lineHeight: '1.7',
      color: '#374151',
      maxWidth: 'none'
    }}>
      <ReactMarkdown
        components={{
          // Custom heading styles for medical reports
          h1: ({ children }) => (
            <h1 style={{
              fontSize: '2rem',
              fontWeight: '700',
              color: '#111827',
              marginBottom: '1.5rem',
              marginTop: '2rem',
              borderBottom: '2px solid #E5E7EB',
              paddingBottom: '0.5rem'
            }}>
              {children}
            </h1>
          ),
          h2: ({ children }) => (
            <h2 style={{
              fontSize: '1.5rem',
              fontWeight: '600',
              color: '#1F2937',
              marginBottom: '1rem',
              marginTop: '2rem'
            }}>
              {children}
            </h2>
          ),
          h3: ({ children }) => (
            <h3 style={{
              fontSize: '1.25rem',
              fontWeight: '600',
              color: '#374151',
              marginBottom: '0.75rem',
              marginTop: '1.5rem'
            }}>
              {children}
            </h3>
          ),
          h4: ({ children }) => (
            <h4 style={{
              fontSize: '1.125rem',
              fontWeight: '600',
              color: '#4B5563',
              marginBottom: '0.5rem',
              marginTop: '1rem'
            }}>
              {children}
            </h4>
          ),
          // Paragraph styling
          p: ({ children }) => (
            <p style={{
              marginBottom: '1rem',
              lineHeight: '1.7'
            }}>
              {children}
            </p>
          ),
          // Strong text styling for medical terms
          strong: ({ children }) => (
            <strong style={{
              fontWeight: '600',
              color: '#1F2937'
            }}>
              {children}
            </strong>
          ),
          // Emphasis styling
          em: ({ children }) => (
            <em style={{
              fontStyle: 'italic',
              color: '#4B5563'
            }}>
              {children}
            </em>
          ),
          // List styling
          ul: ({ children }) => (
            <ul style={{
              marginBottom: '1rem',
              paddingLeft: '1.5rem',
              listStyleType: 'disc'
            }}>
              {children}
            </ul>
          ),
          ol: ({ children }) => (
            <ol style={{
              marginBottom: '1rem',
              paddingLeft: '1.5rem',
              listStyleType: 'decimal'
            }}>
              {children}
            </ol>
          ),
          li: ({ children }) => (
            <li style={{
              marginBottom: '0.5rem',
              lineHeight: '1.6'
            }}>
              {children}
            </li>
          ),
          // Table styling for medical data
          table: ({ children }) => (
            <div style={{ overflowX: 'auto', marginBottom: '1.5rem' }}>
              <table style={{
                width: '100%',
                borderCollapse: 'collapse',
                border: '1px solid #D1D5DB',
                borderRadius: '0.5rem',
                overflow: 'hidden'
              }}>
                {children}
              </table>
            </div>
          ),
          thead: ({ children }) => (
            <thead style={{
              backgroundColor: '#F9FAFB'
            }}>
              {children}
            </thead>
          ),
          th: ({ children }) => (
            <th style={{
              padding: '0.75rem',
              textAlign: 'left',
              fontWeight: '600',
              color: '#374151',
              borderBottom: '1px solid #D1D5DB'
            }}>
              {children}
            </th>
          ),
          td: ({ children }) => (
            <td style={{
              padding: '0.75rem',
              borderBottom: '1px solid #E5E7EB'
            }}>
              {children}
            </td>
          ),
          // Code styling for genetic sequences
          code: ({ children, ...props }) => (
            'inline' in props && props.inline ? (
              <code style={{
                backgroundColor: '#F3F4F6',
                padding: '0.125rem 0.25rem',
                borderRadius: '0.25rem',
                fontSize: '0.875rem',
                fontFamily: 'ui-monospace, SFMono-Regular, "SF Mono", Monaco, Consolas, "Liberation Mono", "Courier New", monospace',
                color: '#1F2937'
              }}>
                {children}
              </code>
            ) : (
              <pre style={{
                backgroundColor: '#F9FAFB',
                padding: '1rem',
                borderRadius: '0.5rem',
                border: '1px solid #E5E7EB',
                overflowX: 'auto',
                marginBottom: '1rem'
              }}>
                <code style={{
                  fontSize: '0.875rem',
                  fontFamily: 'ui-monospace, SFMono-Regular, "SF Mono", Monaco, Consolas, "Liberation Mono", "Courier New", monospace',
                  color: '#374151'
                }}>
                  {children}
                </code>
              </pre>
            )
          ),
          // Blockquote styling for important medical notes
          blockquote: ({ children }) => (
            <blockquote style={{
              borderLeft: '4px solid #3B82F6',
              paddingLeft: '1rem',
              marginLeft: '0',
              marginBottom: '1rem',
              backgroundColor: '#EFF6FF',
              padding: '1rem',
              borderRadius: '0.5rem',
              fontStyle: 'italic'
            }}>
              {children}
            </blockquote>
          ),
          // Horizontal rule styling
          hr: () => (
            <hr style={{
              border: 'none',
              borderTop: '1px solid #E5E7EB',
              margin: '2rem 0'
            }} />
          )
        }}
      >
        {content}
      </ReactMarkdown>
    </div>
  );
};

export default MarkdownRenderer; 