import type { ReactNode, FC } from 'react'
import { Link } from 'react-router-dom'
import { invoke } from '@tauri-apps/api/core'

const GeneKnowLogo = () => (
  <svg className="w-8 h-8" style={{ color: '#0066cc' }} viewBox="0 0 240 240" fill="none" xmlns="http://www.w3.org/2000/svg">
    <defs>
      <linearGradient id="dnaGradient1" x1="0%" y1="0%" x2="0%" y2="100%">
        <stop offset="0%" style={{stopColor:'#0066cc', stopOpacity:1}} />
        <stop offset="100%" style={{stopColor:'#003d7a', stopOpacity:1}} />
      </linearGradient>
      <linearGradient id="dnaGradient2" x1="0%" y1="0%" x2="0%" y2="100%">
        <stop offset="0%" style={{stopColor:'#00a8e6', stopOpacity:1}} />
        <stop offset="100%" style={{stopColor:'#0066cc', stopOpacity:1}} />
      </linearGradient>
    </defs>
    <g transform="translate(120, 120)">
      <g>
        <path d="M -25,-70 C -25,-50 25,-40 25,-20 C 25,0 -25,10 -25,30 C -25,50 25,60 25,80" 
              fill="none" stroke="url(#dnaGradient1)" strokeWidth="5" strokeLinecap="round"/>
        <path d="M 25,-70 C 25,-50 -25,-40 -25,-20 C -25,0 25,10 25,30 C 25,50 -25,60 -25,80" 
              fill="none" stroke="url(#dnaGradient2)" strokeWidth="5" strokeLinecap="round"/>
        <line x1="-25" y1="-60" x2="25" y2="-60" stroke="url(#dnaGradient1)" strokeWidth="3"/>
        <line x1="-15" y1="-45" x2="15" y2="-45" stroke="url(#dnaGradient2)" strokeWidth="3"/>
        <line x1="-25" y1="-30" x2="25" y2="-30" stroke="url(#dnaGradient1)" strokeWidth="3"/>
        <line x1="-15" y1="-15" x2="15" y2="-15" stroke="url(#dnaGradient2)" strokeWidth="3"/>
        <line x1="-25" y1="0" x2="25" y2="0" stroke="url(#dnaGradient1)" strokeWidth="3"/>
        <line x1="-15" y1="15" x2="15" y2="15" stroke="url(#dnaGradient2)" strokeWidth="3"/>
        <line x1="-25" y1="30" x2="25" y2="30" stroke="url(#dnaGradient1)" strokeWidth="3"/>
        <line x1="-15" y1="45" x2="15" y2="45" stroke="url(#dnaGradient2)" strokeWidth="3"/>
        <line x1="-25" y1="60" x2="25" y2="60" stroke="url(#dnaGradient1)" strokeWidth="3"/>
      </g>
    </g>
  </svg>
)

const Header = () => (
  <header style={{
    background: '#FFFFFF',
    backdropFilter: 'blur(10px)',
    position: 'fixed',
    top: 0,
    left: 0,
    right: 0,
    zIndex: 1030,
    borderBottom: '1px solid #E5E7EB',
    padding: '1rem 0',
    boxShadow: '0 1px 2px 0 rgba(0, 0, 0, 0.05)'
  }}>
    <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
      <Link to="/" style={{ 
        textDecoration: 'none', 
        display: 'flex', 
        alignItems: 'center', 
        gap: '0.5rem' 
      }}>
        <GeneKnowLogo />
        <h1 style={{ 
          fontSize: '1.5rem', 
          fontWeight: 'bold', 
          color: '#1F2937', 
          margin: 0 
        }}>GeneKnow</h1>
      </Link>
      <nav style={{ display: 'flex', alignItems: 'center', gap: '2rem' }}>
        <Link to="/how-it-works" style={{ 
          color: '#4B5563', 
          textDecoration: 'none', 
          transition: 'color 200ms ease',
          fontWeight: '500'
        }}
        onMouseEnter={(e) => e.currentTarget.style.color = '#1F2937'}
        onMouseLeave={(e) => e.currentTarget.style.color = '#4B5563'}>
          How It Works
        </Link>
        <Link to="/features" style={{ 
          color: '#4B5563', 
          textDecoration: 'none', 
          transition: 'color 200ms ease',
          fontWeight: '500'
        }}
        onMouseEnter={(e) => e.currentTarget.style.color = '#1F2937'}
        onMouseLeave={(e) => e.currentTarget.style.color = '#4B5563'}>
          Features
        </Link>
        <Link to="/privacy" style={{ 
          color: '#4B5563', 
          textDecoration: 'none', 
          transition: 'color 200ms ease',
          fontWeight: '500'
        }}
        onMouseEnter={(e) => e.currentTarget.style.color = '#1F2937'}
        onMouseLeave={(e) => e.currentTarget.style.color = '#4B5563'}>
          Privacy First
        </Link>
      </nav>
    </div>
  </header>
)

const Footer = () => (
  <footer style={{ background: 'var(--gray-800)', color: '#FFFFFF', padding: '2rem 0' }}>
    <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
      {/* Row with copyright & disclaimer on left and GitHub button on right */}
      <div style={{ 
        display: 'flex', 
        justifyContent: 'space-between', 
        alignItems: 'center',
        flexWrap: 'wrap',
        gap: '1rem'
      }}>
        {/* Copyright & Disclaimer - left side */}
        <div style={{ flex: '1', minWidth: '300px' }}>
          <p style={{ 
            fontSize: '0.875rem', 
            color: 'var(--gray-400)', 
            margin: '0 0 0.5rem 0',
            textAlign: 'left'
          }}>
            &copy; 2025 GeneKnow. All rights reserved.
          </p>
          <p style={{ 
            fontSize: '0.875rem', 
            color: 'var(--gray-400)', 
            margin: 0,
            textAlign: 'left'
          }}>
            Disclaimer: This tool is for informational purposes only and is not a substitute for professional medical advice, diagnosis, or treatment.
          </p>
        </div>
        
        {/* GitHub Link - right side */}
        <div style={{ flexShrink: 0 }}>
          <button 
            onClick={async () => {
              try {
                await invoke('plugin:shell|open', {
                  path: 'https://github.com/astrogirlnim/GeneKnow'
                });
              } catch (error) {
                console.error('Failed to open GitHub link:', error);
              }
            }}
            style={{
              display: 'inline-flex',
              alignItems: 'center',
              gap: '0.5rem',
              color: 'var(--gray-400)',
              background: 'none',
              border: 'none',
              padding: '0',
              cursor: 'pointer',
              fontSize: '0.875rem',
              fontWeight: '400',
              transition: 'color 200ms ease'
            }}
            onMouseEnter={(e) => {
              e.currentTarget.style.color = '#FFFFFF';
            }}
            onMouseLeave={(e) => {
              e.currentTarget.style.color = 'var(--gray-400)';
            }}
          >
            <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
              <path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/>
            </svg>
            View on GitHub
          </button>
        </div>
      </div>
    </div>
  </footer>
)

interface LayoutProps {
  children: ReactNode
}

const Layout: FC<LayoutProps> = ({ children }) => (
  <div style={{ fontFamily: 'var(--font-family-sans)', background: '#FFFFFF', color: 'var(--gray-800)', minHeight: '100vh' }}>
    <Header />
    <main style={{ paddingTop: '4rem', minHeight: 'calc(100vh - 4rem)' }}>
      {children}
    </main>
    <Footer />
  </div>
)

export default Layout 