import type { ReactNode, FC } from 'react'
import { Link } from 'react-router-dom'

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
  <footer style={{ background: 'var(--gray-800)', color: '#FFFFFF', padding: '2rem 0', textAlign: 'center' }}>
    <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
      <p>&copy; 2025 GeneKnow. All rights reserved.</p>
      <p style={{ fontSize: '0.875rem', color: 'var(--gray-400)', marginTop: '0.5rem' }}>
        Disclaimer: This tool is for informational purposes only and is not a substitute for professional medical advice, diagnosis, or treatment.
      </p>
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