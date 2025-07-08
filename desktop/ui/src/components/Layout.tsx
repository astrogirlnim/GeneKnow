import React from 'react'
import { Link } from 'react-router-dom'

const GeneKnowLogo = () => (
  <svg className="w-8 h-8" style={{ color: '#2563EB' }} viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
    <path d="M12 21C16.9706 21 21 16.9706 21 12C21 7.02944 16.9706 3 12 3C7.02944 3 3 7.02944 3 12C3 16.9706 7.02944 21 12 21Z" stroke="currentColor" strokeWidth="1.5"/>
    <g stroke="currentColor" strokeWidth="1.75" strokeLinecap="round" strokeLinejoin="round">
      <path d="M15.18 8.5C12.08 10.33 11.08 13.67 13.18 16.5"/>
      <path d="M8.82 7.5C11.92 9.33 12.92 12.67 10.82 15.5"/>
      <path d="M14.2 10.5H9.8"/>
      <path d="M15 12.5H9"/>
      <path d="M13.2 14.5H10.8"/>
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
  children: React.ReactNode
}

const Layout: React.FC<LayoutProps> = ({ children }) => (
  <div style={{ fontFamily: 'var(--font-family-sans)', background: '#FFFFFF', color: 'var(--gray-800)', minHeight: '100vh' }}>
    <Header />
    <main style={{ paddingTop: '4rem', minHeight: 'calc(100vh - 4rem)' }}>
      {children}
    </main>
    <Footer />
  </div>
)

export default Layout 