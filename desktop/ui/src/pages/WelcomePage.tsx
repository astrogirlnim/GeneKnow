import { useNavigate } from 'react-router-dom'
import Layout from '../components/Layout'

const WelcomePage = () => {
  const navigate = useNavigate()

  const handleTestNow = () => {
    navigate('/upload')
  }

  return (
    <Layout>
      <section style={{
        background: 'radial-gradient(circle at top left, rgba(239, 246, 255, 1) 0%, rgba(255, 255, 255, 1) 50%)',
        padding: '1rem 0 1rem 0',
        height: 'calc(100vh - 10rem)',
        minHeight: '400px'
      }}>
        <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem', textAlign: 'center', height: '100%', display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
          <div style={{ maxWidth: '48rem', margin: '0 auto' }}>
            <span style={{
              color: '#2563EB',
              fontWeight: '600',
              background: '#DBEAFE',
              borderRadius: '9999px',
              padding: '0.4rem 0.8rem',
              fontSize: '0.8rem',
              display: 'inline-block'
            }}>Your Personal Genomic Insights</span>
            <h1 style={{
              marginTop: '0.8rem',
              fontSize: 'clamp(1.8rem, 4vw, 2.5rem)',
              fontWeight: 'bold',
              letterSpacing: '-0.02em',
              color: '#111827',
              lineHeight: '1.1'
            }}>
              Understand Your Genomic Health, Instantly and Privately.
            </h1>
            <p style={{
              marginTop: '0.8rem',
              fontSize: '1rem',
              lineHeight: '1.6',
              color: '#4B5563',
              maxWidth: '42rem',
              margin: '0.8rem auto 0'
            }}>
              Get a private cancer risk assessment by analyzing your genomic data securely on your own device. Start by selecting a file—your data never leaves your computer.
            </p>
            <div style={{ marginTop: '1.5rem' }}>
              <button 
                onClick={handleTestNow}
                style={{
                  backgroundColor: '#2563EB',
                  color: '#FFFFFF',
                  border: 'none',
                  borderRadius: '0.5rem',
                  padding: '0.75rem 2rem',
                  fontSize: '1rem',
                  fontWeight: '600',
                  cursor: 'pointer',
                  transition: 'all 200ms ease',
                  boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)'
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.backgroundColor = '#1D4ED8';
                  e.currentTarget.style.transform = 'translateY(-1px)';
                  e.currentTarget.style.boxShadow = '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.backgroundColor = '#2563EB';
                  e.currentTarget.style.transform = 'translateY(0)';
                  e.currentTarget.style.boxShadow = '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)';
                }}>
                Test Now
              </button>
            </div>
          </div>
        </div>
      </section>
    </Layout>
  )
}

export default WelcomePage 