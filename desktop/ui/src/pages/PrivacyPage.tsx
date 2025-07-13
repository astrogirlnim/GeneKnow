import Layout from '../components/Layout'

const CheckBadgeIcon = () => (
  <div style={{
    width: '16px',
    height: '16px',
    backgroundColor: '#22C55E',
    borderRadius: '50%',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    marginRight: '8px',
    marginTop: '2px',
    flexShrink: 0
  }}>
    <svg width="10" height="8" viewBox="0 0 10 8" fill="none">
      <path d="M1 4L3.5 6.5L9 1" stroke="#FFFFFF" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round"/>
    </svg>
  </div>
)

const PrivacyPage = () => (
  <Layout>
    <section style={{ 
      background: 'radial-gradient(circle at top left, rgba(239, 246, 255, 1) 0%, rgba(255, 255, 255, 1) 50%)', 
      padding: '3rem 0 2rem 0' 
    }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem', textAlign: 'center' }}>
        <h1 style={{
          fontSize: 'clamp(2.5rem, 5vw, 3.75rem)',
          fontWeight: 'bold',
          letterSpacing: '-0.02em',
          color: '#111827',
          lineHeight: '1.1',
          marginBottom: '1rem'
        }}>
          Privacy First
        </h1>
        <p style={{
          fontSize: '1.125rem',
          lineHeight: '1.75',
          color: '#4B5563',
          maxWidth: '42rem',
          margin: '0 auto'
        }}>
          Your genomic data is the most personal information you have. We believe it should stay private and under your complete control.
        </p>
      </div>
    </section>

    <section style={{ background: '#FFFFFF', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '4rem', flexWrap: 'wrap' }}>
          <div style={{ flex: '1', minWidth: '300px' }}>
            <span style={{
              color: '#2563EB',
              fontWeight: '600',
              background: '#DBEAFE',
              borderRadius: '9999px',
              padding: '0.5rem 1rem',
              fontSize: '0.875rem',
              display: 'inline-block'
            }}>Uncompromising Security</span>
            <h2 style={{ marginTop: '1rem', fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: '#111827' }}>
              Your Data is Your Own. Period.
            </h2>
            <p style={{ marginTop: '1.5rem', fontSize: '1.125rem', color: '#4B5563', lineHeight: '1.75' }}>
              In an age of data breaches, we believe your most personal information should remain in your hands. GeneKnow is built on a foundation of privacy-by-design. By processing everything locally, we eliminate the risk of server-side breaches and data misuse.
            </p>
            <ul style={{ marginTop: '1.5rem', listStyle: 'none', padding: 0 }}>
              <li style={{ display: 'flex', alignItems: 'flex-start', marginBottom: '1rem' }}>
                <CheckBadgeIcon />
                <span style={{ color: '#4B5563' }}>
                  <strong>No Cloud Uploads:</strong> Your genomic file is never sent to us or any third party.
                </span>
              </li>
              <li style={{ display: 'flex', alignItems: 'flex-start' }}>
                <CheckBadgeIcon />
                <span style={{ color: '#4B5563' }}>
                  <strong>Secure by Design:</strong> Built with Rust and Tauri for a secure, sandboxed application environment.
                </span>
              </li>
            </ul>
          </div>
          <div style={{ flex: '1', minWidth: '300px' }}>
            <div style={{
              width: '100%',
              height: '300px',
              borderRadius: '0.75rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              overflow: 'hidden',
              boxShadow: '0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)'
            }}>
              <img 
                src="/privacy-first.png" 
                alt="Privacy First - Secure Local Analysis"
                style={{
                  width: '100%',
                  height: '100%',
                  objectFit: 'cover',
                  borderRadius: '0.75rem'
                }}
              />
            </div>
          </div>
        </div>
      </div>
    </section>

    <section style={{ background: '#F9FAFB', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
          <h2 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: '#111827' }}>
            Privacy Guarantees
          </h2>
          <p style={{ marginTop: '1rem', fontSize: '1.125rem', color: '#4B5563', maxWidth: '32rem', margin: '1rem auto 0' }}>
            Our commitment to your privacy goes beyond just keeping your data local.
          </p>
        </div>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', gap: '2rem' }}>
          <div style={{ 
            textAlign: 'center',
            background: '#FFFFFF',
            border: '1px solid #E5E7EB',
            borderRadius: '0.75rem',
            padding: '2rem',
            boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)'
          }}>
            <div style={{
              margin: '0 auto 1.25rem',
              width: '3rem',
              height: '3rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              borderRadius: '50%',
              background: '#D1FAE5',
              fontSize: '1.5rem'
            }}>
              🔒
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
              Zero Data Collection
            </h3>
            <p style={{ color: '#4B5563' }}>
              We don't collect, store, or analyze any of your genomic data. Everything stays on your device.
            </p>
          </div>
          <div style={{ 
            textAlign: 'center',
            background: '#FFFFFF',
            border: '1px solid #E5E7EB',
            borderRadius: '0.75rem',
            padding: '2rem',
            boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)'
          }}>
            <div style={{
              margin: '0 auto 1.25rem',
              width: '3rem',
              height: '3rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              borderRadius: '50%',
              background: '#D1FAE5',
              fontSize: '1.5rem'
            }}>
              🛡️
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
              No Third-Party Sharing
            </h3>
            <p style={{ color: '#4B5563' }}>
              Your data never leaves your computer, so there's no risk of unauthorized sharing or breaches.
            </p>
          </div>
          <div style={{ 
            textAlign: 'center',
            background: '#FFFFFF',
            border: '1px solid #E5E7EB',
            borderRadius: '0.75rem',
            padding: '2rem',
            boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)'
          }}>
            <div style={{
              margin: '0 auto 1.25rem',
              width: '3rem',
              height: '3rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              borderRadius: '50%',
              background: '#D1FAE5',
              fontSize: '1.5rem'
            }}>
              🔐
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
              Offline Processing
            </h3>
            <p style={{ color: '#4B5563' }}>
              Analysis happens completely offline. No internet connection required during processing.
            </p>
          </div>
        </div>
      </div>
    </section>

    <section style={{ background: '#FFFFFF', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
          <h2 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: '#111827' }}>
            Technical Privacy Implementation
          </h2>
        </div>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(350px, 1fr))', gap: '3rem' }}>
          <div>
            <h3 style={{ fontSize: '1.5rem', fontWeight: '600', color: '#111827', marginBottom: '1rem' }}>
              Local-Only Architecture
            </h3>
            <ul style={{ listStyle: 'none', padding: 0 }}>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                💻 All processing happens on your local machine
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                🚫 No network requests during analysis
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                🔒 Data never leaves your device
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                📱 Works completely offline
              </li>
            </ul>
          </div>
          <div>
            <h3 style={{ fontSize: '1.5rem', fontWeight: '600', color: '#111827', marginBottom: '1rem' }}>
              Security Measures
            </h3>
            <ul style={{ listStyle: 'none', padding: 0 }}>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                🔐 Memory-safe Rust implementation
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                🛡️ Sandboxed Tauri environment
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                🔒 Encrypted temporary files
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                🗑️ Automatic cleanup of processing artifacts
              </li>
            </ul>
          </div>
        </div>
      </div>
    </section>
  </Layout>
)

export default PrivacyPage 