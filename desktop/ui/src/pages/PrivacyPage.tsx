import React from 'react'
import Layout from '../components/Layout'

const CheckBadgeIcon = () => (
  <svg className="h-6 w-6" style={{ color: 'var(--success)', marginRight: '12px' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M9 12.75L11.25 15 15 9.75m-3-7.036A11.959 11.959 0 013.598 6 11.99 11.99 0 003 9.749c0 5.592 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.31-.21-2.571-.598-3.751h-.152c-3.196 0-6.1-1.248-8.25-3.286zm0 13.036h.008v.008h-.008v-.008z" />
  </svg>
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
          color: 'var(--gray-900)',
          lineHeight: '1.1',
          marginBottom: '1rem'
        }}>
          Privacy First
        </h1>
        <p style={{
          fontSize: '1.125rem',
          lineHeight: '1.75',
          color: 'var(--gray-600)',
          maxWidth: '42rem',
          margin: '0 auto'
        }}>
          Your genomic data is the most personal information you have. We believe it should stay private and under your complete control.
        </p>
      </div>
    </section>

    <section style={{ background: 'white', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '4rem', flexWrap: 'wrap' }}>
          <div style={{ flex: '1', minWidth: '300px' }}>
            <span style={{
              color: 'var(--primary-blue)',
              fontWeight: '600',
              background: 'rgba(59, 130, 246, 0.1)',
              borderRadius: '9999px',
              padding: '0.5rem 1rem',
              fontSize: '0.875rem'
            }}>Uncompromising Security</span>
            <h2 style={{ marginTop: '1rem', fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: 'var(--gray-900)' }}>
              Your Data is Your Own. Period.
            </h2>
            <p style={{ marginTop: '1.5rem', fontSize: '1.125rem', color: 'var(--gray-600)', lineHeight: '1.75' }}>
              In an age of data breaches, we believe your most personal information should remain in your hands. GeneKnow is built on a foundation of privacy-by-design. By processing everything locally, we eliminate the risk of server-side breaches and data misuse.
            </p>
            <ul style={{ marginTop: '1.5rem', listStyle: 'none', padding: 0 }}>
              <li style={{ display: 'flex', alignItems: 'flex-start', marginBottom: '1rem' }}>
                <CheckBadgeIcon />
                <span style={{ color: 'var(--gray-600)' }}>
                  <strong>No Cloud Uploads:</strong> Your genomic file is never sent to us or any third party.
                </span>
              </li>
              <li style={{ display: 'flex', alignItems: 'flex-start' }}>
                <CheckBadgeIcon />
                <span style={{ color: 'var(--gray-600)' }}>
                  <strong>Secure by Design:</strong> Built with Rust and Tauri for a secure, sandboxed application environment.
                </span>
              </li>
            </ul>
          </div>
          <div style={{ flex: '1', minWidth: '300px' }}>
            <div style={{
              width: '100%',
              height: '300px',
              background: 'linear-gradient(135deg, #EBF4FF 0%, #3B82F6 100%)',
              borderRadius: 'var(--radius-lg)',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              color: 'white',
              fontSize: '1.25rem',
              fontWeight: '600',
              textAlign: 'center',
              boxShadow: 'var(--shadow-xl)'
            }}>
              🔒 Secure Analysis<br />Illustration
            </div>
          </div>
        </div>
      </div>
    </section>

    <section style={{ background: 'var(--gray-50)', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
          <h2 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: 'var(--gray-900)' }}>
            Privacy Guarantees
          </h2>
          <p style={{ marginTop: '1rem', fontSize: '1.125rem', color: 'var(--gray-600)', maxWidth: '32rem', margin: '1rem auto 0' }}>
            Our commitment to your privacy goes beyond just keeping your data local.
          </p>
        </div>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', gap: '2rem' }}>
          <div className="card" style={{ textAlign: 'center' }}>
            <div style={{
              margin: '0 auto 1.25rem',
              width: '3rem',
              height: '3rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              borderRadius: '50%',
              background: 'rgba(16, 185, 129, 0.1)',
              fontSize: '1.5rem'
            }}>
              🔒
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: 'var(--gray-900)', marginBottom: '0.5rem' }}>
              Zero Data Collection
            </h3>
            <p style={{ color: 'var(--gray-600)' }}>
              We don't collect, store, or analyze any of your genomic data. Everything stays on your device.
            </p>
          </div>
          <div className="card" style={{ textAlign: 'center' }}>
            <div style={{
              margin: '0 auto 1.25rem',
              width: '3rem',
              height: '3rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              borderRadius: '50%',
              background: 'rgba(16, 185, 129, 0.1)',
              fontSize: '1.5rem'
            }}>
              🛡️
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: 'var(--gray-900)', marginBottom: '0.5rem' }}>
              No Third-Party Sharing
            </h3>
            <p style={{ color: 'var(--gray-600)' }}>
              Your data never leaves your computer, so there's no risk of unauthorized sharing or breaches.
            </p>
          </div>
          <div className="card" style={{ textAlign: 'center' }}>
            <div style={{
              margin: '0 auto 1.25rem',
              width: '3rem',
              height: '3rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              borderRadius: '50%',
              background: 'rgba(16, 185, 129, 0.1)',
              fontSize: '1.5rem'
            }}>
              🔐
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: 'var(--gray-900)', marginBottom: '0.5rem' }}>
              Offline Processing
            </h3>
            <p style={{ color: 'var(--gray-600)' }}>
              Analysis happens completely offline. No internet connection required during processing.
            </p>
          </div>
        </div>
      </div>
    </section>

    <section style={{ background: 'white', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
          <h2 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: 'var(--gray-900)' }}>
            Technical Privacy Implementation
          </h2>
        </div>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(350px, 1fr))', gap: '3rem' }}>
          <div>
            <h3 style={{ fontSize: '1.5rem', fontWeight: '600', color: 'var(--gray-900)', marginBottom: '1rem' }}>
              Local-Only Architecture
            </h3>
            <ul style={{ listStyle: 'none', padding: 0 }}>
              <li style={{ marginBottom: '0.75rem', color: 'var(--gray-600)' }}>
                💻 All processing happens on your local machine
              </li>
              <li style={{ marginBottom: '0.75rem', color: 'var(--gray-600)' }}>
                🚫 No network requests during analysis
              </li>
              <li style={{ marginBottom: '0.75rem', color: 'var(--gray-600)' }}>
                🔒 Data never leaves your device
              </li>
              <li style={{ marginBottom: '0.75rem', color: 'var(--gray-600)' }}>
                📱 Works completely offline
              </li>
            </ul>
          </div>
          <div>
            <h3 style={{ fontSize: '1.5rem', fontWeight: '600', color: 'var(--gray-900)', marginBottom: '1rem' }}>
              Security Measures
            </h3>
            <ul style={{ listStyle: 'none', padding: 0 }}>
              <li style={{ marginBottom: '0.75rem', color: 'var(--gray-600)' }}>
                🔐 Memory-safe Rust implementation
              </li>
              <li style={{ marginBottom: '0.75rem', color: 'var(--gray-600)' }}>
                🛡️ Sandboxed Tauri environment
              </li>
              <li style={{ marginBottom: '0.75rem', color: 'var(--gray-600)' }}>
                🔒 Encrypted temporary files
              </li>
              <li style={{ marginBottom: '0.75rem', color: 'var(--gray-600)' }}>
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