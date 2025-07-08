import React from 'react'

// --- Helper Components for Icons ---

const GeneKnowLogo = () => (
  <svg className="w-8 h-8" style={{ color: 'var(--primary-blue)' }} viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
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

const CloudArrowUpIcon = () => (
  <svg className="h-6 w-6" style={{ color: 'var(--primary-blue)' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M12 16.5V9.75m0 0l3 3m-3-3l-3 3M6.75 19.5a4.5 4.5 0 01-1.41-8.775 5.25 5.25 0 0110.233-2.33 3 3 0 013.758 3.848A3.752 3.752 0 0118 19.5H6.75z" />
  </svg>
)

const CheckCircleIcon = () => (
  <svg className="h-6 w-6" style={{ color: 'var(--primary-blue)' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M9 12.75L11.25 15 15 9.75M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
  </svg>
)

const DocumentTextIcon = () => (
  <svg className="h-6 w-6" style={{ color: 'var(--primary-blue)' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M10.5 6h9.75M10.5 6a1.5 1.5 0 11-3 0m3 0a1.5 1.5 0 10-3 0M3.75 6H7.5m3 12h9.75m-9.75 0a1.5 1.5 0 01-3 0m3 0a1.5 1.5 0 00-3 0M3.75 18H7.5m9-6h3.75m-3.75 0a1.5 1.5 0 01-3 0m3 0a1.5 1.5 0 00-3 0m-9.75 0h9.75" />
  </svg>
)

const CheckBadgeIcon = () => (
  <svg className="h-6 w-6" style={{ color: 'var(--success)', marginRight: '12px' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M9 12.75L11.25 15 15 9.75m-3-7.036A11.959 11.959 0 013.598 6 11.99 11.99 0 003 9.749c0 5.592 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.31-.21-2.571-.598-3.751h-.152c-3.196 0-6.1-1.248-8.25-3.286zm0 13.036h.008v.008h-.008v-.008z" />
  </svg>
)

// --- Page Section Components ---

const Header = () => (
  <header style={{
    background: 'rgba(255, 255, 255, 0.9)',
    backdropFilter: 'blur(10px)',
    position: 'fixed',
    top: 0,
    left: 0,
    right: 0,
    zIndex: 'var(--z-fixed)',
    borderBottom: '1px solid var(--gray-200)',
    padding: '1rem 0'
  }}>
    <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
      <div className="flex items-center gap-2">
        <GeneKnowLogo />
        <h1 style={{ fontSize: '1.5rem', fontWeight: 'bold', color: 'var(--gray-800)', margin: 0 }}>GeneKnow</h1>
      </div>
      <nav className="hidden" style={{ display: 'flex', alignItems: 'center', gap: '2rem' }}>
        <a href="#features" style={{ color: 'var(--gray-600)', textDecoration: 'none', transition: 'color var(--transition-normal)' }}>Features</a>
        <a href="#privacy" style={{ color: 'var(--gray-600)', textDecoration: 'none', transition: 'color var(--transition-normal)' }}>Privacy First</a>
        <a href="#how-it-works" style={{ color: 'var(--gray-600)', textDecoration: 'none', transition: 'color var(--transition-normal)' }}>How It Works</a>
      </nav>
    </div>
  </header>
)

const HeroSection = () => (
  <section style={{
    background: 'radial-gradient(circle at top left, rgba(239, 246, 255, 1) 0%, rgba(255, 255, 255, 1) 50%)',
    padding: '5rem 0 8rem 0'
  }}>
    <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem', textAlign: 'center' }}>
      <div style={{ maxWidth: '48rem', margin: '0 auto' }}>
        <span style={{
          color: 'var(--primary-blue)',
          fontWeight: '600',
          background: 'rgba(59, 130, 246, 0.1)',
          borderRadius: '9999px',
          padding: '0.5rem 1rem',
          fontSize: '0.875rem'
        }}>Your Personal Genomic Insights</span>
        <h2 style={{
          marginTop: '1rem',
          fontSize: 'clamp(2.5rem, 5vw, 3.75rem)',
          fontWeight: 'bold',
          letterSpacing: '-0.02em',
          color: 'var(--gray-900)',
          lineHeight: '1.1'
        }}>
          Understand Your Genomic Health, Privately.
        </h2>
        <p style={{
          marginTop: '1.5rem',
          fontSize: '1.125rem',
          lineHeight: '1.75',
          color: 'var(--gray-600)',
          maxWidth: '42rem',
          margin: '1.5rem auto 0'
        }}>
          Our desktop application analyzes your genomic data file (`.fastq`) directly on your computer to provide a secure cancer risk assessment. No data uploads, no cloud processing, complete privacy.
        </p>
        <div style={{ marginTop: '2.5rem' }}>
          <button className="btn btn-primary btn-large" style={{ transform: 'none', transition: 'all var(--transition-normal)' }}>
            Test Now
          </button>
        </div>
      </div>
    </div>
  </section>
)

const FeatureCard = ({ icon, title, children }: { icon: React.ReactNode; title: string; children: React.ReactNode }) => (
  <div className="card" style={{ textAlign: 'center', height: '100%' }}>
    <div style={{
      margin: '0 auto 1.25rem',
      width: '3rem',
      height: '3rem',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      borderRadius: '50%',
      background: 'rgba(59, 130, 246, 0.1)'
    }}>
      {icon}
    </div>
    <h4 style={{ marginTop: '1.25rem', fontSize: '1.25rem', fontWeight: '600', color: 'var(--gray-900)' }}>{title}</h4>
    <p style={{ marginTop: '0.5rem', color: 'var(--gray-600)' }}>{children}</p>
  </div>
)

const FeaturesSection = () => (
  <section id="features" style={{ background: 'var(--gray-50)', padding: '5rem 0' }}>
    <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
      <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
        <h3 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: 'var(--gray-900)' }}>
          A New Standard in Personal Genetic Analysis
        </h3>
        <p style={{ marginTop: '1rem', fontSize: '1.125rem', color: 'var(--gray-600)', maxWidth: '32rem', margin: '1rem auto 0' }}>
          We combine powerful technology with an unwavering commitment to your privacy.
        </p>
      </div>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', gap: '2rem' }}>
        <FeatureCard icon={<CloudArrowUpIcon />} title="Local Processing">
          Your `.fastq` file is analyzed on your machine. Your sensitive genetic data never leaves your computer.
        </FeatureCard>
        <FeatureCard icon={<CheckCircleIcon />} title="State-of-the-Art AI">
          Leverages advanced models to provide insights on potential risks based on the latest scientific research.
        </FeatureCard>
        <FeatureCard icon={<DocumentTextIcon />} title="Comprehensive Reports">
          Receive a clear, understandable report generated by Llama 3.1, explaining the findings without jargon.
        </FeatureCard>
      </div>
    </div>
  </section>
)

const PrivacySection = () => (
  <section id="privacy" style={{ background: 'white', padding: '5rem 0' }}>
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
          <h3 style={{ marginTop: '1rem', fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: 'var(--gray-900)' }}>
            Your Data is Your Own. Period.
          </h3>
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
)

const HowItWorksSection = () => (
  <section id="how-it-works" style={{ background: 'var(--gray-50)', padding: '5rem 0' }}>
    <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
      <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
        <h3 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: 'var(--gray-900)' }}>
          Get Your Report in 3 Simple Steps
        </h3>
      </div>
      <div style={{ position: 'relative' }}>
        <div style={{
          position: 'absolute',
          top: '2.5rem',
          left: '0',
          right: '0',
          height: '2px',
          background: 'repeating-linear-gradient(to right, var(--gray-300) 0, var(--gray-300) 10px, transparent 10px, transparent 20px)',
          zIndex: '1'
        }}></div>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))', gap: '3rem', position: 'relative', zIndex: '2' }}>
          <div style={{ textAlign: 'center' }}>
            <div style={{
              margin: '0 auto',
              width: '5rem',
              height: '5rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              borderRadius: '50%',
              background: 'white',
              boxShadow: 'var(--shadow-lg)',
              border: '4px solid var(--primary-blue)',
              color: 'var(--primary-blue)',
              fontSize: '1.5rem',
              fontWeight: 'bold'
            }}>1</div>
            <h4 style={{ marginTop: '1.5rem', fontSize: '1.25rem', fontWeight: '600', color: 'var(--gray-900)' }}>Select Your File</h4>
            <p style={{ marginTop: '0.5rem', color: 'var(--gray-600)' }}>
              Open the app and select your `.fastq` file from your local disk.
            </p>
          </div>
          <div style={{ textAlign: 'center' }}>
            <div style={{
              margin: '0 auto',
              width: '5rem',
              height: '5rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              borderRadius: '50%',
              background: 'white',
              boxShadow: 'var(--shadow-lg)',
              border: '4px solid var(--primary-blue)',
              color: 'var(--primary-blue)',
              fontSize: '1.5rem',
              fontWeight: 'bold'
            }}>2</div>
            <h4 style={{ marginTop: '1.5rem', fontSize: '1.25rem', fontWeight: '600', color: 'var(--gray-900)' }}>Run Analysis</h4>
            <p style={{ marginTop: '0.5rem', color: 'var(--gray-600)' }}>
              Click "Test Now". Our app performs the analysis locally using its built-in AI models.
            </p>
          </div>
          <div style={{ textAlign: 'center' }}>
            <div style={{
              margin: '0 auto',
              width: '5rem',
              height: '5rem',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              borderRadius: '50%',
              background: 'white',
              boxShadow: 'var(--shadow-lg)',
              border: '4px solid var(--primary-blue)',
              color: 'var(--primary-blue)',
              fontSize: '1.5rem',
              fontWeight: 'bold'
            }}>3</div>
            <h4 style={{ marginTop: '1.5rem', fontSize: '1.25rem', fontWeight: '600', color: 'var(--gray-900)' }}>View Your Report</h4>
            <p style={{ marginTop: '0.5rem', color: 'var(--gray-600)' }}>
              Receive and save your private, easy-to-understand health insights report.
            </p>
          </div>
        </div>
      </div>
    </div>
  </section>
)

const Footer = () => (
  <footer style={{ background: 'var(--gray-800)', color: 'white', padding: '2rem 0', textAlign: 'center' }}>
    <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
      <p>&copy; 2025 GeneKnow. All rights reserved.</p>
      <p style={{ fontSize: '0.875rem', color: 'var(--gray-400)', marginTop: '0.5rem' }}>
        Disclaimer: This tool is for informational purposes only and is not a substitute for professional medical advice, diagnosis, or treatment.
      </p>
    </div>
  </footer>
)

// --- Main App Component ---
export default function App() {
  return (
    <div style={{ fontFamily: 'var(--font-family-sans)', background: 'white', color: 'var(--gray-800)', minHeight: '100vh' }}>
      <Header />
      <main style={{ paddingTop: '4rem' }}>
        <HeroSection />
        <FeaturesSection />
        <PrivacySection />
        <HowItWorksSection />
      </main>
      <Footer />
    </div>
  )
}

