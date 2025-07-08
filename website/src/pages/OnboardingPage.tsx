import Layout from '../components/Layout'

const CloudArrowUpIcon = () => (
  <svg className="h-6 w-6" style={{ color: '#2563EB' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M12 16.5V9.75m0 0l3 3m-3-3l-3 3M6.75 19.5a4.5 4.5 0 01-1.41-8.775 5.25 5.25 0 0110.233-2.33 3 3 0 013.758 3.848A3.752 3.752 0 0118 19.5H6.75z" />
  </svg>
)

const CheckCircleIcon = () => (
  <svg className="h-6 w-6" style={{ color: '#2563EB' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M9 12.75L11.25 15 15 9.75M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
  </svg>
)

const DocumentTextIcon = () => (
  <svg className="h-6 w-6" style={{ color: '#2563EB' }} fill="none" viewBox="0 0 24 24" strokeWidth="1.5" stroke="currentColor">
    <path strokeLinecap="round" strokeLinejoin="round" d="M10.5 6h9.75M10.5 6a1.5 1.5 0 11-3 0m3 0a1.5 1.5 0 10-3 0M3.75 6H7.5m3 12h9.75m-9.75 0a1.5 1.5 0 01-3 0m3 0a1.5 1.5 0 00-3 0M3.75 18H7.5m9-6h3.75m-3.75 0a1.5 1.5 0 01-3 0m3 0a1.5 1.5 0 00-3 0m-9.75 0h9.75" />
  </svg>
)

const DownloadButton = ({ os, icon, href, version = "1.0.0", size = "~50MB" }: { 
  os: string; 
  icon: string; 
  href: string; 
  version?: string; 
  size?: string; 
}) => (
  <a 
    href={href}
    style={{
      display: 'block',
      padding: '1.5rem 2rem',
      backgroundColor: '#2563EB',
      color: '#FFFFFF',
      textDecoration: 'none',
      borderRadius: '0.75rem',
      textAlign: 'center',
      minWidth: '200px',
      transition: 'all 200ms ease',
      boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)'
    }}
    onMouseEnter={(e) => {
      e.currentTarget.style.backgroundColor = '#1D4ED8';
      e.currentTarget.style.transform = 'translateY(-2px)';
      e.currentTarget.style.boxShadow = '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)';
    }}
    onMouseLeave={(e) => {
      e.currentTarget.style.backgroundColor = '#2563EB';
      e.currentTarget.style.transform = 'translateY(0)';
      e.currentTarget.style.boxShadow = '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)';
    }}
  >
    <div style={{ fontSize: '2rem', marginBottom: '0.5rem' }}>{icon}</div>
    <div style={{ fontWeight: 'bold', fontSize: '1.1rem' }}>Download for {os}</div>
    <div style={{ fontSize: '0.875rem', opacity: 0.8, marginTop: '0.25rem' }}>
      Version {version} • {size}
    </div>
  </a>
)

const OnboardingPage = () => {
  const downloadLinks = {
    windows: 'https://github.com/yourusername/Geneknow/releases/latest/download/geneknow-windows.exe',
    macos: 'https://github.com/yourusername/Geneknow/releases/latest/download/geneknow-macos.dmg',
    linux: 'https://github.com/yourusername/Geneknow/releases/latest/download/geneknow-linux.AppImage'
  }

  return (
    <Layout>
      {/* Hero Section */}
      <section style={{
        background: 'radial-gradient(circle at top left, rgba(239, 246, 255, 1) 0%, rgba(255, 255, 255, 1) 50%)',
        padding: '4rem 0',
        minHeight: '100vh',
        display: 'flex',
        alignItems: 'center'
      }}>
        <div style={{ maxWidth: '1400px', margin: '0 auto', padding: '0 2rem', width: '100%' }}>
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '4rem', alignItems: 'center' }}>
            
            {/* Left Side - Hero Content */}
            <div>
              <span style={{
                color: '#2563EB',
                fontWeight: '600',
                background: '#DBEAFE',
                borderRadius: '9999px',
                padding: '0.5rem 1rem',
                fontSize: '0.875rem',
                display: 'inline-block'
              }}>Your Personal Genomic Insights</span>
              
              <h1 style={{
                marginTop: '1.5rem',
                fontSize: 'clamp(2.5rem, 5vw, 4rem)',
                fontWeight: 'bold',
                letterSpacing: '-0.02em',
                color: '#111827',
                lineHeight: '1.1'
              }}>
                Understand Your Genomic Health, Privately.
              </h1>
              
              <p style={{
                marginTop: '1.5rem',
                fontSize: '1.25rem',
                lineHeight: '1.75',
                color: '#4B5563'
              }}>
                Analyze your genomic data file (`.fastq`) directly on your computer for secure cancer risk assessment. No uploads, no cloud processing.
              </p>

              {/* Key Features - Inline */}
              <div style={{ 
                display: 'grid', 
                gridTemplateColumns: 'repeat(3, 1fr)', 
                gap: '1.5rem', 
                marginTop: '2.5rem' 
              }}>
                <div style={{ textAlign: 'center' }}>
                  <div style={{ color: '#2563EB', fontSize: '2rem', marginBottom: '0.5rem' }}>🔒</div>
                  <div style={{ fontSize: '0.875rem', fontWeight: '600', color: '#111827' }}>100% Local</div>
                </div>
                <div style={{ textAlign: 'center' }}>
                  <div style={{ color: '#2563EB', fontSize: '2rem', marginBottom: '0.5rem' }}>🤖</div>
                  <div style={{ fontSize: '0.875rem', fontWeight: '600', color: '#111827' }}>AI-Powered</div>
                </div>
                <div style={{ textAlign: 'center' }}>
                  <div style={{ color: '#2563EB', fontSize: '2rem', marginBottom: '0.5rem' }}>📊</div>
                  <div style={{ fontSize: '0.875rem', fontWeight: '600', color: '#111827' }}>Clear Reports</div>
                </div>
              </div>
            </div>

            {/* Right Side - Download Section */}
            <div style={{ 
              background: '#FFFFFF',
              borderRadius: '1rem',
              padding: '3rem',
              boxShadow: '0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)',
              border: '1px solid #E5E7EB'
            }}>
              <h2 style={{
                fontSize: '1.875rem',
                fontWeight: 'bold',
                color: '#111827',
                textAlign: 'center',
                marginBottom: '1rem'
              }}>
                Download GeneKnow
              </h2>
              
              <p style={{
                textAlign: 'center',
                color: '#4B5563',
                marginBottom: '2rem',
                fontSize: '1rem'
              }}>
                Choose your operating system to get started:
              </p>
              
              <div style={{ 
                display: 'grid', 
                gap: '1rem',
                marginBottom: '2rem'
              }}>
                <DownloadButton 
                  os="Windows" 
                  icon="🖥️" 
                  href={downloadLinks.windows}
                />
                <DownloadButton 
                  os="macOS" 
                  icon="🍎" 
                  href={downloadLinks.macos}
                />
                <DownloadButton 
                  os="Linux" 
                  icon="🐧" 
                  href={downloadLinks.linux}
                />
              </div>
              
              <div style={{ 
                textAlign: 'center', 
                fontSize: '0.875rem', 
                color: '#6B7280',
                borderTop: '1px solid #E5E7EB',
                paddingTop: '1rem'
              }}>
                System Requirements: 8GB RAM, 2GB storage
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* How It Works + Privacy - Compact */}
      <section style={{ background: '#F9FAFB', padding: '4rem 0' }}>
        <div style={{ maxWidth: '1400px', margin: '0 auto', padding: '0 2rem' }}>
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '4rem', alignItems: 'start' }}>
            
            {/* How It Works */}
            <div>
              <h2 style={{ 
                fontSize: '1.875rem', 
                fontWeight: 'bold', 
                color: '#111827', 
                marginBottom: '2rem',
                textAlign: 'center'
              }}>
                How It Works
              </h2>
              
              <div style={{ display: 'grid', gap: '1.5rem' }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '1rem' }}>
                  <div style={{
                    width: '3rem',
                    height: '3rem',
                    background: '#2563EB',
                    borderRadius: '50%',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    color: '#FFFFFF',
                    fontWeight: 'bold',
                    fontSize: '1.25rem',
                    flexShrink: 0
                  }}>1</div>
                  <div>
                    <h3 style={{ fontWeight: '600', color: '#111827', marginBottom: '0.25rem' }}>
                      Select Your File
                    </h3>
                    <p style={{ color: '#4B5563', fontSize: '0.875rem' }}>
                      Choose your `.fastq` file from your computer
                    </p>
                  </div>
                </div>
                
                <div style={{ display: 'flex', alignItems: 'center', gap: '1rem' }}>
                  <div style={{
                    width: '3rem',
                    height: '3rem',
                    background: '#2563EB',
                    borderRadius: '50%',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    color: '#FFFFFF',
                    fontWeight: 'bold',
                    fontSize: '1.25rem',
                    flexShrink: 0
                  }}>2</div>
                  <div>
                    <h3 style={{ fontWeight: '600', color: '#111827', marginBottom: '0.25rem' }}>
                      Run Analysis
                    </h3>
                    <p style={{ color: '#4B5563', fontSize: '0.875rem' }}>
                      AI analyzes your data locally on your machine
                    </p>
                  </div>
                </div>
                
                <div style={{ display: 'flex', alignItems: 'center', gap: '1rem' }}>
                  <div style={{
                    width: '3rem',
                    height: '3rem',
                    background: '#2563EB',
                    borderRadius: '50%',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    color: '#FFFFFF',
                    fontWeight: 'bold',
                    fontSize: '1.25rem',
                    flexShrink: 0
                  }}>3</div>
                  <div>
                    <h3 style={{ fontWeight: '600', color: '#111827', marginBottom: '0.25rem' }}>
                      Get Your Report
                    </h3>
                    <p style={{ color: '#4B5563', fontSize: '0.875rem' }}>
                      Receive clear, actionable health insights
                    </p>
                  </div>
                </div>
              </div>
            </div>

            {/* Privacy Highlights */}
            <div>
              <h2 style={{ 
                fontSize: '1.875rem', 
                fontWeight: 'bold', 
                color: '#111827', 
                marginBottom: '2rem',
                textAlign: 'center'
              }}>
                Privacy Guaranteed
              </h2>
              
              <div style={{ display: 'grid', gap: '1.5rem' }}>
                <div style={{ 
                  display: 'flex', 
                  alignItems: 'center', 
                  gap: '1rem',
                  background: '#FFFFFF',
                  padding: '1.5rem',
                  borderRadius: '0.75rem',
                  boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
                }}>
                  <div style={{ fontSize: '1.5rem', flexShrink: 0 }}>🔒</div>
                  <div>
                    <h3 style={{ fontWeight: '600', color: '#111827', marginBottom: '0.25rem' }}>
                      Zero Data Collection
                    </h3>
                    <p style={{ color: '#4B5563', fontSize: '0.875rem' }}>
                      Your data never leaves your computer
                    </p>
                  </div>
                </div>
                
                <div style={{ 
                  display: 'flex', 
                  alignItems: 'center', 
                  gap: '1rem',
                  background: '#FFFFFF',
                  padding: '1.5rem',
                  borderRadius: '0.75rem',
                  boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
                }}>
                  <div style={{ fontSize: '1.5rem', flexShrink: 0 }}>🛡️</div>
                  <div>
                    <h3 style={{ fontWeight: '600', color: '#111827', marginBottom: '0.25rem' }}>
                      Offline Processing
                    </h3>
                    <p style={{ color: '#4B5563', fontSize: '0.875rem' }}>
                      Works completely without internet connection
                    </p>
                  </div>
                </div>
                
                <div style={{ 
                  display: 'flex', 
                  alignItems: 'center', 
                  gap: '1rem',
                  background: '#FFFFFF',
                  padding: '1.5rem',
                  borderRadius: '0.75rem',
                  boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
                }}>
                  <div style={{ fontSize: '1.5rem', flexShrink: 0 }}>🔐</div>
                  <div>
                    <h3 style={{ fontWeight: '600', color: '#111827', marginBottom: '0.25rem' }}>
                      Secure by Design
                    </h3>
                    <p style={{ color: '#4B5563', fontSize: '0.875rem' }}>
                      Built with Rust and Tauri for maximum security
                    </p>
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Features - Compact */}
      <section style={{ background: '#FFFFFF', padding: '4rem 0' }}>
        <div style={{ maxWidth: '1400px', margin: '0 auto', padding: '0 2rem' }}>
          <h2 style={{ 
            fontSize: '1.875rem', 
            fontWeight: 'bold', 
            color: '#111827', 
            textAlign: 'center',
            marginBottom: '3rem'
          }}>
            Why Choose GeneKnow?
          </h2>
          
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: '2rem' }}>
            <div style={{ 
              textAlign: 'center',
              padding: '2rem',
              background: '#F9FAFB',
              borderRadius: '0.75rem',
              border: '1px solid #E5E7EB'
            }}>
              <div style={{
                margin: '0 auto 1rem',
                width: '3rem',
                height: '3rem',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                borderRadius: '50%',
                background: '#DBEAFE'
              }}>
                <CloudArrowUpIcon />
              </div>
              <h3 style={{ fontSize: '1.125rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
                Local Processing
              </h3>
              <p style={{ color: '#4B5563', fontSize: '0.875rem' }}>
                Your sensitive genetic data never leaves your computer
              </p>
            </div>
            
            <div style={{ 
              textAlign: 'center',
              padding: '2rem',
              background: '#F9FAFB',
              borderRadius: '0.75rem',
              border: '1px solid #E5E7EB'
            }}>
              <div style={{
                margin: '0 auto 1rem',
                width: '3rem',
                height: '3rem',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                borderRadius: '50%',
                background: '#DBEAFE'
              }}>
                <CheckCircleIcon />
              </div>
              <h3 style={{ fontSize: '1.125rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
                State-of-the-Art AI
              </h3>
              <p style={{ color: '#4B5563', fontSize: '0.875rem' }}>
                Advanced models provide insights based on latest research
              </p>
            </div>
            
            <div style={{ 
              textAlign: 'center',
              padding: '2rem',
              background: '#F9FAFB',
              borderRadius: '0.75rem',
              border: '1px solid #E5E7EB'
            }}>
              <div style={{
                margin: '0 auto 1rem',
                width: '3rem',
                height: '3rem',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                borderRadius: '50%',
                background: '#DBEAFE'
              }}>
                <DocumentTextIcon />
              </div>
              <h3 style={{ fontSize: '1.125rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
                Clear Reports
              </h3>
              <p style={{ color: '#4B5563', fontSize: '0.875rem' }}>
                Easy-to-understand reports without medical jargon
              </p>
            </div>
          </div>
        </div>
      </section>
    </Layout>
  )
}

export default OnboardingPage 