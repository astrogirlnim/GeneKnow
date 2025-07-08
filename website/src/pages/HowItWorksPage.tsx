import Layout from '../components/Layout'

const HowItWorksPage = () => (
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
          How It Works
        </h1>
        <p style={{
          fontSize: '1.125rem',
          lineHeight: '1.75',
          color: '#4B5563',
          maxWidth: '42rem',
          margin: '0 auto'
        }}>
          GeneKnow makes genomic analysis simple and secure. Follow these three easy steps to get your personalized health insights.
        </p>
      </div>
    </section>

    <section style={{ background: '#F9FAFB', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
          <h2 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: '#111827' }}>
            Get Your Report in 3 Simple Steps
          </h2>
        </div>
        <div style={{ position: 'relative' }}>
          <div style={{
            position: 'absolute',
            top: '2.5rem',
            left: '0',
            right: '0',
            height: '2px',
            background: 'repeating-linear-gradient(to right, #D1D5DB 0, #D1D5DB 10px, transparent 10px, transparent 20px)',
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
                background: '#FFFFFF',
                boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
                border: '4px solid #2563EB',
                color: '#2563EB',
                fontSize: '1.5rem',
                fontWeight: 'bold'
              }}>1</div>
              <h3 style={{ marginTop: '1.5rem', fontSize: '1.25rem', fontWeight: '600', color: '#111827' }}>Select Your File</h3>
              <p style={{ marginTop: '0.5rem', color: '#4B5563' }}>
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
                background: '#FFFFFF',
                boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
                border: '4px solid #2563EB',
                color: '#2563EB',
                fontSize: '1.5rem',
                fontWeight: 'bold'
              }}>2</div>
              <h3 style={{ marginTop: '1.5rem', fontSize: '1.25rem', fontWeight: '600', color: '#111827' }}>Run Analysis</h3>
              <p style={{ marginTop: '0.5rem', color: '#4B5563' }}>
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
                background: '#FFFFFF',
                boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
                border: '4px solid #2563EB',
                color: '#2563EB',
                fontSize: '1.5rem',
                fontWeight: 'bold'
              }}>3</div>
              <h3 style={{ marginTop: '1.5rem', fontSize: '1.25rem', fontWeight: '600', color: '#111827' }}>View Your Report</h3>
              <p style={{ marginTop: '0.5rem', color: '#4B5563' }}>
                Receive and save your private, easy-to-understand health insights report.
              </p>
            </div>
          </div>
        </div>
      </div>
    </section>

    <section style={{ background: '#FFFFFF', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
          <h2 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: '#111827' }}>
            What Happens During Analysis
          </h2>
          <p style={{ marginTop: '1rem', fontSize: '1.125rem', color: '#4B5563', maxWidth: '32rem', margin: '1rem auto 0' }}>
            Here's what GeneKnow does with your genomic data during the analysis process.
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
              background: '#DBEAFE',
              fontSize: '1.5rem'
            }}>
              🧬
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
              Data Preprocessing
            </h3>
            <p style={{ color: '#4B5563' }}>
              Your FASTQ file is securely parsed and prepared for analysis using established bioinformatics protocols.
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
              background: '#DBEAFE',
              fontSize: '1.5rem'
            }}>
              🔬
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
              Variant Analysis
            </h3>
            <p style={{ color: '#4B5563' }}>
              AI models identify genetic variants and compare them against known health-related genetic markers.
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
              background: '#DBEAFE',
              fontSize: '1.5rem'
            }}>
              📊
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
              Risk Assessment
            </h3>
            <p style={{ color: '#4B5563' }}>
              Sophisticated algorithms calculate potential health risks based on your unique genetic profile.
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
              background: '#DBEAFE',
              fontSize: '1.5rem'
            }}>
              📝
            </div>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '0.5rem' }}>
              Report Generation
            </h3>
            <p style={{ color: '#4B5563' }}>
              Llama 3.1 creates a personalized, easy-to-understand report with actionable health insights.
            </p>
          </div>
        </div>
      </div>
    </section>

    <section style={{ background: '#F9FAFB', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
          <h2 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: '#111827' }}>
            Technical Requirements
          </h2>
          <p style={{ marginTop: '1rem', fontSize: '1.125rem', color: '#4B5563', maxWidth: '32rem', margin: '1rem auto 0' }}>
            Make sure your system meets these requirements for optimal performance.
          </p>
        </div>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(350px, 1fr))', gap: '3rem' }}>
          <div>
            <h3 style={{ fontSize: '1.5rem', fontWeight: '600', color: '#111827', marginBottom: '1rem' }}>
              System Requirements
            </h3>
            <ul style={{ listStyle: 'none', padding: 0 }}>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                💻 Windows 10/11, macOS 10.15+, or Linux
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                🧠 8GB RAM minimum (16GB recommended)
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                💾 2GB free storage space
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                ⚡ Modern multi-core processor
              </li>
            </ul>
          </div>
          <div>
            <h3 style={{ fontSize: '1.5rem', fontWeight: '600', color: '#111827', marginBottom: '1rem' }}>
              File Requirements
            </h3>
            <ul style={{ listStyle: 'none', padding: 0 }}>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                📁 FASTQ file format (.fastq or .fq)
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                🔬 Whole genome sequencing data
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                📏 Minimum 30x coverage recommended
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                ✅ Human genomic data (hg38 reference)
              </li>
            </ul>
          </div>
        </div>
      </div>
    </section>

    <section style={{ background: '#FFFFFF', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem', textAlign: 'center' }}>
        <h2 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: '#111827', marginBottom: '1rem' }}>
          Ready to Get Started?
        </h2>
        <p style={{ fontSize: '1.125rem', color: '#4B5563', maxWidth: '32rem', margin: '0 auto 2rem' }}>
          Download GeneKnow and start your personalized genomic analysis journey today.
        </p>
        <button style={{
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
          Start Analysis
        </button>
      </div>
    </section>
  </Layout>
)

export default HowItWorksPage 