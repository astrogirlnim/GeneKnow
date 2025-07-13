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

    <section style={{ background: '#F3F4F6', padding: '5rem 0' }}>
      <div style={{ maxWidth: '1900px', margin: '0 auto', padding: '0 1.5rem' }}>
        <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
          <h2 style={{ fontSize: '1.875rem', fontWeight: 'bold', letterSpacing: '-0.02em', color: '#111827' }}>
            Technical Pipeline Architecture
          </h2>
          <p style={{ marginTop: '1rem', fontSize: '1.125rem', color: '#4B5563', maxWidth: '48rem', margin: '1rem auto 0' }}>
            For healthcare professionals: Our LangGraph orchestration system ensures accurate, reproducible genomic analysis through a sophisticated multi-stage pipeline.
          </p>
        </div>

        {/* Pipeline Visualization */}
        <div style={{
          background: '#FFFFFF',
          borderRadius: '1rem',
          padding: '2rem',
          boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)',
          border: '1px solid #E5E7EB',
          marginBottom: '3rem',
          overflowX: 'auto',
          width: '100%',
          maxWidth: '1800px',
          margin: '0 auto 3rem auto'
        }}>
          <div style={{
            display: 'flex',
            flexDirection: 'column',
            alignItems: 'center',
            gap: '2rem'
          }}>
            {/* Title */}
            <div style={{
              fontSize: '1.5rem',
              fontWeight: '600',
              color: '#1F2937'
            }}>
              LangGraph Orchestration
            </div>

            {/* Horizontal Pipeline Flow */}
            <div style={{ 
              display: 'flex', 
              alignItems: 'center', 
              gap: '0.5rem',
              minWidth: '1050px',
              padding: '0.5rem',
              justifyContent: 'flex-start',
              margin: '0'
            }}>
              {/* Preprocess */}
              <div style={{
                background: '#E3F2FD',
                border: '2px solid #1976D2',
                borderRadius: '0.5rem',
                padding: '0.75rem',
                textAlign: 'center',
                minWidth: '90px',
                width: '90px'
              }}>
                <div style={{ fontWeight: '600', color: '#1976D2', fontSize: '0.7rem' }}>Preprocess<br/>Data</div>
              </div>

              {/* Arrow */}
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <div style={{ width: '15px', height: '2px', background: '#1976D2' }} />
                <div style={{ 
                  width: '0', 
                  height: '0', 
                  borderTop: '4px solid transparent',
                  borderBottom: '4px solid transparent',
                  borderLeft: '6px solid #1976D2'
                }} />
              </div>

              {/* Variant Calling / QC */}
              <div style={{
                background: '#E3F2FD',
                border: '2px solid #1976D2',
                borderRadius: '0.5rem',
                padding: '0.75rem',
                textAlign: 'center',
                minWidth: '90px',
                width: '90px'
              }}>
                <div style={{ fontWeight: '600', color: '#1976D2', fontSize: '0.65rem' }}>Variant Calling<br/>QC Filtering</div>
              </div>

              {/* Arrow */}
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <div style={{ width: '15px', height: '2px', background: '#1976D2' }} />
                <div style={{ 
                  width: '0', 
                  height: '0', 
                  borderTop: '4px solid transparent',
                  borderBottom: '4px solid transparent',
                  borderLeft: '6px solid #1976D2'
                }} />
              </div>

              {/* Population Mapper */}
              <div style={{
                background: '#E3F2FD',
                border: '2px solid #1976D2',
                borderRadius: '0.5rem',
                padding: '0.75rem',
                textAlign: 'center',
                minWidth: '90px',
                width: '90px'
              }}>
                <div style={{ fontWeight: '600', color: '#1976D2', fontSize: '0.65rem' }}>Population<br/>Mapper</div>
              </div>

              {/* Branching Arrow - Adjusted for new box heights */}
              <div style={{ position: 'relative', width: '25px', height: '170px', display: 'flex', alignItems: 'center' }}>
                <svg width="25" height="170" style={{ position: 'absolute', left: '0', top: '50%', transform: 'translateY(-50%)' }}>
                  {/* Main horizontal line from Population Mapper */}
                  <line x1="0" y1="85" x2="12" y2="85" stroke="#1976D2" strokeWidth="2" />
                  {/* Vertical distribution line - only between arrows */}
                  <line x1="12" y1="15" x2="12" y2="155" stroke="#1976D2" strokeWidth="2" />
                  {/* Branches to each parallel node - precisely centered on each node */}
                  <line x1="12" y1="15" x2="20" y2="15" stroke="#1976D2" strokeWidth="2" />
                  <polygon points="20,12 25,15 20,18" fill="#1976D2" />
                  
                  <line x1="12" y1="50" x2="20" y2="50" stroke="#1976D2" strokeWidth="2" />
                  <polygon points="20,47 25,50 20,53" fill="#1976D2" />
                  
                  <line x1="12" y1="85" x2="20" y2="85" stroke="#1976D2" strokeWidth="2" />
                  <polygon points="20,82 25,85 20,88" fill="#1976D2" />
                  
                  <line x1="12" y1="120" x2="20" y2="120" stroke="#1976D2" strokeWidth="2" />
                  <polygon points="20,117 25,120 20,123" fill="#1976D2" />
                  
                  <line x1="12" y1="155" x2="20" y2="155" stroke="#1976D2" strokeWidth="2" />
                  <polygon points="20,152 25,155 20,158" fill="#1976D2" />
                </svg>
              </div>

              {/* Parallel Variant Models - Adjusted spacing */}
              <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
                {/* TCGA Mapper */}
                <div style={{
                  background: '#FFF3E0',
                  border: '2px solid #F57C00',
                  borderRadius: '0.5rem',
                  padding: '0.5rem 0.75rem',
                  textAlign: 'center',
                  minWidth: '105px',
                  width: '105px',
                  height: '30px',
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center'
                }}>
                  <div style={{ fontWeight: '600', color: '#E65100', fontSize: '0.65rem' }}>TCGA Mapper</div>
                </div>

                {/* CADD Scoring */}
                <div style={{
                  background: '#FFF3E0',
                  border: '2px solid #F57C00',
                  borderRadius: '0.5rem',
                  padding: '0.5rem 0.75rem',
                  textAlign: 'center',
                  minWidth: '105px',
                  width: '105px',
                  height: '30px',
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center'
                }}>
                  <div style={{ fontWeight: '600', color: '#E65100', fontSize: '0.65rem' }}>CADD Scoring</div>
                </div>

                {/* ClinVar Annotator */}
                <div style={{
                  background: '#FFF3E0',
                  border: '2px solid #F57C00',
                  borderRadius: '0.5rem',
                  padding: '0.5rem 0.75rem',
                  textAlign: 'center',
                  minWidth: '105px',
                  width: '105px',
                  height: '30px',
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center'
                }}>
                  <div style={{ fontWeight: '600', color: '#E65100', fontSize: '0.6rem' }}>ClinVar<br/>Annotator</div>
                </div>

                {/* PRS Calculator */}
                <div style={{
                  background: '#FFF3E0',
                  border: '2px solid #F57C00',
                  borderRadius: '0.5rem',
                  padding: '0.5rem 0.75rem',
                  textAlign: 'center',
                  minWidth: '105px',
                  width: '105px',
                  height: '30px',
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center'
                }}>
                  <div style={{ fontWeight: '600', color: '#E65100', fontSize: '0.6rem' }}>PRS<br/>Calculator</div>
                </div>

                {/* Pathway Burden */}
                <div style={{
                  background: '#FFF3E0',
                  border: '2px solid #F57C00',
                  borderRadius: '0.5rem',
                  padding: '0.5rem 0.75rem',
                  textAlign: 'center',
                  minWidth: '105px',
                  width: '105px',
                  height: '30px',
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center'
                }}>
                  <div style={{ fontWeight: '600', color: '#E65100', fontSize: '0.6rem' }}>Pathway<br/>Burden</div>
                </div>
              </div>

              {/* Merging Arrow - Adjusted for new box heights */}
              <div style={{ position: 'relative', width: '25px', height: '170px', display: 'flex', alignItems: 'center' }}>
                <svg width="25" height="170" style={{ position: 'absolute', left: '0', top: '50%', transform: 'translateY(-50%)' }}>
                  {/* Branches from each parallel node - precisely centered on each node */}
                  <line x1="0" y1="15" x2="12" y2="15" stroke="#1976D2" strokeWidth="2" />
                  <line x1="0" y1="50" x2="12" y2="50" stroke="#1976D2" strokeWidth="2" />
                  <line x1="0" y1="85" x2="12" y2="85" stroke="#1976D2" strokeWidth="2" />
                  <line x1="0" y1="120" x2="12" y2="120" stroke="#1976D2" strokeWidth="2" />
                  <line x1="0" y1="155" x2="12" y2="155" stroke="#1976D2" strokeWidth="2" />
                  {/* Vertical collection line - only between first and last arrows */}
                  <line x1="12" y1="15" x2="12" y2="155" stroke="#1976D2" strokeWidth="2" />
                  {/* Main horizontal line to ML Risk Model with arrow head */}
                  <line x1="12" y1="85" x2="20" y2="85" stroke="#1976D2" strokeWidth="2" />
                  <polygon points="20,82 25,85 20,88" fill="#1976D2" />
                </svg>
              </div>

              {/* Arrow */}
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <div style={{ width: '15px', height: '2px', background: '#1976D2' }} />
                <div style={{ 
                  width: '0', 
                  height: '0', 
                  borderTop: '4px solid transparent',
                  borderBottom: '4px solid transparent',
                  borderLeft: '6px solid #1976D2'
                }} />
              </div>

              {/* ML Risk Model */}
              <div style={{
                background: '#F3E5F5',
                border: '2px solid #9C27B0',
                borderRadius: '0.5rem',
                padding: '0.75rem',
                textAlign: 'center',
                minWidth: '90px',
                width: '90px'
              }}>
                <div style={{ fontWeight: '600', color: '#7B1FA2', fontSize: '0.65rem' }}>ML Risk<br/>Model</div>
              </div>

              {/* Arrow */}
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <div style={{ width: '15px', height: '2px', background: '#1976D2' }} />
                <div style={{ 
                  width: '0', 
                  height: '0', 
                  borderTop: '4px solid transparent',
                  borderBottom: '4px solid transparent',
                  borderLeft: '6px solid #1976D2'
                }} />
              </div>

              {/* Validation and Metrics */}
              <div style={{
                background: '#E3F2FD',
                border: '2px solid #1976D2',
                borderRadius: '0.5rem',
                padding: '0.75rem',
                textAlign: 'center',
                minWidth: '90px',
                width: '90px'
              }}>
                <div style={{ fontWeight: '600', color: '#1976D2', fontSize: '0.65rem' }}>Validation<br/>and Metrics</div>
              </div>

              {/* Arrow */}
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <div style={{ width: '15px', height: '2px', background: '#1976D2' }} />
                <div style={{ 
                  width: '0', 
                  height: '0', 
                  borderTop: '4px solid transparent',
                  borderBottom: '4px solid transparent',
                  borderLeft: '6px solid #1976D2'
                }} />
              </div>

              {/* Report Writer */}
              <div style={{
                background: '#E3F2FD',
                border: '2px solid #1976D2',
                borderRadius: '0.5rem',
                padding: '0.75rem',
                textAlign: 'center',
                minWidth: '90px',
                width: '90px'
              }}>
                <div style={{ fontWeight: '600', color: '#1976D2', fontSize: '0.65rem' }}>Report<br/>Writer</div>
              </div>

              {/* Arrow */}
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <div style={{ width: '15px', height: '2px', background: '#388E3C' }} />
                <div style={{ 
                  width: '0', 
                  height: '0', 
                  borderTop: '4px solid transparent',
                  borderBottom: '4px solid transparent',
                  borderLeft: '6px solid #388E3C'
                }} />
              </div>

              {/* Visualizations + Frontend */}
              <div style={{
                background: '#E0F2F1',
                border: '2px solid #388E3C',
                borderRadius: '0.5rem',
                padding: '0.75rem',
                textAlign: 'center',
                minWidth: '90px',
                width: '90px'
              }}>
                <div style={{ fontWeight: '600', color: '#388E3C', fontSize: '0.6rem' }}>Visualizations<br/>Frontend</div>
              </div>
            </div>

            {/* Legend */}
            <div style={{ 
              display: 'flex', 
              gap: '1.5rem', 
              marginTop: '1rem',
              fontSize: '0.75rem',
              color: '#6B7280',
              flexWrap: 'wrap',
              justifyContent: 'center'
            }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <div style={{ 
                  width: '16px', 
                  height: '8px', 
                  background: '#E3F2FD',
                  border: '1px solid #1976D2',
                  borderRadius: '2px'
                }} />
                <span>Processing Node</span>
              </div>
              <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <div style={{ 
                  width: '16px', 
                  height: '8px', 
                  background: '#FFF3E0',
                  border: '1px solid #F57C00',
                  borderRadius: '2px'
                }} />
                <span>Parallel Analysis</span>
              </div>
              <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <div style={{ 
                  width: '16px', 
                  height: '8px', 
                  background: '#F3E5F5',
                  border: '1px solid #9C27B0',
                  borderRadius: '2px'
                }} />
                <span>Machine Learning</span>
              </div>
              <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <div style={{ 
                  width: '16px', 
                  height: '8px', 
                  background: '#E0F2F1',
                  border: '1px solid #388E3C',
                  borderRadius: '2px'
                }} />
                <span>Output</span>
              </div>
            </div>
          </div>
        </div>

        {/* Technical Details Grid */}
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', gap: '2rem' }}>
          <div style={{
            background: '#FFFFFF',
            borderRadius: '0.75rem',
            padding: '2rem',
            border: '1px solid #E5E7EB',
            boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
          }}>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '1rem' }}>
              Quality Metrics
            </h3>
            <ul style={{ listStyle: 'none', padding: 0, color: '#4B5563', lineHeight: '1.8' }}>
              <li>• Minimum Quality Score: 30 (Phred)</li>
              <li>• Minimum Read Depth: 10x</li>
              <li>• Allele Frequency Threshold: 1%</li>
              <li>• Population Database: gnomAD v3.1</li>
              <li>• Reference Genome: GRCh38/hg38</li>
            </ul>
          </div>

          <div style={{
            background: '#FFFFFF',
            borderRadius: '0.75rem',
            padding: '2rem',
            border: '1px solid #E5E7EB',
            boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
          }}>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '1rem' }}>
              Clinical Databases
            </h3>
            <ul style={{ listStyle: 'none', padding: 0, color: '#4B5563', lineHeight: '1.8' }}>
              <li>• TCGA: 2,828 patient cohort</li>
              <li>• ClinVar: Latest release</li>
              <li>• CADD: v1.6 annotations</li>
              <li>• PRS Catalog: 500+ traits</li>
              <li>• KEGG: Pathway annotations</li>
            </ul>
          </div>

          <div style={{
            background: '#FFFFFF',
            borderRadius: '0.75rem',
            padding: '2rem',
            border: '1px solid #E5E7EB',
            boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
          }}>
            <h3 style={{ fontSize: '1.25rem', fontWeight: '600', color: '#111827', marginBottom: '1rem' }}>
              Processing Standards
            </h3>
            <ul style={{ listStyle: 'none', padding: 0, color: '#4B5563', lineHeight: '1.8' }}>
              <li>• ACMG/AMP variant classification</li>
              <li>• VCF 4.3 specification</li>
              <li>• HGVS nomenclature</li>
              <li>• CAP/CLIA compliance ready</li>
              <li>• HIPAA-compliant processing</li>
            </ul>
          </div>
        </div>
      </div>
    </section>

    <section style={{ background: '#FFFFFF', padding: '5rem 0' }}>
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
                Windows 10/11, macOS 10.15+, or Linux
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                8GB RAM minimum (16GB recommended)
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                2GB free storage space
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                Modern multi-core processor
              </li>
            </ul>
          </div>
          <div>
            <h3 style={{ fontSize: '1.5rem', fontWeight: '600', color: '#111827', marginBottom: '1rem' }}>
              File Requirements
            </h3>
            <ul style={{ listStyle: 'none', padding: 0 }}>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                FASTQ file format (.fastq or .fq)
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                Whole genome sequencing data
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                Minimum 30x coverage recommended
              </li>
              <li style={{ marginBottom: '0.75rem', color: '#4B5563' }}>
                Human genomic data (hg38 reference)
              </li>
            </ul>
          </div>
        </div>
      </div>
    </section>

    <section style={{ background: '#F9FAFB', padding: '5rem 0' }}>
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