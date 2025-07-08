import { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom'
import Layout from '../components/Layout'

interface Release {
  assets: Asset[];
  tag_name: string;
  name: string;
  published_at: string;
}

interface Asset {
  id: number;
  name: string;
  browser_download_url: string;
  size: number;
}

const WelcomePage = () => {
  const navigate = useNavigate()
  const [release, setRelease] = useState<Release | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const handleTestNow = () => {
    navigate('/upload')
  }

  useEffect(() => {
    const fetchLatestRelease = async () => {
      try {
        setLoading(true);
        // TODO: Replace 'YOUR_GITHUB_USERNAME' with the actual GitHub username or organization name.
        const response = await fetch('https://api.github.com/repos/trist/GeneKnow/releases/latest');
        if (!response.ok) {
          throw new Error(`Failed to fetch latest release: ${response.statusText}`);
        }
        const data = await response.json();
        setRelease(data);
      } catch (e: any) {
        setError(e.message);
        console.error(e);
      } finally {
        setLoading(false);
      }
    };

    fetchLatestRelease();
  }, []);


  return (
    <Layout>
      <section style={{
        background: 'radial-gradient(circle at top left, rgba(239, 246, 255, 1) 0%, rgba(255, 255, 255, 1) 50%)',
        padding: '5rem 0 8rem 0'
      }}>
        <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem', textAlign: 'center' }}>
          <div style={{ maxWidth: '48rem', margin: '0 auto' }}>
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
              marginTop: '1rem',
              fontSize: 'clamp(2.5rem, 5vw, 3.75rem)',
              fontWeight: 'bold',
              letterSpacing: '-0.02em',
              color: '#111827',
              lineHeight: '1.1'
            }}>
              Understand Your Genomic Health, Privately.
            </h1>
            <p style={{
              marginTop: '1.5rem',
              fontSize: '1.125rem',
              lineHeight: '1.75',
              color: '#4B5563',
              maxWidth: '42rem',
              margin: '1.5rem auto 0'
            }}>
              Our desktop application analyzes your genomic data file (`.fastq`) directly on your computer to provide a secure cancer risk assessment. No data uploads, no cloud processing, complete privacy.
            </p>
            <div style={{ marginTop: '2.5rem' }}>
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

      <section style={{ padding: '4rem 1.5rem', backgroundColor: '#F9FAFB' }}>
        <div style={{ maxWidth: '1200px', margin: '0 auto', textAlign: 'center' }}>
          <h2 style={{
            fontSize: '2.25rem',
            fontWeight: 'bold',
            color: '#111827',
            letterSpacing: '-0.02em',
          }}>
            Download Latest Version
          </h2>
          {release && (
            <p style={{
              marginTop: '0.5rem',
              color: '#4B5563'
            }}>
              Version {release.name} ({release.tag_name}) - Released on {new Date(release.published_at).toLocaleDateString()}
            </p>
          )}

          <div style={{
            marginTop: '2.5rem',
            display: 'flex',
            flexWrap: 'wrap',
            gap: '1.5rem',
            justifyContent: 'center'
          }}>
            {loading && <p>Loading download links...</p>}
            {error && <p style={{ color: 'red' }}>{error}</p>}
            {release?.assets && release.assets.map(asset => (
              <a
                key={asset.id}
                href={asset.browser_download_url}
                target="_blank"
                rel="noopener noreferrer"
                style={{
                  display: 'inline-block',
                  backgroundColor: '#FFFFFF',
                  color: '#1F2937',
                  border: '1px solid #E5E7EB',
                  borderRadius: '0.5rem',
                  padding: '1rem 2rem',
                  fontSize: '1rem',
                  fontWeight: '600',
                  cursor: 'pointer',
                  transition: 'all 200ms ease',
                  textDecoration: 'none',
                  boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
                  minWidth: '250px'
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.borderColor = '#D1D5DB';
                  e.currentTarget.style.boxShadow = '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.borderColor = '#E5E7EB';
                  e.currentTarget.style.boxShadow = '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)';
                }}
              >
                {asset.name}
                <span style={{ display: 'block', fontSize: '0.875rem', color: '#6B7280', marginTop: '0.25rem' }}>
                  {(asset.size / 1024 / 1024).toFixed(2)} MB
                </span>
              </a>
            ))}
            {release && release.assets.length === 0 && (
              <p>No downloads available for the latest release.</p>
            )}
          </div>
        </div>
      </section>
    </Layout>
  )
}

export default WelcomePage 