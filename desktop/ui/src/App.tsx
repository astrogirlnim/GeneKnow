import { useState, useEffect } from 'react'
import reactLogo from './assets/react.svg'
import viteLogo from '/vite.svg'
import { useLogger } from './hooks/useLogger'
import { GenomicProcessing } from './components/GenomicProcessing'

function App() {
  const [count, setCount] = useState(0)
  const logger = useLogger()

  useEffect(() => {
    logger.info('GenePredict application initialized', { version: '0.1.0' })
  }, [logger])

  return (
    <div className="app-container">
      <div className="main-content">
        <div className="logo-section">
          <a href="https://vite.dev" target="_blank" className="logo">
            <img src={viteLogo} className="logo" alt="Vite logo" />
          </a>
          <div className="dna-emoji">🧬</div>
          <a href="https://react.dev" target="_blank" className="logo">
            <img src={reactLogo} className="logo react" alt="React logo" />
          </a>
        </div>
        
        <h1 className="main-title">
          <span className="brand-name">Gene</span>Predict
        </h1>
        <p className="subtitle">
          AI-Powered Genomic Risk Assessment Platform
        </p>
        
        <div className="card">
          <button 
            onClick={() => {
              const newCount = count + 1
              setCount(newCount)
              logger.debug('Sample analysis counter incremented', { count: newCount })
            }}
            className="btn btn-primary btn-large"
          >
            Sample Analysis Count: {count}
          </button>
          <p className="mt-4 text-sm" style={{ color: 'var(--gray-500)' }}>
            Edit <code className="font-mono p-1 rounded" style={{ backgroundColor: 'var(--gray-100)' }}>src/App.tsx</code> and save to test hot reload
          </p>
        </div>
        
        <div className="text-sm" style={{ color: 'var(--gray-500)' }}>
          <p>🔬 Secure • 🏠 Local Processing • 🌍 Open Source</p>
        </div>
      </div>
    </div>
  )
}

export default App
