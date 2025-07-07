import { useState, useEffect } from 'react'
import reactLogo from './assets/react.svg'
import viteLogo from '/vite.svg'
import { useLogger } from './hooks/useLogger'

function App() {
  const [count, setCount] = useState(0)
  const logger = useLogger()

  useEffect(() => {
    logger.info('GenePredict application initialized', { version: '0.1.0' })
  }, [logger])

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100 flex items-center justify-center">
      <div className="max-w-2xl mx-auto p-8 text-center">
        <div className="flex justify-center items-center space-x-4 mb-8">
          <a href="https://vite.dev" target="_blank" className="hover:opacity-75 transition-opacity">
            <img src={viteLogo} className="h-16 w-16" alt="Vite logo" />
          </a>
          <div className="text-4xl">🧬</div>
          <a href="https://react.dev" target="_blank" className="hover:opacity-75 transition-opacity">
            <img src={reactLogo} className="h-16 w-16 animate-spin-slow" alt="React logo" />
          </a>
        </div>
        
        <h1 className="text-5xl font-bold text-gray-800 mb-4">
          <span className="text-blue-600">Gene</span>Predict
        </h1>
        <p className="text-lg text-gray-600 mb-8">
          AI-Powered Genomic Risk Assessment Platform
        </p>
        
        <div className="bg-white rounded-lg shadow-lg p-6 mb-8">
          <button 
            onClick={() => {
              const newCount = count + 1
              setCount(newCount)
              logger.debug('Sample analysis counter incremented', { count: newCount })
            }}
            className="bg-blue-600 hover:bg-blue-700 text-white font-semibold py-3 px-6 rounded-lg transition-colors duration-200 shadow-md hover:shadow-lg"
          >
            Sample Analysis Count: {count}
          </button>
          <p className="mt-4 text-gray-500">
            Edit <code className="bg-gray-100 px-2 py-1 rounded">src/App.tsx</code> and save to test hot reload
          </p>
        </div>
        
        <div className="text-sm text-gray-500">
          <p>🔬 Secure • 🏠 Local Processing • 🌍 Open Source</p>
        </div>
      </div>
    </div>
  )
}

export default App
