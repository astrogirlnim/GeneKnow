import React, { useState, useEffect } from 'react'
import { useNavigate } from 'react-router-dom'
import Layout from '../components/Layout'
import { apiConfig } from '../api/apiConfig'

interface ModelConfig {
  backend: 'ollama' | 'huggingface' | 'none'
  model_name: string | null
  temperature: number
  max_tokens: number
  style: 'clinician' | 'technical' | 'patient'
}

interface AvailableModels {
  ollama: string[]
  huggingface: string[]
}

const SettingsPage: React.FC = () => {
  const navigate = useNavigate()
  const [config, setConfig] = useState<ModelConfig>({
    backend: 'ollama',
    model_name: null,
    temperature: 0.3,
    max_tokens: 2000,
    style: 'clinician'
  })
  const [availableModels, setAvailableModels] = useState<AvailableModels>({
    ollama: [],
    huggingface: []
  })
  const [isLoading, setIsLoading] = useState(true)
  const [isSaving, setIsSaving] = useState(false)
  const [saveMessage, setSaveMessage] = useState<string | null>(null)
  const [backendStatus, setBackendStatus] = useState<{
    ollama: boolean
    huggingface: boolean
  }>({
    ollama: false,
    huggingface: false
  })

  // Load current configuration and available models
  useEffect(() => {
    loadConfiguration()
  }, [])

  const loadConfiguration = async () => {
    try {
      setIsLoading(true)
      
      // Load current config
      const response = await fetch(`${apiConfig.getBaseUrl()}/api/report-generator/config`)
      if (response.ok) {
        const data = await response.json()
        setConfig({
          backend: data.backend || 'ollama',
          model_name: data.model_name,
          temperature: data.temperature || 0.3,
          max_tokens: data.max_tokens || 2000,
          style: data.style || 'clinician'
        })
      }

      // Load available models and backend status
      const modelsResponse = await fetch(`${apiConfig.getBaseUrl()}/api/report-generator/available-models`)
      if (modelsResponse.ok) {
        const modelsData = await modelsResponse.json()
        setAvailableModels(modelsData.models || { ollama: [], huggingface: [] })
        setBackendStatus(modelsData.status || { ollama: false, huggingface: false })
      }
    } catch (error) {
      console.error('Failed to load configuration:', error)
    } finally {
      setIsLoading(false)
    }
  }

  const handleSave = async () => {
    try {
      setIsSaving(true)
      setSaveMessage(null)

      const response = await fetch(`${apiConfig.getBaseUrl()}/api/report-generator/config`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify(config)
      })

      if (response.ok) {
        setSaveMessage('Settings saved successfully!')
        setTimeout(() => setSaveMessage(null), 3000)
      } else {
        setSaveMessage('Failed to save settings. Please try again.')
      }
    } catch (error) {
      console.error('Failed to save configuration:', error)
      setSaveMessage('Failed to save settings. Please try again.')
    } finally {
      setIsSaving(false)
    }
  }

  const handleBackendChange = (backend: 'ollama' | 'huggingface' | 'none') => {
    setConfig(prev => ({
      ...prev,
      backend,
      model_name: null // Reset model selection when backend changes
    }))
  }

  const handleModelChange = (modelName: string) => {
    setConfig(prev => ({
      ...prev,
      model_name: modelName === 'auto' ? null : modelName
    }))
  }

  const getRecommendedModels = (backend: 'ollama' | 'huggingface') => {
    const recommendations = {
      ollama: ['llama3.1:8b', 'llama3', 'mistral', 'phi3'],
      huggingface: ['microsoft/phi-2', 'google/flan-t5-base', 'microsoft/DialoGPT-medium']
    }
    return recommendations[backend] || []
  }

    if (isLoading) {
    return (
      <Layout>
        <section style={{
          background: 'radial-gradient(circle at top left, rgba(239, 246, 255, 1) 0%, rgba(255, 255, 255, 1) 50%)',
          padding: '3rem 0',
          minHeight: 'calc(100vh - 4rem)',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center'
        }}>
          <div style={{
            background: 'rgba(255, 255, 255, 0.95)',
            backdropFilter: 'blur(10px)',
            border: '1px solid rgba(255, 255, 255, 0.2)',
            borderRadius: '1rem',
            padding: '2rem',
            boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
            textAlign: 'center'
          }}>
            <div style={{ 
              fontSize: '1.125rem', 
              color: '#111827',
              fontWeight: '600'
            }}>
              Loading settings...
            </div>
          </div>
        </section>
      </Layout>
    )
  }

  return (
    <Layout>
      <section style={{
        background: 'radial-gradient(circle at top left, rgba(239, 246, 255, 1) 0%, rgba(255, 255, 255, 1) 50%)',
        padding: '3rem 0 2rem 0'
      }}>
        <div style={{ maxWidth: '1200px', margin: '0 auto', padding: '0 1.5rem' }}>
          {/* Header */}
          <div style={{ textAlign: 'center', marginBottom: '3rem' }}>
            <button
              onClick={() => navigate(-1)}
              style={{
                display: 'inline-flex',
                alignItems: 'center',
                gap: '0.5rem',
                color: '#6B7280',
                background: 'none',
                border: 'none',
                cursor: 'pointer',
                fontSize: '0.875rem',
                marginBottom: '1.5rem',
                padding: '0.5rem 1rem',
                borderRadius: '0.5rem',
                transition: 'all 200ms ease'
              }}
              onMouseEnter={(e) => {
                e.currentTarget.style.background = 'rgba(255, 255, 255, 0.8)';
                e.currentTarget.style.color = '#111827';
              }}
              onMouseLeave={(e) => {
                e.currentTarget.style.background = 'none';
                e.currentTarget.style.color = '#6B7280';
              }}
            >
              ← Back
            </button>
            
            <span style={{
              color: '#2563EB',
              fontWeight: '600',
              background: '#DBEAFE',
              borderRadius: '9999px',
              padding: '0.5rem 1rem',
              fontSize: '0.875rem',
              display: 'inline-block',
              marginBottom: '1rem'
            }}>
              AI Model Configuration
            </span>
            
            <h1 style={{
              fontSize: 'clamp(2.5rem, 5vw, 3.75rem)',
              fontWeight: 'bold',
              letterSpacing: '-0.02em',
              color: '#111827',
              lineHeight: '1.1',
              marginBottom: '1rem'
            }}>
              Report Generation Settings
            </h1>
            
            <p style={{
              fontSize: '1.125rem',
              lineHeight: '1.75',
              color: '#4B5563',
              maxWidth: '42rem',
              margin: '0 auto'
            }}>
              Configure your preferred AI model for generating genomic reports. Choose from local models or template-based generation.
            </p>
          </div>

          {/* Settings Card */}
          <div style={{ maxWidth: '800px', margin: '0 auto' }}>
            <div style={{
              background: 'rgba(255, 255, 255, 0.95)',
              backdropFilter: 'blur(10px)',
              border: '1px solid rgba(255, 255, 255, 0.2)',
              borderRadius: '1rem',
              padding: '2rem',
              boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
              transition: 'all 200ms ease'
            }}>
              
              {/* Backend Selection */}
              <div style={{ marginBottom: '2rem' }}>
                <label style={{ 
                  display: 'block', 
                  fontSize: '1.125rem', 
                  fontWeight: '700', 
                  color: '#111827',
                  marginBottom: '1rem'
                }}>
                  AI Backend
                </label>
                <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                  {/* Ollama Option */}
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    padding: '1.5rem',
                    border: `2px solid ${config.backend === 'ollama' ? '#2563EB' : '#E5E7EB'}`,
                    borderRadius: '0.75rem',
                    background: config.backend === 'ollama' ? 'linear-gradient(135deg, rgba(59, 130, 246, 0.1) 0%, rgba(139, 92, 246, 0.1) 100%)' : '#FFFFFF',
                    cursor: 'pointer',
                    transition: 'all 200ms ease',
                    boxShadow: config.backend === 'ollama' ? '0 4px 6px -1px rgba(0, 0, 0, 0.1)' : '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
                  }}
                  onClick={() => handleBackendChange('ollama')}
                  onMouseEnter={(e) => {
                    if (config.backend !== 'ollama') {
                      e.currentTarget.style.transform = 'translateY(-2px)';
                      e.currentTarget.style.boxShadow = '0 4px 6px -1px rgba(0, 0, 0, 0.1)';
                    }
                  }}
                  onMouseLeave={(e) => {
                    if (config.backend !== 'ollama') {
                      e.currentTarget.style.transform = 'translateY(0)';
                      e.currentTarget.style.boxShadow = '0 1px 3px 0 rgba(0, 0, 0, 0.1)';
                    }
                  }}>
                    <input
                      type="radio"
                      name="backend"
                      value="ollama"
                      checked={config.backend === 'ollama'}
                      onChange={() => handleBackendChange('ollama')}
                      style={{ 
                        marginRight: '1rem',
                        width: '1.25rem',
                        height: '1.25rem',
                        accentColor: '#2563EB'
                      }}
                    />
                    <div style={{ flex: 1 }}>
                      <div style={{ 
                        display: 'flex', 
                        alignItems: 'center', 
                        gap: '0.75rem',
                        marginBottom: '0.5rem'
                      }}>
                        <span style={{ fontWeight: '700', fontSize: '1.125rem', color: '#111827' }}>
                          Ollama (Local)
                        </span>
                        <span style={{
                          padding: '0.25rem 0.75rem',
                          fontSize: '0.75rem',
                          fontWeight: '600',
                          borderRadius: '9999px',
                          background: backendStatus.ollama ? '#DCFCE7' : '#FEF2F2',
                          color: backendStatus.ollama ? '#166534' : '#991B1B'
                        }}>
                          {backendStatus.ollama ? 'Available' : 'Not Available'}
                        </span>
                      </div>
                      <div style={{ fontSize: '1rem', color: '#6B7280', lineHeight: '1.5' }}>
                        Run models locally with Ollama. Recommended for privacy and performance.
                      </div>
                    </div>
                  </div>

                  {/* Hugging Face Option */}
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    padding: '1.5rem',
                    border: `2px solid ${config.backend === 'huggingface' ? '#2563EB' : '#E5E7EB'}`,
                    borderRadius: '0.75rem',
                    background: config.backend === 'huggingface' ? 'linear-gradient(135deg, rgba(59, 130, 246, 0.1) 0%, rgba(139, 92, 246, 0.1) 100%)' : '#FFFFFF',
                    cursor: 'pointer',
                    transition: 'all 200ms ease',
                    boxShadow: config.backend === 'huggingface' ? '0 4px 6px -1px rgba(0, 0, 0, 0.1)' : '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
                  }}
                  onClick={() => handleBackendChange('huggingface')}
                  onMouseEnter={(e) => {
                    if (config.backend !== 'huggingface') {
                      e.currentTarget.style.transform = 'translateY(-2px)';
                      e.currentTarget.style.boxShadow = '0 4px 6px -1px rgba(0, 0, 0, 0.1)';
                    }
                  }}
                  onMouseLeave={(e) => {
                    if (config.backend !== 'huggingface') {
                      e.currentTarget.style.transform = 'translateY(0)';
                      e.currentTarget.style.boxShadow = '0 1px 3px 0 rgba(0, 0, 0, 0.1)';
                    }
                  }}>
                    <input
                      type="radio"
                      name="backend"
                      value="huggingface"
                      checked={config.backend === 'huggingface'}
                      onChange={() => handleBackendChange('huggingface')}
                      style={{ 
                        marginRight: '1rem',
                        width: '1.25rem',
                        height: '1.25rem',
                        accentColor: '#2563EB'
                      }}
                    />
                    <div style={{ flex: 1 }}>
                      <div style={{ 
                        display: 'flex', 
                        alignItems: 'center', 
                        gap: '0.75rem',
                        marginBottom: '0.5rem'
                      }}>
                        <span style={{ fontWeight: '700', fontSize: '1.125rem', color: '#111827' }}>
                          Hugging Face Transformers
                        </span>
                        <span style={{
                          padding: '0.25rem 0.75rem',
                          fontSize: '0.75rem',
                          fontWeight: '600',
                          borderRadius: '9999px',
                          background: backendStatus.huggingface ? '#DCFCE7' : '#FEF2F2',
                          color: backendStatus.huggingface ? '#166534' : '#991B1B'
                        }}>
                          {backendStatus.huggingface ? 'Available' : 'Not Available'}
                        </span>
                      </div>
                      <div style={{ fontSize: '1rem', color: '#6B7280', lineHeight: '1.5' }}>
                        Use Hugging Face models directly. Requires transformers library.
                      </div>
                    </div>
                  </div>

                  {/* None/Fallback Option */}
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    padding: '1.5rem',
                    border: `2px solid ${config.backend === 'none' ? '#2563EB' : '#E5E7EB'}`,
                    borderRadius: '0.75rem',
                    background: config.backend === 'none' ? 'linear-gradient(135deg, rgba(59, 130, 246, 0.1) 0%, rgba(139, 92, 246, 0.1) 100%)' : '#FFFFFF',
                    cursor: 'pointer',
                    transition: 'all 200ms ease',
                    boxShadow: config.backend === 'none' ? '0 4px 6px -1px rgba(0, 0, 0, 0.1)' : '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
                  }}
                  onClick={() => handleBackendChange('none')}
                  onMouseEnter={(e) => {
                    if (config.backend !== 'none') {
                      e.currentTarget.style.transform = 'translateY(-2px)';
                      e.currentTarget.style.boxShadow = '0 4px 6px -1px rgba(0, 0, 0, 0.1)';
                    }
                  }}
                  onMouseLeave={(e) => {
                    if (config.backend !== 'none') {
                      e.currentTarget.style.transform = 'translateY(0)';
                      e.currentTarget.style.boxShadow = '0 1px 3px 0 rgba(0, 0, 0, 0.1)';
                    }
                  }}>
                    <input
                      type="radio"
                      name="backend"
                      value="none"
                      checked={config.backend === 'none'}
                      onChange={() => handleBackendChange('none')}
                      style={{ 
                        marginRight: '1rem',
                        width: '1.25rem',
                        height: '1.25rem',
                        accentColor: '#2563EB'
                      }}
                    />
                    <div style={{ flex: 1 }}>
                      <div style={{ 
                        display: 'flex', 
                        alignItems: 'center', 
                        gap: '0.75rem',
                        marginBottom: '0.5rem'
                      }}>
                        <span style={{ fontWeight: '700', fontSize: '1.125rem', color: '#111827' }}>
                          Template-Based (No AI)
                        </span>
                        <span style={{
                          padding: '0.25rem 0.75rem',
                          fontSize: '0.75rem',
                          fontWeight: '600',
                          borderRadius: '9999px',
                          background: '#DCFCE7',
                          color: '#166534'
                        }}>
                          Always Available
                        </span>
                      </div>
                      <div style={{ fontSize: '1rem', color: '#6B7280', lineHeight: '1.5' }}>
                        Generate reports using templates without AI assistance.
                      </div>
                    </div>
                  </div>
                </div>
              </div>

              {/* Model Selection */}
              {config.backend !== 'none' && (
                <div style={{ marginBottom: '2rem' }}>
                  <label style={{ 
                    display: 'block', 
                    fontSize: '1.125rem', 
                    fontWeight: '700', 
                    color: '#111827',
                    marginBottom: '1rem'
                  }}>
                    Model Selection
                  </label>
                  <select
                    value={config.model_name || 'auto'}
                    onChange={(e) => handleModelChange(e.target.value)}
                    style={{
                      width: '100%',
                      padding: '1rem',
                      border: '2px solid #E5E7EB',
                      borderRadius: '0.75rem',
                      fontSize: '1rem',
                      background: '#FFFFFF',
                      color: '#111827',
                      transition: 'all 200ms ease',
                      outline: 'none'
                    }}
                    onFocus={(e) => {
                      e.currentTarget.style.borderColor = '#2563EB';
                      e.currentTarget.style.boxShadow = '0 0 0 3px rgba(59, 130, 246, 0.1)';
                    }}
                    onBlur={(e) => {
                      e.currentTarget.style.borderColor = '#E5E7EB';
                      e.currentTarget.style.boxShadow = 'none';
                    }}
                  >
                    <option value="auto">Auto-detect best model</option>
                    {availableModels[config.backend].map(model => (
                      <option key={model} value={model}>
                        {model}
                        {config.backend !== 'none' && getRecommendedModels(config.backend).includes(model) && ' (Recommended)'}
                      </option>
                    ))}
                  </select>
                  <div style={{ fontSize: '0.875rem', color: '#6B7280', marginTop: '0.75rem', fontStyle: 'italic' }}>
                    {config.model_name 
                      ? `Using specific model: ${config.model_name}`
                      : 'Will automatically select the best available model'
                    }
                  </div>
                </div>
              )}

              {/* Report Style */}
              <div style={{ marginBottom: '2rem' }}>
                <label style={{ 
                  display: 'block', 
                  fontSize: '1.125rem', 
                  fontWeight: '700', 
                  color: '#111827',
                  marginBottom: '1rem'
                }}>
                  Report Style
                </label>
                <select
                  value={config.style}
                  onChange={(e) => setConfig(prev => ({ ...prev, style: e.target.value as 'clinician' | 'technical' | 'patient' }))}
                  style={{
                    width: '100%',
                    padding: '1rem',
                    border: '2px solid #E5E7EB',
                    borderRadius: '0.75rem',
                    fontSize: '1rem',
                    background: '#FFFFFF',
                    color: '#111827',
                    transition: 'all 200ms ease',
                    outline: 'none'
                  }}
                  onFocus={(e) => {
                    e.currentTarget.style.borderColor = '#2563EB';
                    e.currentTarget.style.boxShadow = '0 0 0 3px rgba(59, 130, 246, 0.1)';
                  }}
                  onBlur={(e) => {
                    e.currentTarget.style.borderColor = '#E5E7EB';
                    e.currentTarget.style.boxShadow = 'none';
                  }}
                >
                  <option value="clinician">Clinician (Medical professionals)</option>
                  <option value="technical">Technical (Researchers/Scientists)</option>
                  <option value="patient">Patient (General audience)</option>
                </select>
              </div>

              {/* Advanced Settings */}
              {config.backend !== 'none' && (
                <div style={{ marginBottom: '2rem' }}>
                  <h3 style={{ 
                    fontSize: '1.125rem', 
                    fontWeight: '700', 
                    color: '#111827',
                    marginBottom: '1rem'
                  }}>
                    Advanced Settings
                  </h3>
                  
                  {/* Temperature */}
                  <div style={{ marginBottom: '1.5rem' }}>
                    <label style={{ 
                      display: 'block', 
                      fontSize: '1rem', 
                      fontWeight: '600', 
                      color: '#111827',
                      marginBottom: '0.75rem'
                    }}>
                      Temperature: {config.temperature}
                    </label>
                    <input
                      type="range"
                      min="0"
                      max="1"
                      step="0.1"
                      value={config.temperature}
                      onChange={(e) => setConfig(prev => ({ ...prev, temperature: parseFloat(e.target.value) }))}
                      style={{ 
                        width: '100%',
                        height: '0.5rem',
                        borderRadius: '0.25rem',
                        background: '#E5E7EB',
                        outline: 'none',
                        accentColor: '#2563EB'
                      }}
                    />
                    <div style={{ 
                      fontSize: '0.875rem', 
                      color: '#6B7280', 
                      marginTop: '0.5rem',
                      fontStyle: 'italic'
                    }}>
                      Lower values = more focused, Higher values = more creative
                    </div>
                  </div>

                  {/* Max Tokens */}
                  <div style={{ marginBottom: '1.5rem' }}>
                    <label style={{ 
                      display: 'block', 
                      fontSize: '1rem', 
                      fontWeight: '600', 
                      color: '#111827',
                      marginBottom: '0.75rem'
                    }}>
                      Max Tokens
                    </label>
                    <input
                      type="number"
                      min="500"
                      max="4000"
                      value={config.max_tokens}
                      onChange={(e) => setConfig(prev => ({ ...prev, max_tokens: parseInt(e.target.value) }))}
                      style={{
                        width: '100%',
                        padding: '0.75rem',
                        border: '2px solid #E5E7EB',
                        borderRadius: '0.5rem',
                        fontSize: '1rem',
                        background: '#FFFFFF',
                        color: '#111827',
                        transition: 'all 200ms ease',
                        outline: 'none'
                      }}
                      onFocus={(e) => {
                        e.currentTarget.style.borderColor = '#2563EB';
                        e.currentTarget.style.boxShadow = '0 0 0 3px rgba(59, 130, 246, 0.1)';
                      }}
                      onBlur={(e) => {
                        e.currentTarget.style.borderColor = '#E5E7EB';
                        e.currentTarget.style.boxShadow = 'none';
                      }}
                    />
                    <div style={{ 
                      fontSize: '0.875rem', 
                      color: '#6B7280', 
                      marginTop: '0.5rem',
                      fontStyle: 'italic'
                    }}>
                      Maximum length of generated reports
                    </div>
                  </div>
                </div>
              )}

              {/* Save Button */}
              <div style={{ 
                display: 'flex', 
                justifyContent: 'space-between', 
                alignItems: 'center',
                paddingTop: '1.5rem',
                borderTop: '2px solid #F3F4F6'
              }}>
                {saveMessage && (
                  <div style={{
                    padding: '0.75rem 1rem',
                    borderRadius: '0.5rem',
                    background: saveMessage.includes('success') ? '#DCFCE7' : '#FEF2F2',
                    color: saveMessage.includes('success') ? '#166534' : '#991B1B',
                    fontSize: '0.875rem',
                    fontWeight: '600'
                  }}>
                    {saveMessage}
                  </div>
                )}
                <button
                  onClick={handleSave}
                  disabled={isSaving}
                  style={{
                    padding: '0.75rem 2rem',
                    background: isSaving ? '#9CA3AF' : 'linear-gradient(135deg, #2563EB 0%, #1D4ED8 100%)',
                    color: '#FFFFFF',
                    border: 'none',
                    borderRadius: '0.5rem',
                    fontSize: '1rem',
                    fontWeight: '600',
                    cursor: isSaving ? 'not-allowed' : 'pointer',
                    transition: 'all 200ms ease',
                    marginLeft: 'auto',
                    boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)'
                  }}
                  onMouseEnter={(e) => {
                    if (!isSaving) {
                      e.currentTarget.style.background = 'linear-gradient(135deg, #1D4ED8 0%, #1E40AF 100%)';
                      e.currentTarget.style.transform = 'translateY(-1px)';
                      e.currentTarget.style.boxShadow = '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)';
                    }
                  }}
                  onMouseLeave={(e) => {
                    if (!isSaving) {
                      e.currentTarget.style.background = 'linear-gradient(135deg, #2563EB 0%, #1D4ED8 100%)';
                      e.currentTarget.style.transform = 'translateY(0)';
                      e.currentTarget.style.boxShadow = '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)';
                    }
                  }}
                >
                  {isSaving ? 'Saving...' : 'Save Settings'}
                </button>
              </div>
            </div>

            {/* Help Section */}
            <div style={{
              background: 'rgba(249, 250, 251, 0.8)',
              backdropFilter: 'blur(10px)',
              border: '1px solid rgba(229, 231, 235, 0.5)',
              borderRadius: '1rem',
              padding: '2rem',
              marginTop: '2rem'
            }}>
              <h3 style={{ 
                fontSize: '1.25rem', 
                fontWeight: '700', 
                color: '#111827',
                marginBottom: '1rem'
              }}>
                🚀 Getting Started
              </h3>
              <div style={{ fontSize: '1rem', color: '#4B5563', lineHeight: '1.6' }}>
                <div style={{ marginBottom: '1rem' }}>
                  <strong style={{ color: '#111827' }}>Ollama:</strong> Install Ollama and pull models like "llama3" or "mistral" for local AI processing.
                </div>
                <div style={{ marginBottom: '1rem' }}>
                  <strong style={{ color: '#111827' }}>Hugging Face:</strong> Requires the transformers library. Models will be downloaded automatically.
                </div>
                <div>
                  <strong style={{ color: '#111827' }}>Template-Based:</strong> Uses pre-defined templates without AI for basic reports.
                </div>
              </div>
            </div>
          </div>
        </div>
      </section>
    </Layout>
  )
}

export default SettingsPage 