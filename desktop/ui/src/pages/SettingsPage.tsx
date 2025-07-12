import React, { useState, useEffect } from 'react'
import { useNavigate } from 'react-router-dom'
import { apiConfig } from '../api/apiConfig'

interface ReportGeneratorConfig {
  backend: 'ollama' | 'none'
  model_name: string | null
  temperature: number
  style: 'clinical' | 'technical' | 'patient'
  include_recommendations: boolean
  include_glossary: boolean
}

interface AvailableModels {
  ollama: string[]
}

// Add CSS animation for loading spinner
const spinnerStyle = `
  @keyframes spin {
    0% { transform: rotate(0deg); }
    100% { transform: rotate(360deg); }
  }
`

// Inject the CSS
if (typeof document !== 'undefined') {
  const style = document.createElement('style')
  style.textContent = spinnerStyle
  document.head.appendChild(style)
}

// Custom Dropdown Component
interface CustomDropdownProps {
  value: string
  onChange: (value: string) => void
  options: { value: string; label: string; description?: string }[]
  placeholder?: string
  disabled?: boolean
}

const CustomDropdown: React.FC<CustomDropdownProps> = ({ 
  value, 
  onChange, 
  options, 
  placeholder = 'Select an option',
  disabled = false 
}) => {
  const [isOpen, setIsOpen] = useState(false)
  const [selectedOption, setSelectedOption] = useState(
    options.find(opt => opt.value === value) || null
  )

  useEffect(() => {
    setSelectedOption(options.find(opt => opt.value === value) || null)
  }, [value, options])

  const handleSelect = (option: { value: string; label: string; description?: string }) => {
    setSelectedOption(option)
    onChange(option.value)
    setIsOpen(false)
  }

  if (disabled) {
    return (
      <div style={{
        width: '100%',
        padding: '12px 16px',
        border: '1px solid #E5E7EB',
        borderRadius: '8px',
        fontSize: '14px',
        background: '#F9FAFB',
        color: '#9CA3AF',
        cursor: 'not-allowed'
      }}>
        {selectedOption?.label || placeholder}
      </div>
    )
  }

  return (
    <div style={{ position: 'relative', width: '100%' }}>
      {/* Dropdown Button */}
      <button
        type="button"
        onClick={() => setIsOpen(!isOpen)}
        style={{
          width: '100%',
          padding: '12px 16px',
          border: `1px solid ${isOpen ? '#2563EB' : '#E5E7EB'}`,
          borderRadius: '8px',
          fontSize: '14px',
          background: '#FFFFFF',
          color: '#111827',
          cursor: 'pointer',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          transition: 'all 200ms ease',
          outline: 'none',
          boxShadow: isOpen ? '0 0 0 3px rgba(37, 99, 235, 0.1)' : 'none'
        }}
        onMouseEnter={(e) => {
          if (!isOpen) {
            e.currentTarget.style.borderColor = '#D1D5DB';
          }
        }}
        onMouseLeave={(e) => {
          if (!isOpen) {
            e.currentTarget.style.borderColor = '#E5E7EB';
          }
        }}
      >
        <span>{selectedOption?.label || placeholder}</span>
        <svg
          style={{
            width: '20px',
            height: '20px',
            transform: isOpen ? 'rotate(180deg)' : 'rotate(0deg)',
            transition: 'transform 200ms ease',
            color: '#6B7280'
          }}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Dropdown Menu */}
      {isOpen && (
        <div style={{
          position: 'absolute',
          top: '100%',
          left: 0,
          right: 0,
          marginTop: '4px',
          background: '#FFFFFF',
          border: '1px solid #E5E7EB',
          borderRadius: '8px',
          boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
          zIndex: 50,
          maxHeight: '300px',
          overflowY: 'auto'
        }}>
          {options.map((option, index) => (
            <button
              key={option.value}
              type="button"
              onClick={() => handleSelect(option)}
              style={{
                width: '100%',
                padding: '12px 16px',
                textAlign: 'left',
                border: 'none',
                background: selectedOption?.value === option.value 
                  ? '#F3F4F6'
                  : 'transparent',
                color: '#111827',
                cursor: 'pointer',
                transition: 'all 150ms ease',
                borderRadius: index === 0 ? '8px 8px 0 0' : index === options.length - 1 ? '0 0 8px 8px' : '0',
                fontSize: '14px',
                outline: 'none'
              }}
              onMouseEnter={(e) => {
                if (selectedOption?.value !== option.value) {
                  e.currentTarget.style.background = '#F9FAFB';
                }
              }}
              onMouseLeave={(e) => {
                if (selectedOption?.value !== option.value) {
                  e.currentTarget.style.background = 'transparent';
                }
              }}
            >
              <div style={{ fontWeight: selectedOption?.value === option.value ? '600' : '400', marginBottom: option.description ? '4px' : '0' }}>
                {option.label}
              </div>
              {option.description && (
                <div style={{ fontSize: '12px', color: '#6B7280', lineHeight: '1.4' }}>
                  {option.description}
                </div>
              )}
            </button>
          ))}
        </div>
      )}
    </div>
  )
}

export const SettingsPage: React.FC = () => {
  const navigate = useNavigate()
  const [config, setConfig] = useState<ReportGeneratorConfig>({
    backend: 'ollama',
    model_name: null,
    temperature: 0.3,
    style: 'clinical',
    include_recommendations: true,
    include_glossary: true
  })
  const [availableModels, setAvailableModels] = useState<AvailableModels>({
    ollama: []
  })
  const [isLoading, setIsLoading] = useState(true)
  const [isSaving, setIsSaving] = useState(false)
  const [saveMessage, setSaveMessage] = useState<string | null>(null)
  const [backendStatus, setBackendStatus] = useState<{
    ollama: boolean
  }>({
    ollama: false
  })
  const [modelWarmingStatus, setModelWarmingStatus] = useState<{
    isWarming: boolean
    warmingModel: string | null
    warmingProgress: string
    cachedModels: string[]
  }>({
    isWarming: false,
    warmingModel: null,
    warmingProgress: '',
    cachedModels: []
  })

  // Load current configuration and available models
  useEffect(() => {
    loadConfiguration()
    loadModelStatus()
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
          style: data.style || 'clinical',
          include_recommendations: data.include_recommendations || true,
          include_glossary: data.include_glossary || true
        })
      }

      // Load available models and backend status
      const modelsResponse = await fetch(`${apiConfig.getBaseUrl()}/api/report-generator/available-models`)
      if (modelsResponse.ok) {
        const modelsData = await modelsResponse.json()
        setAvailableModels(modelsData.models || { ollama: [] })
        setBackendStatus(modelsData.status || { ollama: false })
      }
    } catch (error) {
      console.error('Failed to load configuration:', error)
    } finally {
      setIsLoading(false)
    }
  }

  const loadModelStatus = async () => {
    try {
      const response = await fetch(`${apiConfig.getBaseUrl()}/api/report-generator/model-status`)
      if (response.ok) {
        const data = await response.json()
        setModelWarmingStatus(prev => ({
          ...prev,
          cachedModels: data.cached_models || []
        }))
      }
    } catch (error) {
      console.error('Failed to load model status:', error)
    }
  }

  const warmUpModel = async (modelName: string, backend: string) => {
    if (backend !== 'ollama' || !modelName || modelName === 'auto') {
      return
    }

    // Check if model is already cached
    if (modelWarmingStatus.cachedModels.includes(modelName)) {
      console.log(`Model ${modelName} is already cached`)
      return
    }

    setModelWarmingStatus(prev => ({
      ...prev,
      isWarming: true,
      warmingModel: modelName,
      warmingProgress: 'Initializing model...'
    }))

    try {
      const response = await fetch(`${apiConfig.getBaseUrl()}/api/report-generator/warm-model`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({
          model_name: modelName,
          backend: backend
        })
      })

      if (response.ok) {
        const data = await response.json()
        setModelWarmingStatus(prev => ({
          ...prev,
          isWarming: false,
          warmingModel: null,
          warmingProgress: 'Model ready! 🎉',
          cachedModels: [...prev.cachedModels, modelName]
        }))
        
        // Show a browser notification if possible (for when user has navigated away)
        if ('Notification' in window && Notification.permission === 'granted') {
          new Notification('GeneKnow Model Ready', {
            body: `${modelName} is now ready for report generation`,
            icon: '/favicon.ico'
          })
        }
        
        // Clear the success message after a few seconds
        setTimeout(() => {
          setModelWarmingStatus(prev => ({
            ...prev,
            warmingProgress: ''
          }))
        }, 5000)
      } else {
        throw new Error('Failed to warm up model')
      }
    } catch (error) {
      console.error('Failed to warm up model:', error)
      setModelWarmingStatus(prev => ({
        ...prev,
        isWarming: false,
        warmingModel: null,
        warmingProgress: 'Failed to warm up model ❌'
      }))
      
      // Clear the error message after a few seconds
      setTimeout(() => {
        setModelWarmingStatus(prev => ({
          ...prev,
          warmingProgress: ''
        }))
      }, 5000)
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

  const handleBackendChange = (backend: 'ollama' | 'none') => {
    const previousBackend = config.backend
    const previousModel = config.model_name
    
    setConfig(prev => ({
      ...prev,
      backend,
      model_name: null // Reset model selection when backend changes
    }))
    
    // If switching to Ollama and there was a previously selected model, 
    // we could potentially warm it up, but since we're resetting model_name,
    // the user will need to select a new model which will trigger warming
  }

  const handleModelChange = (modelName: string) => {
    const actualModelName = modelName === 'auto' ? null : modelName
    setConfig(prev => ({
      ...prev,
      model_name: actualModelName
    }))

    // Request notification permission for background loading updates
    if (config.backend === 'ollama' && actualModelName && 'Notification' in window) {
      if (Notification.permission === 'default') {
        Notification.requestPermission().then(permission => {
          if (permission === 'granted') {
            console.log('Notification permission granted')
          }
        })
      }
    }

    // Trigger model warming for Ollama models
    if (config.backend === 'ollama' && actualModelName) {
      warmUpModel(actualModelName, config.backend)
    }
  }

  const getRecommendedModels = (backend: 'ollama' | 'none') => {
    const models = {
      ollama: ['llama3', 'mistral', 'codellama', 'phi', 'neural-chat']
    }
    return backend === 'none' ? [] : models[backend] || []
  }

  const getModelDescription = (modelName: string) => {
    const descriptions: Record<string, string> = {
      // Ollama models
      'llama3': 'Latest Llama model, excellent for medical reports',
      'mistral': 'Fast and efficient, good balance of speed and quality',
      'codellama': 'Specialized for technical content',
      'phi': 'Smaller model, faster generation',
      'neural-chat': 'Optimized for conversational style',
      'llama2': 'Previous generation, still very capable',
      'vicuna': 'Good general-purpose model',
      'orca-mini': 'Compact model for quick generation'
    }
    
    return descriptions[modelName] || 'General purpose model'
  }

  const getModelOptions = () => {
    if (config.backend === 'none') return []
    
    const models = availableModels.ollama || []
    const recommended = getRecommendedModels(config.backend)
    
    return [
      { value: 'auto', label: 'Auto-detect best model', description: 'Automatically select the best available model' },
      ...models.map(model => {
        const isRecommended = recommended.includes(model)
        const description = getModelDescription(model)
        
        return {
          value: model,
          label: model + (isRecommended ? ' (Recommended)' : ''),
          description: description
        }
      })
    ]
  }

  const getStyleOptions = () => [
    { 
      value: 'clinical', 
      label: 'Clinical', 
      description: 'Medical professionals - detailed clinical language' 
    },
    { 
      value: 'technical', 
      label: 'Technical', 
      description: 'Researchers/Scientists - comprehensive technical details' 
    },
    { 
      value: 'patient', 
      label: 'Patient', 
      description: 'General audience - accessible language and explanations' 
    }
  ]

  if (isLoading) {
    return (
      <div style={{
        minHeight: '100vh',
        background: '#F8FAFC',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center'
      }}>
        <div style={{
          background: '#FFFFFF',
          border: '1px solid #E5E7EB',
          borderRadius: '12px',
          padding: '24px',
          boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
          textAlign: 'center'
        }}>
          <div style={{ fontSize: '16px', fontWeight: '500', color: '#111827' }}>
            Loading Settings...
          </div>
        </div>
      </div>
    )
  }

  return (
    <div style={{
      minHeight: '100vh',
      background: '#F8FAFC'
    }}>
      {/* Header */}
      <div style={{
        background: '#FFFFFF',
        borderBottom: '1px solid #E5E7EB',
        padding: '16px 24px',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between'
      }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '16px' }}>
          <button
            onClick={() => navigate(-1)}
            style={{
              background: 'none',
              border: 'none',
              color: '#6B7280',
              fontSize: '14px',
              cursor: 'pointer',
              display: 'flex',
              alignItems: 'center',
              gap: '8px',
              padding: '8px 12px',
              borderRadius: '6px',
              transition: 'all 200ms ease'
            }}
            onMouseEnter={(e) => {
              e.currentTarget.style.background = '#F3F4F6';
              e.currentTarget.style.color = '#374151';
            }}
            onMouseLeave={(e) => {
              e.currentTarget.style.background = 'none';
              e.currentTarget.style.color = '#6B7280';
            }}
          >
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
            </svg>
            Back
          </button>
          
          <div style={{
            background: '#DBEAFE',
            color: '#2563EB',
            padding: '4px 12px',
            borderRadius: '16px',
            fontSize: '12px',
            fontWeight: '600'
          }}>
            AI Model Configuration
          </div>
        </div>
      </div>

      {/* Main Content */}
      <div style={{ padding: '24px', maxWidth: '1200px', margin: '0 auto' }}>
        {/* Page Title */}
        <div style={{ marginBottom: '32px' }}>
          <h1 style={{
            fontSize: '32px',
            fontWeight: '700',
            color: '#111827',
            marginBottom: '8px'
          }}>
            Settings
          </h1>
          <p style={{
            fontSize: '16px',
            color: '#6B7280',
            lineHeight: '1.5'
          }}>
            Configure your preferred AI model for generating genomic reports.
          </p>
        </div>

        {/* Settings Card */}
        <div style={{
          background: '#FFFFFF',
          border: '1px solid #E5E7EB',
          borderRadius: '12px',
          boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
          overflow: 'hidden'
        }}>
          
          {/* AI Backend Section */}
          <div style={{ padding: '24px', borderBottom: '1px solid #F3F4F6' }}>
            <h2 style={{ 
              fontSize: '18px', 
              fontWeight: '600', 
              color: '#111827',
              marginBottom: '16px'
            }}>
              AI Backend
            </h2>
            
            <div style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
              {/* Ollama Option */}
              <div style={{
                display: 'flex',
                alignItems: 'center',
                padding: '16px',
                border: `2px solid ${config.backend === 'ollama' ? '#2563EB' : '#E5E7EB'}`,
                borderRadius: '8px',
                background: config.backend === 'ollama' ? '#F0F9FF' : '#FFFFFF',
                cursor: 'pointer',
                transition: 'all 200ms ease'
              }}
              onClick={() => handleBackendChange('ollama')}
              onMouseEnter={(e) => {
                if (config.backend !== 'ollama') {
                  e.currentTarget.style.borderColor = '#D1D5DB';
                }
              }}
              onMouseLeave={(e) => {
                if (config.backend !== 'ollama') {
                  e.currentTarget.style.borderColor = '#E5E7EB';
                }
              }}>
                <input
                  type="radio"
                  name="backend"
                  value="ollama"
                  checked={config.backend === 'ollama'}
                  onChange={() => handleBackendChange('ollama')}
                  style={{ 
                    marginRight: '12px',
                    width: '16px',
                    height: '16px',
                    accentColor: '#2563EB'
                  }}
                />
                <div style={{ flex: 1 }}>
                  <div style={{ 
                    display: 'flex', 
                    alignItems: 'center', 
                    gap: '8px',
                    marginBottom: '4px'
                  }}>
                    <span style={{ fontWeight: '600', fontSize: '16px', color: '#111827' }}>
                      Ollama (Local)
                    </span>
                    <span style={{
                      padding: '2px 8px',
                      fontSize: '12px',
                      fontWeight: '500',
                      borderRadius: '12px',
                      background: backendStatus.ollama ? '#DCFCE7' : '#FEF2F2',
                      color: backendStatus.ollama ? '#166534' : '#991B1B'
                    }}>
                      {backendStatus.ollama ? 'Available' : 'Not Available'}
                    </span>
                    <a
                      href="https://ollama.com"
                      target="_blank"
                      rel="noopener noreferrer"
                      onClick={(e) => e.stopPropagation()}
                      style={{
                        padding: '2px 8px',
                        fontSize: '12px',
                        fontWeight: '500',
                        borderRadius: '12px',
                        background: '#E0E7FF',
                        color: '#3730A3',
                        textDecoration: 'none',
                        transition: 'all 200ms ease'
                      }}
                      onMouseEnter={(e) => {
                        e.currentTarget.style.background = '#C7D2FE';
                      }}
                      onMouseLeave={(e) => {
                        e.currentTarget.style.background = '#E0E7FF';
                      }}
                    >
                      Download Ollama
                    </a>
                  </div>
                  <div style={{ fontSize: '14px', color: '#6B7280', lineHeight: '1.4' }}>
                    Run models locally with Ollama. Recommended for privacy and performance.
                  </div>
                </div>
              </div>

              {/* Template-Based Option */}
              <div style={{
                display: 'flex',
                alignItems: 'center',
                padding: '16px',
                border: `2px solid ${config.backend === 'none' ? '#2563EB' : '#E5E7EB'}`,
                borderRadius: '8px',
                background: config.backend === 'none' ? '#F0F9FF' : '#FFFFFF',
                cursor: 'pointer',
                transition: 'all 200ms ease'
              }}
              onClick={() => handleBackendChange('none')}
              onMouseEnter={(e) => {
                if (config.backend !== 'none') {
                  e.currentTarget.style.borderColor = '#D1D5DB';
                }
              }}
              onMouseLeave={(e) => {
                if (config.backend !== 'none') {
                  e.currentTarget.style.borderColor = '#E5E7EB';
                }
              }}>
                <input
                  type="radio"
                  name="backend"
                  value="none"
                  checked={config.backend === 'none'}
                  onChange={() => handleBackendChange('none')}
                  style={{ 
                    marginRight: '12px',
                    width: '16px',
                    height: '16px',
                    accentColor: '#2563EB'
                  }}
                />
                <div style={{ flex: 1 }}>
                  <div style={{ 
                    display: 'flex', 
                    alignItems: 'center', 
                    gap: '8px',
                    marginBottom: '4px'
                  }}>
                    <span style={{ fontWeight: '600', fontSize: '16px', color: '#111827' }}>
                      Template-Based (No AI)
                    </span>
                    <span style={{
                      padding: '2px 8px',
                      fontSize: '12px',
                      fontWeight: '500',
                      borderRadius: '12px',
                      background: '#DCFCE7',
                      color: '#166534'
                    }}>
                      Always Available
                    </span>
                    <span style={{
                      padding: '2px 8px',
                      fontSize: '12px',
                      fontWeight: '500',
                      borderRadius: '12px',
                      background: '#E0E7FF',
                      color: '#3730A3'
                    }}>
                      Instant
                    </span>
                  </div>
                  <div style={{ fontSize: '14px', color: '#6B7280', lineHeight: '1.4' }}>
                    Generate reports using templates without AI assistance. Fast and reliable.
                  </div>
                </div>
              </div>
            </div>
            
            {/* Installation Help Section */}
            <div style={{ 
              marginTop: '16px', 
              padding: '16px', 
              background: '#F8FAFC', 
              borderRadius: '8px',
              border: '1px solid #E2E8F0'
            }}>
              <h3 style={{ fontSize: '14px', fontWeight: '600', color: '#111827', marginBottom: '8px' }}>
                🚀 Quick Setup Guide
              </h3>
              <div style={{ fontSize: '12px', color: '#4B5563', lineHeight: '1.4' }}>
                <p style={{ margin: '0' }}>
                  <strong>Ollama:</strong> Download and install from{' '}
                  <a href="https://ollama.com" target="_blank" rel="noopener noreferrer" style={{ color: '#2563EB' }}>
                    ollama.com
                  </a>
                  , then run <code style={{ background: '#E5E7EB', padding: '2px 4px', borderRadius: '3px' }}>ollama pull llama3</code> in terminal.
                </p>
              </div>
            </div>
          </div>

          {/* Model Selection */}
          {config.backend !== 'none' && (
            <div style={{ padding: '24px', borderBottom: '1px solid #F3F4F6' }}>
              <h2 style={{ 
                fontSize: '18px', 
                fontWeight: '600', 
                color: '#111827',
                marginBottom: '8px'
              }}>
                Model Selection
              </h2>
              <p style={{ fontSize: '14px', color: '#6B7280', marginBottom: '16px' }}>
                Choose which specific model to use for report generation.
              </p>
              
              <CustomDropdown
                value={config.model_name || 'auto'}
                onChange={handleModelChange}
                options={getModelOptions()}
                placeholder="Select a model"
                disabled={!backendStatus.ollama}
              />
              
              <div style={{ fontSize: '12px', color: '#6B7280', marginTop: '8px' }}>
                {config.model_name 
                  ? `Using specific model: ${config.model_name}`
                  : 'Will automatically select the best available model'
                }
              </div>
              
              {/* Model Warming Status */}
              {config.backend === 'ollama' && modelWarmingStatus.isWarming && (
                <div style={{ marginTop: '12px' }}>
                  <div style={{
                    padding: '12px',
                    background: '#FEF3C7',
                    borderRadius: '8px',
                    fontSize: '14px',
                    color: '#92400E'
                  }}>
                    <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                      <span>🔄</span>
                      <span>Loading model: {modelWarmingStatus.warmingModel}</span>
                    </div>
                    <div style={{ fontSize: '12px', marginTop: '4px', color: '#B45309' }}>
                      This may take 1-2 minutes for first-time loading...
                    </div>
                  </div>
                </div>
              )}
              
              {modelWarmingStatus.warmingProgress && !modelWarmingStatus.isWarming && (
                <div style={{ 
                  padding: '8px 12px', 
                  background: modelWarmingStatus.warmingProgress.includes('ready') ? '#F0FDF4' : '#FEF2F2', 
                  borderRadius: '6px',
                  border: `1px solid ${modelWarmingStatus.warmingProgress.includes('ready') ? '#BBF7D0' : '#FECACA'}`,
                  display: 'flex',
                  alignItems: 'center',
                  gap: '8px'
                }}>
                  <div style={{ fontSize: '12px', color: modelWarmingStatus.warmingProgress.includes('ready') ? '#166534' : '#991B1B', fontWeight: '500' }}>
                    {modelWarmingStatus.warmingProgress}
                  </div>
                </div>
              )}
              
              {modelWarmingStatus.cachedModels.length > 0 && (
                <div style={{ 
                  marginTop: '6px',
                  fontSize: '11px', 
                  color: '#059669',
                  fontWeight: '500'
                }}>
                  ✅ Ready models: {modelWarmingStatus.cachedModels.join(', ')}
                </div>
              )}
            </div>
          )}

          {/* Report Style */}
          <div style={{ padding: '24px', borderBottom: '1px solid #F3F4F6' }}>
            <h2 style={{ 
              fontSize: '18px', 
              fontWeight: '600', 
              color: '#111827',
              marginBottom: '8px'
            }}>
              Report Style
            </h2>
            <p style={{ fontSize: '14px', color: '#6B7280', marginBottom: '16px' }}>
              Choose the writing style and target audience for your reports.
            </p>
            
            <CustomDropdown
              value={config.style}
              onChange={(value) => setConfig(prev => ({ ...prev, style: value as 'clinical' | 'technical' | 'patient' }))}
              options={getStyleOptions()}
              placeholder="Select report style"
            />
          </div>

          {/* Advanced Settings Section */}
          <div style={{ 
            padding: '20px',
            background: '#F9FAFB',
            borderRadius: '8px',
            marginTop: '24px'
          }}>
            <h3 style={{ fontSize: '16px', fontWeight: '600', color: '#111827', marginBottom: '16px' }}>
              Advanced Settings
            </h3>
            
            <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
              {/* Temperature Slider */}
              <div>
                <label style={{ 
                  display: 'flex', 
                  alignItems: 'center',
                  fontSize: '14px', 
                  fontWeight: '500', 
                  color: '#374151',
                  marginBottom: '8px'
                }}>
                  Temperature
                  <span style={{ 
                    marginLeft: '8px',
                    padding: '2px 8px',
                    fontSize: '12px',
                    background: '#E5E7EB',
                    borderRadius: '12px',
                    color: '#6B7280'
                  }}>
                    {config.temperature}
                  </span>
                </label>
                <input
                  type="range"
                  min="0.1"
                  max="1.0"
                  step="0.1"
                  value={config.temperature}
                  onChange={(e) => setConfig(prev => ({ ...prev, temperature: parseFloat(e.target.value) }))}
                  style={{ width: '100%' }}
                />
                <div style={{ 
                  display: 'flex', 
                  justifyContent: 'space-between',
                  fontSize: '12px',
                  color: '#9CA3AF',
                  marginTop: '4px'
                }}>
                  <span>More Focused</span>
                  <span>More Creative</span>
                </div>
              </div>

              {/* Include Recommendations */}
              <label style={{ 
                display: 'flex', 
                alignItems: 'center',
                fontSize: '14px',
                cursor: 'pointer'
              }}>
                <input
                  type="checkbox"
                  checked={config.include_recommendations}
                  onChange={(e) => setConfig(prev => ({ ...prev, include_recommendations: e.target.checked }))}
                  style={{ marginRight: '8px' }}
                />
                <span style={{ color: '#374151', fontWeight: '500' }}>Include Clinical Recommendations</span>
              </label>

              {/* Include Glossary */}
              <label style={{ 
                display: 'flex', 
                alignItems: 'center',
                fontSize: '14px',
                cursor: 'pointer'
              }}>
                <input
                  type="checkbox"
                  checked={config.include_glossary}
                  onChange={(e) => setConfig(prev => ({ ...prev, include_glossary: e.target.checked }))}
                  style={{ marginRight: '8px' }}
                />
                <span style={{ color: '#374151', fontWeight: '500' }}>Include Medical Glossary</span>
              </label>
            </div>
          </div>

          {/* Save Section */}
          <div style={{ padding: '24px', display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
            <div>
              {saveMessage && (
                <div style={{
                  padding: '8px 12px',
                  borderRadius: '6px',
                  fontSize: '14px',
                  fontWeight: '500',
                  background: saveMessage.includes('successfully') ? '#DCFCE7' : '#FEF2F2',
                  color: saveMessage.includes('successfully') ? '#166534' : '#991B1B'
                }}>
                  {saveMessage}
                </div>
              )}
            </div>
            
            <button
              onClick={handleSave}
              disabled={isSaving}
              style={{
                background: isSaving ? '#9CA3AF' : '#2563EB',
                color: '#FFFFFF',
                border: 'none',
                borderRadius: '8px',
                padding: '12px 24px',
                fontSize: '14px',
                fontWeight: '600',
                cursor: isSaving ? 'not-allowed' : 'pointer',
                transition: 'all 200ms ease',
                outline: 'none'
              }}
              onMouseEnter={(e) => {
                if (!isSaving) {
                  e.currentTarget.style.background = '#1D4ED8';
                }
              }}
              onMouseLeave={(e) => {
                if (!isSaving) {
                  e.currentTarget.style.background = '#2563EB';
                }
              }}
            >
              {isSaving ? 'Saving...' : 'Save Settings'}
            </button>
          </div>
        </div>
      </div>
    </div>
  )
}

export default SettingsPage 