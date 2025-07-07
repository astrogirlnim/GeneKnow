import React, { useState, useRef } from 'react'
import { invoke } from '@tauri-apps/api/core'
import { 
  Upload, 
  FileText, 
  Activity, 
  Shield, 
  Zap, 
  Database,
  ChevronRight,
  AlertCircle,
  CheckCircle,
  Clock,
  Dna
} from 'lucide-react'

interface FileUploadStatus {
  file: File | null
  status: 'idle' | 'uploading' | 'processing' | 'completed' | 'error'
  progress: number
  error?: string
}

interface RiskAssessment {
  overallRisk: 'low' | 'medium' | 'high'
  confidence: number
  variants: number
  genes: string[]
  riskScore: number
}

function App() {
  const [uploadStatus, setUploadStatus] = useState<FileUploadStatus>({
    file: null,
    status: 'idle',
    progress: 0
  })
  const [riskAssessment, setRiskAssessment] = useState<RiskAssessment | null>(null)
  const [isDragOver, setIsDragOver] = useState(false)
  const fileInputRef = useRef<HTMLInputElement>(null)

  console.log('🧬 App component rendered - GenePredict v0.1.0')
  console.log('📊 Current upload status:', uploadStatus.status)
  console.log('🔬 Risk assessment data:', riskAssessment)

  const handleFileUpload = async (file: File) => {
    console.log('📁 File upload initiated:', file.name, file.size, 'bytes')
    
    setUploadStatus({ file, status: 'uploading', progress: 0 })
    
    try {
      // Simulate file validation
      if (!file.name.match(/\.(vcf|bam|fastq|fq)$/i)) {
        throw new Error('Invalid file type. Please upload VCF, BAM, or FASTQ files.')
      }
      
      console.log('✅ File validation passed')
      
      // Simulate upload progress
      for (let i = 0; i <= 100; i += 10) {
        await new Promise(resolve => setTimeout(resolve, 100))
        setUploadStatus(prev => ({ ...prev, progress: i }))
      }
      
      console.log('🔄 Processing genomic data...')
      setUploadStatus(prev => ({ ...prev, status: 'processing' }))
      
      // Simulate processing time
      await new Promise(resolve => setTimeout(resolve, 2000))
      
      // Mock risk assessment results
      const mockRiskAssessment: RiskAssessment = {
        overallRisk: 'medium',
        confidence: 87,
        variants: 1247,
        genes: ['BRCA1', 'BRCA2', 'TP53', 'PALB2', 'ATM'],
        riskScore: 6.3
      }
      
      console.log('🎯 Risk assessment completed:', mockRiskAssessment)
      setRiskAssessment(mockRiskAssessment)
      setUploadStatus(prev => ({ ...prev, status: 'completed' }))
      
    } catch (error) {
      console.error('❌ File processing error:', error)
      setUploadStatus(prev => ({ 
        ...prev, 
        status: 'error', 
        error: error instanceof Error ? error.message : 'Unknown error'
      }))
    }
  }

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault()
    setIsDragOver(true)
  }

  const handleDragLeave = (e: React.DragEvent) => {
    e.preventDefault()
    setIsDragOver(false)
  }

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault()
    setIsDragOver(false)
    
    const files = Array.from(e.dataTransfer.files)
    if (files.length > 0) {
      handleFileUpload(files[0])
    }
  }

  const handleFileSelect = (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = e.target.files
    if (files && files.length > 0) {
      handleFileUpload(files[0])
    }
  }

  const resetUpload = () => {
    console.log('🔄 Resetting upload state')
    setUploadStatus({ file: null, status: 'idle', progress: 0 })
    setRiskAssessment(null)
    if (fileInputRef.current) {
      fileInputRef.current.value = ''
    }
  }

  const StatusIcon = ({ status }: { status: FileUploadStatus['status'] }) => {
    switch (status) {
      case 'uploading':
      case 'processing':
        return <Clock className="w-5 h-5 text-blue-500 animate-spin" />
      case 'completed':
        return <CheckCircle className="w-5 h-5 text-green-500" />
      case 'error':
        return <AlertCircle className="w-5 h-5 text-red-500" />
      default:
        return <Upload className="w-5 h-5 text-gray-400" />
    }
  }

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100">
      {/* Header */}
      <header className="bg-white shadow-sm border-b border-gray-200">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex items-center justify-between h-16">
            <div className="flex items-center space-x-3">
                             <div className="p-2 bg-primary-100 rounded-lg">
                 <Dna className="w-6 h-6 text-primary-600" />
               </div>
              <div>
                <h1 className="text-xl font-bold text-gray-900">GenePredict</h1>
                <p className="text-sm text-gray-500">AI for Genomic Risk Assessment</p>
              </div>
            </div>
            <div className="flex items-center space-x-4">
              <div className="flex items-center space-x-2 px-3 py-1 bg-green-100 rounded-full">
                <Shield className="w-4 h-4 text-green-600" />
                <span className="text-sm font-medium text-green-800">Local Only</span>
              </div>
              <div className="flex items-center space-x-2 px-3 py-1 bg-blue-100 rounded-full">
                <Zap className="w-4 h-4 text-blue-600" />
                <span className="text-sm font-medium text-blue-800">Offline</span>
              </div>
            </div>
          </div>
        </div>
      </header>

      {/* Main Content */}
      <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {/* Upload Section */}
          <div className="lg:col-span-2">
            <div className="card">
              <div className="mb-6">
                <h2 className="text-xl font-semibold text-gray-900 mb-2">Upload Genomic Data</h2>
                <p className="text-gray-600">
                  Upload your VCF, BAM, or FASTQ files to begin risk assessment analysis.
                </p>
              </div>

              {uploadStatus.status === 'idle' && (
                <div
                  className={`file-upload-zone ${isDragOver ? 'dragover' : ''}`}
                  onDragOver={handleDragOver}
                  onDragLeave={handleDragLeave}
                  onDrop={handleDrop}
                  onClick={() => fileInputRef.current?.click()}
                >
                  <Upload className="w-12 h-12 text-gray-400 mx-auto mb-4" />
                  <p className="text-lg font-medium text-gray-900 mb-2">
                    Drag & drop your genomic files here
                  </p>
                  <p className="text-gray-500 mb-4">
                    or click to browse your computer
                  </p>
                  <div className="flex flex-wrap justify-center gap-2 text-sm text-gray-400">
                    <span className="px-2 py-1 bg-gray-100 rounded">VCF</span>
                    <span className="px-2 py-1 bg-gray-100 rounded">BAM</span>
                    <span className="px-2 py-1 bg-gray-100 rounded">FASTQ</span>
                  </div>
                  <input
                    ref={fileInputRef}
                    type="file"
                    accept=".vcf,.bam,.fastq,.fq"
                    onChange={handleFileSelect}
                    className="hidden"
                  />
                </div>
              )}

              {uploadStatus.status !== 'idle' && (
                <div className="space-y-4">
                  <div className="flex items-center space-x-3 p-4 bg-gray-50 rounded-lg">
                    <StatusIcon status={uploadStatus.status} />
                    <div className="flex-1">
                      <p className="font-medium text-gray-900">
                        {uploadStatus.file?.name}
                      </p>
                      <p className="text-sm text-gray-500">
                        {uploadStatus.status === 'uploading' && 'Uploading...'}
                        {uploadStatus.status === 'processing' && 'Processing genomic data...'}
                        {uploadStatus.status === 'completed' && 'Analysis complete'}
                        {uploadStatus.status === 'error' && 'Error occurred'}
                      </p>
                    </div>
                  </div>

                  {uploadStatus.status === 'uploading' && (
                    <div className="w-full bg-gray-200 rounded-full h-2">
                      <div 
                        className="bg-primary-600 h-2 rounded-full transition-all duration-300"
                        style={{ width: `${uploadStatus.progress}%` }}
                      />
                    </div>
                  )}

                  {uploadStatus.status === 'error' && (
                    <div className="p-3 bg-red-50 border border-red-200 rounded-lg">
                      <p className="text-sm text-red-800">{uploadStatus.error}</p>
                    </div>
                  )}

                  {uploadStatus.status !== 'processing' && (
                    <button
                      onClick={resetUpload}
                      className="btn-outline"
                    >
                      Upload Another File
                    </button>
                  )}
                </div>
              )}
            </div>
          </div>

          {/* Info Panel */}
          <div className="space-y-6">
            <div className="card">
              <h3 className="text-lg font-semibold text-gray-900 mb-4">Features</h3>
              <div className="space-y-3">
                <div className="flex items-start space-x-3">
                  <Activity className="w-5 h-5 text-primary-600 mt-0.5" />
                  <div>
                    <p className="font-medium text-gray-900">AI Risk Assessment</p>
                    <p className="text-sm text-gray-500">TensorFlow-powered genomic analysis</p>
                  </div>
                </div>
                <div className="flex items-start space-x-3">
                  <Shield className="w-5 h-5 text-green-600 mt-0.5" />
                  <div>
                    <p className="font-medium text-gray-900">Privacy First</p>
                    <p className="text-sm text-gray-500">All processing happens locally</p>
                  </div>
                </div>
                <div className="flex items-start space-x-3">
                  <FileText className="w-5 h-5 text-blue-600 mt-0.5" />
                  <div>
                    <p className="font-medium text-gray-900">Interpretable Reports</p>
                    <p className="text-sm text-gray-500">Clear, actionable insights</p>
                  </div>
                </div>
                <div className="flex items-start space-x-3">
                  <Database className="w-5 h-5 text-purple-600 mt-0.5" />
                  <div>
                    <p className="font-medium text-gray-900">Reference Data</p>
                    <p className="text-sm text-gray-500">1000 Genomes & ClinVar</p>
                  </div>
                </div>
              </div>
            </div>

            {riskAssessment && (
              <div className="card">
                <h3 className="text-lg font-semibold text-gray-900 mb-4">Risk Assessment</h3>
                <div className="space-y-4">
                  <div className="p-4 bg-gradient-to-r from-blue-50 to-indigo-50 rounded-lg">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-sm font-medium text-gray-700">Overall Risk</span>
                      <span className={`px-2 py-1 text-xs font-medium rounded-full ${
                        riskAssessment.overallRisk === 'low' ? 'bg-green-100 text-green-800' :
                        riskAssessment.overallRisk === 'medium' ? 'bg-yellow-100 text-yellow-800' :
                        'bg-red-100 text-red-800'
                      }`}>
                        {riskAssessment.overallRisk.toUpperCase()}
                      </span>
                    </div>
                    <div className="text-2xl font-bold text-gray-900 mb-1">
                      {riskAssessment.riskScore}/10
                    </div>
                    <div className="text-sm text-gray-500">
                      {riskAssessment.confidence}% confidence
                    </div>
                  </div>
                  
                  <div className="grid grid-cols-2 gap-4">
                    <div className="text-center p-3 bg-gray-50 rounded-lg">
                      <div className="text-xl font-bold text-gray-900">
                        {riskAssessment.variants.toLocaleString()}
                      </div>
                      <div className="text-sm text-gray-500">Variants</div>
                    </div>
                    <div className="text-center p-3 bg-gray-50 rounded-lg">
                      <div className="text-xl font-bold text-gray-900">
                        {riskAssessment.genes.length}
                      </div>
                      <div className="text-sm text-gray-500">Risk Genes</div>
                    </div>
                  </div>
                  
                  <div>
                    <h4 className="font-medium text-gray-900 mb-2">Key Genes</h4>
                    <div className="flex flex-wrap gap-2">
                      {riskAssessment.genes.map(gene => (
                        <span key={gene} className="px-2 py-1 bg-blue-100 text-blue-800 text-xs font-medium rounded">
                          {gene}
                        </span>
                      ))}
                    </div>
                  </div>
                  
                  <button className="btn-primary w-full mt-4">
                    Generate Full Report
                    <ChevronRight className="w-4 h-4 ml-1" />
                  </button>
                </div>
              </div>
            )}
          </div>
        </div>
      </main>
    </div>
  )
}

export default App
