/**
 * Example component demonstrating GeneKnow Pipeline API usage
 * Shows how the API auto-starts without manual Python commands
 */
import { useState, useEffect } from 'react';
import { invoke } from '@tauri-apps/api/core';
// TODO: For Tauri v2, dialog API might need different import or plugin
// For now, using a file input as fallback
import { useGeneKnowPipeline } from '../api/geneknowPipeline';

export function ApiTestExample() {
  const { client } = useGeneKnowPipeline();
  const [apiStatus, setApiStatus] = useState<'checking' | 'running' | 'stopped' | 'starting'>('checking');
  const [processing, setProcessing] = useState(false);
  const [progress, setProgress] = useState(0);
  const [currentStep, setCurrentStep] = useState('');
  const [results, setResults] = useState<any>(null);
  const [error, setError] = useState<string>('');

  // Check API status on mount
  useEffect(() => {
    checkApiStatus();
  }, []);

  const checkApiStatus = async () => {
    try {
      const status = await invoke<any>('get_api_server_status');
      setApiStatus(status.status === 'healthy' ? 'running' : 'stopped');
    } catch (e) {
      setApiStatus('stopped');
    }
  };

  const startApiServer = async () => {
    try {
      setApiStatus('starting');
      setError('');
      const started = await invoke<boolean>('start_api_server');
      
      if (started) {
        setApiStatus('running');
        // Connect WebSocket after server starts
        client.connectWebSocket();
      } else {
        throw new Error('Failed to start API server');
      }
    } catch (e: any) {
      setError(e.toString());
      setApiStatus('stopped');
    }
  };

  const selectAndProcessFile = async () => {
    try {
      setError('');
      
      // Ensure API is running (Tauri will auto-start if needed)
      if (apiStatus !== 'running') {
        await startApiServer();
      }

      // For demo purposes, using a hardcoded test file path
      // In a real app, you'd use Tauri's file dialog plugin
      // Example: const selected = await open({ ... });
      const testFilePath = '/path/to/test/file.vcf';
      
      // In production, replace the above with actual file selection:
      // const selected = await invoke<string>('select_file_dialog');
      // if (!selected) return;

      setProcessing(true);
      setProgress(0);
      setResults(null);

      // Process file through Tauri command (which auto-starts API if needed)
      const jobId = await invoke<string>('process_genomic_file', {
        filePath: testFilePath,
        preferences: {
          language: 'en',
          include_technical: true,
          patient_data: {
            age: 45,
            sex: 'F'
          }
        }
      });

      console.log('Job started:', jobId);

      // Subscribe to real-time progress updates
      client.subscribeToJobProgress(jobId, (update) => {
        setProgress(update.progress);
        setCurrentStep(update.current_step || '');
        console.log(`Progress: ${update.progress}% - ${update.current_step}`);
      });

      // Wait for completion
      let completed = false;
      while (!completed) {
        const status = await invoke<any>('get_job_status', { jobId });
        
        if (status.status === 'completed') {
          completed = true;
          const jobResults = await invoke<any>('get_job_results', { jobId });
          setResults(jobResults);
          setProcessing(false);
          client.unsubscribeFromJobProgress(jobId);
        } else if (status.status === 'failed') {
          throw new Error(status.error || 'Processing failed');
        }

        await new Promise(resolve => setTimeout(resolve, 1000));
      }

    } catch (error: any) {
      console.error('Processing error:', error);
      setProcessing(false);
      setError(error.toString());
    }
  };

  return (
    <div className="p-6 max-w-4xl mx-auto">
      <h2 className="text-2xl font-bold mb-4">GeneKnow Pipeline API Test</h2>
      
      {/* API Status */}
      <div className="mb-6 p-4 bg-gray-100 rounded">
        <h3 className="font-semibold mb-2">API Server Status</h3>
        <div className="flex items-center gap-4">
          <span className="flex items-center gap-2">
            <div className={`w-3 h-3 rounded-full ${
              apiStatus === 'running' ? 'bg-green-500' : 
              apiStatus === 'starting' ? 'bg-yellow-500 animate-pulse' : 
              'bg-red-500'
            }`} />
            {apiStatus === 'checking' && 'Checking...'}
            {apiStatus === 'running' && 'API Running'}
            {apiStatus === 'starting' && 'Starting API...'}
            {apiStatus === 'stopped' && 'API Stopped'}
          </span>
          
          {apiStatus === 'stopped' && (
            <button
              onClick={startApiServer}
              className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600"
            >
              Start API Server
            </button>
          )}
          
          <button
            onClick={checkApiStatus}
            className="px-4 py-2 bg-gray-500 text-white rounded hover:bg-gray-600"
          >
            Refresh Status
          </button>
        </div>
      </div>

      {/* File Processing */}
      <div className="mb-6">
        <button
          onClick={selectAndProcessFile}
          disabled={processing || apiStatus !== 'running'}
          className={`px-6 py-3 rounded font-semibold ${
            processing || apiStatus !== 'running'
              ? 'bg-gray-300 text-gray-500 cursor-not-allowed'
              : 'bg-blue-500 text-white hover:bg-blue-600'
          }`}
        >
          {processing ? 'Processing...' : 'Select & Process Genomic File'}
        </button>
      </div>

      {/* Progress */}
      {processing && (
        <div className="mb-6 p-4 bg-blue-50 rounded">
          <h3 className="font-semibold mb-2">Processing Progress</h3>
          <div className="mb-2">
            <div className="w-full bg-gray-200 rounded-full h-4">
              <div 
                className="bg-blue-500 h-4 rounded-full transition-all duration-300"
                style={{ width: `${progress}%` }}
              />
            </div>
          </div>
          <p className="text-sm text-gray-600">
            {progress}% - {currentStep || 'Initializing...'}
          </p>
        </div>
      )}

      {/* Results */}
      {results && (
        <div className="mb-6 p-4 bg-green-50 rounded">
          <h3 className="font-semibold mb-2">Results</h3>
          <div className="space-y-2 text-sm">
            <p><strong>Variants found:</strong> {results.variant_count || 0}</p>
            <p><strong>Processing time:</strong> {results.processing_time_seconds?.toFixed(2)}s</p>
            
            {results.risk_scores && Object.keys(results.risk_scores).length > 0 && (
              <div>
                <strong>Risk Scores:</strong>
                <ul className="ml-4 mt-1">
                  {Object.entries(results.risk_scores).map(([type, score]) => (
                    <li key={type}>
                      {type}: {((score as number) * 100).toFixed(1)}%
                    </li>
                  ))}
                </ul>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Error Display */}
      {error && (
        <div className="p-4 bg-red-50 text-red-700 rounded">
          <strong>Error:</strong> {error}
        </div>
      )}

      {/* Instructions */}
      <div className="mt-8 p-4 bg-gray-50 rounded text-sm text-gray-600">
        <h4 className="font-semibold mb-2">How it works:</h4>
        <ul className="list-disc ml-5 space-y-1">
          <li>The API server auto-starts when needed (no manual Python commands!)</li>
          <li>In production builds, it starts automatically when the app launches</li>
          <li>The server runs locally on port 5001</li>
          <li>WebSocket provides real-time progress updates</li>
          <li>The server stops automatically when the app closes</li>
        </ul>
      </div>
    </div>
  );
} 