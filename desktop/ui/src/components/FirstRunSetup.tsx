import React, { useState, useEffect, useCallback } from 'react';
import { invoke } from '@tauri-apps/api/core';

interface SetupStep {
  id: string;
  title: string;
  description: string;
  status: 'pending' | 'running' | 'completed' | 'error';
  error?: string;
}

interface FirstRunSetupProps {
  onComplete: () => void;
}

export const FirstRunSetup: React.FC<FirstRunSetupProps> = ({ onComplete }) => {
  const [steps, setSteps] = useState<SetupStep[]>([
    {
      id: 'check-env',
      title: 'Checking Environment',
      description: 'Verifying Python runtime and dependencies',
      status: 'pending'
    },
    {
      id: 'init-db',
      title: 'Creating Database',
      description: 'Building population variants database (this may take a few minutes)',
      status: 'pending'
    },
    {
      id: 'start-server',
      title: 'Starting Pipeline Server',
      description: 'Initializing the genomic analysis pipeline',
      status: 'pending'
    },
    {
      id: 'verify',
      title: 'Verifying Setup',
      description: 'Testing all components are working correctly',
      status: 'pending'
    }
  ]);

  // Track current step for debugging if needed
  // const [currentStep, setCurrentStep] = useState(0);
  const [progress, setProgress] = useState(0);

  const updateStep = useCallback((stepId: string, status: SetupStep['status'], error?: string) => {
    setSteps(prevSteps => 
      prevSteps.map(step => 
        step.id === stepId 
          ? { ...step, status, error } 
          : step
      )
    );
  }, []);

  const runSetup = useCallback(async () => {
    try {
      // Step 1: Check environment
      updateStep('check-env', 'running');
      
      const envCheck = await invoke<boolean>('check_environment');
      if (!envCheck) {
        throw new Error('Environment check failed');
      }
      
      updateStep('check-env', 'completed');
      setProgress(25);

      // Step 2: Initialize database if needed
      updateStep('init-db', 'running');
      
      const dbExists = await invoke<boolean>('check_database_exists');
      if (!dbExists) {
        console.log('Database not found, creating...');
        
        // Initialize database with progress updates
        await invoke('initialize_database', {
          cancerOnly: true,
          onProgress: (percent: number) => {
            setProgress(25 + (percent * 0.4)); // 25-65% for database
          }
        });
      }
      
      updateStep('init-db', 'completed');
      setProgress(65);

      // Step 3: Start API server
      updateStep('start-server', 'running');
      
      const serverStarted = await invoke<boolean>('start_api_server');
      if (!serverStarted) {
        throw new Error('Failed to start pipeline server');
      }
      
      // Wait for server to be ready
      let retries = 30; // 30 seconds timeout
      let serverReady = false;
      
      while (retries > 0 && !serverReady) {
        try {
          serverReady = await invoke<boolean>('check_api_health');
          if (!serverReady) {
            await new Promise(resolve => setTimeout(resolve, 1000));
            retries--;
          }
        } catch {
          await new Promise(resolve => setTimeout(resolve, 1000));
          retries--;
        }
        setProgress(65 + ((30 - retries) * 0.8)); // 65-89% for server startup
      }
      
      if (!serverReady) {
        throw new Error('Pipeline server failed to start');
      }
      
      updateStep('start-server', 'completed');
      setProgress(90);

      // Step 4: Verify everything works
      updateStep('verify', 'running');
      
      // Test a simple pipeline operation
      const testResult = await invoke<boolean>('test_pipeline_connectivity');
      if (!testResult) {
        throw new Error('Pipeline connectivity test failed');
      }
      
      updateStep('verify', 'completed');
      setProgress(100);

      // Complete setup after a brief delay
      setTimeout(() => {
        onComplete();
      }, 1000);

    } catch (error) {
      console.error('Setup failed:', error);
      const errorMessage = error instanceof Error ? error.message : 'Unknown error';
      
      // Update the failed step
      const failedStep = steps.find(s => s.status === 'running');
      if (failedStep) {
        updateStep(failedStep.id, 'error', errorMessage);
      }
    }
  }, [steps, updateStep, onComplete]);

  useEffect(() => {
    runSetup();
  }, [runSetup]);

  const getStepIcon = (status: SetupStep['status']) => {
    switch (status) {
      case 'completed':
        return (
          <svg className="w-6 h-6 text-green-500" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
          </svg>
        );
      case 'running':
        return (
          <svg className="animate-spin h-6 w-6 text-blue-500" fill="none" viewBox="0 0 24 24">
            <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
            <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
          </svg>
        );
      case 'error':
        return (
          <svg className="w-6 h-6 text-red-500" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
          </svg>
        );
      default:
        return (
          <div className="w-6 h-6 rounded-full border-2 border-gray-300 bg-white"></div>
        );
    }
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-indigo-50 via-white to-blue-50 flex items-center justify-center p-4">
      <div className="max-w-2xl w-full">
        <div className="bg-white rounded-2xl shadow-xl p-8">
          <div className="text-center mb-8">
            <h1 className="text-3xl font-bold text-gray-900 mb-2">
              Welcome to GeneKnow
            </h1>
            <p className="text-gray-600">
              Setting up your genomic analysis environment for the first time
            </p>
          </div>

          <div className="space-y-6 mb-8">
            {steps.map((step) => (
              <div key={step.id} className="flex items-start space-x-4">
                <div className="flex-shrink-0 mt-1">
                  {getStepIcon(step.status)}
                </div>
                <div className="flex-1">
                  <h3 className={`text-lg font-semibold ${
                    step.status === 'error' ? 'text-red-700' : 
                    step.status === 'completed' ? 'text-green-700' : 
                    step.status === 'running' ? 'text-blue-700' : 
                    'text-gray-700'
                  }`}>
                    {step.title}
                  </h3>
                  <p className="text-sm text-gray-600 mt-1">
                    {step.description}
                  </p>
                  {step.error && (
                    <p className="text-sm text-red-600 mt-2">
                      Error: {step.error}
                    </p>
                  )}
                </div>
              </div>
            ))}
          </div>

          <div className="mb-6">
            <div className="flex justify-between text-sm text-gray-600 mb-2">
              <span>Overall Progress</span>
              <span>{Math.round(progress)}%</span>
            </div>
            <div className="w-full bg-gray-200 rounded-full h-3">
              <div 
                className="bg-gradient-to-r from-blue-500 to-indigo-600 h-3 rounded-full transition-all duration-500"
                style={{ width: `${progress}%` }}
              ></div>
            </div>
          </div>

          {steps.some(s => s.status === 'error') && (
            <div className="bg-red-50 border border-red-200 rounded-lg p-4">
              <p className="text-red-800 text-sm">
                Setup encountered an error. Please check the logs or contact support.
              </p>
              <button 
                onClick={runSetup}
                className="mt-3 px-4 py-2 bg-red-600 text-white rounded-lg hover:bg-red-700 transition-colors"
              >
                Retry Setup
              </button>
            </div>
          )}

          <div className="text-center text-sm text-gray-500 mt-8">
            This is a one-time setup process. Your data never leaves your device.
          </div>
        </div>
      </div>
    </div>
  );
}; 