/**
 * GeneKnow Pipeline Tauri Integration
 * Provides integration with the Tauri backend for genomic processing
 */

import { invoke } from '@tauri-apps/api/core';
import type { 
  Job, 
  PipelineResult, 
  UserPreferences, 
  JobProgress,
  HealthStatus
} from './geneknowPipeline';

export class GeneKnowTauriClient {
  private progressCallbacks: Map<string, (progress: JobProgress) => void> = new Map();
  private pollIntervals: Map<string, NodeJS.Timeout> = new Map();

  // Ensure API server is running
  async ensureApiServerRunning(): Promise<boolean> {
    try {
      console.log('Checking API server health...');
      const isHealthy = await invoke<boolean>('check_api_health');
      console.log('API health check result:', isHealthy);
      
      if (!isHealthy) {
        console.log('API server not healthy, attempting to start...');
        const startResult = await invoke<boolean>('start_api_server');
        console.log('API server start result:', startResult);
        
        if (!startResult) {
          throw new Error('Failed to start API server');
        }
        
        return startResult;
      }
      return true;
    } catch (error) {
      console.error('Failed to ensure API server is running:', error);
      console.error('Error type:', typeof error);
      console.error('Error message:', error instanceof Error ? error.message : String(error));
      throw new Error('Failed to start GeneKnow API server');
    }
  }

  // Get API server status
  async getApiServerStatus(): Promise<HealthStatus> {
    try {
      return await invoke('get_api_server_status');
    } catch (error) {
      console.error('Failed to get API server status:', error);
      throw error;
    }
  }

  // Process a genomic file
  async processGenomicFile(
    filePath: string, 
    preferences?: UserPreferences
  ): Promise<string> {
    try {
      console.log('Processing genomic file:', filePath);
      console.log('Preferences:', preferences);
      
      // Ensure server is running
      console.log('Ensuring API server is running...');
      await this.ensureApiServerRunning();
      console.log('API server is running');
      
      // Process the file
      console.log('Calling Tauri process_genomic_file command...');
      const jobId = await invoke<string>('process_genomic_file', {
        filePath,
        preferences: preferences || {}
      });
      
      console.log('Tauri command returned:', jobId);
      console.log('Job ID type:', typeof jobId);
      
      if (!jobId || typeof jobId !== 'string' || jobId.trim() === '') {
        throw new Error('Invalid job ID returned from Tauri command');
      }
      
      return jobId;
    } catch (error) {
      console.error('Failed to process genomic file:', error);
      console.error('Error details:', error);
      throw error;
    }
  }

  // Get job status
  async getJobStatus(jobId: string): Promise<Job> {
    try {
      const status = await invoke<Job>('get_job_status', { jobId });
      return status;
    } catch (error) {
      console.error('Failed to get job status:', error);
      throw error;
    }
  }

  // Get job results
  async getJobResults(jobId: string): Promise<PipelineResult> {
    try {
      console.log('Getting job results for job ID:', jobId);
      const results = await invoke<PipelineResult>('get_job_results', { jobId });
      console.log('Raw results from invoke:', results);
      console.log('Results keys:', results ? Object.keys(results) : 'no results');
      
      if (!results) {
        throw new Error('No results returned from API');
      }
      
      return results;
    } catch (error) {
      console.error('Failed to get job results:', error);
      console.error('Error details:', error);
      throw error;
    }
  }

  // Subscribe to job progress with polling
  subscribeToJobProgress(jobId: string, callback: (progress: JobProgress) => void): void {
    // Store callback
    this.progressCallbacks.set(jobId, callback);
    
    // Start polling for status
    const pollInterval = setInterval(async () => {
      try {
        const job = await this.getJobStatus(jobId);
        const progress: JobProgress = {
          job_id: jobId,
          status: job.status,
          progress: job.progress,
          current_step: job.current_step
        };
        
        callback(progress);
        
        // Stop polling if job is complete
        if (job.status === 'completed' || job.status === 'failed' || job.status === 'cancelled') {
          this.unsubscribeFromJobProgress(jobId);
        }
      } catch (error) {
        console.error('Error polling job status:', error);
      }
    }, 1000); // Poll every second
    
    this.pollIntervals.set(jobId, pollInterval);
  }

  // Unsubscribe from job progress
  unsubscribeFromJobProgress(jobId: string): void {
    // Clear callback
    this.progressCallbacks.delete(jobId);
    
    // Clear poll interval
    const interval = this.pollIntervals.get(jobId);
    if (interval) {
      clearInterval(interval);
      this.pollIntervals.delete(jobId);
    }
  }

  // Wait for job completion
  async waitForJobCompletion(
    jobId: string,
    onProgress?: (progress: JobProgress) => void
  ): Promise<Job> {
    return new Promise((resolve, reject) => {
      console.log('Waiting for job completion:', jobId);
      
      // Subscribe to progress updates
      this.subscribeToJobProgress(jobId, (progress) => {
        console.log('Job progress update:', progress);
        
        if (onProgress) {
          onProgress(progress);
        }
        
        // Check if job is complete
        if (progress.status === 'completed') {
          console.log('Job completed, getting final status');
          this.unsubscribeFromJobProgress(jobId);
          this.getJobStatus(jobId).then(finalStatus => {
            console.log('Final job status:', finalStatus);
            resolve(finalStatus);
          }).catch(reject);
        } else if (progress.status === 'failed' || progress.status === 'cancelled') {
          console.log('Job failed or cancelled:', progress.status);
          this.unsubscribeFromJobProgress(jobId);
          this.getJobStatus(jobId).then(job => {
            console.log('Failed job details:', job);
            reject(new Error(job.error || `Job ${job.status}`));
          }).catch(reject);
        }
      });
    });
  }

  // Process file and wait for results
  async processFileAndWaitForResults(
    filePath: string,
    preferences?: UserPreferences,
    onProgress?: (progress: JobProgress) => void
  ): Promise<PipelineResult> {
    try {
      console.log('Starting processFileAndWaitForResults for:', filePath);
      
      // Start processing
      console.log('Step 1: Starting file processing...');
      const jobId = await this.processGenomicFile(filePath, preferences);
      console.log('Step 1 Complete - Got job ID:', jobId);
      
      if (!jobId) {
        throw new Error('No job ID returned from file processing');
      }
      
      // Wait for completion
      console.log('Step 2: Waiting for job completion...');
      const job = await this.waitForJobCompletion(jobId, onProgress);
      console.log('Step 2 Complete - Job completed:', job);
      
      if (!job) {
        throw new Error('No job data returned from completion wait');
      }
      
      // Get and return results
      console.log('Step 3: Getting job results...');
      const results = await this.getJobResults(jobId);
      console.log('Step 3 Complete - Got results:', results);
      console.log('Results type:', typeof results);
      console.log('Results structure:', results ? Object.keys(results) : 'no results');
      
      if (!results) {
        throw new Error('No results returned from job');
      }
      
      return results;
    } catch (error) {
      console.error('Error in processFileAndWaitForResults:', error);
      console.error('Error type:', typeof error);
      console.error('Error message:', error instanceof Error ? error.message : String(error));
      console.error('Error stack:', error instanceof Error ? error.stack : 'No stack trace');
      throw error;
    }
  }

  // Stop API server (useful for cleanup)
  async stopApiServer(): Promise<void> {
    try {
      await invoke('stop_api_server');
    } catch (error) {
      console.error('Failed to stop API server:', error);
    }
  }
}

// Default client instance
export const geneKnowTauriClient = new GeneKnowTauriClient();

// React Hook for Tauri integration
export function useGeneKnowTauri() {
  const client = geneKnowTauriClient;

  return {
    client,
    ensureApiRunning: () => client.ensureApiServerRunning(),
    getApiStatus: () => client.getApiServerStatus(),
    processFile: (filePath: string, preferences?: UserPreferences) => 
      client.processGenomicFile(filePath, preferences),
    getJobStatus: (jobId: string) => client.getJobStatus(jobId),
    getJobResults: (jobId: string) => client.getJobResults(jobId),
    subscribeToProgress: (jobId: string, callback: (progress: JobProgress) => void) => 
      client.subscribeToJobProgress(jobId, callback),
    processAndWait: (filePath: string, preferences?: UserPreferences, onProgress?: (progress: JobProgress) => void) =>
      client.processFileAndWaitForResults(filePath, preferences, onProgress),
    stopApiServer: () => client.stopApiServer()
  };
} 