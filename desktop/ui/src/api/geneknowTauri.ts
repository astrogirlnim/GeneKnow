/**
 * GeneKnow Pipeline Tauri Integration
 * Provides integration with the Tauri backend for genomic processing
 */

import { invoke } from '@tauri-apps/api/core';
import type { 
  HealthStatus, 
  Job, 
  PipelineResult, 
  UserPreferences, 
  JobProgress,
  PipelineInfo
} from './geneknowPipeline';

export class GeneKnowTauriClient {
  private progressCallbacks: Map<string, (progress: JobProgress) => void> = new Map();
  private pollIntervals: Map<string, NodeJS.Timeout> = new Map();

  // Ensure API server is running
  async ensureApiServerRunning(): Promise<boolean> {
    try {
      const isHealthy = await invoke<boolean>('check_api_health');
      if (!isHealthy) {
        return await invoke<boolean>('start_api_server');
      }
      return true;
    } catch (error) {
      console.error('Failed to ensure API server is running:', error);
      throw new Error('Failed to start GeneKnow API server');
    }
  }

  // Get API server status
  async getApiServerStatus(): Promise<any> {
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
      // Ensure server is running
      await this.ensureApiServerRunning();
      
      // Process the file
      const jobId = await invoke<string>('process_genomic_file', {
        filePath,
        preferences: preferences || {}
      });
      
      return jobId;
    } catch (error) {
      console.error('Failed to process genomic file:', error);
      throw error;
    }
  }

  // Get job status
  async getJobStatus(jobId: string): Promise<Job> {
    try {
      const status = await invoke<any>('get_job_status', { jobId });
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
      const results = await invoke<any>('get_job_results', { jobId });
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
      const jobId = await this.processGenomicFile(filePath, preferences);
      console.log('Got job ID:', jobId);
      
      // Wait for completion
      const job = await this.waitForJobCompletion(jobId, onProgress);
      console.log('Job completed:', job);
      
      // Get and return results
      const results = await this.getJobResults(jobId);
      console.log('Got results:', results);
      console.log('Results type:', typeof results);
      
      return results;
    } catch (error) {
      console.error('Error in processFileAndWaitForResults:', error);
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