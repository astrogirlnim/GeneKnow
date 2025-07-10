/**
 * GeneKnow Pipeline API Client
 * Provides type-safe access to the LangGraph genomic processing pipeline
 */

import { io, Socket } from 'socket.io-client';
import { apiConfig, ensureApiConfig } from './apiConfig';
import React from 'react';

// Types and Interfaces

export interface FileFormat {
  extension: string;
  description: string;
  compressed: string[];
  paired_end_support: boolean;
}

export interface PipelineNode {
  id: string;
  name: string;
  description: string;
}

export interface PipelineCapabilities {
  supported_formats: string[];
  max_file_size_gb: number;
  pipeline_nodes: PipelineNode[];
  cancer_types: string[];
  output_formats: string[];
  languages: string[];
}

export interface PipelineInfo {
  name: string;
  version: string;
  capabilities: PipelineCapabilities;
}

export interface UserPreferences {
  language?: string;
  include_technical?: boolean;
  patient_data?: {
    age?: number;
    sex?: 'M' | 'F';
    family_history?: boolean;
  };
}

export interface Job {
  id: string;
  status: 'pending' | 'processing' | 'completed' | 'failed' | 'cancelled';
  file_type: string;
  preferences: UserPreferences;
  created_at: string;
  started_at?: string;
  completed_at?: string;
  progress: number;
  current_step?: string;
  result?: JobResult;
  error?: string;
}

// Define proper types for report sections
export interface ReportSection {
  title: string;
  content: string;
  severity?: 'low' | 'medium' | 'high';
  technical_details?: string;
}

export interface JobResult {
  variant_count: number;
  risk_scores: Record<string, number>;
  report_sections: Record<string, ReportSection>;
  processing_time: number;
}

export interface JobProgress {
  job_id: string;
  status: string;
  progress: number;
  current_step?: string;
}

export interface HealthStatus {
  status: string;
  timestamp: string;
  service: string;
  version: string;
  jobs_active: number;
}

// Define proper types for pipeline results
export interface PipelineResult {
  pipeline_status: 'completed' | 'failed' | 'cancelled';
  variant_count: number;
  risk_scores: Record<string, number>;
  report_sections: Record<string, ReportSection>;
  processing_time_seconds: number;
  cadd_stats?: {
    variants_scored: number;
    mean_phred: number;
    max_phred: number;
    variants_gt20: number;
    variants_in_cancer_genes: number;
  };
  structured_json?: {
    cadd_summary?: {
      enabled: boolean;
      variants_scored: number;
      mean_phred_score: number;
      max_phred_score: number;
      high_impact_variants: number;
      cancer_gene_variants: number;
      description: string;
    };
    metrics?: {
      timestamp: string;
      pipeline_version: string;
      confidence_metrics: {
        mean_model_confidence: number;
        risk_score_range: number;
        ml_fusion_confidence: number;
      };
      variant_metrics: {
        total_variants: number;
        pathogenic_variants: number;
        benign_variants: number;
        uncertain_variants: number;
        high_cadd_variants: number;
        genes_affected: number;
        mean_cadd_score: number;
        max_cadd_score: number;
      };
      prediction_metrics: {
        max_risk_score: number;
        high_risk_cancer_count: number;
        ml_risk_category: string;
        aggregate_confidence: number;
      };
      prs_metrics: {
        high_prs_cancer_count: number;
        max_prs_percentile: number;
        mean_prs_score: number;
        overall_confidence: string;
      };
      pathway_metrics: {
        high_burden_pathway_count: number;
        mean_pathway_burden: number;
        pathway_risk_level: string;
      };
      overall_assessment: {
        high_risk_cancers: string[];
        max_risk_score: number;
        clinical_action_needed: boolean;
        risk_category: string;
      };
      performance_indicators: {
        variant_coverage: boolean;
        model_confidence_adequate: boolean;
        sufficient_evidence: boolean;
      };
      validation_metrics: {
        ground_truth_available: boolean;
        validation_ready: boolean;
        auc_roc?: number;
        sensitivity?: number;
        specificity?: number;
        f1_score?: number;
        mcc?: number;
      };
    };
    metrics_summary?: {
      highest_risk_cancer: string;
      highest_risk_score: number;
      pathogenic_variant_count: number;
      confidence_level: string;
      validation_ready: boolean;
    };
    [key: string]: unknown;
  };
  variants?: Array<{
    gene: string;
    position: string;
    type: string;
    impact: string;
  }>;
  errors?: string[];
}

// Define WebSocket error type
export interface WebSocketError {
  message: string;
  code?: string;
  type?: string;
}

// API Client Class
export class GeneKnowPipelineClient {
  private socket: Socket | null = null;
  private progressCallbacks: Map<string, (progress: JobProgress) => void> = new Map();
  private baseUrl: string = '';
  private initialized: boolean = false;

  constructor() {
    // URL will be set dynamically
  }

  private async ensureInitialized(): Promise<void> {
    if (!this.initialized) {
      await ensureApiConfig();
      this.baseUrl = apiConfig.getBaseUrl();
      this.initialized = true;
    }
  }

  // Health Check
  async checkHealth(): Promise<HealthStatus> {
    await this.ensureInitialized();
    const response = await fetch(`${this.baseUrl}/api/health`);
    if (!response.ok) throw new Error('Health check failed');
    return response.json();
  }

  // Pipeline Information
  async getPipelineInfo(): Promise<PipelineInfo> {
    await this.ensureInitialized();
    const response = await fetch(`${this.baseUrl}/api/pipeline-info`);
    if (!response.ok) throw new Error('Failed to get pipeline info');
    return response.json();
  }

  // Supported Formats
  async getSupportedFormats(): Promise<{ formats: FileFormat[] }> {
    const response = await fetch(`${this.baseUrl}/api/supported-formats`);
    if (!response.ok) throw new Error('Failed to get supported formats');
    return response.json();
  }

  // Upload and Process File
  async uploadFile(file: File, preferences?: UserPreferences): Promise<{
    job_id: string;
    filename: string;
    file_size: number;
    file_type: string;
    status: string;
    message: string;
  }> {
    await this.ensureInitialized();
    const formData = new FormData();
    formData.append('file', file);
    if (preferences) {
      formData.append('preferences', JSON.stringify(preferences));
    }

    const response = await fetch(`${this.baseUrl}/api/upload`, {
      method: 'POST',
      body: formData,
    });

    if (!response.ok) {
      const error = await response.json();
      throw new Error(error.error || 'Upload failed');
    }

    return response.json();
  }

  // Process Local File (for Tauri integration)
  async processLocalFile(filePath: string, preferences?: UserPreferences): Promise<{
    job_id: string;
    status: string;
    message: string;
  }> {
    const response = await fetch(`${this.baseUrl}/api/process`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        file_path: filePath,
        preferences: preferences || {},
      }),
    });

    if (!response.ok) {
      const error = await response.json();
      throw new Error(error.error || 'Processing failed');
    }

    return response.json();
  }

  // Get Job Status
  async getJobStatus(jobId: string): Promise<Job> {
    const response = await fetch(`${this.baseUrl}/api/status/${jobId}`);
    if (!response.ok) {
      if (response.status === 404) throw new Error('Job not found');
      throw new Error('Failed to get job status');
    }
    return response.json();
  }

  // Get Job Results
  async getJobResults(jobId: string): Promise<PipelineResult> {
    const response = await fetch(`${this.baseUrl}/api/results/${jobId}`);
    if (!response.ok) {
      const error = await response.json();
      throw new Error(error.error || 'Failed to get results');
    }
    return response.json();
  }

  // Download Results
  async downloadResults(jobId: string): Promise<Blob> {
    const response = await fetch(`${this.baseUrl}/api/results/${jobId}/download`);
    if (!response.ok) {
      throw new Error('Failed to download results');
    }
    return response.blob();
  }

  // Cancel Job
  async cancelJob(jobId: string): Promise<{ message: string; job_id: string }> {
    const response = await fetch(`${this.baseUrl}/api/cancel/${jobId}`, {
      method: 'POST',
    });
    if (!response.ok) {
      const error = await response.json();
      throw new Error(error.error || 'Failed to cancel job');
    }
    return response.json();
  }

  // List Jobs
  async listJobs(status?: string, limit?: number): Promise<{
    total: number;
    jobs: Job[];
  }> {
    const params = new URLSearchParams();
    if (status) params.append('status', status);
    if (limit) params.append('limit', limit.toString());

    const response = await fetch(`${this.baseUrl}/api/jobs?${params}`);
    if (!response.ok) throw new Error('Failed to list jobs');
    return response.json();
  }

  // WebSocket Connection Management
  async connectWebSocket(): Promise<void> {
    if (this.socket?.connected) return;

    await this.ensureInitialized();
    const wsUrl = apiConfig.getWsUrl();

    this.socket = io(wsUrl, {
      transports: ['websocket', 'polling'],
    });

    this.socket.on('connect', () => {
      console.log('Connected to GeneKnow Pipeline WebSocket');
    });

    this.socket.on('disconnect', () => {
      console.log('Disconnected from GeneKnow Pipeline WebSocket');
    });

    this.socket.on('job_progress', (data: JobProgress) => {
      const callback = this.progressCallbacks.get(data.job_id);
      if (callback) callback(data);
    });

    this.socket.on('error', (error: WebSocketError) => {
      console.error('WebSocket error:', error);
    });
  }

  disconnectWebSocket(): void {
    if (this.socket) {
      this.socket.disconnect();
      this.socket = null;
    }
  }

  // Subscribe to Job Progress
  subscribeToJobProgress(jobId: string, callback: (progress: JobProgress) => void): void {
    if (!this.socket?.connected) {
      this.connectWebSocket();
    }

    this.progressCallbacks.set(jobId, callback);
    this.socket?.emit('subscribe_job', { job_id: jobId });
  }

  // Unsubscribe from Job Progress
  unsubscribeFromJobProgress(jobId: string): void {
    this.progressCallbacks.delete(jobId);
    this.socket?.emit('unsubscribe_job', { job_id: jobId });
  }

  // Utility: Wait for Job Completion
  async waitForJobCompletion(
    jobId: string,
    onProgress?: (progress: JobProgress) => void,
    pollInterval: number = 1000
  ): Promise<Job> {
    return new Promise((resolve, reject) => {
      // Use WebSocket if available
      if (onProgress) {
        this.subscribeToJobProgress(jobId, onProgress);
      }

      // Polling fallback
      const checkStatus = async () => {
        try {
          const job = await this.getJobStatus(jobId);
          
          if (job.status === 'completed') {
            if (onProgress) this.unsubscribeFromJobProgress(jobId);
            resolve(job);
          } else if (job.status === 'failed' || job.status === 'cancelled') {
            if (onProgress) this.unsubscribeFromJobProgress(jobId);
            reject(new Error(job.error || `Job ${job.status}`));
          } else {
            setTimeout(checkStatus, pollInterval);
          }
        } catch (error) {
          if (onProgress) this.unsubscribeFromJobProgress(jobId);
          reject(error);
        }
      };

      checkStatus();
    });
  }

  // Utility: Process File and Wait for Results
  async processFileAndWaitForResults(
    file: File | string,
    preferences?: UserPreferences,
    onProgress?: (progress: JobProgress) => void
  ): Promise<PipelineResult> {
    let jobInfo;
    
    if (typeof file === 'string') {
      // Local file path (Tauri)
      jobInfo = await this.processLocalFile(file, preferences);
    } else {
      // File upload
      jobInfo = await this.uploadFile(file, preferences);
    }

    // Wait for completion
    await this.waitForJobCompletion(jobInfo.job_id, onProgress);

    // Get results
    return this.getJobResults(jobInfo.job_id);
  }
}

// Default client instance
export const geneKnowClient = new GeneKnowPipelineClient();

// React Hook for easy integration
export function useGeneKnowPipeline() {
  const client = geneKnowClient;

  // Connect WebSocket on mount
  React.useEffect(() => {
    if (typeof window !== 'undefined') {
      client.connectWebSocket().catch(console.error);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []); // client is a stable singleton, so we don't need it in deps

  return {
    client,
    checkHealth: () => client.checkHealth(),
    getPipelineInfo: () => client.getPipelineInfo(),
    getSupportedFormats: () => client.getSupportedFormats(),
    uploadFile: (file: File, preferences?: UserPreferences) => 
      client.uploadFile(file, preferences),
    processLocalFile: (path: string, preferences?: UserPreferences) => 
      client.processLocalFile(path, preferences),
    getJobStatus: (jobId: string) => client.getJobStatus(jobId),
    getJobResults: (jobId: string) => client.getJobResults(jobId),
    cancelJob: (jobId: string) => client.cancelJob(jobId),
    subscribeToProgress: (jobId: string, callback: (progress: JobProgress) => void) => 
      client.subscribeToJobProgress(jobId, callback),
    processAndWait: (file: File | string, preferences?: UserPreferences, onProgress?: (progress: JobProgress) => void) =>
      client.processFileAndWaitForResults(file, preferences, onProgress),
  };
} 