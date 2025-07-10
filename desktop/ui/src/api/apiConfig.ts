/**
 * API Configuration for GeneKnow Pipeline
 * Handles dynamic port allocation for the local API server
 */

import { invoke } from '@tauri-apps/api/core';

// Declare Tauri on window object
declare global {
  interface Window {
    __TAURI__?: any;
  }
}

class ApiConfig {
  private static instance: ApiConfig;
  private port: number | null = null;
  private baseUrl: string | null = null;
  private wsUrl: string | null = null;

  private constructor() {}

  static getInstance(): ApiConfig {
    if (!ApiConfig.instance) {
      ApiConfig.instance = new ApiConfig();
    }
    return ApiConfig.instance;
  }

  async initialize(): Promise<void> {
    try {
      // In Tauri app, get the dynamic port
      if (window.__TAURI__) {
        console.log('Getting API port from Tauri backend...');
        this.port = await invoke<number>('get_api_port');
        console.log('API server running on port:', this.port);
        this.baseUrl = `http://localhost:${this.port}`;
        this.wsUrl = `http://localhost:${this.port}`;
      } else {
        // In development, use environment variables or defaults
        this.baseUrl = import.meta.env.VITE_API_URL || 'http://localhost:5001';
        this.wsUrl = import.meta.env.VITE_WS_URL || 'http://localhost:5001';
        // Extract port from URL
        if (this.baseUrl) {
          const url = new URL(this.baseUrl);
          this.port = parseInt(url.port) || 5001;
        } else {
          this.port = 5001;
        }
      }
    } catch (error) {
      console.error('Failed to get API port, using defaults:', error);
      // Fallback to defaults
      this.port = 5001;
      this.baseUrl = 'http://localhost:5001';
      this.wsUrl = 'http://localhost:5001';
    }
  }

  getPort(): number {
    if (this.port === null) {
      throw new Error('API configuration not initialized. Call initialize() first.');
    }
    return this.port;
  }

  getBaseUrl(): string {
    if (this.baseUrl === null) {
      throw new Error('API configuration not initialized. Call initialize() first.');
    }
    return this.baseUrl;
  }

  getWsUrl(): string {
    if (this.wsUrl === null) {
      throw new Error('API configuration not initialized. Call initialize() first.');
    }
    return this.wsUrl;
  }

  // Reset configuration (useful for testing or reconnection)
  reset(): void {
    this.port = null;
    this.baseUrl = null;
    this.wsUrl = null;
  }
}

export const apiConfig = ApiConfig.getInstance();

// Helper function to ensure API is configured
export async function ensureApiConfig(): Promise<void> {
  try {
    apiConfig.getBaseUrl(); // This will throw if not initialized
  } catch {
    await apiConfig.initialize();
  }
} 