import { invoke } from '@tauri-apps/api/core';

// Types for the API
export interface FastqToVcfOptions {
  reference: string;
  fastq1: string;
  fastq2?: string;
  outputPrefix: string;
  threads?: number;
  aligner?: 'bwa' | 'minimap2' | 'bowtie2';
}

export interface FastqToVcfResult {
  success: boolean;
  vcfFile?: string;
  bamFile?: string;
  variantCount?: number;
  error?: string;
  executionTime?: number;
}

export interface ExtractRegionOptions {
  bedFile: string;
  outputDir: string;
  maxProcesses?: number;
}

export interface ExtractRegionResult {
  success: boolean;
  totalRegions: number;
  successfulExtractions: number;
  totalVariants: number;
  totalSize: number;
  executionTime: number;
  results: string[];
  error?: string;
}

// API wrapper for FASTQ to VCF conversion
export async function convertFastqToVcf(options: FastqToVcfOptions): Promise<FastqToVcfResult> {
  try {
    // Call the Tauri command that will execute the Python script
    const result = await invoke<FastqToVcfResult>('convert_fastq_to_vcf', {
      options
    });
    return result;
  } catch (error) {
    return {
      success: false,
      error: error instanceof Error ? error.message : 'Unknown error occurred'
    };
  }
}

// API wrapper for genomic region extraction
export async function extractGenomicRegions(options: ExtractRegionOptions): Promise<ExtractRegionResult> {
  try {
    // Call the Tauri command that will execute the Python script
    const result = await invoke<ExtractRegionResult>('extract_genomic_regions', {
      options
    });
    return result;
  } catch (error) {
    return {
      success: false,
      totalRegions: 0,
      successfulExtractions: 0,
      totalVariants: 0,
      totalSize: 0,
      executionTime: 0,
      results: [],
      error: error instanceof Error ? error.message : 'Unknown error occurred'
    };
  }
}

// Check if required dependencies are installed
export async function checkDependencies(): Promise<{
  samtools: boolean;
  bcftools: boolean;
  bwa: boolean;
  minimap2: boolean;
  bowtie2: boolean;
}> {
  try {
    const result = await invoke<{
      samtools: boolean;
      bcftools: boolean;
      bwa: boolean;
      minimap2: boolean;
      bowtie2: boolean;
    }>('check_genomic_dependencies');
    return result;
  } catch (error) {
    console.error('Failed to check dependencies:', error);
    return {
      samtools: false,
      bcftools: false,
      bwa: false,
      minimap2: false,
      bowtie2: false
    };
  }
}

// Get available VCF files for extraction
export async function getAvailableVcfFiles(): Promise<Record<string, string>> {
  try {
    const result = await invoke<Record<string, string>>('get_available_vcf_files');
    return result;
  } catch (error) {
    console.error('Failed to get VCF files:', error);
    return {};
  }
}

// Generate test FASTQ files
export async function generateTestFastqFiles(outputDir: string, numReads: number = 1000): Promise<{
  success: boolean;
  file1?: string;
  file2?: string;
  error?: string;
}> {
  try {
    const result = await invoke<{
      success: boolean;
      file1?: string;
      file2?: string;
      error?: string;
    }>('generate_test_fastq', {
      outputDir,
      numReads
    });
    return result;
  } catch (error) {
    return {
      success: false,
      error: error instanceof Error ? error.message : 'Unknown error occurred'
    };
  }
} 