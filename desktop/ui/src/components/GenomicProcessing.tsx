import React, { useState } from 'react';
import { 
  convertFastqToVcf, 
  extractGenomicRegions, 
  checkDependencies,
  generateTestFastqFiles
} from '../api/genomicProcessing';
import type {
  FastqToVcfOptions,
  ExtractRegionOptions 
} from '../api/genomicProcessing';

export const GenomicProcessing: React.FC = () => {
  const [isProcessing, setIsProcessing] = useState(false);
  const [results, setResults] = useState<string>('');
  const [dependencies, setDependencies] = useState<Record<string, boolean>>({});

  // Check dependencies on component mount
  React.useEffect(() => {
    checkDependencies().then(setDependencies);
  }, []);

  const handleFastqToVcf = async () => {
    setIsProcessing(true);
    setResults('Processing FASTQ to VCF conversion...');
    
    try {
      // First generate test files
      const testFiles = await generateTestFastqFiles('/tmp/test_fastq', 1000);
      
      if (!testFiles.success) {
        throw new Error(testFiles.error || 'Failed to generate test files');
      }

      const options: FastqToVcfOptions = {
        reference: 'path/to/reference.fa',
        fastq1: testFiles.file1!,
        fastq2: testFiles.file2!,
        outputPrefix: '/tmp/test_output',
        threads: 4,
        aligner: 'bwa'
      };

      const result = await convertFastqToVcf(options);
      
      if (result.success) {
        setResults(`
          ✅ Conversion successful!
          VCF File: ${result.vcfFile}
          BAM File: ${result.bamFile}
          Variants Found: ${result.variantCount || 0}
          Execution Time: ${result.executionTime?.toFixed(2)}s
        `);
      } else {
        setResults(`❌ Conversion failed: ${result.error}`);
      }
    } catch (error) {
      setResults(`❌ Error: ${error}`);
    } finally {
      setIsProcessing(false);
    }
  };

  const handleExtractRegions = async () => {
    setIsProcessing(true);
    setResults('Extracting genomic regions...');
    
    try {
      const options: ExtractRegionOptions = {
        bedFile: 'path/to/regions.bed',
        outputDir: '/tmp/extracted_regions',
        maxProcesses: 4
      };

      const result = await extractGenomicRegions(options);
      
      if (result.success) {
        setResults(`
          ✅ Extraction successful!
          Total Regions: ${result.totalRegions}
          Successful Extractions: ${result.successfulExtractions}
          Total Variants: ${result.totalVariants}
          Total Size: ${(result.totalSize / 1024 / 1024).toFixed(2)} MB
          Execution Time: ${result.executionTime.toFixed(2)}s
          
          Results:
          ${result.results.join('\n')}
        `);
      } else {
        setResults(`❌ Extraction failed: ${result.error}`);
      }
    } catch (error) {
      setResults(`❌ Error: ${error}`);
    } finally {
      setIsProcessing(false);
    }
  };

  return (
    <div className="p-6 max-w-4xl mx-auto">
      <h1 className="text-2xl font-bold mb-6">Genomic Processing Tools</h1>
      
      {/* Dependencies Status */}
      <div className="mb-6 p-4 bg-gray-100 rounded">
        <h2 className="text-lg font-semibold mb-2">Dependencies Status</h2>
        <div className="grid grid-cols-2 gap-2">
          {Object.entries(dependencies).map(([tool, installed]) => (
            <div key={tool} className="flex items-center">
              <span className={installed ? 'text-green-600' : 'text-red-600'}>
                {installed ? '✅' : '❌'}
              </span>
              <span className="ml-2">{tool}</span>
            </div>
          ))}
        </div>
      </div>

      {/* Action Buttons */}
      <div className="space-y-4">
        <button
          onClick={handleFastqToVcf}
          disabled={isProcessing}
          className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600 disabled:bg-gray-400"
        >
          {isProcessing ? 'Processing...' : 'Convert FASTQ to VCF'}
        </button>

        <button
          onClick={handleExtractRegions}
          disabled={isProcessing}
          className="px-4 py-2 bg-green-500 text-white rounded hover:bg-green-600 disabled:bg-gray-400 ml-4"
        >
          {isProcessing ? 'Processing...' : 'Extract Genomic Regions'}
        </button>
      </div>

      {/* Results Display */}
      {results && (
        <div className="mt-6 p-4 bg-gray-50 rounded">
          <h2 className="text-lg font-semibold mb-2">Results</h2>
          <pre className="whitespace-pre-wrap">{results}</pre>
        </div>
      )}
    </div>
  );
}; 