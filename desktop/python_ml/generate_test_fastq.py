#!/usr/bin/env python3
"""
Generate Test FASTQ Files
Creates synthetic FASTQ files for testing genomic pipelines
"""

import random
import gzip
import os
import json
import argparse
import sys
from pathlib import Path

def generate_fastq(output_dir, filename, num_reads=1000, read_length=150):
    """Generate a synthetic FASTQ file with random sequences"""
    bases = ['A', 'T', 'G', 'C']
    filepath = os.path.join(output_dir, filename)
    
    with gzip.open(filepath, 'wt') as f:
        for i in range(num_reads):
            # Generate random sequence
            seq = ''.join(random.choice(bases) for _ in range(read_length))
            # Use high quality scores (I = 40 in Phred+33)
            qual = 'I' * read_length
            
            # Write FASTQ record
            f.write(f'@read_{i}\n{seq}\n+\n{qual}\n')
    
    return filepath

def main():
    parser = argparse.ArgumentParser(description='Generate synthetic FASTQ files for testing')
    parser.add_argument('--output-dir', required=True, help='Output directory for FASTQ files')
    parser.add_argument('--num-reads', type=int, default=1000, help='Number of reads to generate')
    parser.add_argument('--read-length', type=int, default=150, help='Length of each read')
    parser.add_argument('--json', action='store_true', help='Output results in JSON format')
    
    args = parser.parse_args()
    
    try:
        # Create output directory if it doesn't exist
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Generate paired-end FASTQ files
        file1 = generate_fastq(args.output_dir, 'test_R1.fastq.gz', args.num_reads, args.read_length)
        file2 = generate_fastq(args.output_dir, 'test_R2.fastq.gz', args.num_reads, args.read_length)
        
        if args.json:
            result = {
                'success': True,
                'file1': file1,
                'file2': file2,
                'num_reads': args.num_reads,
                'read_length': args.read_length,
                'output_dir': args.output_dir,
                'error': None
            }
            print(json.dumps(result))
        else:
            print(f"✅ Successfully generated test FASTQ files!")
            print(f"   Output directory: {args.output_dir}")
            print(f"   File 1: {file1}")
            print(f"   File 2: {file2}")
            print(f"   Number of reads: {args.num_reads}")
            print(f"   Read length: {args.read_length}")
    
    except Exception as e:
        error_msg = f"Failed to generate test FASTQ files: {e}"
        if args.json:
            result = {
                'success': False,
                'error': error_msg,
                'file1': None,
                'file2': None,
                'num_reads': None,
                'read_length': None,
                'output_dir': args.output_dir
            }
            print(json.dumps(result))
        else:
            print(f"❌ {error_msg}")
        sys.exit(1)

if __name__ == "__main__":
    main() 