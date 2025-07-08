#!/usr/bin/env python3
"""
FASTQ to VCF Pipeline
Converts raw sequencing reads (FASTQ) to variant calls (VCF)
"""

import os
import subprocess
import argparse
import sys
import json
from pathlib import Path
import time

class FastqToVcfPipeline:
    def __init__(self, reference_genome, threads=4, json_output=False):
        self.reference = reference_genome
        self.threads = threads
        self.json_output = json_output
        self.check_dependencies()
    
    def log(self, message):
        """Log a message (only if not in JSON mode)"""
        if not self.json_output:
            print(message)
    
    def check_dependencies(self):
        """Check if required tools are installed"""
        required_tools = {
            'samtools': 'SAM/BAM file manipulation',
            'bcftools': 'Variant calling and VCF manipulation',
        }
        
        optional_tools = {
            'bwa': 'BWA aligner (most common)',
            'minimap2': 'Minimap2 aligner (faster, good for long reads)',
            'bowtie2': 'Bowtie2 aligner (alternative)'
        }
        
        if not self.json_output:
            self.log("🔍 Checking dependencies...")
        
        # Check required tools
        missing_required = []
        for tool, desc in required_tools.items():
            if not self._tool_exists(tool):
                missing_required.append(f"{tool} - {desc}")
        
        if missing_required:
            error_msg = f"Missing required tools: {', '.join(missing_required)}"
            if self.json_output:
                result = {
                    "success": False,
                    "error": error_msg,
                    "missing_tools": missing_required
                }
                print(json.dumps(result))
                sys.exit(1)
            else:
                self.log(f"\n❌ {error_msg}")
                self.log("Install with: brew install samtools bcftools")
                sys.exit(1)
        
        # Check aligners
        self.available_aligners = []
        if not self.json_output:
            self.log("\n📊 Available aligners:")
        
        for tool, desc in optional_tools.items():
            if self._tool_exists(tool):
                self.available_aligners.append(tool)
                if not self.json_output:
                    self.log(f"   ✅ {tool} - {desc}")
            else:
                if not self.json_output:
                    self.log(f"   ❌ {tool} - {desc}")
        
        if not self.available_aligners:
            error_msg = "No aligners found! Install one with: brew install bwa"
            if self.json_output:
                result = {
                    "success": False,
                    "error": error_msg,
                    "available_aligners": []
                }
                print(json.dumps(result))
                sys.exit(1)
            else:
                self.log(f"\n⚠️  {error_msg}")
                sys.exit(1)
    
    def _tool_exists(self, tool):
        """Check if a tool is in PATH"""
        try:
            subprocess.run(['which', tool], capture_output=True, check=True)
            return True
        except:
            return False
    
    def index_reference(self, aligner='bwa'):
        """Index reference genome for alignment"""
        self.log(f"\n🧬 Indexing reference genome with {aligner}...")
        
        if aligner == 'bwa':
            cmd = ['bwa', 'index', self.reference]
        elif aligner == 'minimap2':
            # Minimap2 doesn't need pre-indexing for small genomes
            self.log("   Minimap2 will index on-the-fly")
            return
        elif aligner == 'bowtie2':
            base = self.reference.replace('.fa', '').replace('.fasta', '')
            cmd = ['bowtie2-build', self.reference, base]
        
        subprocess.run(cmd, check=True)
    
    def align_reads(self, fastq1, fastq2, output_bam, aligner='bwa'):
        """Align FASTQ reads to reference genome"""
        self.log(f"\n🧬 Aligning reads with {aligner}...")
        
        sam_file = output_bam.replace('.bam', '.sam')
        
        if aligner == 'bwa':
            # BWA mem algorithm (recommended for Illumina reads > 70bp)
            cmd = ['bwa', 'mem', '-t', str(self.threads), self.reference, fastq1]
            if fastq2:  # Paired-end
                cmd.append(fastq2)
        
        elif aligner == 'minimap2':
            # Minimap2 (good for long reads or quick alignment)
            cmd = ['minimap2', '-ax', 'sr', '-t', str(self.threads), self.reference, fastq1]
            if fastq2:
                cmd.append(fastq2)
        
        elif aligner == 'bowtie2':
            # Bowtie2
            base = self.reference.replace('.fa', '').replace('.fasta', '')
            cmd = ['bowtie2', '-x', base, '-p', str(self.threads), '-1', fastq1]
            if fastq2:
                cmd.extend(['-2', fastq2])
            else:
                cmd = ['bowtie2', '-x', base, '-p', str(self.threads), '-U', fastq1]
        
        # Run alignment and pipe to SAM file
        self.log(f"   Running: {' '.join(cmd[:4])}...")
        with open(sam_file, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        
        # Convert SAM to sorted BAM
        self.log("   Converting to sorted BAM...")
        subprocess.run([
            'samtools', 'sort', '-@', str(self.threads), 
            '-o', output_bam, sam_file
        ], check=True)
        
        # Index BAM
        subprocess.run(['samtools', 'index', output_bam], check=True)
        
        # Clean up SAM file
        os.remove(sam_file)
        
        return output_bam
    
    def call_variants(self, bam_file, output_vcf, method='bcftools'):
        """Call variants from aligned BAM file"""
        self.log(f"\n🧬 Calling variants with {method}...")
        
        if method == 'bcftools':
            # BCFtools mpileup + call pipeline
            bcf_file = output_vcf.replace('.vcf', '.bcf')
            
            # Generate genotype likelihoods
            cmd1 = [
                'bcftools', 'mpileup', '-Ou',
                '-f', self.reference,
                '--threads', str(self.threads),
                bam_file
            ]
            
            # Call variants
            cmd2 = [
                'bcftools', 'call', '-mv',
                '-Ob', '-o', bcf_file
            ]
            
            self.log("   Generating pileup and calling variants...")
            p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
            p1.stdout.close()
            p2.communicate()
            
            # Convert BCF to VCF
            self.log("   Converting to VCF...")
            subprocess.run([
                'bcftools', 'view', bcf_file, 
                '-o', output_vcf
            ], check=True)
            
            # Clean up
            os.remove(bcf_file)
        
        # Could add GATK, FreeBayes, etc. here
        
        return output_vcf
    
    def run_pipeline(self, fastq1, fastq2=None, output_prefix='output', aligner=None):
        """Run complete FASTQ to VCF pipeline"""
        if not self.json_output:
            self.log("\n🚀 Starting FASTQ to VCF pipeline...")
        
        start_time = time.time()
        
        # Select aligner
        if not aligner:
            if self.available_aligners:
                aligner = self.available_aligners[0]
                self.log(f"   Using {aligner} aligner")
            else:
                error_msg = "No aligner available!"
                if self.json_output:
                    result = {
                        "success": False,
                        "error": error_msg,
                        "available_aligners": self.available_aligners
                    }
                    print(json.dumps(result))
                    return None
                else:
                    self.log(f"❌ {error_msg}")
                    return None
        
        # Define output files
        bam_file = f"{output_prefix}.bam"
        vcf_file = f"{output_prefix}.vcf"
        
        # Run pipeline steps
        try:
            # 1. Index reference (if needed)
            if aligner in ['bwa', 'bowtie2']:
                self.index_reference(aligner)
            
            # 2. Align reads
            self.align_reads(fastq1, fastq2, bam_file, aligner)
            
            # 3. Call variants
            self.call_variants(bam_file, vcf_file)
            
            # 4. Generate stats
            if not self.json_output:
                self.log("\n📊 Pipeline complete! Generating statistics...")
                subprocess.run(['samtools', 'flagstat', bam_file])
            
            # Count variants
            var_count = subprocess.run(
                ['bcftools', 'view', '-H', vcf_file], 
                capture_output=True, text=True
            )
            num_variants = len([line for line in var_count.stdout.strip().split('\n') if line.strip()])
            
            elapsed = time.time() - start_time
            
            if self.json_output:
                # Return structured JSON result
                result = {
                    "success": True,
                    "vcf_file": vcf_file,
                    "bam_file": bam_file,
                    "variant_count": num_variants,
                    "execution_time": elapsed,
                    "aligner_used": aligner,
                    "error": None
                }
                print(json.dumps(result))
            else:
                self.log(f"\n✅ Success!")
                self.log(f"   Time: {elapsed:.1f} seconds")
                self.log(f"   Output BAM: {bam_file}")
                self.log(f"   Output VCF: {vcf_file}")
                self.log(f"   Variants found: {num_variants}")
            
            return vcf_file
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Pipeline failed: {e}"
            if self.json_output:
                result = {
                    "success": False,
                    "error": error_msg,
                    "vcf_file": None,
                    "bam_file": None,
                    "variant_count": None,
                    "execution_time": time.time() - start_time
                }
                print(json.dumps(result))
            else:
                self.log(f"\n❌ {error_msg}")
            return None
        except Exception as e:
            error_msg = f"Unexpected error: {e}"
            if self.json_output:
                result = {
                    "success": False,
                    "error": error_msg,
                    "vcf_file": None,
                    "bam_file": None,
                    "variant_count": None,
                    "execution_time": time.time() - start_time
                }
                print(json.dumps(result))
            else:
                self.log(f"\n❌ {error_msg}")
            return None


def download_test_data():
    """Download small test FASTQ files"""
    print("\n📥 Test Data Sources:\n")
    
    print("1. Generate Synthetic Test Data (Simplest & Most Reliable):")
    print("   # Create synthetic FASTQ files")
    print("   cat > generate_test_fastq.py << 'EOF'")
    print("import random")
    print("import gzip")
    print("")
    print("def generate_fastq(filename, num_reads=1000, read_length=150):")
    print("    bases = ['A', 'T', 'G', 'C']")
    print("    with gzip.open(filename, 'wt') as f:")
    print("        for i in range(num_reads):")
    print("            seq = ''.join(random.choice(bases) for _ in range(read_length))")
    print("            qual = 'I' * read_length")
    print("            f.write(f'@read_{i}\\n{seq}\\n+\\n{qual}\\n')")
    print("")
    print("generate_fastq('test_R1.fastq.gz')")
    print("generate_fastq('test_R2.fastq.gz')")
    print("print('Generated test FASTQ files!')")
    print("EOF")
    print("")
    print("   python3 generate_test_fastq.py")
    print()
    
    print("2. Download from SRA (Most Reliable for Real Data):")
    print("   # Install SRA toolkit if needed:")
    print("     brew install sra-tools")
    print()
    print("   # Download E. coli dataset (limited to 100k reads):")
    print("     fastq-dump --split-files --gzip SRR390728 --maxSpotId 100000")
    print()
    
    print("3. Reference Genomes (Verified Working):")
    print("   # E. coli reference (small, ~4.6MB):")
    print("     curl -o ecoli_ref.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz")
    print("     gunzip ecoli_ref.fa.gz")
    print("     mv GCF_000005845.2_ASM584v2_genomic.fna ecoli_ref.fa")
    print()
    print("   # Human chr20 reference (GRCh37/hg19):")
    print("     wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz")
    print("     gunzip chr20.fa.gz")
    print()
    
    print("4. Quick Test Example:")
    print("   # Generate test data and run pipeline:")
    print("     python3 generate_test_fastq.py")
    print("     python3 fastq_to_vcf_pipeline.py -r ecoli_ref.fa -1 test_R1.fastq.gz -2 test_R2.fastq.gz -o test_output")


def main():
    parser = argparse.ArgumentParser(description='Convert FASTQ files to VCF')
    parser.add_argument('--download-test-data', action='store_true',
                      help='Show where to download test FASTQ files')
    parser.add_argument('-r', '--reference', help='Reference genome (FASTA)')
    parser.add_argument('-1', '--fastq1', help='FASTQ file (or R1 for paired-end)')
    parser.add_argument('-2', '--fastq2', help='R2 FASTQ file (for paired-end)')
    parser.add_argument('-o', '--output', default='output', help='Output prefix')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('-a', '--aligner', choices=['bwa', 'minimap2', 'bowtie2'],
                      help='Aligner to use (auto-detect if not specified)')
    parser.add_argument('--json', action='store_true',
                      help='Output results in JSON format')
    
    args = parser.parse_args()
    
    if args.download_test_data:
        download_test_data()
        return
    
    if not args.reference or not args.fastq1:
        if args.json:
            result = {
                "success": False,
                "error": "Reference genome and at least one FASTQ file required!",
                "vcf_file": None,
                "bam_file": None,
                "variant_count": None
            }
            print(json.dumps(result))
        else:
            print("❌ Error: Reference genome and at least one FASTQ file required!")
            print("\nUsage examples:")
            print("  # Show test data sources:")
            print("  python3 fastq_to_vcf_pipeline.py --download-test-data")
            print()
            print("  # Run pipeline:")
            print("  python3 fastq_to_vcf_pipeline.py -r reference.fa -1 reads.fastq -o output")
            print("  python3 fastq_to_vcf_pipeline.py -r reference.fa -1 R1.fastq -2 R2.fastq -o output")
        return
    
    # Run pipeline
    pipeline = FastqToVcfPipeline(args.reference, args.threads, args.json)
    pipeline.run_pipeline(args.fastq1, args.fastq2, args.output, args.aligner)


if __name__ == "__main__":
    main() 