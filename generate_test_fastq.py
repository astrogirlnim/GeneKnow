import random
import gzip

def generate_fastq(filename, num_reads=1000, read_length=150):
    """Generate synthetic FASTQ file for testing"""
    bases = ['A', 'T', 'G', 'C']
    with gzip.open(filename, 'wt') as f:
        for i in range(num_reads):
            # Generate random sequence
            seq = ''.join(random.choice(bases) for _ in range(read_length))
            # Simple quality scores (all high quality 'I')
            qual = 'I' * read_length
            # Write FASTQ format: @id, sequence, +, quality
            f.write(f'@read_{i}\n{seq}\n+\n{qual}\n')

# Generate paired-end reads
generate_fastq('test_R1.fastq.gz', num_reads=500)
generate_fastq('test_R2.fastq.gz', num_reads=500)
print('✅ Generated test FASTQ files!')
print('   test_R1.fastq.gz - 500 reads')
print('   test_R2.fastq.gz - 500 reads') 