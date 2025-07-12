"""
Create a minimal test reference genome for development.
This is NOT a real reference - just for testing the pipeline.
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def create_test_reference():
    """Create a minimal test reference with a few chromosomes."""

    # Create references directory
    ref_dir = "test_reference"
    os.makedirs(ref_dir, exist_ok=True)

    # Create a few small "chromosomes" with random sequences
    chromosomes = []

    # Add a few genes we're testing (simplified)
    # Chr17 for BRCA1 and TP53
    chr17_seq = "ATCGATCGATCG" * 1000  # 12kb mini chromosome
    chromosomes.append(
        SeqRecord(Seq(chr17_seq), id="chr17", description="Test chromosome 17")
    )

    # Chr5 for APC
    chr5_seq = "GCTAGCTAGCTA" * 1000  # 12kb mini chromosome
    chromosomes.append(
        SeqRecord(Seq(chr5_seq), id="chr5", description="Test chromosome 5")
    )

    # Write FASTA file
    ref_file = os.path.join(ref_dir, "test_genome.fa")
    SeqIO.write(chromosomes, ref_file, "fasta")

    print(f"✅ Created test reference: {ref_file}")
    print("   Note: This is a minimal test reference, NOT suitable for real analysis!")

    # Create BWA index
    print("\n🔨 Creating BWA index...")
    os.system(f"bwa index {ref_file} 2>/dev/null")

    # Create samtools index
    print("🔨 Creating samtools index...")
    os.system(f"samtools faidx {ref_file} 2>/dev/null")

    print("\n✅ Test reference ready for pipeline testing!")
    return ref_file


if __name__ == "__main__":
    create_test_reference()
