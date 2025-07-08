"""
Create a test VCF file with cancer-related variants.
This simulates what DeepVariant would output.
"""
import os
from datetime import datetime

def create_test_vcf():
    """Create a test VCF file with known cancer variants."""
    
    vcf_dir = "test_data"
    os.makedirs(vcf_dir, exist_ok=True)
    
    vcf_content = """##fileformat=VCFv4.2
##fileDate={date}
##source=DeepVariant-1.5.0
##reference=test_genome.fa
##contig=<ID=chr17,length=12000>
##contig=<ID=chr5,length=12000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr17\t41223094\t.\tA\tG\t99.0\tPASS\tDP=45;AF=0.48\tGT:GQ:DP\t0/1:99:45
chr17\t41244936\t.\tG\tA\t87.5\tPASS\tDP=38;AF=0.42\tGT:GQ:DP\t0/1:87:38
chr17\t41256885\t.\tC\tT\t78.2\tPASS\tDP=32;AF=0.38\tGT:GQ:DP\t0/1:78:32
chr17\t7577121\t.\tG\tA\t92.1\tPASS\tDP=40;AF=0.45\tGT:GQ:DP\t0/1:92:40
chr17\t7578190\t.\tC\tT\t85.3\tPASS\tDP=35;AF=0.40\tGT:GQ:DP\t0/1:85:35
chr5\t112173917\t.\tC\tT\t65.2\tPASS\tDP=28;AF=0.36\tGT:GQ:DP\t0/1:65:28
chr5\t112175639\t.\tA\tG\t72.8\tPASS\tDP=30;AF=0.40\tGT:GQ:DP\t0/1:72:30
""".format(date=datetime.now().strftime("%Y%m%d"))
    
    # Write VCF file
    vcf_path = os.path.join(vcf_dir, "test_variants.vcf")
    with open(vcf_path, 'w') as f:
        f.write(vcf_content.strip())
    
    # Also create a variant annotation file (simplified)
    annotation_content = """CHROM\tPOS\tREF\tALT\tGENE\tCONSEQUENCE\tHGVS_C\tHGVS_P\tIMPACT
chr17\t41223094\tA\tG\tBRCA1\tmissense_variant\tc.5266A>G\tp.Ile1756Val\tMODERATE
chr17\t41244936\tG\tA\tBRCA1\tmissense_variant\tc.4689G>A\tp.Met1563Ile\tMODERATE
chr17\t41256885\tC\tT\tBRCA1\tstop_gained\tc.4327C>T\tp.Arg1443*\tHIGH
chr17\t7577121\tG\tA\tTP53\tmissense_variant\tc.743G>A\tp.Arg248Gln\tMODERATE
chr17\t7578190\tC\tT\tTP53\tmissense_variant\tc.488C>T\tp.Ala163Val\tMODERATE
chr5\t112173917\tC\tT\tAPC\tstop_gained\tc.4348C>T\tp.Arg1450*\tHIGH
chr5\t112175639\tA\tG\tAPC\tmissense_variant\tc.4479A>G\tp.Ile1493Met\tMODERATE
"""
    
    annotation_path = os.path.join(vcf_dir, "test_annotations.tsv")
    with open(annotation_path, 'w') as f:
        f.write(annotation_content.strip())
    
    print(f"✅ Created test VCF file: {vcf_path}")
    print(f"✅ Created annotation file: {annotation_path}")
    print(f"   Contains {len(vcf_content.strip().split('\\n')) - 11} variants")
    print("\nVariants include:")
    print("  - BRCA1: 3 variants (2 missense, 1 stop-gained)")
    print("  - TP53: 2 variants (missense)")
    print("  - APC: 2 variants (1 stop-gained, 1 missense)")
    
    return vcf_path, annotation_path

if __name__ == "__main__":
    create_test_vcf() 