"""
Test script to verify LangGraph pipeline parallelization and node functionality.
"""
import time
import json
from datetime import datetime
from graph import run_pipeline, create_genomic_pipeline
import logging
import sys

# Configure logging to capture execution order
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s.%(msecs)03d - %(name)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

logger = logging.getLogger(__name__)

def test_fastq_path():
    """Test the FASTQ → Variant Calling path"""
    print("\n" + "="*60)
    print("TEST 1: FASTQ Input (Variant Calling Path)")
    print("="*60)
    
    start_time = time.time()
    result = run_pipeline(
        "../test_R1.fastq.gz",
        {"patient_data": {"age": 45, "sex": "F"}}
    )
    end_time = time.time()
    
    print(f"\n✅ Pipeline Status: {result['pipeline_status']}")
    print(f"⏱️  Total Time: {end_time - start_time:.3f}s")
    print(f"📋 Completed Nodes: {result['completed_nodes']}")
    print(f"🧬 Variants Found: {len(result.get('raw_variants', []))}")
    print(f"✓ Filtered Variants: {result.get('variant_count', 0)}")
    
    # Check execution path
    if "variant_calling" in result['completed_nodes']:
        print("✅ Variant calling executed (correct path for FASTQ)")
    else:
        print("❌ Variant calling NOT executed (incorrect!)")
    
    return result


def test_vcf_path():
    """Test the VCF → QC Filter path (skips variant calling)"""
    print("\n" + "="*60)
    print("TEST 2: VCF Input (Direct QC Filter Path)")
    print("="*60)
    
    # First, create a test VCF file with proper formatting
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr17	41223094	.	A	G	99	PASS	DP=50;AF=0.5	GT:DP:AF	0/1:50:0.5
chr17	7577121	.	G	A	85	PASS	DP=40;AF=0.45	GT:DP:AF	0/1:40:0.45
chr5	112173917	.	C	T	35	PASS	DP=30;AF=0.4	GT:DP:AF	0/1:30:0.4
"""
    
    with open("test_variants.vcf", "w") as f:
        f.write(vcf_content)
    
    start_time = time.time()
    result = run_pipeline(
        "test_variants.vcf",
        {"patient_data": {"age": 45, "sex": "F"}}
    )
    end_time = time.time()
    
    print(f"\n✅ Pipeline Status: {result['pipeline_status']}")
    print(f"⏱️  Total Time: {end_time - start_time:.3f}s")
    print(f"📋 Completed Nodes: {result['completed_nodes']}")
    print(f"🧬 Variants Loaded from VCF: {len(result.get('raw_variants', []))}")
    print(f"✓ Filtered Variants: {result.get('variant_count', 0)}")
    
    # Check execution path
    if "variant_calling" in result['completed_nodes']:
        print("❌ Variant calling executed (should skip for VCF!)")
    else:
        print("✅ Variant calling skipped (correct path for VCF)")
    
    # Clean up
    import os
    os.remove("test_variants.vcf")
    
    return result


def test_parallel_execution_timing():
    """Test to show parallel execution by adding delays to nodes"""
    print("\n" + "="*60)
    print("TEST 3: Parallel Execution Timing")
    print("="*60)
    
    # Monkey patch nodes to add delays and logging
    from nodes import variant_calling, qc_filter
    
    original_vc_process = variant_calling.process
    original_qc_process = qc_filter.process
    
    def delayed_vc_process(state):
        logger.info("🧬 VARIANT CALLING: Starting (2s delay)")
        time.sleep(2)  # Simulate slow variant calling
        result = original_vc_process(state)
        logger.info("🧬 VARIANT CALLING: Completed")
        return result
    
    def delayed_qc_process(state):
        logger.info("🔍 QC FILTER: Starting (1s delay)")
        time.sleep(1)  # Simulate QC processing
        result = original_qc_process(state)
        logger.info("🔍 QC FILTER: Completed")
        return result
    
    # Apply patches
    variant_calling.process = delayed_vc_process
    qc_filter.process = delayed_qc_process
    
    print("\nIf parallel execution works, QC Filter should start before Variant Calling completes.")
    print("Watch the timestamps...\n")
    
    start_time = time.time()
    
    # This would show parallel execution if it was truly parallel
    # Note: LangGraph's conditional_edges doesn't actually run nodes in parallel
    # by default - it's more about graph structure than concurrent execution
    
    result = run_pipeline(
        "../test_R1.fastq.gz",
        {"patient_data": {"age": 45, "sex": "F"}}
    )
    
    end_time = time.time()
    
    print(f"\n⏱️  Total Time: {end_time - start_time:.3f}s")
    print("Note: True parallel execution would complete faster than sequential (3s total)")
    
    # Restore original functions
    variant_calling.process = original_vc_process
    qc_filter.process = original_qc_process
    
    return result


def test_node_functionality():
    """Test individual node outputs"""
    print("\n" + "="*60)
    print("TEST 4: Node Functionality Verification")
    print("="*60)
    
    result = run_pipeline(
        "../test_R1.fastq.gz",
        {"patient_data": {"age": 45, "sex": "F"}}
    )
    
    # Check each node's output
    print("\n📊 Node Output Verification:")
    
    # 1. File Input
    if result.get('file_metadata'):
        print("✅ file_input: Metadata extracted")
        print(f"   - Format: {result['file_metadata'].get('format')}")
        print(f"   - Read count: {result['file_metadata'].get('estimated_read_count')}")
    
    # 2. Preprocess
    if result.get('aligned_bam_path') or result.get('raw_variants'):
        print("✅ preprocess: Data processed")
        if result.get('aligned_bam_path'):
            print(f"   - BAM created: {result['aligned_bam_path']}")
    
    # 3. Variant Calling
    if result.get('raw_variants') is not None:
        print(f"✅ variant_calling: {len(result['raw_variants'])} variants found")
    
    # 4. QC Filter
    if result.get('filtered_variants') is not None:
        print(f"✅ qc_filter: {len(result['filtered_variants'])} variants passed QC")
        if result['file_metadata'].get('qc_stats'):
            stats = result['file_metadata']['qc_stats']
            print(f"   - Filter rate: {stats['filter_rate_percent']:.1f}%")
    
    # 5. TCGA Mapper
    if result.get('tcga_matches') is not None:
        print(f"✅ tcga_mapper: Matched against {len(result['tcga_matches'])} cancer types")
    
    # 6. Risk Model
    if result.get('risk_scores'):
        print("✅ risk_model: Risk scores calculated")
        for cancer, score in result['risk_scores'].items():
            print(f"   - {cancer}: {score:.1f}%")
    
    # 7. Report Writer
    if result.get('report_sections'):
        print("✅ report_writer: Report sections generated")
        print(f"   - Sections: {list(result['report_sections'].keys())}")
    
    return result


def test_error_handling():
    """Test error handling in parallel paths"""
    print("\n" + "="*60)
    print("TEST 5: Error Handling")
    print("="*60)
    
    # Test with non-existent file
    result = run_pipeline(
        "non_existent_file.fastq",
        {"patient_data": {"age": 45, "sex": "F"}}
    )
    
    print(f"Pipeline Status: {result['pipeline_status']}")
    print(f"Errors: {len(result.get('errors', []))}")
    for error in result.get('errors', []):
        print(f"  - {error['node']}: {error['error']}")
    
    return result


def visualize_pipeline():
    """Show the pipeline structure"""
    print("\n" + "="*60)
    print("PIPELINE STRUCTURE")
    print("="*60)
    
    pipeline = create_genomic_pipeline()
    
    print("\n🔷 Nodes:")
    for i, node in enumerate(pipeline.nodes):
        if node not in ['__start__', '__end__']:
            print(f"  {i}. {node}")
    
    print("\n🔗 Execution Flow:")
    print("  file_input")
    print("      ↓")
    print("  preprocess")
    print("      ↓")
    print("  [Conditional Routing]")
    print("    ├─→ variant_calling (if FASTQ/BAM)")
    print("    └─→ qc_filter (if VCF with variants)")
    print("      ↓")
    print("  merge_parallel")
    print("      ↓")
    print("  tcga_mapper → risk_model → formatter → report_writer")


if __name__ == "__main__":
    print("🧬 LangGraph Pipeline Parallelization Tests")
    print("==========================================")
    
    # Show pipeline structure
    visualize_pipeline()
    
    # Run all tests
    test_fastq_path()
    test_vcf_path()
    test_parallel_execution_timing()
    test_node_functionality()
    test_error_handling()
    
    print("\n✅ All tests completed!")
    print("\nNOTE: LangGraph's conditional routing creates a parallel *structure*")
    print("but doesn't necessarily execute nodes concurrently by default.")
    print("For true parallel execution, you'd need to use async nodes or")
    print("configure the graph with a parallel executor.") 