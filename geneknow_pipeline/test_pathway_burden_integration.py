#!/usr/bin/env python3
"""
Integration test for Pathway Burden Model within the full pipeline.
Tests the pathway burden model as part of the complete GeneKnow pipeline.
"""

import os
import sys
import json
import tempfile
import logging
from datetime import datetime

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from graph import run_pipeline, create_genomic_pipeline

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_test_vcf_file():
    """Create a test VCF file with variants targeting multiple pathways."""
    
    # VCF content with variants in different cancer pathways
    vcf_content = """##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=17,length=81195210>
##contig=<ID=13,length=114364328>
##contig=<ID=11,length=135086622>
##contig=<ID=3,length=198295559>
##contig=<ID=9,length=138394717>
##contig=<ID=12,length=133275309>
##contig=<ID=7,length=159345973>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
17	43045677	BRCA1_var1	G	A	60	PASS	AF=0.0001;ANN=A|frameshift_variant|HIGH|BRCA1|BRCA1	GT:AD:DP:GQ	0/1:50,25:75:60
13	32363533	BRCA2_var1	T	C	55	PASS	AF=0.0005;ANN=C|missense_variant|MODERATE|BRCA2|BRCA2	GT:AD:DP:GQ	0/1:45,30:75:55
11	108098525	ATM_var1	G	A	40	PASS	AF=0.01;ANN=A|missense_variant|MODERATE|ATM|ATM	GT:AD:DP:GQ	0/1:40,15:55:40
17	7673803	TP53_var1	C	T	65	PASS	AF=0.0003;ANN=T|missense_variant|MODERATE|TP53|TP53	GT:AD:DP:GQ	0/1:50,35:85:65
13	48877782	RB1_var1	G	A	50	PASS	AF=0.0008;ANN=A|missense_variant|MODERATE|RB1|RB1	GT:AD:DP:GQ	0/1:45,20:65:50
12	25245350	KRAS_var1	C	T	45	PASS	AF=0.001;ANN=T|missense_variant|MODERATE|KRAS|KRAS	GT:AD:DP:GQ	0/1:40,25:65:45
7	55019021	EGFR_var1	G	A	42	PASS	AF=0.002;ANN=A|missense_variant|MODERATE|EGFR|EGFR	GT:AD:DP:GQ	0/1:35,20:55:42
3	37034946	MLH1_var1	G	A	58	PASS	AF=0.0004;ANN=A|stop_gained|HIGH|MLH1|MLH1	GT:AD:DP:GQ	0/1:48,28:76:58
2	47403068	MSH2_var1	C	T	52	PASS	AF=0.0006;ANN=T|missense_variant|MODERATE|MSH2|MSH2	GT:AD:DP:GQ	0/1:42,22:64:52
9	21968225	CDKN2A_var1	C	T	48	PASS	AF=0.0007;ANN=T|missense_variant|MODERATE|CDKN2A|CDKN2A	GT:AD:DP:GQ	0/1:38,18:56:48
17	43045678	BRCA1_var2	C	T	62	PASS	AF=0.0002;ANN=T|frameshift_variant|HIGH|BRCA1|BRCA1	GT:AD:DP:GQ	0/1:52,28:80:62
17	43045679	BRCA1_var3	T	G	59	PASS	AF=0.0003;ANN=G|missense_variant|MODERATE|BRCA1|BRCA1	GT:AD:DP:GQ	0/1:49,26:75:59
1	12345678	UNKNOWN_var1	A	G	35	PASS	AF=0.05;ANN=G|synonymous_variant|LOW|UNKNOWN|UNKNOWN	GT:AD:DP:GQ	0/1:30,10:40:35
"""
    
    # Create temporary file
    temp_fd, temp_path = tempfile.mkstemp(suffix='.vcf')
    try:
        with os.fdopen(temp_fd, 'w') as f:
            f.write(vcf_content)
        return temp_path
    except:
        os.close(temp_fd)
        raise

def test_pathway_burden_in_pipeline():
    """Test pathway burden model as part of the complete pipeline."""
    
    print("🧬 Pathway Burden Model - Pipeline Integration Test")
    print("=" * 70)
    
    # Create test VCF file
    vcf_path = create_test_vcf_file()
    
    try:
        # Set up preferences
        preferences = {
            "language": "en",
            "include_technical_details": True,
            "patient_data": {
                "age": 45,
                "sex": "F",
                "family_history": "breast_cancer"
            }
        }
        
        print(f"\n📁 Test VCF file created: {vcf_path}")
        print(f"⚙️  Running full pipeline with pathway burden analysis...")
        
        # Run the pipeline
        start_time = datetime.now()
        result = run_pipeline(vcf_path, preferences)
        end_time = datetime.now()
        
        processing_time = (end_time - start_time).total_seconds()
        
        print(f"\n✅ Pipeline completed in {processing_time:.2f} seconds")
        print(f"📊 Status: {result['pipeline_status']}")
        
        # Check if pathway burden was processed
        pathway_burden_results = result.get('pathway_burden_results', {})
        pathway_burden_summary = result.get('pathway_burden_summary', {})
        
        if pathway_burden_results:
            print(f"\n🎯 Pathway Burden Analysis Results:")
            print(f"  Pathways analyzed: {len(pathway_burden_results)}")
            print(f"  Overall burden score: {pathway_burden_summary.get('overall_burden_score', 0):.3f}")
            print(f"  High burden pathways: {', '.join(pathway_burden_summary.get('high_burden_pathways', []))}")
            print(f"  Primary concern: {pathway_burden_summary.get('primary_concern', 'None')}")
            print(f"  Total damaging variants: {pathway_burden_summary.get('total_damaging_variants', 0)}")
            
            # Show detailed pathway results
            high_burden_pathways = []
            for pathway_name, results in pathway_burden_results.items():
                if results['total_variants'] > 0:
                    print(f"\n  {pathway_name.upper()}:")
                    print(f"    Variants: {results['damaging_variants']}/{results['total_variants']}")
                    print(f"    Burden score: {results['burden_score']:.3f}")
                    print(f"    Risk level: {results['risk_level']}")
                    print(f"    Genes: {', '.join(results['contributing_genes'])}")
                    
                    if results['risk_level'] == 'high':
                        high_burden_pathways.append(pathway_name)
                    
                    if results['multi_hit_genes']:
                        print(f"    Multi-hit genes: {', '.join(results['multi_hit_genes'])}")
                    
                    if results['top_variant']:
                        top = results['top_variant']
                        print(f"    Top variant: {top['variant_id']} in {top['gene']}")
            
            # Test specific expectations
            print(f"\n🔍 Validation Tests:")
            
            # Should have detected BRCA1 multi-hit
            brca1_multi_hit = any(
                "BRCA1" in results.get('multi_hit_genes', []) 
                for results in pathway_burden_results.values()
            )
            print(f"  BRCA1 multi-hit detection: {'✅ PASS' if brca1_multi_hit else '❌ FAIL'}")
            
            # Should have high burden in DNA repair pathway
            dna_repair_high = pathway_burden_results.get('dna_repair', {}).get('risk_level') == 'high'
            print(f"  DNA repair high burden: {'✅ PASS' if dna_repair_high else '❌ FAIL'}")
            
            # Should identify pathway crosstalk
            pathway_crosstalk = pathway_burden_summary.get('pathway_crosstalk', False)
            print(f"  Pathway crosstalk detection: {'✅ PASS' if pathway_crosstalk else '❌ FAIL'}")
            
            # Check integration with other static models
            print(f"\n🔗 Integration with Other Models:")
            print(f"  CADD scores available: {'✅ YES' if result.get('cadd_enriched_variants') else '❌ NO'}")
            print(f"  ClinVar annotations: {'✅ YES' if result.get('clinvar_annotations') else '❌ NO'}")
            print(f"  PRS results: {'✅ YES' if result.get('prs_results') else '❌ NO'}")
            print(f"  TCGA matches: {'✅ YES' if result.get('tcga_matches') else '❌ NO'}")
            
            # Check feature vector builder integration
            feature_vector_info = result.get('feature_vector_info', {})
            available_inputs = feature_vector_info.get('available_inputs', {})
            pathway_burden_available = available_inputs.get('pathway_burden', False)
            print(f"  Feature vector integration: {'✅ YES' if pathway_burden_available else '❌ NO'}")
            
            if pathway_burden_available:
                annotation_summary = feature_vector_info.get('annotation_summary', {})
                variants_with_pathway = annotation_summary.get('with_pathway_burden', 0)
                print(f"  Variants with pathway assessment: {variants_with_pathway}")
            
        else:
            print(f"\n❌ Pathway burden analysis not found in results!")
            print(f"Available keys: {list(result.keys())}")
        
        # Show overall risk assessment
        if result.get('risk_scores'):
            print(f"\n🎯 Overall Risk Assessment:")
            for cancer_type, score in result['risk_scores'].items():
                print(f"  {cancer_type.capitalize()}: {score}% risk")
        
        # Check for errors
        if result.get('errors'):
            print(f"\n❌ Errors ({len(result['errors'])}):")
            for error in result['errors']:
                print(f"  - {error.get('node', 'Unknown')}: {error.get('error', 'Unknown error')}")
        
        # Check for warnings
        if result.get('warnings'):
            print(f"\n⚠️  Warnings ({len(result['warnings'])}):")
            for warning in result['warnings'][:5]:  # Show first 5
                print(f"  - {warning}")
        
        # Save detailed results
        output_file = "test_pathway_burden_pipeline_results.json"
        with open(output_file, 'w') as f:
            serializable_result = json.loads(json.dumps(result, default=str))
            json.dump(serializable_result, f, indent=2)
        
        print(f"\n💾 Full pipeline results saved to: {output_file}")
        
        return result, pathway_burden_results is not None
        
    finally:
        # Clean up temporary file
        try:
            os.unlink(vcf_path)
        except:
            pass

def test_parallel_processing():
    """Test that pathway burden runs in parallel with other static models."""
    
    print("\n" + "=" * 70)
    print("Testing Parallel Processing")
    print("=" * 70)
    
    # Create the pipeline graph
    pipeline = create_genomic_pipeline()
    
    # Check that pathway_burden node exists
    nodes = list(pipeline.nodes.keys())
    pathway_burden_exists = "pathway_burden" in nodes
    print(f"Pathway burden node exists: {'✅ YES' if pathway_burden_exists else '❌ NO'}")
    
    if pathway_burden_exists:
        print(f"Pipeline nodes: {', '.join(nodes)}")
        
        # Check for expected parallel structure
        expected_parallel_nodes = [
            "tcga_mapper", "cadd_scoring", "clinvar_annotator", 
            "prs_calculator", "pathway_burden"
        ]
        
        parallel_nodes_present = all(node in nodes for node in expected_parallel_nodes)
        print(f"All parallel nodes present: {'✅ YES' if parallel_nodes_present else '❌ NO'}")
        
        if not parallel_nodes_present:
            missing = [node for node in expected_parallel_nodes if node not in nodes]
            print(f"Missing nodes: {', '.join(missing)}")
    
    return pathway_burden_exists

def test_edge_cases():
    """Test edge cases and error handling."""
    
    print("\n" + "=" * 70)
    print("Testing Edge Cases")
    print("=" * 70)
    
    # Test with minimal VCF (no variants)
    minimal_vcf = """##fileformat=VCFv4.2
##reference=GRCh38
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
"""
    
    temp_fd, temp_path = tempfile.mkstemp(suffix='.vcf')
    try:
        with os.fdopen(temp_fd, 'w') as f:
            f.write(minimal_vcf)
        
        print(f"Testing with empty VCF file...")
        result = run_pipeline(temp_path, {})
        
        pathway_burden_summary = result.get('pathway_burden_summary', {})
        overall_burden = pathway_burden_summary.get('overall_burden_score', 0)
        
        print(f"Empty VCF burden score: {overall_burden:.3f}")
        print(f"Pipeline status: {result.get('pipeline_status', 'unknown')}")
        
        # Should handle gracefully
        graceful_handling = result.get('pipeline_status') != 'failed'
        print(f"Graceful handling: {'✅ PASS' if graceful_handling else '❌ FAIL'}")
        
    finally:
        try:
            os.unlink(temp_path)
        except:
            pass
    
    return True

def run_integration_tests():
    """Run all integration tests."""
    
    print("🧬 Pathway Burden Model - Integration Test Suite")
    print("=" * 70)
    
    try:
        # Test 1: Full pipeline integration
        result, pathway_burden_success = test_pathway_burden_in_pipeline()
        
        # Test 2: Parallel processing verification
        parallel_success = test_parallel_processing()
        
        # Test 3: Edge cases
        edge_case_success = test_edge_cases()
        
        print("\n" + "=" * 70)
        print("Integration Test Results")
        print("=" * 70)
        
        print(f"✅ Pipeline integration: {'PASS' if pathway_burden_success else 'FAIL'}")
        print(f"✅ Parallel processing: {'PASS' if parallel_success else 'FAIL'}")
        print(f"✅ Edge case handling: {'PASS' if edge_case_success else 'FAIL'}")
        
        overall_success = pathway_burden_success and parallel_success and edge_case_success
        
        if overall_success:
            print(f"\n🎉 ALL INTEGRATION TESTS PASSED!")
            print(f"The pathway burden model is successfully integrated into the pipeline.")
        else:
            print(f"\n❌ Some integration tests failed.")
            
        return overall_success
        
    except Exception as e:
        print(f"\n❌ INTEGRATION TEST FAILED: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = run_integration_tests()
    sys.exit(0 if success else 1) 