#!/usr/bin/env python3
"""
Test script to verify parallel execution of TCGA, CADD, ClinVar, PRS, and pathway burden nodes.
"""
import logging
from graph import run_pipeline
import json

# Set up logging to see detailed information
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def test_parallel_execution():
    """Test that all 5 parallel nodes execute correctly."""
    print("\n🧪 Testing Parallel Execution of Pipeline Nodes\n")
    
    # Use a test VCF file with variants
    test_file = "test_data/test_variants.vcf"
    
    # Run the pipeline
    result = run_pipeline(
        test_file,
        {
            "language": "en",
            "include_technical_details": True,
            "patient_data": {"age": 45, "sex": "F"}
        }
    )
    
    # Check pipeline completion
    print(f"Pipeline Status: {result.get('pipeline_status')}")
    print(f"Completed Nodes: {result.get('completed_nodes', [])}")
    
    # Check each parallel node's results
    parallel_nodes = {
        'tcga_mapper': ('tcga_matches', 'tcga_summary'),
        'cadd_scoring': ('cadd_enriched_variants', 'cadd_stats'),
        'clinvar_annotator': ('clinvar_annotations', 'clinvar_stats'),
        'prs_calculator': ('prs_results', 'prs_summary'),
        'pathway_burden': ('pathway_burden_results', 'pathway_burden_summary')
    }
    
    print("\n📊 Parallel Node Results:")
    for node_name, (result_key, summary_key) in parallel_nodes.items():
        has_results = bool(result.get(result_key))
        has_summary = bool(result.get(summary_key))
        status = "✓" if (has_results or has_summary) else "✗"
        
        print(f"\n{status} {node_name}:")
        print(f"  - Has {result_key}: {has_results}")
        print(f"  - Has {summary_key}: {has_summary}")
        
        if has_results:
            results = result.get(result_key)
            if isinstance(results, dict):
                print(f"  - Result count: {len(results)}")
            elif isinstance(results, list):
                print(f"  - Result count: {len(results)}")
        
        if has_summary:
            summary = result.get(summary_key)
            if isinstance(summary, dict):
                print(f"  - Summary keys: {list(summary.keys())[:5]}...")
    
    # Check merge node execution count
    merge_calls = result.get('_merge_static_calls', 0)
    print(f"\n🔄 Merge Static Models Calls: {merge_calls}")
    
    # Check for pathway enriched variants
    has_pathway_enriched = bool(result.get('pathway_enriched_variants'))
    print(f"\n📦 Has pathway_enriched_variants: {has_pathway_enriched}")
    
    # Save detailed results for debugging
    output_file = "test_parallel_execution_results.json"
    with open(output_file, 'w') as f:
        # Filter out large data for readability
        filtered_result = {
            k: v for k, v in result.items() 
            if k not in ['raw_variants', 'filtered_variants', 'cadd_enriched_variants', 'pathway_enriched_variants']
        }
        json.dump(filtered_result, f, indent=2, default=str)
    
    print(f"\n💾 Detailed results saved to: {output_file}")
    
    # Determine success
    all_complete = all(
        node_name in result.get('completed_nodes', [])
        for node_name in parallel_nodes.keys()
    )
    
    if all_complete:
        print("\n✅ All parallel nodes completed successfully!")
    else:
        print("\n❌ Some parallel nodes did not complete")
        missing = [
            node for node in parallel_nodes.keys()
            if node not in result.get('completed_nodes', [])
        ]
        print(f"Missing nodes: {missing}")
    
    return all_complete


if __name__ == "__main__":
    success = test_parallel_execution()
    exit(0 if success else 1) 