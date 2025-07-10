#!/usr/bin/env python3
"""
Test script to verify parallel execution results after fixes.
"""
import logging
import json
from graph import run_pipeline
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(message)s'
)


class DateTimeEncoder(json.JSONEncoder):
    """Custom JSON encoder to handle datetime objects."""
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


def test_parallel_results():
    """Test that all parallel nodes execute correctly with proper results."""
    print("\n🧪 Testing Fixed Parallel Execution Results\n")
    
    # Use a test VCF file
    test_file = "test_data/test_variants.vcf"
    
    # Run the pipeline
    result = run_pipeline(
        test_file,
        {
            "language": "en",
            "include_technical_details": True,
            "patient_data": {"age": 45, "sex": "male"}
        }
    )
    
    # Check pipeline status
    print(f"✅ Pipeline Status: {result.get('pipeline_status')}\n")
    
    # Check completed nodes
    completed_nodes = result.get("completed_nodes", [])
    unique_nodes = list(set(completed_nodes))
    print(f"📊 Completed Nodes: {len(unique_nodes)}")
    
    # Check for duplicates
    if len(completed_nodes) == len(unique_nodes):
        print("✅ No duplicate nodes - state management working correctly!")
    else:
        print(f"❌ Found duplicate nodes: {len(completed_nodes)} total, {len(unique_nodes)} unique")
    
    # Check parallel node status
    parallel_status = result.get("parallel_nodes_status", {})
    print(f"\n🔄 Parallel Node Status:")
    for node, completed in parallel_status.items():
        status = "✅" if completed else "❌"
        print(f"  {status} {node}")
    
    # Check if all parallel nodes completed
    if all(parallel_status.values()):
        print("✅ All parallel nodes completed successfully!")
    else:
        missing = [node for node, completed in parallel_status.items() if not completed]
        print(f"❌ Missing nodes: {missing}")
    
    # Check parallel node results
    print(f"\n📝 Parallel Node Results:")
    print(f"  {'✅' if result.get('tcga_matches') else '❌'} TCGA Matches: {'Found' if result.get('tcga_matches') else 'Missing'}")
    print(f"  {'✅' if result.get('cadd_stats') else '❌'} CADD Stats: {'Found' if result.get('cadd_stats') else 'Missing'}")
    print(f"  {'✅' if result.get('clinvar_annotations') else '❌'} ClinVar Annotations: {'Found' if result.get('clinvar_annotations') else 'Missing'}")
    print(f"  {'✅' if result.get('prs_results') else '❌'} PRS Results: {'Found' if result.get('prs_results') else 'Missing'}")
    print(f"  {'✅' if result.get('pathway_burden_results') else '❌'} Pathway Burden Results: {'Found' if result.get('pathway_burden_results') else 'Missing'}")
    
    # Check merge tracking
    print(f"\n🔀 Merge Node Tracking:")
    print(f"  Parallel merge count: {result.get('parallel_merge_count', 0)}")
    
    # Check for errors
    errors = result.get("errors", [])
    if errors:
        print(f"\n❌ Errors detected: {len(errors)}")
        for error in errors:
            print(f"  - {error}")
    else:
        print(f"\n✅ No errors detected")
    
    # Save detailed results
    with open("test_parallel_results.json", "w") as f:
        # Create a serializable version of the results
        serializable_result = {
            "pipeline_status": result.get("pipeline_status"),
            "completed_nodes": unique_nodes,
            "node_execution_count": len(completed_nodes),
            "unique_node_count": len(unique_nodes),
            "parallel_nodes_status": parallel_status,
            "parallel_merge_count": result.get("parallel_merge_count", 0),
            "errors": errors,
            "has_tcga": bool(result.get("tcga_matches")),
            "has_cadd": bool(result.get("cadd_stats")),
            "has_clinvar": bool(result.get("clinvar_annotations")),
            "has_prs": bool(result.get("prs_results")),
            "has_pathway": bool(result.get("pathway_burden_results")),
            "clinvar_annotation_count": len(result.get("clinvar_annotations", {})),
            "tcga_match_count": sum(len(matches) for matches in result.get("tcga_matches", {}).values()),
            "prs_cancer_types": len(result.get("prs_results", {})),
            "pathway_count": len(result.get("pathway_burden_results", {}))
        }
        json.dump(serializable_result, f, indent=2, cls=DateTimeEncoder)
    
    print(f"\n✅ Test results saved to test_parallel_results.json")
    
    return result


if __name__ == "__main__":
    test_parallel_results() 