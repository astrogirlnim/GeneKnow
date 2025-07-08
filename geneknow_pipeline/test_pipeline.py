"""
Test script for LangGraph genomic pipeline.
Runs the pipeline with mock data to verify end-to-end flow.
"""
import json
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from graph import run_pipeline


def test_pipeline():
    """Run pipeline with test data and display results."""
    
    print("🧬 GeneKnow Pipeline Test - Running with MOCK data")
    print("=" * 60)
    
    # Test with FASTQ file
    test_file = "test-data/sample.fastq.gz"  # ⚠️ MOCK file path
    
    # User preferences
    preferences = {
        "language": "en",
        "include_technical_details": True,
        "risk_threshold_percentage": 30.0  # Lower threshold to see more results
    }
    
    print(f"\n📁 Input file: {test_file}")
    print(f"⚙️  Preferences: {json.dumps(preferences, indent=2)}")
    print("\n🚀 Starting pipeline...\n")
    
    # Run the pipeline
    result = run_pipeline(test_file, preferences)
    
    # Display results
    print("\n✅ Pipeline Results:")
    print("=" * 60)
    
    print(f"\n📊 Status: {result['pipeline_status']}")
    processing_time = result.get('processing_time_seconds')
    if processing_time is not None:
        print(f"⏱️  Processing time: {processing_time:.2f} seconds")
    print(f"✓ Completed nodes: {', '.join(result.get('completed_nodes', []))}")
    
    if result['errors']:
        print(f"\n❌ Errors:")
        for error in result['errors']:
            print(f"  - {error['node']}: {error['error']}")
    
    if result['warnings']:
        print(f"\n⚠️  Warnings ({len(result['warnings'])}):")
        for warning in result['warnings'][:3]:  # Show first 3
            print(f"  - {warning['node']}: {warning['warning']}")
    
    # Show risk assessment
    if 'risk_scores' in result:
        print(f"\n🎯 Risk Assessment:")
        for cancer_type, score in result['risk_scores'].items():
            genes = result['risk_genes'].get(cancer_type, [])
            print(f"  - {cancer_type.capitalize()}: {score}% risk")
            if genes:
                print(f"    Affected genes: {', '.join(genes)}")
    
    # Show variant summary
    if 'filtered_variants' in result:
        print(f"\n🧬 Variants Summary:")
        print(f"  - Total found: {result.get('variant_count', 0)}")
        print(f"  - Passed QC: {len(result['filtered_variants'])}")
        print(f"\n  Top variants:")
        for variant in result['filtered_variants'][:3]:
            print(f"    • {variant['gene']} - {variant['variant_id']}")
    
    # Show report preview
    if 'report_markdown' in result and result['report_markdown']:
        print(f"\n📄 Report Preview (first 500 chars):")
        print("-" * 40)
        print(result['report_markdown'][:500] + "...")
    
    # Save full results
    output_file = "test_pipeline_results.json"
    with open(output_file, 'w') as f:
        # Convert datetime objects to strings for JSON serialization
        serializable_result = json.loads(
            json.dumps(result, default=str)
        )
        json.dump(serializable_result, f, indent=2)
    
    print(f"\n💾 Full results saved to: {output_file}")
    
    return result


if __name__ == "__main__":
    test_pipeline() 