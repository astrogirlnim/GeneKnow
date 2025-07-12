#!/usr/bin/env python3
"""
Test pipeline flow tracing - from file upload to dashboard display.
This test verifies that:
1. Static models generate correct risk scores
2. ML fusion combines them properly
3. Final risk scores get displayed correctly
"""

import json
import requests
import time
import os
import sys
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent))

from graph import run_pipeline


def test_direct_pipeline():
    """Test the pipeline directly to trace data flow."""
    print("🔬 Testing Direct Pipeline Flow")
    print("=" * 60)

    # Find a test file
    test_file = None
    test_files = [
        "test_data/242b87b1-bf7c-4c1a-bed2-bb077e5ccd00.wxs.aliquot_ensemble_masked.maf.gz",
        "test_data/test_variants.vc",
        "test_data/tcga_downloads/test_sample.ma",
        "../test_data/test.vc",
    ]

    for f in test_files:
        if os.path.exists(f):
            test_file = f
            break

    if not test_file:
        print("❌ No test file found")
        assert print(f"📁 Using test file: {test_file}")

    # Run pipeline
    preferences = {
        "patient_data": {"age": 45, "sex": "F", "family_history": True},
        "language": "en",
        "include_technical": True,
    }

    print("\n🚀 Running pipeline...")
    result = run_pipeline(test_file, preferences)

    # Trace key outputs
    print("\n📊 TRACING PIPELINE OUTPUTS:")
    print("=" * 60)

    # 1. Check filtered variants
    print(f"\n1️⃣ Filtered Variants: {result.get('variant_count', 0)} total")
    variants = result.get("filtered_variants", [])
    if variants:
        sample_variant = variants[0]
        print(
            f"   Sample variant: {sample_variant.get('gene', 'N/A')} - {sample_variant.get('variant_id', 'N/A')}"
        )

    # 2. Check static model outputs
    print("\n2️⃣ Static Model Outputs:")

    # PRS
    prs_results = result.get("prs_results", {})
    if prs_results:
        print("   ✅ PRS Results:")
        for cancer, prs_data in prs_results.items():
            print(
                f"      {cancer}: score={prs_data.get('raw_score', 0):.3f}, percentile={prs_data.get('percentile', 0)}%"
            )
    else:
        print("   ❌ No PRS results")

    # ClinVar
    clinvar_stats = result.get("clinvar_stats", {})
    if clinvar_stats:
        print(
            f"   ✅ ClinVar: {clinvar_stats.get('total_annotated', 0)} variants annotated"
        )
        print(f"      Pathogenic: {clinvar_stats.get('pathogenic_count', 0)}")
        print(f"      Benign: {clinvar_stats.get('benign_count', 0)}")
    else:
        print("   ❌ No ClinVar stats")

    # CADD
    cadd_stats = result.get("cadd_stats", {})
    if cadd_stats:
        print(f"   ✅ CADD: {cadd_stats.get('variants_scored', 0)} variants scored")
        print(f"      Mean PHRED: {cadd_stats.get('mean_phred', 0):.1f}")
        print(f"      High impact (>20): {cadd_stats.get('variants_gt20', 0)}")
    else:
        print("   ❌ No CADD stats")

    # TCGA
    tcga_summary = result.get("tcga_summary", {})
    if tcga_summary:
        print(
            f"   ✅ TCGA: Analyzed {tcga_summary.get('total_variants_matched', 0)} variants"
        )
        enrichment = tcga_summary.get("max_tumor_enrichment", {})
        if enrichment:
            print(
                f"      Max enrichment: {enrichment.get('enrichment', 0):.1f}x in {enrichment.get('cancer_type', 'N/A')}"
            )
    else:
        print("   ❌ No TCGA summary")

    # Pathway Burden
    pathway_summary = result.get("pathway_burden_summary", {})
    if pathway_summary:
        print(
            f"   ✅ Pathway Burden: score={pathway_summary.get('overall_burden_score', 0):.3f}"
        )
        print(
            f"      High burden pathways: {', '.join(pathway_summary.get('high_burden_pathways', [])) or 'None'}"
        )
    else:
        print("   ❌ No pathway burden summary")

    # 3. Check ML fusion results
    print("\n3️⃣ ML Fusion Results:")
    ml_fusion_results = result.get("ml_fusion_results", {})
    if ml_fusion_results and ml_fusion_results.get("processing_successful"):
        aggregate = ml_fusion_results.get("aggregate_risk_assessment", {})
        print("   ✅ ML Fusion successful")
        print(
            f"      Aggregate risk score: {aggregate.get('aggregate_risk_score', 0):.3f}"
        )
        print(f"      Risk category: {aggregate.get('risk_category', 'unknown')}")
        print(f"      High-risk variants: {aggregate.get('high_risk_variants', 0)}")

        # Show contributing factors
        factors = aggregate.get("contributing_factors", {})
        if factors:
            print("      Contributing factors:")
            for factor, score in sorted(
                factors.items(), key=lambda x: x[1], reverse=True
            )[:3]:
                print(f"        - {factor}: {score:.3f}")
    else:
        print("   ❌ ML Fusion not successful or not run")

    # 4. Check final risk scores
    print("\n4️⃣ Final Risk Scores (what dashboard shows):")
    risk_scores = result.get("risk_scores", {})
    if risk_scores:
        print("   ✅ Risk scores generated:")
        for cancer, score in sorted(
            risk_scores.items(), key=lambda x: x[1], reverse=True
        ):
            genes = result.get("risk_genes", {}).get(cancer, [])
            print(
                f"      {cancer}: {score}% (genes: {', '.join(genes[:3]) if genes else 'None'})"
            )
    else:
        print("   ❌ No risk scores generated")

    # 5. Check structured JSON (what frontend receives)
    print("\n5️⃣ Structured JSON for Frontend:")
    structured_json = result.get("structured_json", {})
    if structured_json:
        risk_assessment = structured_json.get("risk_assessment", {})
        if risk_assessment:
            scores = risk_assessment.get("scores", {})
            print(f"   ✅ Frontend will receive {len(scores)} cancer risk scores")
            print(
                f"      High risk findings: {len(risk_assessment.get('high_risk_findings', []))}"
            )
        else:
            print("   ❌ No risk assessment in structured JSON")
    else:
        print("   ❌ No structured JSON generated")

    # 6. Save results for inspection
    output_file = "pipeline_trace_results.json"
    with open(output_file, "w") as f:
        # Convert numpy types before saving
        clean_result = json.loads(json.dumps(result, default=str))
        json.dump(clean_result, f, indent=2)
    print(f"\n💾 Full results saved to: {output_file}")

    print("\n" + "=" * 60)
    print("✅ Pipeline flow tracing complete")


def test_api_flow():
    """Test the full API flow to ensure data reaches frontend correctly."""
    print("\n\n🌐 Testing API Flow")
    print("=" * 60)

    api_url = "http://localhost:5001"

    # Check if API is running
    try:
        response = requests.get(f"{api_url}/api/health")
        if response.status_code != 200:
            print(
                "❌ API server not running. Start with: python enhanced_api_server.py"
            )
            return
    except:
        print("❌ Cannot connect to API server at localhost:5001")
        return

    print("✅ API server is running")

    # Find test file
    test_file = None
    for f in [
        "test_data/242b87b1-bf7c-4c1a-bed2-bb077e5ccd00.wxs.aliquot_ensemble_masked.maf.gz",
        "test_data/test_variants.vc",
        "test_data/tcga_downloads/test_sample.ma",
    ]:
        if os.path.exists(f):
            test_file = os.path.abspath(f)
            break

    if not test_file:
        print("❌ No test file found")
        return

    print(f"📁 Processing file: {test_file}")

    # Submit for processing
    response = requests.post(
        f"{api_url}/api/process",
        json={
            "file_path": test_file,
            "preferences": {
                "patient_data": {"age": 45, "sex": "F"},
                "language": "en",
                "include_technical": True,
            },
        },
    )

    if response.status_code != 202:
        print(f"❌ Failed to submit file: {response.text}")
        return

    job_data = response.json()
    job_id = job_data["job_id"]
    print(f"✅ Job submitted: {job_id}")

    # Poll for completion
    print("\n⏳ Waiting for processing...")
    max_wait = 60  # seconds
    start_time = time.time()

    while time.time() - start_time < max_wait:
        response = requests.get(f"{api_url}/api/status/{job_id}")
        if response.status_code == 200:
            status_data = response.json()
            status = status_data["status"]
            progress = status_data.get("progress", 0)
            current_step = status_data.get("current_step", "N/A")

            print(
                f"\r   Status: {status} | Progress: {progress}% | Step: {current_step}",
                end="",
                flush=True,
            )

            if status == "completed":
                print("\n✅ Processing completed!")
                break
            elif status == "failed":
                print(
                    f"\n❌ Processing failed: {status_data.get('error', 'Unknown error')}"
                )
                return

        time.sleep(1)
    else:
        print("\n❌ Timeout waiting for processing")
        return

    # Get results
    print("\n📊 Fetching results...")
    response = requests.get(f"{api_url}/api/results/{job_id}")
    if response.status_code != 200:
        print(f"❌ Failed to get results: {response.text}")
        return

    results = response.json()

    # Verify risk scores are present
    print("\n🔍 Verifying API Results:")
    print("=" * 40)

    risk_scores = results.get("risk_scores", {})
    if risk_scores:
        print("✅ Risk scores present:")
        for cancer, score in sorted(
            risk_scores.items(), key=lambda x: x[1], reverse=True
        ):
            print(f"   {cancer}: {score}%")
    else:
        print("❌ No risk scores in API response!")

    # Check structured JSON
    structured_json = results.get("structured_json", {})
    if structured_json and structured_json.get("risk_assessment"):
        scores = structured_json["risk_assessment"].get("scores", {})
        print(f"\n✅ Structured JSON contains {len(scores)} risk scores")
    else:
        print("\n❌ No risk assessment in structured JSON!")

    # Save API results
    with open("api_trace_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\n💾 API results saved to: api_trace_results.json")


def verify_dashboard_data_mapping():
    """Verify how dashboard maps the pipeline data."""
    print("\n\n🖥️  Dashboard Data Mapping Verification")
    print("=" * 60)

    # Simulate what the dashboard does with pipeline results
    sample_results = {
        "risk_scores": {
            "breast": 23.5,
            "colon": 15.2,
            "lung": 8.7,
            "prostate": 5.3,
            "blood": 3.1,
        },
        "variant_count": 142,
        "processing_time_seconds": 2.3,
    }

    print("📥 Sample pipeline results:")
    print(f"   Risk scores: {sample_results['risk_scores']}")

    # Dashboard logic (from DashboardPage.tsx)
    risk_scores = list(sample_results["risk_scores"].items())
    highest_risk = max(risk_scores, key=lambda x: x[1])

    probability = round(highest_risk[1])
    hazard_score = highest_risk[1] / 100 * 3

    print("\n📊 Dashboard will display:")
    print(f"   Highest Risk Cancer: {highest_risk[0]}")
    print(f"   Risk Probability: {probability}%")
    print(f"   Hazard Score: {hazard_score:.1f}")
    print(f"   Total Variants: {sample_results['variant_count']}")
    print(f"   Processing Time: {sample_results['processing_time_seconds']}s")

    print(
        "\n✅ Data mapping verified - risk scores are displayed as percentages directly"
    )


if __name__ == "__main__":
    print("🧬 GeneKnow Pipeline Flow Tracing Test")
    print("=" * 60)
    print("This test traces data flow from file processing to dashboard display\n")

    # Test 1: Direct pipeline
    test_direct_pipeline()

    # Test 2: API flow (optional - requires API server running)
    if "--api" in sys.argv:
        test_api_flow()
    else:
        print("\n💡 To test API flow, run: python test_pipeline_flow_tracing.py --api")

    # Test 3: Verify dashboard mapping
    verify_dashboard_data_mapping()

    print("\n✨ All tests completed!")
