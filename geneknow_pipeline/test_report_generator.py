#!/usr/bin/env python3
"""
Test script for the report generator module (Phase 1).
Tests LLM detection, fallback mode, and basic report generation.
"""

import sys
from pathlib import Path

# Add the parent directory to the path so we can import the module
sys.path.insert(0, str(Path(__file__).parent))

from nodes.report_generator import process as generate_report
from nodes.report_generator.config import ReportConfig, LLMBackend, ReportStyle
from nodes.report_generator.model_interface import ModelInterface


def create_sample_state():
    """Create a sample state dict that mimics what the pipeline would pass."""
    return {
        "case_id": "TEST_CASE_001",
        "structured_json": {
            "case_id": "TEST_CASE_001",
            "processing_time_seconds": 45.2,
            "variant_count": 1250,
            "risk_scores": {
                "breast": 82.5,
                "ovarian": 45.2,
                "colon": 15.3,
                "lung": 8.7,
                "prostate": 3.2,
                "pancreatic": 2.1,
            },
            "variants": [
                {
                    "gene": "BRCA1",
                    "position": "chr17:41246747",
                    "type": "frameshift",
                    "impact": "HIGH",
                    "clinvar_significance": "Pathogenic",
                    "cadd_score": 35.0,
                    "tcga_cancer_relevance": 0.95,
                    "pathway_damage_assessment": {
                        "pathways_affected": ["DNA repair", "Cell cycle"],
                        "damage_score": 0.89,
                    },
                },
                {
                    "gene": "TP53",
                    "position": "chr17:7577121",
                    "type": "missense",
                    "impact": "MODERATE",
                    "clinvar_significance": "Likely pathogenic",
                    "cadd_score": 28.5,
                    "tcga_cancer_relevance": 0.75,
                    "pathway_damage_assessment": {
                        "pathways_affected": ["p53 signaling"],
                        "damage_score": 0.72,
                    },
                },
            ],
            "metrics": {
                "variant_metrics": {
                    "total_variants": 1250,
                    "high_impact": 45,
                    "moderate_impact": 230,
                    "low_impact": 975,
                },
                "confidence_metrics": {
                    "overall_confidence": 0.92,
                    "data_quality": 0.98,
                },
            },
            "report_sections": {
                "overview": {
                    "title": "Analysis Overview",
                    "content": "Genomic analysis identified significant cancer risk factors.",
                    "severity": "high",
                },
                "key_findings": {
                    "title": "Key Findings",
                    "content": "BRCA1 pathogenic variant detected with high cancer relevance.",
                    "severity": "high",
                },
            },
        },
    }


def test_llm_detection():
    """Test LLM backend detection."""
    print("\n=== Testing LLM Detection ===")

    config = ReportConfig()
    interface = ModelInterface(config)
    print(f"Detected backend: {interface.backend}")
    print(f"Available: {interface.is_available()}")

    if interface.backend == LLMBackend.OLLAMA:
        print("Ollama is available!")
        # List models if Ollama is running
        try:
            import requests

            response = requests.get("http://localhost:11434/api/tags")
            if response.status_code == 200:
                models = response.json().get("models", [])
                print(f"Available Ollama models: {[m['name'] for m in models]}")
        except:
            pass
    elif interface.backend == LLMBackend.HUGGINGFACE:
        print("HuggingFace is available!")
    else:
        print("No LLM backend available - will use fallback mode")


def test_report_generation():
    """Test the full report generation process."""
    print("\n=== Testing Report Generation ===")

    # Create sample state
    state = create_sample_state()

    # Test with different configurations
    # Note: We need to save configs to file for the process function to pick them up
    from nodes.report_generator.config import save_config, load_config

    configs = [
        ReportConfig(),  # Default config
        ReportConfig(style=ReportStyle.TECHNICAL),
        ReportConfig(style=ReportStyle.PATIENT),
    ]

    for i, config in enumerate(configs):
        print(f"\n--- Test {i+1}: {config.style.value} style ---")

        # Only save config if it's different from what's already there
        # Don't override backend settings
        if i == 0:  # Only save the first config to avoid overriding backend
            current_config = load_config()
            config.backend = current_config.backend  # Preserve backend setting
            save_config(config)

        # Generate report
        result = generate_report(state)

        # Check results
        print(f"Original sections preserved: {'report_sections' in result}")
        print(f"Enhanced reports generated: {'enhanced_report_paths' in result}")
        print(f"Report info available: {'report_generator_info' in result}")

        if "enhanced_report_paths" in result:
            for fmt, path in result["enhanced_report_paths"].items():
                print(f"  {fmt}: {path}")

                # Read and display first 500 chars of markdown
                if fmt == "markdown" and Path(path).exists():
                    with open(path, "r") as f:
                        content = f.read()
                        print(f"\n  First 500 chars of {fmt}:")
                        print("  " + "-" * 50)
                        print(content[:500] + "..." if len(content) > 500 else content)
                        print("  " + "-" * 50)

                        # Check for key indicators
                        if "NON-LLM MODE" in content:
                            print("  ✓ Fallback mode indicator found")
                        if "Medical Glossary" in content:
                            print("  ✓ Glossary section found")
                        if "BRCA1" in content:
                            print("  ✓ High-risk variant mentioned")

        if "report_generator_info" in result:
            info = result["report_generator_info"]
            print("\nReport Generator Info:")
            print(f"  Backend: {info.get('backend', 'unknown')}")
            print(f"  Model: {info.get('model', 'none')}")
            print(f"  Generation time: {info.get('generation_time_seconds', 0):.2f}s")


def test_streaming():
    """Test streaming capability."""
    print("\n=== Testing Streaming Support ===")

    config = ReportConfig()
    interface = ModelInterface(config)
    if interface.is_available():
        print("Testing streaming with LLM...")
        chunks = list(
            interface.generate(
                "Generate a brief test report about BRCA1 mutation.", stream=True
            )
        )
        print(f"Received {len(chunks)} chunks")
        if chunks:
            print(f"First chunk: {chunks[0][:50]}...")
    else:
        print("No LLM available - streaming test skipped")


def main():
    """Run all tests."""
    print("GeneKnow Report Generator - Phase 1 Test Suite")
    print("=" * 50)

    # Test 1: LLM Detection
    test_llm_detection()

    # Test 2: Report Generation
    test_report_generation()

    # Test 3: Streaming
    test_streaming()

    print("\n" + "=" * 50)
    print("Phase 1 testing complete!")
    print("\nNext steps:")
    print("1. If no LLM was detected, install Ollama (https://ollama.ai)")
    print("2. Pull a model: ollama pull llama2 (or mistral/codellama)")
    print("3. Re-run this test to see LLM-enhanced reports")


if __name__ == "__main__":
    main()
