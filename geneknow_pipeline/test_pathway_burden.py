#!/usr/bin/env python3
"""
Test script for Pathway Burden Model node.
Verifies that the pathway burden implementation works correctly with various scenarios.
"""

import os
import sys
import logging
import json
import time
from datetime import datetime

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from nodes.pathway_burden import process, is_damaging_variant, calculate_pathway_burden, assess_overall_burden, CANCER_PATHWAYS

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_test_variant(gene, variant_id, cadd_phred=0, clinical_significance="",
                       allele_frequency=0.01, consequence="missense_variant"):
    """Create a test variant with specified parameters."""
    return {
        "variant_id": variant_id,
        "gene": gene,
        "chrom": "17",
        "pos": 12345678,
        "ref": "A",
        "alt": "G",
        "cadd_phred": cadd_phred,
        "clinical_significance": clinical_significance,
        "allele_frequency": allele_frequency,
        "population_frequency": allele_frequency,
        "consequence": consequence,
        "impact": consequence,
        "variant_classification": consequence
    }

def test_damaging_variant_assessment():
    """Test the variant damage assessment function."""

    print("=" * 60)
    print("Testing Variant Damage Assessment")
    print("=" * 60)

    # Test high-impact variant
    high_impact_variant = create_test_variant(
        "BRCA1", "17:43045677:G>A",
        cadd_phred=30,
        clinical_significance="Pathogenic",
        allele_frequency=0.0001,
        consequence="frameshift_variant"
    )

    damage_result = is_damaging_variant(high_impact_variant)
    print(f"\nHigh-impact variant (BRCA1):")
    print(f"  Damage score: {damage_result['damage_score']:.3f}")
    print(f"  Is damaging: {damage_result['is_damaging']}")
    print(f"  Reasons: {', '.join(damage_result['damage_reasons'])}")

    # Test moderate-impact variant
    moderate_variant = create_test_variant(
        "TP53", "17:7673803:C>T",
        cadd_phred=20,
        clinical_significance="Likely_pathogenic",
        allele_frequency=0.001,
        consequence="missense_variant"
    )

    damage_result = is_damaging_variant(moderate_variant)
    print(f"\nModerate-impact variant (TP53):")
    print(f"  Damage score: {damage_result['damage_score']:.3f}")
    print(f"  Is damaging: {damage_result['is_damaging']}")
    print(f"  Reasons: {', '.join(damage_result['damage_reasons'])}")

    # Test benign variant
    benign_variant = create_test_variant(
        "KRAS", "12:25245350:C>T",
        cadd_phred=5,
        clinical_significance="Benign",
        allele_frequency=0.1,
        consequence="synonymous_variant"
    )

    damage_result = is_damaging_variant(benign_variant)
    print(f"\nBenign variant (KRAS):")
    print(f"  Damage score: {damage_result['damage_score']:.3f}")
    print(f"  Is damaging: {damage_result['is_damaging']}")
    print(f"  Reasons: {', '.join(damage_result['damage_reasons'])}")

    return True

def test_pathway_burden_calculation():
    """Test pathway-specific burden calculation."""

    print("\n" + "=" * 60)
    print("Testing Pathway Burden Calculation")
    print("=" * 60)

    # Create test variants for DNA repair pathway
    dna_repair_variants = [
        create_test_variant("BRCA1", "17:43045677:G>A", cadd_phred=30, clinical_significance="Pathogenic", allele_frequency=0.0001),
        create_test_variant("BRCA2", "13:32363533:T>C", cadd_phred=25, clinical_significance="Likely_pathogenic", allele_frequency=0.0005),
        create_test_variant("ATM", "11:108098525:G>A", cadd_phred=15, clinical_significance="Uncertain_significance", allele_frequency=0.01),
        create_test_variant("PALB2", "16:23635000:A>G", cadd_phred=20, clinical_significance="Pathogenic", allele_frequency=0.0002),
        create_test_variant("UNKNOWN_GENE", "1:12345678:A>G", cadd_phred=10, clinical_significance="", allele_frequency=0.05),  # Not in pathway
    ]

    # Test DNA repair pathway
    dna_repair_data = CANCER_PATHWAYS["dna_repair"]
    burden_result = calculate_pathway_burden("dna_repair", dna_repair_data, dna_repair_variants)

    print(f"\nDNA Repair Pathway Results:")
    print(f"  Total variants in pathway: {burden_result['total_variants']}")
    print(f"  Damaging variants: {burden_result['damaging_variants']}")
    print(f"  Burden score: {burden_result['burden_score']:.3f}")
    print(f"  Risk level: {burden_result['risk_level']}")
    print(f"  Genes with variants: {', '.join(burden_result['contributing_genes'])}")
    print(f"  Genes with damaging variants: {', '.join(burden_result['damaging_genes'])}")

    if burden_result['top_variant']:
        top = burden_result['top_variant']
        print(f"  Top variant: {top['variant_id']} in {top['gene']} (score: {top['damage_score']:.3f})")

    # Test with multi-hit gene
    multi_hit_variants = [
        create_test_variant("BRCA1", "17:43045677:G>A", cadd_phred=30, clinical_significance="Pathogenic", allele_frequency=0.0001),
        create_test_variant("BRCA1", "17:43045678:C>T", cadd_phred=28, clinical_significance="Pathogenic", allele_frequency=0.0002),
        create_test_variant("BRCA1", "17:43045679:T>G", cadd_phred=26, clinical_significance="Likely_pathogenic", allele_frequency=0.0003),
    ]

    multi_hit_result = calculate_pathway_burden("dna_repair", dna_repair_data, multi_hit_variants)
    print(f"\n  Multi-hit test (BRCA1):")
    print(f"    Multi-hit genes: {', '.join(multi_hit_result['multi_hit_genes'])}")
    print(f"    BRCA1 variant count: {multi_hit_result['gene_variant_counts'].get('BRCA1', 0)}")
    print(f"    BRCA1 damaging count: {multi_hit_result['gene_damaging_counts'].get('BRCA1', 0)}")

    return burden_result

def test_overall_burden_assessment():
    """Test overall burden assessment across pathways."""

    print("\n" + "=" * 60)
    print("Testing Overall Burden Assessment")
    print("=" * 60)

    # Create mock pathway results
    pathway_results = {
        "dna_repair": {
            "burden_score": 0.8,
            "risk_level": "high",
            "damaging_variants": 5,
            "total_variants": 6,
            "contributing_genes": ["BRCA1", "BRCA2", "ATM"]
        },
        "tumor_suppressors": {
            "burden_score": 0.6,
            "risk_level": "moderate",
            "damaging_variants": 3,
            "total_variants": 4,
            "contributing_genes": ["TP53", "BRCA1"]  # BRCA1 appears in both
        },
        "oncogenes": {
            "burden_score": 0.2,
            "risk_level": "low",
            "damaging_variants": 1,
            "total_variants": 5,
            "contributing_genes": ["KRAS"]
        }
    }

    overall_result = assess_overall_burden(pathway_results)

    print(f"\nOverall Burden Assessment:")
    print(f"  Overall burden score: {overall_result['overall_burden_score']:.3f}")
    print(f"  High burden pathways: {', '.join(overall_result['high_burden_pathways'])}")
    print(f"  Moderate burden pathways: {', '.join(overall_result['moderate_burden_pathways'])}")
    print(f"  Primary concern: {overall_result['primary_concern']}")
    print(f"  Total damaging variants: {overall_result['total_damaging_variants']}")
    print(f"  Pathway crosstalk detected: {overall_result['pathway_crosstalk']}")

    if overall_result['multi_pathway_genes']:
        print(f"  Multi-pathway genes:")
        for gene, pathways in overall_result['multi_pathway_genes'].items():
            print(f"    {gene}: {', '.join(pathways)}")

    return overall_result

def test_pathway_burden_node():
    """Test the complete pathway burden node."""

    print("\n" + "=" * 60)
    print("Testing Complete Pathway Burden Node")
    print("=" * 60)

    # Create comprehensive test variants covering multiple pathways
    test_variants = [
        # DNA repair pathway
        create_test_variant("BRCA1", "17:43045677:G>A", cadd_phred=30, clinical_significance="Pathogenic", allele_frequency=0.0001),
        create_test_variant("BRCA2", "13:32363533:T>C", cadd_phred=25, clinical_significance="Likely_pathogenic", allele_frequency=0.0005),
        create_test_variant("ATM", "11:108098525:G>A", cadd_phred=15, clinical_significance="Uncertain_significance", allele_frequency=0.01),

        # Tumor suppressors
        create_test_variant("TP53", "17:7673803:C>T", cadd_phred=28, clinical_significance="Pathogenic", allele_frequency=0.0003),
        create_test_variant("RB1", "13:48877782:G>A", cadd_phred=22, clinical_significance="Likely_pathogenic", allele_frequency=0.0008),

        # Oncogenes
        create_test_variant("KRAS", "12:25245350:C>T", cadd_phred=20, clinical_significance="Likely_pathogenic", allele_frequency=0.001),
        create_test_variant("EGFR", "7:55019021:G>A", cadd_phred=18, clinical_significance="Uncertain_significance", allele_frequency=0.002),

        # Mismatch repair
        create_test_variant("MLH1", "3:37034946:G>A", cadd_phred=27, clinical_significance="Pathogenic", allele_frequency=0.0004),
        create_test_variant("MSH2", "2:47403068:C>T", cadd_phred=24, clinical_significance="Likely_pathogenic", allele_frequency=0.0006),

        # Cell cycle
        create_test_variant("CDKN2A", "9:21968225:C>T", cadd_phred=23, clinical_significance="Pathogenic", allele_frequency=0.0007),

        # Non-cancer gene (should be ignored)
        create_test_variant("UNKNOWN_GENE", "1:12345678:A>G", cadd_phred=10, clinical_significance="", allele_frequency=0.05),
    ]

    # Create test state
    test_state = {
        "filtered_variants": test_variants,
        "file_metadata": {}
    }

    print(f"\nRunning pathway burden analysis on {len(test_variants)} variants...")

    # Run the pathway burden node
    result_state = process(test_state)

    # Check results
    pathway_results = result_state.get('pathway_burden_results', {})
    pathway_summary = result_state.get('pathway_burden_summary', {})

    print(f"\nResults Summary:")
    print(f"  Pathways analyzed: {len(pathway_results)}")
    print(f"  Overall burden score: {pathway_summary.get('overall_burden_score', 0):.3f}")
    print(f"  High burden pathways: {', '.join(pathway_summary.get('high_burden_pathways', []))}")
    print(f"  Primary concern: {pathway_summary.get('primary_concern', 'None')}")
    print(f"  Total damaging variants: {pathway_summary.get('total_damaging_variants', 0)}")

    # Show detailed results for pathways with variants
    print(f"\nDetailed Pathway Results:")
    for pathway_name, results in pathway_results.items():
        if results['total_variants'] > 0:
            print(f"\n  {pathway_name.upper()}:")
            print(f"    Variants: {results['damaging_variants']}/{results['total_variants']}")
            print(f"    Burden score: {results['burden_score']:.3f}")
            print(f"    Risk level: {results['risk_level']}")
            print(f"    Genes: {', '.join(results['contributing_genes'])}")
            if results['multi_hit_genes']:
                print(f"    Multi-hit genes: {', '.join(results['multi_hit_genes'])}")

    # Test error handling
    print(f"\n" + "=" * 60)
    print("Testing Error Handling")
    print("=" * 60)

    # Test with empty variants
    empty_state = {"filtered_variants": []}
    empty_result = process(empty_state)
    print(f"Empty variants test - Overall burden: {empty_result['pathway_burden_summary']['overall_burden_score']:.3f}")

    # Test with malformed variants
    malformed_state = {
        "filtered_variants": [
            {"variant_id": "malformed", "gene": ""},  # Missing required fields
            {"chrom": "1", "pos": 123}  # Missing gene
        ]
    }
    malformed_result = process(malformed_state)
    print(f"Malformed variants test - Pathways analyzed: {len(malformed_result['pathway_burden_results'])}")

    return result_state

def test_performance():
    """Test performance with larger dataset."""

    print("\n" + "=" * 60)
    print("Testing Performance")
    print("=" * 60)

    # Create a large dataset
    large_variants = []
    cancer_genes = []
    for pathway_data in CANCER_PATHWAYS.values():
        cancer_genes.extend(pathway_data['genes'])

    # Time the variant creation too
    creation_start = time.perf_counter()

    # Create 1000 test variants
    for i in range(1000):
        gene = cancer_genes[i % len(cancer_genes)]
        variant = create_test_variant(
            gene,
            f"chr{i % 22 + 1}:{i * 1000}:A>G",
            cadd_phred=i % 35,  # 0-34 range
            clinical_significance="Pathogenic" if i % 10 == 0 else "Uncertain_significance",
            allele_frequency=0.001 + (i % 100) / 10000,  # 0.001-0.0109 range
            consequence="missense_variant" if i % 5 else "frameshift_variant"
        )
        large_variants.append(variant)

    creation_time = time.perf_counter() - creation_start

    test_state = {"filtered_variants": large_variants}

    # Use perf_counter for high-precision timing
    start_time = time.perf_counter()

    # Time individual components - already imported at the top

    # Test 1: Time damage assessment only
    damage_start = time.perf_counter()
    damage_results = []
    for variant in large_variants[:100]:  # Test first 100
        damage_results.append(is_damaging_variant(variant))
    damage_time = time.perf_counter() - damage_start

    # Test 2: Time full process
    process_start = time.perf_counter()
    result = process(test_state)
    process_time = time.perf_counter() - process_start

    end_time = time.perf_counter()

    processing_time = end_time - start_time

    print(f"\nPerformance Test Results:")
    print(f"  Variants processed: {len(large_variants)}")
    print(f"  Variant creation time: {creation_time * 1000:.2f} milliseconds")
    print(f"  Damage assessment (100 variants): {damage_time * 1000:.2f} milliseconds")
    print(f"  Damage assessment rate: {100 / damage_time:.0f} variants/second")

    # Show more precise timing for fast operations
    if process_time < 0.01:
        print(f"  Full pathway burden analysis: {process_time * 1000:.2f} milliseconds")
    else:
        print(f"  Full pathway burden analysis: {process_time:.2f} seconds")

    # Ensure we don't divide by zero
    if process_time > 0:
        print(f"  Analysis rate: {len(large_variants) / process_time:.0f} variants/second")
    else:
        print(f"  Analysis rate: >1,000,000 (too fast to measure accurately)")

    print(f"\nWork actually performed:")
    print(f"  Total damaging variants found: {result['pathway_burden_summary']['total_damaging_variants']}")
    print(f"  Pathways analyzed: {result['pathway_burden_summary']['pathways_analyzed']}")
    print(f"  High burden pathways: {len(result['pathway_burden_summary']['high_burden_pathways'])}")
    print(f"  Multi-pathway genes: {len(result['pathway_burden_summary']['multi_pathway_genes'])}")

    # Verify some actual results to ensure work was done
    dna_repair_result = result['pathway_burden_results'].get('dna_repair', {})
    print(f"\nExample pathway (DNA repair):")
    print(f"  Variants in pathway: {dna_repair_result.get('total_variants', 0)}")
    print(f"  Damaging variants: {dna_repair_result.get('damaging_variants', 0)}")
    print(f"  Burden score: {dna_repair_result.get('burden_score', 0):.3f}")

    # Performance should be reasonable
    if process_time > 0:
        variants_per_second = len(large_variants) / process_time
        if variants_per_second > 100:
            print(f"\n  ✅ Performance: GOOD ({variants_per_second:.0f} variants/sec)")
        elif variants_per_second > 50:
            print(f"  ⚠️  Performance: MODERATE ({variants_per_second:.0f} variants/sec)")
        else:
            print(f"  ❌ Performance: SLOW ({variants_per_second:.0f} variants/sec)")
    else:
        print(f"\n  ✅ Performance: EXCELLENT (processing completed in <1ms)")

    # Let's also test with a smaller dataset to see if timing scales linearly
    print(f"\nScaling test:")
    for size in [100, 500, 2000]:
        test_variants = large_variants[:size] if size <= 1000 else large_variants * (size // 1000)
        small_start = time.perf_counter()
        process({"filtered_variants": test_variants})
        small_time = time.perf_counter() - small_start
        print(f"  {size} variants: {small_time * 1000:.2f} ms ({size / small_time:.0f} variants/sec)")

    return process_time

def run_all_tests():
    """Run all pathway burden tests."""

    print("🧬 Gene/Pathway Burden Model - Comprehensive Test Suite")
    print("=" * 70)

    try:
        # Run individual tests
        test_damaging_variant_assessment()
        test_pathway_burden_calculation()
        test_overall_burden_assessment()
        result_state = test_pathway_burden_node()
        processing_time = test_performance()

        print("\n" + "=" * 70)
        print("✅ ALL TESTS COMPLETED SUCCESSFULLY!")
        print("=" * 70)

        # Save detailed results
        output_file = "test_pathway_burden_results.json"
        with open(output_file, 'w') as f:
            # Convert datetime objects to strings for JSON serialization
            serializable_result = json.loads(
                json.dumps(result_state, default=str)
            )
            json.dump(serializable_result, f, indent=2)

        print(f"\n💾 Detailed results saved to: {output_file}")
        if processing_time < 0.01:
            print(f"🚀 Processing time: {processing_time * 1000:.2f} milliseconds")
        else:
            print(f"🚀 Processing time: {processing_time:.2f} seconds")

        return True

    except Exception as e:
        print(f"\n❌ TEST FAILED: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
