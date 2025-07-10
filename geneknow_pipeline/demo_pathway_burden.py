#!/usr/bin/env python3
"""
Demonstration of the Pathway Burden Model.
Shows practical usage with real-world examples.
"""

import os
import sys
import json
from datetime import datetime

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from nodes.pathway_burden import process, CANCER_PATHWAYS, is_damaging_variant, calculate_pathway_burden

def demo_cancer_pathways():
    """Demonstrate the cancer pathways defined in the model."""
    
    print("🧬 Cancer Pathways in the Burden Model")
    print("=" * 70)
    
    for pathway_name, pathway_data in CANCER_PATHWAYS.items():
        print(f"\n🔹 {pathway_name.upper().replace('_', ' ')}")
        print(f"  Description: {pathway_data['description']}")
        print(f"  Cancer relevance: {pathway_data['cancer_relevance']:.0%}")
        print(f"  Weight: {pathway_data['weight']:.1f}")
        print(f"  Genes ({len(pathway_data['genes'])}): {', '.join(pathway_data['genes'][:8])}")
        if len(pathway_data['genes']) > 8:
            print(f"    ... and {len(pathway_data['genes']) - 8} more")
    
    print(f"\n📊 Total pathways: {len(CANCER_PATHWAYS)}")
    total_genes = sum(len(pathway_data['genes']) for pathway_data in CANCER_PATHWAYS.values())
    unique_genes = len(set(gene for pathway_data in CANCER_PATHWAYS.values() for gene in pathway_data['genes']))
    print(f"📊 Total gene entries: {total_genes}")
    print(f"📊 Unique genes: {unique_genes}")
    print(f"📊 Average genes per pathway: {total_genes / len(CANCER_PATHWAYS):.1f}")

def demo_variant_assessment():
    """Demonstrate variant damage assessment."""
    
    print("\n\n🔍 Variant Damage Assessment Examples")
    print("=" * 70)
    
    # Example variants with different damage profiles
    variants = [
        {
            "name": "High-impact BRCA1 variant",
            "data": {
                "variant_id": "17:43045677:G>A",
                "gene": "BRCA1",
                "cadd_phred": 32,
                "clinical_significance": "Pathogenic",
                "allele_frequency": 0.0001,
                "consequence": "frameshift_variant"
            }
        },
        {
            "name": "Moderate-impact TP53 variant",
            "data": {
                "variant_id": "17:7673803:C>T",
                "gene": "TP53",
                "cadd_phred": 22,
                "clinical_significance": "Likely_pathogenic",
                "allele_frequency": 0.002,
                "consequence": "missense_variant"
            }
        },
        {
            "name": "Low-impact KRAS variant",
            "data": {
                "variant_id": "12:25245350:C>T",
                "gene": "KRAS",
                "cadd_phred": 8,
                "clinical_significance": "Uncertain_significance",
                "allele_frequency": 0.05,
                "consequence": "synonymous_variant"
            }
        },
        {
            "name": "Benign common variant",
            "data": {
                "variant_id": "1:12345678:A>G",
                "gene": "RANDOM_GENE",
                "cadd_phred": 3,
                "clinical_significance": "Benign",
                "allele_frequency": 0.15,
                "consequence": "synonymous_variant"
            }
        }
    ]
    
    for variant in variants:
        print(f"\n🔸 {variant['name']}")
        print(f"  Variant: {variant['data']['variant_id']} ({variant['data']['gene']})")
        print(f"  CADD score: {variant['data']['cadd_phred']}")
        print(f"  Clinical significance: {variant['data']['clinical_significance']}")
        print(f"  Allele frequency: {variant['data']['allele_frequency']:.6f}")
        print(f"  Consequence: {variant['data']['consequence']}")
        
        # Assess damage
        damage_result = is_damaging_variant(variant['data'])
        print(f"  → Damage score: {damage_result['damage_score']:.3f}")
        print(f"  → Is damaging: {damage_result['is_damaging']}")
        if damage_result['damage_reasons']:
            print(f"  → Reasons: {', '.join(damage_result['damage_reasons'])}")

def demo_pathway_analysis():
    """Demonstrate pathway-specific burden analysis."""
    
    print("\n\n🎯 Pathway-Specific Burden Analysis")
    print("=" * 70)
    
    # Create a realistic set of variants
    demo_variants = [
        # DNA repair pathway - high burden
        {"variant_id": "17:43045677:G>A", "gene": "BRCA1", "cadd_phred": 32, "clinical_significance": "Pathogenic", "allele_frequency": 0.0001, "consequence": "frameshift_variant"},
        {"variant_id": "17:43045678:C>T", "gene": "BRCA1", "cadd_phred": 28, "clinical_significance": "Pathogenic", "allele_frequency": 0.0002, "consequence": "stop_gained"},
        {"variant_id": "13:32363533:T>C", "gene": "BRCA2", "cadd_phred": 26, "clinical_significance": "Likely_pathogenic", "allele_frequency": 0.0005, "consequence": "missense_variant"},
        {"variant_id": "11:108098525:G>A", "gene": "ATM", "cadd_phred": 24, "clinical_significance": "Pathogenic", "allele_frequency": 0.0008, "consequence": "missense_variant"},
        
        # Tumor suppressors - moderate burden
        {"variant_id": "17:7673803:C>T", "gene": "TP53", "cadd_phred": 22, "clinical_significance": "Likely_pathogenic", "allele_frequency": 0.002, "consequence": "missense_variant"},
        {"variant_id": "13:48877782:G>A", "gene": "RB1", "cadd_phred": 18, "clinical_significance": "Uncertain_significance", "allele_frequency": 0.005, "consequence": "missense_variant"},
        
        # Oncogenes - low burden
        {"variant_id": "12:25245350:C>T", "gene": "KRAS", "cadd_phred": 15, "clinical_significance": "Uncertain_significance", "allele_frequency": 0.01, "consequence": "missense_variant"},
        
        # Mismatch repair - high burden
        {"variant_id": "3:37034946:G>A", "gene": "MLH1", "cadd_phred": 29, "clinical_significance": "Pathogenic", "allele_frequency": 0.0003, "consequence": "stop_gained"},
        {"variant_id": "2:47403068:C>T", "gene": "MSH2", "cadd_phred": 25, "clinical_significance": "Likely_pathogenic", "allele_frequency": 0.0006, "consequence": "missense_variant"},
        
        # Non-cancer genes (should not contribute)
        {"variant_id": "1:12345678:A>G", "gene": "RANDOM_GENE", "cadd_phred": 10, "clinical_significance": "", "allele_frequency": 0.05, "consequence": "synonymous_variant"},
    ]
    
    # Analyze specific pathways
    pathways_to_analyze = ["dna_repair", "tumor_suppressors", "oncogenes", "mismatch_repair"]
    
    for pathway_name in pathways_to_analyze:
        pathway_data = CANCER_PATHWAYS[pathway_name]
        burden_result = calculate_pathway_burden(pathway_name, pathway_data, demo_variants)
        
        print(f"\n🔸 {pathway_name.upper().replace('_', ' ')} PATHWAY")
        print(f"  Total variants in pathway: {burden_result['total_variants']}")
        print(f"  Damaging variants: {burden_result['damaging_variants']}")
        print(f"  Burden score: {burden_result['burden_score']:.3f}")
        print(f"  Risk level: {burden_result['risk_level'].upper()}")
        
        if burden_result['contributing_genes']:
            print(f"  Contributing genes: {', '.join(burden_result['contributing_genes'])}")
        
        if burden_result['multi_hit_genes']:
            print(f"  Multi-hit genes: {', '.join(burden_result['multi_hit_genes'])}")
        
        if burden_result['top_variant']:
            top = burden_result['top_variant']
            print(f"  Top variant: {top['variant_id']} in {top['gene']} (damage: {top['damage_score']:.3f})")

def demo_full_analysis():
    """Demonstrate full pathway burden analysis."""
    
    print("\n\n🚀 Complete Pathway Burden Analysis")
    print("=" * 70)
    
    # Create a comprehensive variant set
    demo_variants = [
        # High-impact variants
        {"variant_id": "17:43045677:G>A", "gene": "BRCA1", "cadd_phred": 32, "clinical_significance": "Pathogenic", "allele_frequency": 0.0001, "consequence": "frameshift_variant"},
        {"variant_id": "17:43045678:C>T", "gene": "BRCA1", "cadd_phred": 28, "clinical_significance": "Pathogenic", "allele_frequency": 0.0002, "consequence": "stop_gained"},
        {"variant_id": "13:32363533:T>C", "gene": "BRCA2", "cadd_phred": 26, "clinical_significance": "Likely_pathogenic", "allele_frequency": 0.0005, "consequence": "missense_variant"},
        {"variant_id": "17:7673803:C>T", "gene": "TP53", "cadd_phred": 30, "clinical_significance": "Pathogenic", "allele_frequency": 0.0003, "consequence": "missense_variant"},
        {"variant_id": "3:37034946:G>A", "gene": "MLH1", "cadd_phred": 29, "clinical_significance": "Pathogenic", "allele_frequency": 0.0004, "consequence": "stop_gained"},
        
        # Moderate-impact variants
        {"variant_id": "11:108098525:G>A", "gene": "ATM", "cadd_phred": 20, "clinical_significance": "Likely_pathogenic", "allele_frequency": 0.002, "consequence": "missense_variant"},
        {"variant_id": "2:47403068:C>T", "gene": "MSH2", "cadd_phred": 22, "clinical_significance": "Likely_pathogenic", "allele_frequency": 0.0015, "consequence": "missense_variant"},
        {"variant_id": "12:25245350:C>T", "gene": "KRAS", "cadd_phred": 18, "clinical_significance": "Uncertain_significance", "allele_frequency": 0.008, "consequence": "missense_variant"},
        
        # Low-impact variants
        {"variant_id": "7:55019021:G>A", "gene": "EGFR", "cadd_phred": 12, "clinical_significance": "Uncertain_significance", "allele_frequency": 0.02, "consequence": "missense_variant"},
        {"variant_id": "9:21968225:C>T", "gene": "CDKN2A", "cadd_phred": 15, "clinical_significance": "Uncertain_significance", "allele_frequency": 0.015, "consequence": "missense_variant"},
    ]
    
    # Run full analysis
    test_state = {"filtered_variants": demo_variants}
    result = process(test_state)
    
    # Display results
    pathway_results = result.get('pathway_burden_results', {})
    pathway_summary = result.get('pathway_burden_summary', {})
    
    print(f"\n📊 OVERALL RESULTS")
    print(f"  Overall burden score: {pathway_summary.get('overall_burden_score', 0):.3f}")
    print(f"  Total damaging variants: {pathway_summary.get('total_damaging_variants', 0)}")
    print(f"  High burden pathways: {', '.join(pathway_summary.get('high_burden_pathways', []))}")
    print(f"  Moderate burden pathways: {', '.join(pathway_summary.get('moderate_burden_pathways', []))}")
    print(f"  Primary concern: {pathway_summary.get('primary_concern', 'None')}")
    print(f"  Pathway crosstalk: {pathway_summary.get('pathway_crosstalk', False)}")
    
    # Show pathway details
    print(f"\n📋 PATHWAY DETAILS")
    for pathway_name, results in pathway_results.items():
        if results['total_variants'] > 0:
            print(f"\n  {pathway_name.upper().replace('_', ' ')}:")
            print(f"    Risk level: {results['risk_level'].upper()}")
            print(f"    Burden score: {results['burden_score']:.3f}")
            print(f"    Variants: {results['damaging_variants']}/{results['total_variants']}")
            print(f"    Genes: {', '.join(results['contributing_genes'])}")
            
            if results['multi_hit_genes']:
                print(f"    Multi-hit genes: {', '.join(results['multi_hit_genes'])}")
            
            if results['top_variant']:
                top = results['top_variant']
                print(f"    Top variant: {top['variant_id']} in {top['gene']}")
    
    # Show multi-pathway genes
    if pathway_summary.get('multi_pathway_genes'):
        print(f"\n🔗 MULTI-PATHWAY GENES (Crosstalk)")
        for gene, pathways in pathway_summary['multi_pathway_genes'].items():
            print(f"  {gene}: {', '.join(pathways)}")
    
    # Save results
    output_file = "demo_pathway_burden_results.json"
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    
    print(f"\n💾 Detailed results saved to: {output_file}")

def demo_clinical_interpretation():
    """Demonstrate clinical interpretation of results."""
    
    print("\n\n🏥 Clinical Interpretation Guidelines")
    print("=" * 70)
    
    print("""
🔹 BURDEN SCORE INTERPRETATION:
  • 0.0 - 0.1: Minimal burden (low concern)
  • 0.1 - 0.3: Low burden (monitor)
  • 0.3 - 0.6: Moderate burden (consider genetic counseling)
  • 0.6 - 1.0: High burden (recommend genetic counseling)

🔹 RISK LEVEL INTERPRETATION:
  • Minimal: Variants unlikely to contribute to cancer risk
  • Low: Some variants of concern, routine monitoring
  • Moderate: Multiple concerning variants, consider counseling
  • High: Significant variant burden, genetic counseling recommended

🔹 MULTI-HIT GENES:
  • Genes with multiple damaging variants
  • May indicate compound heterozygosity
  • Higher confidence in pathogenicity
  • Priority for clinical follow-up

🔹 PATHWAY CROSSTALK:
  • Genes affecting multiple pathways
  • May have broader impact on cancer risk
  • Consider cumulative effects across pathways
  • Important for comprehensive risk assessment

🔹 CLINICAL RECOMMENDATIONS:
  • High burden (>0.6): Genetic counseling, family screening
  • Moderate burden (0.3-0.6): Consider counseling, monitoring
  • Multi-hit genes: Functional studies, family history
  • Pathway crosstalk: Comprehensive risk assessment
    """)

def main():
    """Run the demonstration."""
    
    print("🧬 Pathway Burden Model - Live Demonstration")
    print("=" * 80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    try:
        # Run demonstrations
        demo_cancer_pathways()
        demo_variant_assessment()
        demo_pathway_analysis()
        demo_full_analysis()
        demo_clinical_interpretation()
        
        print("\n" + "=" * 80)
        print("🎉 DEMONSTRATION COMPLETED SUCCESSFULLY!")
        print("=" * 80)
        print("""
The Pathway Burden Model is now ready for use. Key features:

✅ 10 cancer-relevant pathways analyzed
✅ Multi-criteria variant damage assessment
✅ Pathway-specific burden scoring
✅ Multi-hit gene detection
✅ Pathway crosstalk analysis
✅ Clinical interpretation guidelines

Next steps:
1. Run the full test suite: python run_pathway_burden_tests.py
2. Integrate with your genomic data pipeline
3. Customize pathways for specific cancer types
4. Add clinical decision support rules
        """)
        
    except Exception as e:
        print(f"\n❌ DEMONSTRATION FAILED: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 