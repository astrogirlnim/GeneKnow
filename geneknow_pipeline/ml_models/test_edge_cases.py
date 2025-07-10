#!/usr/bin/env python3
"""
Test Edge Cases for the Fusion Layer

This script tests challenging scenarios to validate the fusion layer's robustness.
"""

from test_external_variant import test_variant

def main():
    print('🔬 FUSION LAYER EDGE CASE TESTING')
    print('='*60)
    
    # Edge case 1: Conflicting evidence
    print("\n🤔 CONFLICTING EVIDENCE TESTS")
    test_variant(
        prs_score=0.9,                    # Very high inherited risk
        clinvar_classification='benign',      # But variant is benign
        cadd_score=35.0,                  # High predicted impact (conflict!)
        tcga_enrichment=0.8,              # No cancer enrichment
        gene_burden_score=0.0,            # No other variants
        variant_name='High PRS + CADD vs Benign ClinVar'
    )
    
    test_variant(
        prs_score=0.1,                    # Low inherited risk
        clinvar_classification='pathogenic',  # But ClinVar says pathogenic
        cadd_score=5.0,                   # Low predicted impact
        tcga_enrichment=15.0,             # High cancer enrichment (conflict!)
        gene_burden_score=0.0,
        variant_name='Pathogenic ClinVar + High TCGA vs Low Scores'
    )
    
    # Edge case 2: Minimal/Average evidence
    print("\n📊 MINIMAL EVIDENCE TESTS")
    test_variant(
        prs_score=0.5,                    # Average everything
        clinvar_classification='not_found',   
        cadd_score=12.0,                  # Moderate impact
        tcga_enrichment=1.0,              # No enrichment
        gene_burden_score=1.0,            # Some variants
        variant_name='Average Across All Scores'
    )
    
    test_variant(
        prs_score=0.1,                    # Very low inherited risk
        clinvar_classification='benign',      
        cadd_score=0.5,                   # Very low impact
        tcga_enrichment=0.2,              # Low enrichment
        gene_burden_score=0.0,            
        variant_name='Ultra Low Risk Everything'
    )
    
    # Edge case 3: Extreme high scores
    print("\n🚨 EXTREME HIGH RISK TESTS")
    test_variant(
        prs_score=1.0,                    # Maximum inherited risk
        clinvar_classification='pathogenic',  # Pathogenic
        cadd_score=50.0,                  # Maximum CADD
        tcga_enrichment=20.0,             # Maximum enrichment
        gene_burden_score=10.0,           # Maximum burden
        variant_name='Maximum Risk Across All Factors'
    )
    
    # Edge case 4: Single strong signal
    print("\n⚡ SINGLE STRONG SIGNAL TESTS")
    test_variant(
        prs_score=0.1,                    # Low everything else
        clinvar_classification='pathogenic',  # Only ClinVar is strong
        cadd_score=2.0,                   
        tcga_enrichment=0.5,              
        gene_burden_score=0.0,            
        variant_name='Only ClinVar Pathogenic Signal'
    )
    
    test_variant(
        prs_score=0.1,                    # Low everything else
        clinvar_classification='not_found',   # Only CADD is strong
        cadd_score=45.0,                  
        tcga_enrichment=0.8,              
        gene_burden_score=0.0,            
        variant_name='Only High CADD Signal'
    )
    
    test_variant(
        prs_score=0.1,                    # Low everything else
        clinvar_classification='not_found',   # Only TCGA is strong
        cadd_score=3.0,                  
        tcga_enrichment=18.0,             
        gene_burden_score=0.0,            
        variant_name='Only High TCGA Enrichment Signal'
    )
    
    print("\n" + "="*60)
    print("🧠 KEY INSIGHTS FROM EDGE CASE TESTING:")
    print("• How does the model handle conflicting evidence?")
    print("• Does it over-rely on any single feature?")
    print("• Can it make reasonable decisions with minimal data?")
    print("• Does it behave sensibly at extreme values?")

if __name__ == "__main__":
    main() 