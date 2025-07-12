#!/usr/bin/env python3
"""
Test External Variants with the Fusion Layer

Use this script to test real-world variants that weren't in our training data.
"""

from fusion_layer import FusionLayer, StaticModelInputs


def test_variant(
    prs_score: float,
    clinvar_classification: str,
    cadd_score: float,
    tcga_enrichment: float,
    gene_burden_score: float,
    variant_name: str = "Unknown",
):
    """
    Test a single variant with the fusion layer.

    Args:
        prs_score: Polygenic risk score (0.0-1.0)
        clinvar_classification: 'pathogenic', 'benign', 'uncertain', 'not_found'
        cadd_score: CADD deleteriousness score (0.0-50.0)
        tcga_enrichment: Fold enrichment in tumors (0.1-20.0)
        gene_burden_score: Number of damaging variants in key genes (0.0-10.0)
        variant_name: Descriptive name for the variant
    """
    # Load trained model
    fusion = FusionLayer()
    fusion.load_model("best_fusion_model.pkl")

    # Create input
    test_input = StaticModelInputs(
        prs_score=prs_score,
        clinvar_classification=clinvar_classification,
        cadd_score=cadd_score,
        tcga_enrichment=tcga_enrichment,
        gene_burden_score=gene_burden_score,
    )

    # Make prediction
    prediction = fusion.predict(test_input)

    # Print results
    print(f"\n=== {variant_name} ===")
    print(f"Risk Score: {prediction.risk_score:.3f}")
    print(f"Risk Category: {prediction.risk_category}")
    print(f"Confidence: {prediction.confidence:.3f}")

    # Show top contributing factors
    print("Top Contributing Factors:")
    sorted_factors = sorted(prediction.contributing_factors.items(), key=lambda x: x[1], reverse=True)
    for factor, contribution in sorted_factors[:5]:
        print(f"  {factor}: {contribution:.3f}")

    return prediction


def main():
    """Test some example external variants."""
    print("🧬 Testing External Variants with Fusion Layer")
    print("=" * 60)

    # Example 1: Real BRCA2 pathogenic variant (c.5946delT)
    test_variant(
        prs_score=0.5,  # Average population risk
        clinvar_classification="pathogenic",  # Known pathogenic in ClinVar
        cadd_score=36.0,  # Very high CADD score
        tcga_enrichment=12.0,  # Strong enrichment in breast cancers
        gene_burden_score=1.0,  # Single high-impact variant
        variant_name="BRCA2 c.5946delT (Known Pathogenic)",
    )

    # Example 2: TP53 missense variant of uncertain significance
    test_variant(
        prs_score=0.4,  # Below average inherited risk
        clinvar_classification="uncertain",  # VUS in ClinVar
        cadd_score=24.5,  # High predicted impact
        tcga_enrichment=8.5,  # Strong cancer association
        gene_burden_score=0.0,  # Only variant in pathway
        variant_name="TP53 R248W (VUS with High Evidence)",
    )

    # Example 3: Rare variant in cancer gene, not in ClinVar
    test_variant(
        prs_score=0.8,  # High polygenic risk
        clinvar_classification="not_found",  # Novel variant
        cadd_score=31.0,  # Very high predicted impact
        tcga_enrichment=6.8,  # Moderate cancer enrichment
        gene_burden_score=2.0,  # Additional variants in pathway
        variant_name="ATM Novel Missense (Not in ClinVar)",
    )

    # Example 4: Common population variant
    test_variant(
        prs_score=0.3,  # Low inherited risk
        clinvar_classification="benign",  # Known benign
        cadd_score=1.2,  # Very low impact
        tcga_enrichment=0.9,  # No enrichment
        gene_burden_score=0.0,  # No other variants
        variant_name="Common Benign SNP",
    )

    print("\n" + "=" * 60)
    print("💡 How to interpret results:")
    print("  • Risk Score 0.0-0.25: Low risk")
    print("  • Risk Score 0.25-0.5: Moderate risk")
    print("  • Risk Score 0.5-0.75: High risk")
    print("  • Risk Score 0.75-1.0: Very high risk")
    print("\n✨ The model combines evidence from multiple sources!")


if __name__ == "__main__":
    main()

    print("\n🔬 To test your own variant, use:")
    print("test_variant(")
    print("    prs_score=0.6,")
    print("    clinvar_classification='uncertain',")
    print("    cadd_score=25.0,")
    print("    tcga_enrichment=3.0,")
    print("    gene_burden_score=1.0,")
    print("    variant_name='Your Variant Name'")
    print(")")
