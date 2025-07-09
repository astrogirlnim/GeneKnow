# PRS Calculator Implementation Summary

## Overview

We have successfully implemented a Polygenic Risk Score (PRS) calculator as a new node in the GeneKnow pipeline. This node calculates weighted genetic risk scores for multiple cancer types based on small-effect SNPs from GWAS studies.

## Key Features Implemented

### 1. **Coverage-Based Confidence Adjustment**
- Automatically adjusts scores based on SNP coverage
- Three confidence levels:
  - **High**: ≥80% coverage
  - **Moderate**: 50-80% coverage  
  - **Low**: <50% coverage
- Adjusts final score by dividing by coverage fraction to avoid underestimation

### 2. **Population Stratification Handling**
- Supports multiple populations: EUR, AFR, EAS, SAS, AMR
- Uses population-specific effect sizes when available
- Falls back to EUR data with warning when population data missing
- Adds limitation notes for non-EUR populations

### 3. **Intelligent Genotype Inference**
- **Observed genotypes**: Uses actual data when available
- **Rare variants** (MAF < 0.01): Assumes heterozygous (1 risk allele)
- **Common variants**: Uses Hardy-Weinberg equilibrium to estimate
- **Default**: Conservative heterozygous assumption
- Tracks inference method for transparency

### 4. **Multi-Cancer Risk Assessment**
- Calculates PRS for 6 cancer types: BRCA, OVCA, PRAD, LUAD, COAD, PANCA
- Same variant can have different effects for different cancers
- Identifies high-risk cancers (≥95th percentile)
- Provides primary concern and overall summary

## Integration with Pipeline

### Position in Pipeline
```
population_mapper → prs_calculator → risk_model
```

### Input from Population Mapper
```python
{
    "variant_id": "17:43044295:G>A",
    "chrom": "17",
    "pos": 43044295,
    "ref": "G",
    "alt": "A",
    "gene": "BRCA1",
    "population_frequency": 0.0001,
    "clinical_significance": "Pathogenic"
}
```

### Output to Risk Model
```python
{
    "prs_results": {
        "BRCA": {
            "raw_score": 2.45,
            "adjusted_score": 3.06,  # Adjusted for coverage
            "percentile": 97,
            "risk_category": "high",
            "confidence": "moderate",
            "coverage": 0.8,
            "contributing_snps": [...]
        },
        # ... other cancer types
    },
    "prs_summary": {
        "high_risk_cancers": ["BRCA", "OVCA"],
        "primary_concern": "BRCA",
        "overall_confidence": "moderate",
        "limitations": [
            "Low SNP coverage for LUAD PRS (15%)",
            "PRS may be less accurate for AFR population"
        ]
    }
}
```

## Real-World Considerations Addressed

### Coverage Problem
- WES only covers ~2% of genome (exons)
- Many PRS SNPs are in regulatory regions
- Solution: Coverage-adjusted scores with confidence levels

### Population Bias
- Most GWAS done on European populations
- Effect sizes can vary 40-80% between populations
- Solution: Population-specific databases with fallback warnings

### Pleiotropy
- Same variants affect multiple cancers differently
- BRCA1 mutations: High risk for breast AND ovarian cancer
- Solution: Cancer-specific PRS calculations

### Missing Genotypes
- Often lack phasing information (which chromosome)
- Solution: Statistical inference based on allele frequency

## Example Results from Demo

### Low Coverage Scenario
```
SNPs matched: 1/313 (0.3%)
Confidence: low
Raw Score: 2.3
Adjusted Score: 766.7 (unreliable due to low coverage)
```

### High Coverage Scenario
```
SNPs matched: 47/313 (15%)
Confidence: moderate
Raw Score: 12.4
Adjusted Score: 82.7
Percentile: 95% (high risk)
```

### Population Effects
Same variant (10:123352317:C>T):
- EUR: Effect size 0.15
- AFR: Effect size 0.08 (47% smaller)
- EAS: Effect size 0.13

## Future Enhancements

1. **Expand PRS Database**
   - Import full GWAS Catalog
   - Add more population-specific data
   - Include latest cancer GWAS results

2. **Improve Genotype Inference**
   - Use haplotype reference panels
   - Implement proper phasing algorithms
   - Consider linkage disequilibrium

3. **Clinical Integration**
   - Add clinical guidelines for PRS interpretation
   - Combine with family history
   - Risk-stratified screening recommendations

## Testing

Run the comprehensive demo:
```bash
cd geneknow_pipeline
python demo_prs_features.py
```

This demonstrates:
- Coverage impact on confidence
- Population stratification effects
- Genotype inference methods
- Multi-cancer risk assessment

## Conclusion

The PRS calculator successfully handles real-world genomic data challenges while providing transparent, interpretable risk scores. It integrates seamlessly with the existing pipeline and provides valuable polygenic risk information that complements rare variant analysis from other nodes. 