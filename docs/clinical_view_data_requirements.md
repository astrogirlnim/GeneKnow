# Clinical View Data Requirements

## Overview
This document lists all the backend data that needs to be provided to eliminate mock data usage in the Clinical View page.

## Current Backend Response Structure

The backend currently provides:
```json
{
  "variant_count": 123,
  "risk_scores": {
    "breast": 23.5,
    "colon": 15.2,
    "lung": 8.7,
    "prostate": 12.1,
    "blood": 5.3
  },
  "structured_json": {
    "summary": {
      "total_variants_found": 123,
      "variants_passed_qc": 110,
      "high_risk_findings": 3
    },
    "variant_details": [
      {
        "gene": "BRCA1",
        "variant": "chr17:41223094:A>G",
        "consequence": "missense_variant",
        "quality_metrics": {...},
        "tcga_matches": {...}
      }
    ]
  }
}
```

## Missing Backend Data by Subtab

### 1. Genomic Analysis Tab

#### Mutation Type Distribution
**Required Format:**
```json
"mutation_types": {
  "snv": 66,        // Single Nucleotide Variants
  "indel": 39,      // Insertions/Deletions
  "cnv": 19,        // Copy Number Variants
  "structural": 6   // Structural Variants
}
```

#### Mutational Signatures
**Required Format:**
```json
"mutation_signatures": [
  {
    "signature": "SBS1",
    "name": "Spontaneous deamination",
    "contribution": 0.35,
    "description": "C>T transitions at CpG sites",
    "etiology": "Aging-related mutations"
  },
  {
    "signature": "SBS3",
    "name": "BRCA1/BRCA2 deficiency",
    "contribution": 0.28,
    "description": "Homologous recombination deficiency",
    "etiology": "Inherited DNA repair defects"
  }
]
```

### 2. Variant Heatmap Tab

#### Variant Transformations
**Required Format:**
```json
"variant_details": [
  {
    "gene": "BRCA1",
    "variant_id": "chr17:41223094:A>G",
    "consequence": "missense_variant",
    "mutation_type": "snv",
    "protein_change": "p.Ile1756Val",
    "functional_impact": "Loss of function",
    "transformation": {
      "original": "ATT",
      "mutated": "GTT",
      "amino_acid_change": "Ile→Val",
      "effect": "Structural disruption of DNA binding domain"
    }
  }
]
```

#### Structural Variants
**Required Format:**
```json
"structural_variants": [
  {
    "type": "deletion",
    "chromosome": "chr17",
    "start": 41196312,
    "end": 41277500,
    "size": 81188,
    "genes_affected": ["BRCA1"],
    "clinical_significance": "pathogenic",
    "functional_impact": "Complete gene deletion"
  }
]
```

#### Copy Number Variants
**Required Format:**
```json
"copy_number_variants": [
  {
    "gene": "ERBB2",
    "chromosome": "chr17",
    "copy_number": 6,
    "normal_copy_number": 2,
    "fold_change": 3.0,
    "clinical_significance": "pathogenic",
    "cancer_relevance": "HER2 amplification in breast cancer"
  }
]
```

### 3. Pathway Analysis Tab

#### Pathway Disruption Analysis
**Required Format:**
```json
"pathway_analysis": {
  "disrupted_pathways": [
    {
      "name": "DNA Repair",
      "pathway_id": "DNA_REPAIR",
      "significance": 95,
      "affected_genes": ["BRCA1", "TP53"],
      "mutations": [
        {
          "gene": "BRCA1",
          "type": "frameshift",
          "effect": "Complete loss of homologous recombination"
        }
      ]
    }
  ],
  "cancer_pathway_associations": {
    "breast": ["DNA_REPAIR", "CELL_CYCLE"],
    "colon": ["WNT_SIGNALING", "APOPTOSIS"]
  }
}
```

#### Gene Interaction Data
**Required Format:**
```json
"gene_interactions": [
  {
    "gene1": "BRCA1",
    "gene2": "TP53",
    "interaction_type": "protein-protein",
    "pathway": "DNA_REPAIR"
  }
]
```

### 4. Clinical Report Tab

#### Survival Analysis Data
**Required Format:**
```json
"survival_analysis": {
  "patient_profile": {
    "risk_category": "high",
    "estimated_survival": [
      {"age": 30, "probability": 0.98},
      {"age": 40, "probability": 0.92},
      {"age": 50, "probability": 0.85},
      {"age": 60, "probability": 0.75},
      {"age": 70, "probability": 0.60},
      {"age": 80, "probability": 0.40}
    ]
  },
  "population_average": [
    {"age": 30, "probability": 0.99},
    {"age": 40, "probability": 0.97},
    {"age": 50, "probability": 0.94},
    {"age": 60, "probability": 0.88},
    {"age": 70, "probability": 0.75},
    {"age": 80, "probability": 0.50}
  ]
}
```

#### Clinical Recommendations
**Required Format:**
```json
"clinical_recommendations": [
  {
    "cancer_type": "breast",
    "risk_level": "high",
    "risk_percentage": 72.5,
    "recommendation": "Enhanced breast cancer screening recommended",
    "screening_protocol": {
      "test": "MRI + Mammography",
      "frequency": "Annual",
      "start_age": "25 or 10 years before youngest diagnosis in family"
    },
    "prevention_options": [
      "Risk-reducing mastectomy",
      "Chemoprevention with tamoxifen",
      "Lifestyle modifications"
    ]
  }
]
```

### 5. Family Analysis Tab

#### Family Risk Assessment
**Required Format:**
```json
"family_analysis": {
  "hereditary_risk": {
    "pattern": "autosomal_dominant",
    "penetrance": 0.85,
    "inheritance_probability": {
      "children": 0.50,
      "siblings": 0.50,
      "parents": 0.50
    }
  },
  "cascade_testing_recommendations": [
    {
      "relationship": "first_degree_relatives",
      "testing_priority": "high",
      "estimated_carrier_probability": 0.50
    }
  ],
  "family_history_impact": {
    "maternal_line": "high_risk",
    "paternal_line": "average_risk"
  }
}
```

## Backend Implementation Priority

1. **High Priority** (Needed for basic functionality):
   - Mutation type distribution
   - Variant transformations
   - Clinical recommendations

2. **Medium Priority** (Enhances analysis):
   - Mutational signatures
   - Pathway disruption analysis
   - Copy number variants

3. **Low Priority** (Advanced features):
   - Survival analysis curves
   - Family risk assessment
   - Gene interaction networks

## API Response Update

The enhanced API response should include all the above data structures in the `structured_json` field or as separate top-level fields in the response. 