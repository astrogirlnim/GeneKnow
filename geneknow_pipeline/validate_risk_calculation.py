"""
Validate and explain risk calculations for a MAF file.
Shows which genes are contributing to each cancer risk.
"""
import sys
import gzip
import pandas as pd
from collections import defaultdict
import json

def analyze_maf_file(maf_path):
    """Analyze a MAF file and show what's driving the risk scores."""
    
    print(f"\n🔍 Analyzing: {maf_path}")
    print("=" * 80)
    
    # Read MAF file
    if maf_path.endswith('.gz'):
        df = pd.read_csv(maf_path, sep='\t', compression='gzip', low_memory=False)
    else:
        df = pd.read_csv(maf_path, sep='\t', low_memory=False)
    
    print(f"\n📊 File Statistics:")
    print(f"   Total variants: {len(df)}")
    print(f"   Unique genes: {df['Hugo_Symbol'].nunique()}")
    print(f"   Unique samples: {df['Tumor_Sample_Barcode'].nunique()}")
    
    # Load our cancer gene lists
    with open('models/model_config.json', 'r') as f:
        config = json.load(f)
    
    cancer_genes = config['cancer_genes']
    
    # Expanded gene risk weights (from risk_model.py)
    gene_risks = {
        "BRCA1": {"breast": 60, "prostate": 20},
        "BRCA2": {"breast": 45, "prostate": 25},
        "TP53": {"breast": 30, "colon": 25, "lung": 30, "bone": 40, "blood": 35},
        "APC": {"colon": 50},
        "KRAS": {"colon": 30, "lung": 35},
        "JAK2": {"blood": 40},
        "FLT3": {"blood": 35},
        "EGFR": {"lung": 40},
        "ALK": {"lung": 35},
        "NPM1": {"blood": 35},
        "DNMT3A": {"blood": 30},
        "IDH1": {"blood": 30, "bone": 35},
        "IDH2": {"blood": 30, "bone": 35},
        "RUNX1": {"blood": 25},
        "TET2": {"blood": 25},
        "ASXL1": {"blood": 25},
        "CEBPA": {"blood": 20},
        "KIT": {"blood": 20},
        "NRAS": {"blood": 20},
        "CBL": {"blood": 15},
        "EZH2": {"blood": 15},
        "STAT5B": {"blood": 20},
        "RB1": {"bone": 50},
        "EWSR1": {"bone": 45},
        "FLI1": {"bone": 40},
        "CDKN2A": {"bone": 35},
        "MDM2": {"bone": 30},
        "PIK3CA": {"breast": 25},
        "PALB2": {"breast": 35},
        "ATM": {"breast": 30},
        "CHEK2": {"breast": 25},
        "SMAD4": {"colon": 30},
        "BRAF": {"colon": 35},
        "MSH2": {"colon": 40},
        "MLH1": {"colon": 40},
        "ROS1": {"lung": 35},
        "MET": {"lung": 30},
        "AR": {"prostate": 40},
        "PTEN": {"prostate": 35},
        "FOXA1": {"prostate": 25},
        "SPOP": {"prostate": 20}
    }
    
    # Find which cancer genes are in this file
    variant_genes = set(df['Hugo_Symbol'].unique())
    
    print(f"\n🧬 Cancer Gene Analysis:")
    print(f"   Total genes with variants: {len(variant_genes)}")
    
    # Calculate risk for each cancer type
    cancer_results = {}
    
    for cancer_type, model_genes in cancer_genes.items():
        affected_genes = []
        risk_contribution = []
        base_risk = {"breast": 3.0, "colon": 1.7, "lung": 2.4, 
                    "prostate": 2.9, "blood": 1.7, "bone": 1.9}.get(cancer_type, 2.0)
        
        total_risk = base_risk
        
        for gene in variant_genes:
            if gene in gene_risks and cancer_type in gene_risks[gene]:
                affected_genes.append(gene)
                contribution = gene_risks[gene][cancer_type]
                
                # Apply dampening for high variant counts
                if len(df) > 100:
                    dampening = min(100 / len(df), 1.0)
                    contribution *= dampening
                
                risk_contribution.append((gene, contribution))
                total_risk += contribution
        
        total_risk = min(total_risk, 95.0)
        cancer_results[cancer_type] = {
            "risk": total_risk,
            "affected_genes": affected_genes,
            "contributions": sorted(risk_contribution, key=lambda x: x[1], reverse=True)
        }
    
    # Sort by risk
    sorted_cancers = sorted(cancer_results.items(), key=lambda x: x[1]["risk"], reverse=True)
    
    print("\n📈 Risk Calculation Breakdown:")
    print("-" * 80)
    
    for cancer, data in sorted_cancers:
        print(f"\n{cancer.upper()} CANCER: {data['risk']:.1f}%")
        base_risks = {'breast': 3.0, 'colon': 1.7, 'lung': 2.4, 'prostate': 2.9, 'blood': 1.7, 'bone': 1.9}
        print(f"   Base risk: {base_risks.get(cancer, 2.0)}%")
        print(f"   Affected genes: {len(data['affected_genes'])}")
        
        if data['affected_genes']:
            print(f"   Gene contributions (with dampening factor {min(100/len(df), 1.0):.3f}):")
            for gene, contrib in data['contributions'][:10]:  # Show top 10
                print(f"      - {gene}: +{contrib:.1f}%")
            if len(data['contributions']) > 10:
                print(f"      ... and {len(data['contributions']) - 10} more genes")
    
    # Show variant classification breakdown
    print(f"\n🔬 Variant Classifications:")
    variant_classes = df['Variant_Classification'].value_counts()
    for vclass, count in variant_classes.head(10).items():
        print(f"   {vclass}: {count}")
    
    # Check for highly mutated genes
    print(f"\n🎯 Most Mutated Genes (might indicate hypermutation):")
    gene_counts = df['Hugo_Symbol'].value_counts()
    for gene, count in gene_counts.head(10).items():
        in_cancer_genes = any(gene in cancer_genes[ct] for ct in cancer_genes)
        cancer_marker = "⚠️  CANCER GENE" if gene in gene_risks else ""
        print(f"   {gene}: {count} variants {cancer_marker}")
    
    # Validation suggestions
    print(f"\n✅ Validation Suggestions:")
    print(f"1. With {len(df)} variants, this sample may be hypermutated")
    print(f"2. Normal samples typically have <100 somatic variants")
    print(f"3. Tumor mutational burden (TMB):")
    print(f"   - Low: <5 mutations/Mb")
    print(f"   - Intermediate: 5-20 mutations/Mb")  
    print(f"   - High: >20 mutations/Mb")
    print(f"4. Current variant density suggests {'HIGH' if len(df) > 1000 else 'MODERATE'} TMB")
    
    return cancer_results


def compare_with_population():
    """Show population statistics for context."""
    print(f"\n📊 Population Cancer Risk Baselines (for comparison):")
    print("   Breast: 12.9% (lifetime risk for women)")
    print("   Prostate: 12.5% (lifetime risk for men)")
    print("   Lung: 6.3% (lifetime risk)")
    print("   Colon: 4.3% (lifetime risk)")
    print("   Blood (all types): 1.8% (lifetime risk)")
    print("   Bone: 0.1% (lifetime risk)")
    print("\n   Note: Your calculated risks should be interpreted relative to these baselines")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python validate_risk_calculation.py <path_to_maf_file>")
        sys.exit(1)
    
    maf_file = sys.argv[1]
    analyze_maf_file(maf_file)
    compare_with_population() 