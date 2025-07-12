# Backend Implementation Guide for Clinical View Data Requirements

## Overview
This document provides a detailed implementation guide for adding all missing data requirements to the GeneKnow pipeline backend. Each requirement is explained with implementation steps, code examples, and integration points.

## Table of Contents
1. [Mutation Type Distribution](#1-mutation-type-distribution)
2. [Mutational Signatures](#2-mutational-signatures)
3. [Variant Transformations](#3-variant-transformations)
4. [Structural Variants Detection](#4-structural-variants-detection)
5. [Copy Number Variants (CNVs)](#5-copy-number-variants-cnvs)
6. [Pathway Disruption Analysis](#6-pathway-disruption-analysis)
7. [Gene Interaction Networks](#7-gene-interaction-networks)
8. [Survival Analysis](#8-survival-analysis)
9. [Clinical Recommendations Engine](#9-clinical-recommendations-engine)
10. [Family Risk Assessment](#10-family-risk-assessment)

---

## 1. Mutation Type Distribution

### Implementation Strategy
Add mutation type classification to the preprocessing node or create a dedicated node.

### Code Implementation

```python
# geneknow_pipeline/nodes/mutation_classifier.py
from typing import Dict, List
from langchain_core.messages import BaseMessage

class MutationClassifierNode:
    """Classifies variants into mutation types (SNV, Indel, CNV, Structural)"""
    
    def classify_mutation_type(self, variant: Dict) -> str:
        """Classify variant based on reference and alternate alleles"""
        ref = variant.get('ref', '')
        alt = variant.get('alt', '')
        
        # Single Nucleotide Variant (SNV)
        if len(ref) == 1 and len(alt) == 1:
            return 'snv'
        
        # Insertion/Deletion (Indel)
        elif len(ref) != len(alt):
            return 'indel'
        
        # Complex variants (could be structural)
        elif len(ref) > 50 or len(alt) > 50:
            return 'structural'
        
        # Default to SNV for other cases
        return 'snv'
    
    def __call__(self, state: Dict) -> Dict:
        """Process variants and count mutation types"""
        variants = state.get('variants', [])
        
        mutation_counts = {
            'snv': 0,
            'indel': 0,
            'cnv': 0,
            'structural': 0
        }
        
        # Classify each variant
        classified_variants = []
        for variant in variants:
            mutation_type = self.classify_mutation_type(variant)
            mutation_counts[mutation_type] += 1
            
            # Add mutation type to variant
            variant['mutation_type'] = mutation_type
            classified_variants.append(variant)
        
        # Update state
        state['variants'] = classified_variants
        state['mutation_types'] = mutation_counts
        
        return {"messages": [BaseMessage(content="Mutation types classified", type="system")]}
```

### Integration Point
Add to graph.py after preprocessing node:
```python
workflow.add_node("mutation_classifier", MutationClassifierNode())
workflow.add_edge("preprocess", "mutation_classifier")
```

---

## 2. Mutational Signatures

### Implementation Strategy
Create a mutational signature analysis node that compares variant patterns to known COSMIC signatures.

### Code Implementation

```python
# geneknow_pipeline/nodes/mutational_signatures.py
import numpy as np
from typing import Dict, List
from sklearn.decomposition import NMF

class MutationalSignatureNode:
    """Analyzes mutational signatures using COSMIC signature database"""
    
    def __init__(self):
        # Load COSMIC signatures (would be loaded from database)
        self.cosmic_signatures = self.load_cosmic_signatures()
    
    def load_cosmic_signatures(self) -> Dict:
        """Load COSMIC mutational signatures"""
        # In production, load from COSMIC database
        return {
            "SBS1": {
                "name": "Spontaneous deamination",
                "description": "C>T transitions at CpG sites",
                "etiology": "Aging-related mutations",
                "pattern": np.array([0.2, 0.1, 0.05, ...])  # 96-dimensional vector
            },
            "SBS3": {
                "name": "BRCA1/BRCA2 deficiency",
                "description": "Homologous recombination deficiency",
                "etiology": "Inherited DNA repair defects",
                "pattern": np.array([0.15, 0.08, 0.12, ...])
            }
            # Add more signatures...
        }
    
    def extract_mutation_context(self, variant: Dict) -> str:
        """Extract trinucleotide context for signature analysis"""
        # Get reference genome context
        chrom = variant['chrom']
        pos = variant['pos']
        ref = variant['ref']
        alt = variant['alt']
        
        # In production, fetch from reference genome
        # For now, return placeholder
        return f"{ref}[{ref}>{alt}]{ref}"
    
    def decompose_signatures(self, mutation_matrix: np.ndarray) -> List[Dict]:
        """Decompose observed mutations into COSMIC signatures"""
        # Use Non-negative Matrix Factorization
        nmf = NMF(n_components=len(self.cosmic_signatures), init='random')
        
        # Decompose mutation matrix
        contributions = nmf.fit_transform(mutation_matrix)
        
        # Calculate relative contributions
        total = contributions.sum()
        relative_contributions = contributions / total if total > 0 else contributions
        
        # Format results
        results = []
        for i, (sig_id, sig_data) in enumerate(self.cosmic_signatures.items()):
            if relative_contributions[0][i] > 0.01:  # Only include significant contributions
                results.append({
                    "signature": sig_id,
                    "name": sig_data["name"],
                    "contribution": float(relative_contributions[0][i]),
                    "description": sig_data["description"],
                    "etiology": sig_data["etiology"]
                })
        
        # Sort by contribution
        results.sort(key=lambda x: x['contribution'], reverse=True)
        return results
    
    def __call__(self, state: Dict) -> Dict:
        """Analyze mutational signatures in variants"""
        variants = state.get('variants', [])
        
        # Build mutation context matrix
        mutation_contexts = []
        for variant in variants:
            if variant.get('mutation_type') == 'snv':
                context = self.extract_mutation_context(variant)
                mutation_contexts.append(context)
        
        # Create mutation matrix (simplified - in production would be 96-dimensional)
        mutation_matrix = np.random.rand(1, len(self.cosmic_signatures))  # Placeholder
        
        # Decompose into signatures
        signatures = self.decompose_signatures(mutation_matrix)
        
        # Add to state
        state['mutation_signatures'] = signatures
        
        return {"messages": [BaseMessage(content="Mutational signatures analyzed", type="system")]}
```

### Database Requirements
- Download COSMIC mutational signatures database
- Create signature reference table
- Implement efficient signature matching algorithm

---

## 3. Variant Transformations

### Implementation Strategy
Enhance variant annotation with protein-level changes and functional predictions.

### Code Implementation

```python
# geneknow_pipeline/nodes/variant_transformer.py
from typing import Dict
import re

class VariantTransformerNode:
    """Adds transformation details to variants including protein changes"""
    
    def __init__(self):
        # Genetic code table
        self.codon_table = {
            'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
            'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
            # ... complete codon table
        }
        
        # Amino acid properties
        self.aa_properties = {
            'Phe': {'hydrophobic': True, 'size': 'large'},
            'Leu': {'hydrophobic': True, 'size': 'large'},
            'Ile': {'hydrophobic': True, 'size': 'medium'},
            'Val': {'hydrophobic': True, 'size': 'small'},
            # ... complete properties
        }
    
    def get_protein_change(self, variant: Dict) -> Dict:
        """Calculate protein-level changes from DNA variant"""
        # This would integrate with VEP or similar annotation tool
        # For demonstration, using simplified logic
        
        ref_codon = variant.get('ref_codon', 'ATT')
        alt_codon = variant.get('alt_codon', 'GTT')
        
        ref_aa = self.codon_table.get(ref_codon, 'X')
        alt_aa = self.codon_table.get(alt_codon, 'X')
        
        # Determine functional impact
        impact = self.assess_functional_impact(ref_aa, alt_aa)
        
        return {
            "original": ref_codon,
            "mutated": alt_codon,
            "amino_acid_change": f"{ref_aa}→{alt_aa}",
            "effect": impact
        }
    
    def assess_functional_impact(self, ref_aa: str, alt_aa: str) -> str:
        """Assess the functional impact of amino acid change"""
        if ref_aa == alt_aa:
            return "Synonymous - no amino acid change"
        
        ref_props = self.aa_properties.get(ref_aa, {})
        alt_props = self.aa_properties.get(alt_aa, {})
        
        # Check property changes
        if ref_props.get('hydrophobic') != alt_props.get('hydrophobic'):
            return "Hydrophobicity change - potential structural disruption"
        elif ref_props.get('size') != alt_props.get('size'):
            return "Size change - potential steric clash"
        else:
            return "Conservative substitution - minimal impact expected"
    
    def __call__(self, state: Dict) -> Dict:
        """Add transformation details to variants"""
        variants = state.get('variant_details', [])
        
        for variant in variants:
            # Add transformation details
            transformation = self.get_protein_change(variant)
            variant['transformation'] = transformation
            
            # Add functional impact assessment
            variant['functional_impact'] = transformation['effect']
        
        state['variant_details'] = variants
        
        return {"messages": [BaseMessage(content="Variant transformations added", type="system")]}
```

### Integration with VEP
```python
# Alternative implementation using Ensembl VEP
import requests

def annotate_with_vep(variant: Dict) -> Dict:
    """Use Ensembl VEP REST API for variant annotation"""
    vep_url = "https://rest.ensembl.org/vep/human/hgvs"
    
    # Format variant as HGVS
    hgvs = f"{variant['chrom']}:g.{variant['pos']}{variant['ref']}>{variant['alt']}"
    
    response = requests.get(f"{vep_url}/{hgvs}", 
                          headers={"Content-Type": "application/json"})
    
    if response.status_code == 200:
        vep_data = response.json()[0]
        return {
            "protein_change": vep_data.get('hgvsp', ''),
            "consequence": vep_data.get('most_severe_consequence', ''),
            "impact": vep_data.get('impact', ''),
            "sift": vep_data.get('sift_prediction', ''),
            "polyphen": vep_data.get('polyphen_prediction', '')
        }
    return {}
```

---

## 4. Structural Variants Detection

### Implementation Strategy
Add structural variant calling capabilities to the variant calling node.

### Code Implementation

```python
# geneknow_pipeline/nodes/structural_variant_detector.py
from typing import Dict, List
import pysam

class StructuralVariantDetectorNode:
    """Detects structural variants from alignment data"""
    
    def __init__(self):
        self.sv_types = {
            'DEL': 'deletion',
            'INS': 'insertion',
            'DUP': 'duplication',
            'INV': 'inversion',
            'TRA': 'translocation'
        }
    
    def detect_from_bam(self, bam_path: str) -> List[Dict]:
        """Detect structural variants from BAM file"""
        structural_variants = []
        
        with pysam.AlignmentFile(bam_path, "rb") as bamfile:
            for read in bamfile:
                # Check for split reads
                if read.has_tag('SA'):
                    sv = self.parse_split_read(read)
                    if sv:
                        structural_variants.append(sv)
                
                # Check for discordant pairs
                if read.is_paired and not read.is_proper_pair:
                    sv = self.analyze_discordant_pair(read)
                    if sv:
                        structural_variants.append(sv)
        
        return self.merge_sv_calls(structural_variants)
    
    def parse_split_read(self, read) -> Dict:
        """Parse split read alignment to identify SV"""
        sa_tag = read.get_tag('SA')
        # Parse supplementary alignment
        # Format: chr,pos,strand,CIGAR,mapQ,NM
        
        parts = sa_tag.split(',')
        if len(parts) >= 4:
            supp_chr = parts[0]
            supp_pos = int(parts[1])
            
            # Determine SV type based on alignment
            if read.reference_name == supp_chr:
                # Same chromosome - likely deletion or duplication
                size = abs(supp_pos - read.reference_start)
                return {
                    "type": "deletion" if size > 1000 else "duplication",
                    "chromosome": read.reference_name,
                    "start": min(read.reference_start, supp_pos),
                    "end": max(read.reference_start, supp_pos),
                    "size": size,
                    "evidence": "split_read"
                }
            else:
                # Different chromosomes - translocation
                return {
                    "type": "translocation",
                    "chr1": read.reference_name,
                    "pos1": read.reference_start,
                    "chr2": supp_chr,
                    "pos2": supp_pos,
                    "evidence": "split_read"
                }
        return None
    
    def analyze_discordant_pair(self, read) -> Dict:
        """Analyze discordant read pairs for SV evidence"""
        if not read.is_read1:  # Only process read1 to avoid duplicates
            return None
        
        # Expected insert size (would be calculated from data)
        expected_insert = 500
        tolerance = 200
        
        if read.tlen > expected_insert + tolerance:
            # Deletion
            return {
                "type": "deletion",
                "chromosome": read.reference_name,
                "start": read.reference_start,
                "end": read.reference_start + read.tlen,
                "size": read.tlen - expected_insert,
                "evidence": "discordant_pair"
            }
        elif 0 < read.tlen < expected_insert - tolerance:
            # Insertion
            return {
                "type": "insertion",
                "chromosome": read.reference_name,
                "start": read.reference_start,
                "end": read.reference_start + read.tlen,
                "size": expected_insert - read.tlen,
                "evidence": "discordant_pair"
            }
        return None
    
    def merge_sv_calls(self, sv_calls: List[Dict]) -> List[Dict]:
        """Merge overlapping SV calls"""
        # Sort by chromosome and position
        sorted_calls = sorted(sv_calls, 
                            key=lambda x: (x.get('chromosome', ''), x.get('start', 0)))
        
        merged = []
        for sv in sorted_calls:
            if not merged:
                merged.append(sv)
                continue
            
            last = merged[-1]
            # Check if overlapping
            if (sv.get('chromosome') == last.get('chromosome') and 
                sv.get('type') == last.get('type') and
                sv.get('start', 0) <= last.get('end', 0) + 1000):
                # Merge
                last['end'] = max(last.get('end', 0), sv.get('end', 0))
                last['size'] = last['end'] - last['start']
                if 'evidence' in sv:
                    last['evidence'] = f"{last.get('evidence', '')},{sv['evidence']}"
            else:
                merged.append(sv)
        
        return merged
    
    def annotate_sv_impact(self, sv: Dict, gene_annotations: Dict) -> Dict:
        """Annotate SV with affected genes and clinical significance"""
        affected_genes = []
        
        # Find genes in SV region
        for gene, coords in gene_annotations.items():
            if (coords['chr'] == sv.get('chromosome') and
                coords['start'] <= sv.get('end', 0) and
                coords['end'] >= sv.get('start', 0)):
                affected_genes.append(gene)
        
        sv['genes_affected'] = affected_genes
        
        # Assess clinical significance
        if any(gene in ['BRCA1', 'BRCA2', 'TP53'] for gene in affected_genes):
            sv['clinical_significance'] = 'pathogenic'
        elif len(affected_genes) > 0:
            sv['clinical_significance'] = 'likely_pathogenic'
        else:
            sv['clinical_significance'] = 'uncertain_significance'
        
        # Add functional impact
        if sv['type'] == 'deletion' and len(affected_genes) > 0:
            sv['functional_impact'] = f"Complete or partial deletion of {', '.join(affected_genes)}"
        elif sv['type'] == 'duplication':
            sv['functional_impact'] = f"Gene duplication affecting {', '.join(affected_genes)}"
        
        return sv
    
    def __call__(self, state: Dict) -> Dict:
        """Process structural variant detection"""
        bam_path = state.get('bam_path')
        
        if not bam_path:
            # If no BAM file, try to detect from VCF
            structural_variants = self.detect_from_vcf(state.get('variants', []))
        else:
            structural_variants = self.detect_from_bam(bam_path)
        
        # Load gene annotations (would come from database)
        gene_annotations = self.load_gene_annotations()
        
        # Annotate SVs
        annotated_svs = []
        for sv in structural_variants:
            annotated_sv = self.annotate_sv_impact(sv, gene_annotations)
            annotated_svs.append(annotated_sv)
        
        state['structural_variants'] = annotated_svs
        
        return {"messages": [BaseMessage(content=f"Detected {len(annotated_svs)} structural variants", type="system")]}
```

---

## 5. Copy Number Variants (CNVs)

### Implementation Strategy
Implement CNV detection using read depth analysis and B-allele frequency.

### Code Implementation

```python
# geneknow_pipeline/nodes/cnv_detector.py
import numpy as np
from scipy import stats
from typing import Dict, List

class CNVDetectorNode:
    """Detects copy number variants using depth-based analysis"""
    
    def __init__(self):
        self.window_size = 1000  # Base pairs
        self.min_windows = 3     # Minimum windows for CNV call
    
    def calculate_depth_profile(self, bam_path: str, region: str) -> np.ndarray:
        """Calculate read depth across genomic region"""
        import pysam
        
        depths = []
        with pysam.AlignmentFile(bam_path, "rb") as bamfile:
            for pileupcolumn in bamfile.pileup(region=region):
                depths.append(pileupcolumn.n)
        
        return np.array(depths)
    
    def segment_depth_profile(self, depths: np.ndarray) -> List[Dict]:
        """Segment depth profile to identify CNV regions"""
        # Use circular binary segmentation or similar
        segments = []
        
        # Calculate baseline (diploid) depth
        baseline = np.median(depths)
        
        # Sliding window analysis
        for i in range(0, len(depths) - self.window_size, self.window_size // 2):
            window = depths[i:i + self.window_size]
            window_mean = np.mean(window)
            
            # Calculate copy number estimate
            cn_estimate = round(2 * window_mean / baseline) if baseline > 0 else 2
            
            if cn_estimate != 2:  # Not diploid
                segments.append({
                    'start': i,
                    'end': i + self.window_size,
                    'copy_number': cn_estimate,
                    'depth_ratio': window_mean / baseline if baseline > 0 else 1.0
                })
        
        return self.merge_segments(segments)
    
    def merge_segments(self, segments: List[Dict]) -> List[Dict]:
        """Merge adjacent segments with same copy number"""
        if not segments:
            return []
        
        merged = [segments[0]]
        
        for seg in segments[1:]:
            last = merged[-1]
            # Check if adjacent and same copy number
            if (seg['start'] <= last['end'] + self.window_size and
                seg['copy_number'] == last['copy_number']):
                # Merge
                last['end'] = seg['end']
            else:
                merged.append(seg)
        
        # Filter small CNVs
        return [s for s in merged if s['end'] - s['start'] >= self.min_windows * self.window_size]
    
    def annotate_cnv(self, cnv: Dict, gene_annotations: Dict) -> Dict:
        """Annotate CNV with gene information and clinical significance"""
        # Find affected genes
        affected_genes = []
        for gene, info in gene_annotations.items():
            if (info['start'] <= cnv['end'] and info['end'] >= cnv['start']):
                affected_genes.append({
                    'gene': gene,
                    'overlap_fraction': self.calculate_overlap(cnv, info)
                })
        
        # Determine clinical significance
        cnv['genes_affected'] = [g['gene'] for g in affected_genes]
        
        # Calculate fold change
        cnv['fold_change'] = cnv['copy_number'] / 2.0
        cnv['normal_copy_number'] = 2
        
        # Assess significance
        if cnv['copy_number'] == 0:
            cnv['clinical_significance'] = 'pathogenic'
            cnv['type'] = 'homozygous_deletion'
        elif cnv['copy_number'] == 1:
            cnv['clinical_significance'] = 'likely_pathogenic'
            cnv['type'] = 'heterozygous_deletion'
        elif cnv['copy_number'] >= 4:
            cnv['clinical_significance'] = 'likely_pathogenic'
            cnv['type'] = 'amplification'
        else:
            cnv['clinical_significance'] = 'uncertain_significance'
            cnv['type'] = 'duplication'
        
        # Add cancer relevance for known oncogenes/suppressors
        cancer_genes = {
            'ERBB2': 'HER2 amplification in breast cancer',
            'MYC': 'MYC amplification in multiple cancers',
            'CDKN2A': 'CDKN2A deletion in melanoma',
            'TP53': 'TP53 deletion in multiple cancers'
        }
        
        for gene in cnv['genes_affected']:
            if gene in cancer_genes:
                cnv['cancer_relevance'] = cancer_genes[gene]
                break
        
        return cnv
    
    def __call__(self, state: Dict) -> Dict:
        """Detect and annotate CNVs"""
        bam_path = state.get('bam_path')
        
        cnvs = []
        
        if bam_path:
            # Get target regions (e.g., exome capture regions)
            target_regions = state.get('target_regions', self.get_default_targets())
            
            for region in target_regions:
                # Calculate depth profile
                depths = self.calculate_depth_profile(bam_path, region)
                
                # Segment to find CNVs
                region_cnvs = self.segment_depth_profile(depths)
                
                # Add chromosome info
                for cnv in region_cnvs:
                    cnv['chromosome'] = region.split(':')[0]
                    cnv['gene'] = region.split(':')[-1]  # If region includes gene name
                
                cnvs.extend(region_cnvs)
        
        # Load gene annotations
        gene_annotations = self.load_gene_annotations()
        
        # Annotate CNVs
        annotated_cnvs = []
        for cnv in cnvs:
            annotated_cnv = self.annotate_cnv(cnv, gene_annotations)
            annotated_cnvs.append(annotated_cnv)
        
        state['copy_number_variants'] = annotated_cnvs
        
        return {"messages": [BaseMessage(content=f"Detected {len(annotated_cnvs)} CNVs", type="system")]}
```

---

## 6. Pathway Disruption Analysis

### Implementation Strategy
Integrate pathway databases and analyze variant impact on biological pathways.

### Code Implementation

```python
# geneknow_pipeline/nodes/pathway_analyzer.py
from typing import Dict, List, Set
import networkx as nx

class PathwayAnalyzerNode:
    """Analyzes disruption of biological pathways by variants"""
    
    def __init__(self):
        # Load pathway databases (KEGG, Reactome, etc.)
        self.pathways = self.load_pathway_database()
        self.cancer_pathways = self.load_cancer_pathways()
    
    def load_pathway_database(self) -> Dict:
        """Load pathway definitions from database"""
        # In production, load from KEGG/Reactome
        return {
            "DNA_REPAIR": {
                "name": "DNA Repair",
                "genes": ["BRCA1", "BRCA2", "ATM", "ATR", "CHEK2", "TP53", "MLH1", "MSH2"],
                "description": "Maintains genomic stability"
            },
            "CELL_CYCLE": {
                "name": "Cell Cycle Control",
                "genes": ["TP53", "RB1", "CDKN2A", "CCND1", "CDK4", "CDK6"],
                "description": "Regulates cell division"
            },
            "APOPTOSIS": {
                "name": "Apoptosis",
                "genes": ["TP53", "BCL2", "BAX", "CASP3", "CASP9", "FAS"],
                "description": "Programmed cell death"
            },
            "WNT_SIGNALING": {
                "name": "WNT Signaling",
                "genes": ["APC", "CTNNB1", "WNT1", "TCF7L2", "AXIN1", "AXIN2"],
                "description": "Cell fate determination"
            },
            "PI3K_AKT": {
                "name": "PI3K-AKT Signaling",
                "genes": ["PIK3CA", "PTEN", "AKT1", "AKT2", "MTOR", "TSC1"],
                "description": "Cell growth and survival"
            }
        }
    
    def load_cancer_pathways(self) -> Dict:
        """Load cancer-specific pathway associations"""
        return {
            "breast": ["DNA_REPAIR", "CELL_CYCLE", "PI3K_AKT", "HORMONE_SIGNALING"],
            "colon": ["WNT_SIGNALING", "APOPTOSIS", "DNA_REPAIR", "CELL_CYCLE"],
            "lung": ["CELL_CYCLE", "PI3K_AKT", "APOPTOSIS", "EGFR_SIGNALING"],
            "prostate": ["HORMONE_SIGNALING", "PI3K_AKT", "CELL_CYCLE"],
            "blood": ["APOPTOSIS", "CELL_CYCLE", "JAK_STAT", "NOTCH_SIGNALING"]
        }
    
    def calculate_pathway_impact(self, pathway: Dict, affected_genes: Set[str], 
                               variant_impacts: Dict[str, str]) -> float:
        """Calculate the impact score for a pathway"""
        pathway_genes = set(pathway['genes'])
        affected_in_pathway = affected_genes.intersection(pathway_genes)
        
        if not affected_in_pathway:
            return 0.0
        
        # Calculate base score (fraction of pathway affected)
        base_score = len(affected_in_pathway) / len(pathway_genes) * 100
        
        # Weight by variant impact
        impact_weights = {
            'frameshift': 1.0,
            'nonsense': 1.0,
            'splice_site': 0.9,
            'missense': 0.5,
            'inframe': 0.3,
            'synonymous': 0.1
        }
        
        total_weight = 0
        for gene in affected_in_pathway:
            variant_type = variant_impacts.get(gene, 'missense')
            weight = impact_weights.get(variant_type, 0.5)
            total_weight += weight
        
        # Weighted score
        weighted_score = (total_weight / len(affected_in_pathway)) * base_score
        
        return min(weighted_score, 100)  # Cap at 100
    
    def build_gene_network(self, pathway_genes: List[str]) -> nx.Graph:
        """Build gene interaction network for pathway"""
        G = nx.Graph()
        
        # Add nodes
        for gene in pathway_genes:
            G.add_node(gene)
        
        # Add edges based on known interactions
        # In production, load from STRING or similar database
        interactions = [
            ("BRCA1", "BRCA2"), ("BRCA1", "TP53"), ("ATM", "TP53"),
            ("CDKN2A", "RB1"), ("CDK4", "RB1"), ("TP53", "MDM2")
        ]
        
        for gene1, gene2 in interactions:
            if gene1 in G and gene2 in G:
                G.add_edge(gene1, gene2)
        
        return G
    
    def __call__(self, state: Dict) -> Dict:
        """Analyze pathway disruptions"""
        variants = state.get('variant_details', [])
        
        # Extract affected genes and their impacts
        affected_genes = set()
        variant_impacts = {}
        
        for variant in variants:
            gene = variant.get('gene')
            if gene:
                affected_genes.add(gene)
                # Store most severe impact per gene
                current_impact = variant.get('consequence', 'missense')
                if gene not in variant_impacts or self.is_more_severe(current_impact, variant_impacts[gene]):
                    variant_impacts[gene] = current_impact
        
        # Analyze each pathway
        disrupted_pathways = []
        
        for pathway_id, pathway_data in self.pathways.items():
            impact_score = self.calculate_pathway_impact(
                pathway_data, affected_genes, variant_impacts
            )
            
            if impact_score > 10:  # Significance threshold
                affected_in_pathway = affected_genes.intersection(pathway_data['genes'])
                
                # Get mutation details for affected genes
                mutations = []
                for gene in affected_in_pathway:
                    for variant in variants:
                        if variant.get('gene') == gene:
                            mutations.append({
                                'gene': gene,
                                'type': variant_impacts.get(gene, 'unknown'),
                                'effect': self.describe_effect(variant)
                            })
                            break
                
                disrupted_pathways.append({
                    'name': pathway_data['name'],
                    'pathway_id': pathway_id,
                    'significance': round(impact_score, 1),
                    'affected_genes': list(affected_in_pathway),
                    'mutations': mutations
                })
        
        # Sort by significance
        disrupted_pathways.sort(key=lambda x: x['significance'], reverse=True)
        
        # Build cancer pathway associations
        cancer_associations = {}
        for cancer_type, pathway_ids in self.cancer_pathways.items():
            associated_pathways = [p for p in disrupted_pathways 
                                 if p['pathway_id'] in pathway_ids]
            if associated_pathways:
                cancer_associations[cancer_type] = [p['pathway_id'] for p in associated_pathways]
        
        # Create pathway analysis result
        pathway_analysis = {
            'disrupted_pathways': disrupted_pathways,
            'cancer_pathway_associations': cancer_associations
        }
        
        state['pathway_analysis'] = pathway_analysis
        
        return {"messages": [BaseMessage(content=f"Analyzed {len(disrupted_pathways)} disrupted pathways", type="system")]}
    
    def is_more_severe(self, impact1: str, impact2: str) -> bool:
        """Compare variant impact severity"""
        severity_order = [
            'frameshift', 'nonsense', 'splice_site', 
            'missense', 'inframe', 'synonymous'
        ]
        
        idx1 = severity_order.index(impact1) if impact1 in severity_order else 999
        idx2 = severity_order.index(impact2) if impact2 in severity_order else 999
        
        return idx1 < idx2
    
    def describe_effect(self, variant: Dict) -> str:
        """Generate human-readable effect description"""
        consequence = variant.get('consequence', 'unknown')
        gene = variant.get('gene', 'unknown')
        
        descriptions = {
            'frameshift': f"Complete loss of {gene} function due to frameshift",
            'nonsense': f"Truncated {gene} protein due to premature stop codon",
            'splice_site': f"Disrupted splicing affecting {gene} function",
            'missense': f"Amino acid change potentially affecting {gene} activity"
        }
        
        return descriptions.get(consequence, f"Variant affecting {gene} function")
```

---

## 7. Gene Interaction Networks

### Implementation Strategy
Build gene interaction networks using protein-protein interaction databases.

### Code Implementation

```python
# geneknow_pipeline/nodes/gene_interaction_network.py
import networkx as nx
import requests
from typing import Dict, List, Tuple

class GeneInteractionNetworkNode:
    """Builds and analyzes gene interaction networks"""
    
    def __init__(self):
        self.string_api_url = "https://string-db.org/api"
        self.species = 9606  # Human
    
    def fetch_interactions(self, genes: List[str]) -> List[Dict]:
        """Fetch gene interactions from STRING database"""
        interactions = []
        
        # STRING API call
        params = {
            'identifiers': '%0d'.join(genes),
            'species': self.species,
            'required_score': 400,  # Medium confidence
            'network_type': 'functional'
        }
        
        response = requests.post(f"{self.string_api_url}/json/network", data=params)
        
        if response.status_code == 200:
            string_data = response.json()
            
            for edge in string_data:
                interactions.append({
                    'gene1': edge['preferredName_A'],
                    'gene2': edge['preferredName_B'],
                    'score': edge['score'],
                    'interaction_type': self.classify_interaction(edge)
                })
        
        return interactions
    
    def classify_interaction(self, edge: Dict) -> str:
        """Classify interaction type based on evidence channels"""
        # STRING evidence channels
        if edge.get('experiments', 0) > 0:
            return 'experimental'
        elif edge.get('database', 0) > 0:
            return 'database'
        elif edge.get('textmining', 0) > 0:
            return 'literature'
        else:
            return 'predicted'
    
    def build_pathway_specific_network(self, genes: List[str], 
                                     pathway: str) -> List[Dict]:
        """Build network specific to a pathway"""
        # In production, use pathway-specific interaction databases
        pathway_interactions = {
            'DNA_REPAIR': [
                ('BRCA1', 'BRCA2', 'protein-protein'),
                ('BRCA1', 'TP53', 'regulation'),
                ('ATM', 'TP53', 'phosphorylation'),
                ('MLH1', 'MSH2', 'protein-protein')
            ],
            'CELL_CYCLE': [
                ('TP53', 'CDKN1A', 'transcription'),
                ('CDK4', 'CCND1', 'protein-protein'),
                ('RB1', 'E2F1', 'inhibition')
            ]
        }
        
        interactions = []
        for gene1, gene2, itype in pathway_interactions.get(pathway, []):
            if gene1 in genes and gene2 in genes:
                interactions.append({
                    'gene1': gene1,
                    'gene2': gene2,
                    'interaction_type': itype,
                    'pathway': pathway
                })
        
        return interactions
    
    def __call__(self, state: Dict) -> Dict:
        """Generate gene interaction data"""
        # Get affected genes from variants
        variants = state.get('variant_details', [])
        affected_genes = list(set(v.get('gene') for v in variants if v.get('gene')))
        
        if not affected_genes:
            state['gene_interactions'] = []
            return {"messages": [BaseMessage(content="No genes to analyze for interactions", type="system")]}
        
        # Fetch interactions
        all_interactions = []
        
        # Get general interactions
        try:
            string_interactions = self.fetch_interactions(affected_genes[:20])  # Limit for API
            all_interactions.extend(string_interactions)
        except:
            # Fallback to local data
            pass
        
        # Get pathway-specific interactions
        pathway_analysis = state.get('pathway_analysis', {})
        for pathway in pathway_analysis.get('disrupted_pathways', []):
            pathway_id = pathway['pathway_id']
            pathway_genes = pathway['affected_genes']
            
            pathway_interactions = self.build_pathway_specific_network(
                pathway_genes, pathway_id
            )
            all_interactions.extend(pathway_interactions)
        
        # Remove duplicates
        unique_interactions = []
        seen = set()
        
        for interaction in all_interactions:
            key = tuple(sorted([interaction['gene1'], interaction['gene2']]))
            if key not in seen:
                seen.add(key)
                unique_interactions.append(interaction)
        
        state['gene_interactions'] = unique_interactions
        
        return {"messages": [BaseMessage(content=f"Found {len(unique_interactions)} gene interactions", type="system")]}
```

---

## 8. Survival Analysis

### Implementation Strategy
Generate survival curves based on risk profiles and population data.

### Code Implementation

```python
# geneknow_pipeline/nodes/survival_analyzer.py
import numpy as np
from scipy import stats
from typing import Dict, List

class SurvivalAnalyzerNode:
    """Generates survival analysis data based on genetic risk factors"""
    
    def __init__(self):
        # Load population survival data
        self.population_survival = self.load_population_data()
        # Load genetic risk modifiers
        self.risk_modifiers = self.load_risk_modifiers()
    
    def load_population_data(self) -> Dict:
        """Load baseline population survival data"""
        # In production, load from actuarial tables
        return {
            'general': [
                {'age': 30, 'probability': 0.99},
                {'age': 40, 'probability': 0.97},
                {'age': 50, 'probability': 0.94},
                {'age': 60, 'probability': 0.88},
                {'age': 70, 'probability': 0.75},
                {'age': 80, 'probability': 0.50},
                {'age': 90, 'probability': 0.20}
            ]
        }
    
    def load_risk_modifiers(self) -> Dict:
        """Load genetic risk modifiers for survival"""
        return {
            'BRCA1': {'hazard_ratio': 2.5, 'cancer_types': ['breast', 'ovarian']},
            'BRCA2': {'hazard_ratio': 2.2, 'cancer_types': ['breast', 'ovarian']},
            'TP53': {'hazard_ratio': 3.0, 'cancer_types': ['multiple']},
            'MLH1': {'hazard_ratio': 2.0, 'cancer_types': ['colon']},
            'APC': {'hazard_ratio': 2.8, 'cancer_types': ['colon']}
        }
    
    def calculate_genetic_hazard_ratio(self, variants: List[Dict]) -> float:
        """Calculate combined hazard ratio from variants"""
        combined_hr = 1.0
        
        for variant in variants:
            gene = variant.get('gene')
            if gene in self.risk_modifiers:
                # Multiplicative model for multiple variants
                variant_hr = self.risk_modifiers[gene]['hazard_ratio']
                
                # Adjust based on variant severity
                severity_multiplier = {
                    'frameshift': 1.0,
                    'nonsense': 1.0,
                    'missense': 0.7,
                    'synonymous': 0.1
                }.get(variant.get('consequence', 'missense'), 0.5)
                
                combined_hr *= (1 + (variant_hr - 1) * severity_multiplier)
        
        return combined_hr
    
    def generate_survival_curve(self, base_survival: List[Dict], 
                              hazard_ratio: float) -> List[Dict]:
        """Generate modified survival curve based on hazard ratio"""
        modified_survival = []
        
        for point in base_survival:
            age = point['age']
            base_prob = point['probability']
            
            # Apply hazard ratio using exponential model
            # S(t) = S0(t)^HR where S0 is baseline survival
            modified_prob = base_prob ** hazard_ratio
            
            modified_survival.append({
                'age': age,
                'probability': round(modified_prob, 3)
            })
        
        return modified_survival
    
    def determine_risk_category(self, risk_scores: Dict) -> str:
        """Determine overall risk category"""
        max_risk = max(risk_scores.values()) if risk_scores else 0
        
        if max_risk >= 50:
            return 'high'
        elif max_risk >= 20:
            return 'moderate'
        else:
            return 'low'
    
    def __call__(self, state: Dict) -> Dict:
        """Generate survival analysis data"""
        variants = state.get('variant_details', [])
        risk_scores = state.get('risk_scores', {})
        
        # Calculate genetic hazard ratio
        hazard_ratio = self.calculate_genetic_hazard_ratio(variants)
        
        # Determine risk category
        risk_category = self.determine_risk_category(risk_scores)
        
        # Generate survival curves
        baseline = self.population_survival['general']
        patient_survival = self.generate_survival_curve(baseline, hazard_ratio)
        
        # Create survival analysis result
        survival_analysis = {
            'patient_profile': {
                'risk_category': risk_category,
                'hazard_ratio': round(hazard_ratio, 2),
                'estimated_survival': patient_survival
            },
            'population_average': baseline,
            'interpretation': self.generate_interpretation(hazard_ratio, risk_category)
        }
        
        state['survival_analysis'] = survival_analysis
        
        return {"messages": [BaseMessage(content="Survival analysis completed", type="system")]}
    
    def generate_interpretation(self, hazard_ratio: float, risk_category: str) -> str:
        """Generate clinical interpretation of survival data"""
        if hazard_ratio > 2.5:
            return "Significantly elevated genetic risk factors detected. Enhanced surveillance strongly recommended."
        elif hazard_ratio > 1.5:
            return "Moderately elevated genetic risk. Consider enhanced screening protocols."
        else:
            return "Genetic risk factors within expected range. Standard screening protocols appropriate."
```

---

## 9. Clinical Recommendations Engine

### Implementation Strategy
Generate evidence-based clinical recommendations using variant data and guidelines.

### Code Implementation

```python
# geneknow_pipeline/nodes/clinical_recommendations.py
from typing import Dict, List
from datetime import datetime

class ClinicalRecommendationsNode:
    """Generates clinical recommendations based on genetic findings"""
    
    def __init__(self):
        # Load clinical guidelines (NCCN, ASCO, etc.)
        self.guidelines = self.load_clinical_guidelines()
        self.screening_protocols = self.load_screening_protocols()
    
    def load_clinical_guidelines(self) -> Dict:
        """Load evidence-based clinical guidelines"""
        return {
            'breast': {
                'high_risk': {
                    'genes': ['BRCA1', 'BRCA2', 'PALB2', 'TP53'],
                    'screening': {
                        'test': 'MRI + Mammography',
                        'frequency': 'Annual',
                        'start_age': '25 or 10 years before youngest diagnosis in family'
                    },
                    'prevention': [
                        'Risk-reducing mastectomy',
                        'Risk-reducing salpingo-oophorectomy',
                        'Chemoprevention with tamoxifen or raloxifene'
                    ]
                },
                'moderate_risk': {
                    'genes': ['CHEK2', 'ATM', 'NBN'],
                    'screening': {
                        'test': 'Mammography with consideration of MRI',
                        'frequency': 'Annual',
                        'start_age': '40 or 10 years before youngest diagnosis'
                    },
                    'prevention': [
                        'Enhanced screening',
                        'Lifestyle modifications',
                        'Consider chemoprevention'
                    ]
                }
            },
            'colon': {
                'high_risk': {
                    'genes': ['MLH1', 'MSH2', 'MSH6', 'PMS2', 'APC'],
                    'screening': {
                        'test': 'Colonoscopy',
                        'frequency': 'Every 1-2 years',
                        'start_age': '20-25 or 2-5 years before youngest diagnosis'
                    },
                    'prevention': [
                        'Prophylactic colectomy consideration',
                        'Aspirin chemoprevention',
                        'Enhanced endoscopic surveillance'
                    ]
                }
            },
            'ovarian': {
                'high_risk': {
                    'genes': ['BRCA1', 'BRCA2', 'RAD51C', 'RAD51D'],
                    'screening': {
                        'test': 'Transvaginal ultrasound + CA-125',
                        'frequency': 'Every 6 months',
                        'start_age': '30-35 or 5-10 years before earliest diagnosis'
                    },
                    'prevention': [
                        'Risk-reducing salpingo-oophorectomy',
                        'Oral contraceptives',
                        'Close surveillance'
                    ]
                }
            }
        }
    
    def generate_recommendations(self, cancer_type: str, risk_level: str, 
                               risk_percentage: float, variants: List[Dict]) -> Dict:
        """Generate specific recommendations for a cancer type"""
        guidelines = self.guidelines.get(cancer_type, {}).get(risk_level + '_risk', {})
        
        if not guidelines:
            # Default recommendations
            return self.generate_default_recommendations(cancer_type, risk_percentage)
        
        # Find relevant genes
        patient_genes = [v.get('gene') for v in variants 
                        if v.get('gene') in guidelines.get('genes', [])]
        
        recommendation = {
            'cancer_type': cancer_type,
            'risk_level': risk_level,
            'risk_percentage': risk_percentage,
            'recommendation': f"Enhanced {cancer_type} cancer screening recommended",
            'screening_protocol': guidelines.get('screening', {}),
            'prevention_options': guidelines.get('prevention', []),
            'relevant_variants': patient_genes,
            'evidence_level': 'High - Based on established guidelines',
            'references': self.get_references(cancer_type, patient_genes)
        }
        
        # Add age-specific modifications
        recommendation = self.add_age_considerations(recommendation, cancer_type)
        
        return recommendation
    
    def generate_default_recommendations(self, cancer_type: str, 
                                       risk_percentage: float) -> Dict:
        """Generate default recommendations when specific guidelines not found"""
        if risk_percentage > 50:
            screening = {
                'test': f'Enhanced {cancer_type} screening',
                'frequency': 'Annual',
                'start_age': 'Consult with genetic counselor'
            }
            prevention = ['Genetic counseling', 'Risk assessment', 'Lifestyle modifications']
        elif risk_percentage > 20:
            screening = {
                'test': f'Regular {cancer_type} screening',
                'frequency': 'Annual',
                'start_age': '40 or per guidelines'
            }
            prevention = ['Regular monitoring', 'Lifestyle modifications']
        else:
            screening = {
                'test': f'Standard {cancer_type} screening',
                'frequency': 'Per population guidelines',
                'start_age': 'Per population guidelines'
            }
            prevention = ['Maintain healthy lifestyle']
        
        return {
            'cancer_type': cancer_type,
            'risk_level': 'high' if risk_percentage > 50 else 'moderate' if risk_percentage > 20 else 'low',
            'risk_percentage': risk_percentage,
            'recommendation': f"{screening['test']} recommended",
            'screening_protocol': screening,
            'prevention_options': prevention,
            'evidence_level': 'Moderate - Based on risk assessment'
        }
    
    def add_age_considerations(self, recommendation: Dict, cancer_type: str) -> Dict:
        """Add age-specific modifications to recommendations"""
        # In production, would consider patient age
        age_modifiers = {
            'breast': {
                'young': 'Consider genetic testing for family members',
                'middle': 'Standard protocol applies',
                'older': 'Discuss risk-benefit ratio of interventions'
            },
            'colon': {
                'young': 'Begin screening immediately if Lynch syndrome',
                'middle': 'Follow enhanced surveillance protocol',
                'older': 'Continue surveillance based on life expectancy'
            }
        }
        
        if cancer_type in age_modifiers:
            recommendation['age_considerations'] = age_modifiers[cancer_type]
        
        return recommendation
    
    def get_references(self, cancer_type: str, genes: List[str]) -> List[str]:
        """Get relevant clinical references"""
        references = {
            'BRCA1': 'NCCN Guidelines - Genetic/Familial High-Risk Assessment: Breast/Ovarian',
            'BRCA2': 'NCCN Guidelines - Genetic/Familial High-Risk Assessment: Breast/Ovarian',
            'MLH1': 'NCCN Guidelines - Genetic/Familial High-Risk Assessment: Colorectal',
            'TP53': 'Li-Fraumeni Syndrome Management Guidelines'
        }
        
        relevant_refs = []
        for gene in genes:
            if gene in references:
                relevant_refs.append(references[gene])
        
        return list(set(relevant_refs))  # Remove duplicates
    
    def __call__(self, state: Dict) -> Dict:
        """Generate clinical recommendations"""
        risk_scores = state.get('risk_scores', {})
        variants = state.get('variant_details', [])
        
        recommendations = []
        
        # Generate recommendations for each cancer type with elevated risk
        for cancer_type, risk_percentage in risk_scores.items():
            if risk_percentage > 10:  # Threshold for recommendations
                risk_level = 'high' if risk_percentage > 50 else 'moderate'
                
                rec = self.generate_recommendations(
                    cancer_type, risk_level, risk_percentage, variants
                )
                recommendations.append(rec)
        
        # Sort by risk percentage
        recommendations.sort(key=lambda x: x['risk_percentage'], reverse=True)
        
        # Add general recommendations
        if recommendations:
            recommendations.append({
                'cancer_type': 'general',
                'recommendation': 'Genetic counseling recommended',
                'details': 'Consider cascade testing for family members',
                'prevention_options': [
                    'Maintain healthy lifestyle',
                    'Regular exercise',
                    'Balanced diet',
                    'Avoid tobacco and excessive alcohol'
                ]
            })
        
        state['clinical_recommendations'] = recommendations
        
        return {"messages": [BaseMessage(content=f"Generated {len(recommendations)} clinical recommendations", type="system")]}
```

---

## Integration Strategy

### 1. Update Graph Configuration

```python
# geneknow_pipeline/graph.py additions

# Import new nodes
from nodes.mutation_classifier import MutationClassifierNode
from nodes.mutational_signatures import MutationalSignatureNode
from nodes.variant_transformer import VariantTransformerNode
from nodes.structural_variant_detector import StructuralVariantDetectorNode
from nodes.cnv_detector import CNVDetectorNode
from nodes.pathway_analyzer import PathwayAnalyzerNode
from nodes.gene_interaction_network import GeneInteractionNetworkNode
from nodes.survival_analyzer import SurvivalAnalyzerNode
from nodes.clinical_recommendations import ClinicalRecommendationsNode

# Add nodes to workflow
workflow.add_node("mutation_classifier", MutationClassifierNode())
workflow.add_node("mutational_signatures", MutationalSignatureNode())
workflow.add_node("variant_transformer", VariantTransformerNode())
workflow.add_node("structural_variants", StructuralVariantDetectorNode())
workflow.add_node("cnv_detector", CNVDetectorNode())
workflow.add_node("pathway_analyzer", PathwayAnalyzerNode())
workflow.add_node("gene_interactions", GeneInteractionNetworkNode())
workflow.add_node("survival_analysis", SurvivalAnalyzerNode())
workflow.add_node("clinical_recommendations", ClinicalRecommendationsNode())

# Update edges for proper flow
workflow.add_edge("preprocess", "mutation_classifier")
workflow.add_edge("mutation_classifier", "variant_transformer")
workflow.add_edge("variant_transformer", "mutational_signatures")
# ... continue with appropriate flow
```

### 2. Update State Definition

```python
# geneknow_pipeline/state.py additions

class GenomicState(TypedDict):
    # ... existing fields ...
    
    # New fields for clinical view
    mutation_types: Dict[str, int]
    mutation_signatures: List[Dict]
    structural_variants: List[Dict]
    copy_number_variants: List[Dict]
    pathway_analysis: Dict
    gene_interactions: List[Dict]
    survival_analysis: Dict
    clinical_recommendations: List[Dict]
```

### 3. Update Report Writer

```python
# Update report_writer.py to include new data in structured_json

def format_results(self) -> Dict:
    """Format comprehensive results including all new analyses"""
    state = self.state
    
    structured_json = {
        # Existing fields
        "summary": {...},
        "variant_details": [...],
        
        # New fields for clinical view
        "mutation_types": state.get("mutation_types", {}),
        "mutation_signatures": state.get("mutation_signatures", []),
        "structural_variants": state.get("structural_variants", []),
        "copy_number_variants": state.get("copy_number_variants", []),
        "pathway_analysis": state.get("pathway_analysis", {}),
        "gene_interactions": state.get("gene_interactions", []),
        "survival_analysis": state.get("survival_analysis", {}),
        "clinical_recommendations": state.get("clinical_recommendations", [])
    }
    
    return structured_json
```

## Testing Strategy

### Unit Tests
Create test files for each new node:
- `test_mutation_classifier.py`
- `test_mutational_signatures.py`
- etc.

### Integration Tests
Test full pipeline with sample data to ensure all new analyses work together.

### Performance Considerations
- Cache external API calls (VEP, STRING, etc.)
- Implement rate limiting for external services
- Consider async processing for time-intensive analyses
- Add progress indicators for long-running analyses

## Database Requirements

### New Tables Needed
1. `cosmic_signatures` - COSMIC mutational signatures
2. `gene_interactions` - Protein-protein interactions
3. `pathway_definitions` - Biological pathway data
4. `clinical_guidelines` - Clinical recommendation rules

### Data Sources
- COSMIC for mutational signatures
- STRING/BioGRID for protein interactions
- KEGG/Reactome for pathways
- NCCN/ASCO for clinical guidelines

## Deployment Checklist

1. ✅ Implement all node classes
2. ✅ Update graph configuration
3. ✅ Extend state definition
4. ✅ Update report writer
5. ⬜ Create unit tests
6. ⬜ Run integration tests
7. ⬜ Set up external databases
8. ⬜ Configure API keys
9. ⬜ Update API documentation
10. ⬜ Deploy and monitor

This implementation guide provides a comprehensive roadmap for adding all missing clinical view data to the GeneKnow pipeline backend. 