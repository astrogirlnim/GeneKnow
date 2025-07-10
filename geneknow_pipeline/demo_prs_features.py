#!/usr/bin/env python3
"""
Comprehensive demonstration of PRS calculator features.
Shows how the node handles all special considerations.
"""
import sys
import os
import logging
import importlib.util

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

# Direct import of PRS calculator
spec = importlib.util.spec_from_file_location('prs_calculator', 'nodes/prs_calculator.py')
if spec and spec.loader:
    prs_calculator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prs_calculator)
else:
    raise ImportError("Could not load prs_calculator module")

def demo_coverage_impact():
    """Demonstrate how coverage affects PRS confidence."""
    logger.info("\n" + "="*60)
    logger.info("DEMO 1: Coverage Impact on PRS Confidence")
    logger.info("="*60)
    
    # Test with different numbers of matched SNPs
    base_variant = {
        'chrom': '17',
        'pos': 43044295,
        'ref': 'G',
        'alt': 'A',
        'population_frequency': 0.0001
    }
    
    # Scenario 1: Low coverage (1 SNP found)
    state_low = {
        'filtered_variants': [{
            'variant_id': '17:43044295:G>A',
            **base_variant
        }],
        'completed_nodes': [],
        'errors': [],
        'patient_population': 'EUR'
    }
    
    result_low = prs_calculator.process(state_low)
    brca_low = result_low['prs_results']['BRCA']
    
    logger.info(f"\nLow Coverage Scenario:")
    logger.info(f"  SNPs matched: {brca_low['matched_snps']}/{brca_low['total_snps']}")
    logger.info(f"  Coverage: {brca_low['coverage']:.1%}")
    logger.info(f"  Confidence: {brca_low['confidence']}")
    logger.info(f"  Raw Score: {brca_low['raw_score']:.3f}")
    logger.info(f"  Adjusted Score: {brca_low['adjusted_score']:.3f}")
    
    # Scenario 2: High coverage (multiple SNPs)
    state_high = {
        'filtered_variants': [
            {'variant_id': '17:43044295:G>A', **base_variant},
            {'variant_id': '13:32911888:A>G', 'chrom': '13', 'pos': 32911888, 
             'ref': 'A', 'alt': 'G', 'population_frequency': 0.0002},
            {'variant_id': '10:123352317:C>T', 'chrom': '10', 'pos': 123352317,
             'ref': 'C', 'alt': 'T', 'population_frequency': 0.23},
            {'variant_id': '5:56031884:T>C', 'chrom': '5', 'pos': 56031884,
             'ref': 'T', 'alt': 'C', 'population_frequency': 0.31}
        ],
        'completed_nodes': [],
        'errors': [],
        'patient_population': 'EUR'
    }
    
    result_high = prs_calculator.process(state_high)
    brca_high = result_high['prs_results']['BRCA']
    
    logger.info(f"\nHigh Coverage Scenario:")
    logger.info(f"  SNPs matched: {brca_high['matched_snps']}/{brca_high['total_snps']}")
    logger.info(f"  Coverage: {brca_high['coverage']:.1%}")
    logger.info(f"  Confidence: {brca_high['confidence']}")
    logger.info(f"  Raw Score: {brca_high['raw_score']:.3f}")
    logger.info(f"  Adjusted Score: {brca_high['adjusted_score']:.3f}")


def demo_population_stratification():
    """Demonstrate population-specific PRS calculations."""
    logger.info("\n" + "="*60)
    logger.info("DEMO 2: Population Stratification Effects")
    logger.info("="*60)
    
    # Same variant, different populations
    test_variant = {
        'variant_id': '10:123352317:C>T',
        'chrom': '10',
        'pos': 123352317,
        'ref': 'C',
        'alt': 'T',
        'population_frequency': 0.23
    }
    
    populations = ['EUR', 'AFR', 'EAS']
    results = {}
    
    for pop in populations:
        state = {
            'filtered_variants': [test_variant],
            'completed_nodes': [],
            'errors': [],
            'patient_population': pop
        }
        
        result = prs_calculator.process(state)
        if 'BRCA' in result['prs_results']:
            results[pop] = result['prs_results']['BRCA']
    
    logger.info("\nSame variant across populations:")
    logger.info(f"  Variant: {test_variant['variant_id']}")
    
    for pop, data in results.items():
        logger.info(f"\n  {pop} population:")
        logger.info(f"    Effect size: {data['contributing_snps'][0]['effect_size'] if data['contributing_snps'] else 'N/A'}")
        logger.info(f"    Score: {data['raw_score']:.3f}")
        logger.info(f"    Population used: {data['population']}")


def demo_genotype_inference():
    """Demonstrate genotype inference methods."""
    logger.info("\n" + "="*60)
    logger.info("DEMO 3: Genotype Inference Methods")
    logger.info("="*60)
    
    # Test different scenarios
    variants = [
        # Rare variant without genotype
        {
            'variant_id': '17:43044295:G>A',
            'chrom': '17',
            'pos': 43044295,
            'ref': 'G',
            'alt': 'A',
            'population_frequency': 0.0001  # Rare
        },
        # Common variant without genotype
        {
            'variant_id': '10:123352317:C>T',
            'chrom': '10',
            'pos': 123352317,
            'ref': 'C',
            'alt': 'T',
            'population_frequency': 0.23  # Common
        },
        # Variant with explicit genotype
        {
            'variant_id': '5:56031884:T>C',
            'chrom': '5',
            'pos': 56031884,
            'ref': 'T',
            'alt': 'C',
            'population_frequency': 0.31,
            'genotype': '1/1'  # Homozygous alt
        }
    ]
    
    state = {
        'filtered_variants': variants,
        'completed_nodes': [],
        'errors': [],
        'patient_population': 'EUR'
    }
    
    result = prs_calculator.process(state)
    
    # Show inference methods used
    logger.info("\nGenotype Inference Results:")
    
    for cancer_type in ['BRCA', 'OVCA']:
        if cancer_type in result['prs_results']:
            logger.info(f"\n{cancer_type} contributing SNPs:")
            for snp in result['prs_results'][cancer_type]['contributing_snps']:
                logger.info(f"  {snp['variant_id']}:")
                logger.info(f"    Risk alleles: {snp['risk_alleles']}")
                logger.info(f"    Inference method: {snp['inference_method']}")
                logger.info(f"    Contribution: {snp['contribution']:.3f}")


def demo_multi_cancer_risk():
    """Demonstrate multi-cancer risk assessment."""
    logger.info("\n" + "="*60)
    logger.info("DEMO 4: Multi-Cancer Risk Assessment")
    logger.info("="*60)
    
    # BRCA1 mutation affects multiple cancers differently
    state = {
        'filtered_variants': [{
            'variant_id': '17:43044295:G>A',
            'chrom': '17',
            'pos': 43044295,
            'ref': 'G',
            'alt': 'A',
            'population_frequency': 0.0001,
            'gene': 'BRCA1'
        }],
        'completed_nodes': [],
        'errors': [],
        'patient_population': 'EUR'
    }
    
    result = prs_calculator.process(state)
    
    logger.info("\nBRCA1 variant impact across cancer types:")
    
    # Sort by score
    cancer_scores = []
    for cancer_type, data in result['prs_results'].items():
        if data['raw_score'] > 0:
            cancer_scores.append((cancer_type, data))
    
    cancer_scores.sort(key=lambda x: x[1]['raw_score'], reverse=True)
    
    for cancer_type, data in cancer_scores:
        logger.info(f"\n  {cancer_type}:")
        logger.info(f"    Score: {data['raw_score']:.3f}")
        logger.info(f"    Percentile: {data['percentile']}%")
        logger.info(f"    Risk Category: {data['risk_category']}")
    
    # Show summary
    summary = result['prs_summary']
    logger.info(f"\nHigh-risk cancers: {', '.join(summary['high_risk_cancers']) or 'None'}")
    logger.info(f"Primary concern: {summary['primary_concern'] or 'None'}")


if __name__ == "__main__":
    # Create PRS database
    prs_calculator.create_prs_database()
    
    # Run all demos
    demo_coverage_impact()
    demo_population_stratification()
    demo_genotype_inference()
    demo_multi_cancer_risk()
    
    logger.info("\n" + "="*60)
    logger.info("Demo complete! PRS calculator handles:")
    logger.info("  ✓ Coverage-based confidence adjustment")
    logger.info("  ✓ Population-specific effect sizes")
    logger.info("  ✓ Intelligent genotype inference")
    logger.info("  ✓ Multi-cancer risk assessment")
    logger.info("="*60) 