#!/usr/bin/env python3
"""Test ML fusion integration in the pipeline."""

import sys
import json
import os

# Add to path
sys.path.append('geneknow_pipeline')

from graph import run_pipeline

# Create a test MAF file with a specific variant
test_maf_content = """Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\tdbSNP_Val_Status\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Validation_Allele1\tTumor_Validation_Allele2\tMatch_Norm_Validation_Allele1\tMatch_Norm_Validation_Allele2\tVerification_Status\tValidation_Status\tMutation_Status\tSequencing_Phase\tSequence_Source\tValidation_Method\tScore\tBAM_File\tSequencer\tTumor_Sample_UUID\tMatched_Norm_Sample_UUID\tHGVSc\tHGVSp\tHGVSp_Short\tTranscript_ID\tExon_Number\tCODONS\tProtein_position\tCDS_position\tCDS_position\ttrans_pep_id\tAmino_acids\tCodons\tExisting_variation\tALLELE_NUM\tDISTANCE\tSTRAND_VEP\tSYMBOL\tSYMBOL_SOURCE\tHGNC_ID\tBIOTYPE\tCANONICAL\tCCDS\tENSP\tSWISSPROT\tTREMBL\tUNIPARC\tRefSeq\tSIFT\tPolyPhen\tEXON\tINTRON\tDOMAINS\tAF\tAFR_AF\tAMR_AF\tASN_AF\tEAS_AF\tEUR_AF\tSAS_AF\tAA_AF\tEA_AF\tCLIN_SIG\tSOMATIC\tPUBMED\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE\tIMPACT\tPICK\tVARIANT_CLASS\tTSL\tHGVS_OFFSET\tPHASTCONS\tMINIMISED\tExAC_AF\tExAC_AF_Adj\tExAC_AF_AFR\tExAC_AF_AMR\tExAC_AF_EAS\tExAC_AF_FIN\tExAC_AF_NFE\tExAC_AF_OTH\tExAC_AF_SAS\tGENE_PHENO\tFILTER\tCOSMIC\tCENTERS\tCONTEXT\tDBGAP_MAF\tEA_MAF\tSOURCE\tNCBI_VAL_CITED\tNCBI_VAL_TESTED\tt_alt_count\tt_ref_count\tn_alt_count\tn_ref_count\tGCC_VALIDATION_Tumor_Seq_Allele1\tGCC_VALIDATION_Tumor_Seq_Allele2\tGCC_VALIDATION_Match_Norm_Seq_Allele1\tGCC_VALIDATION_Match_Norm_Seq_Allele2\tGCC_VALIDATION_Verification_Status\tGCC_VALIDATION_Validation_Status\tGCC_VALIDATION_Mutation_Status\tHOTSPOT\tONCOGENE\tTSG\tNOTE\tKRAS_amplicons\tKRAS_hotspots\tvac_annotation\tEFFECT\tONCOGENIC\tCDS_STRAND\tALT\tREF\tPOS\tVAR_ID\tREF_CONTEXT\tPAN_COUNT\tPAN_FREQ\tQUERY_PAN_COUNT\tQUERY_PAN_FREQ\tQUERY_PAN_STATUS\tVAR_UUID\tREVIEW_STATUS\n
BRAF\t673\t.\tGRCh38\tchr7\t140787560\t140787560\t+\tMissense_Mutation\tSNP\tG\tA\tA\t.\t.\tTEST-SAMPLE\t.\t.\t.\t.\t.\t.\t.\t.\t.\tSomatic\t.\t.\t.\t.\t.\tIllumina HiSeq\t.\t.\tc.1798G>A\tp.V600M\tp.V600M\tENST00000288602\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tBRAF\t.\t.\tprotein_coding\tYES\t.\t.\t.\t.\t.\t.\ttolerated(0.23)\tbenign(0.152)\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tMODERATE\t.\tSNV\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t100\t50\t0\t100\t.\t.\t.\t.\t.\t.\t.\t.\t.\tY\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t."""

# Write test MAF file
with open('test_ml_fusion_demo.maf', 'w') as f:
    f.write(test_maf_content)

print("=" * 80)
print("Testing ML Fusion Integration in GeneKnow Pipeline")
print("=" * 80)

# Run pipeline
result = run_pipeline('test_ml_fusion_demo.maf')

# Check if ML fusion ran
print("\n1. Checking if ML fusion node ran:")
if result.get('ml_fusion_results'):
    ml_results = result['ml_fusion_results']
    print(f"   ✅ ML fusion completed: {ml_results.get('processing_successful', False)}")
    
    if ml_results.get('aggregate_risk_assessment'):
        aggregate = ml_results['aggregate_risk_assessment']
        print(f"   - Aggregate risk score: {aggregate.get('aggregate_risk_score', 0):.3f}")
        print(f"   - Risk category: {aggregate.get('risk_category', 'unknown')}")
        print(f"   - Confidence: {aggregate.get('confidence', 0):.3f}")
        print(f"   - High-risk variants: {aggregate.get('high_risk_variants', 0)}")
else:
    print("   ❌ ML fusion results not found in state")

# Check if risk model used ML fusion
print("\n2. Checking if risk model used ML fusion:")
if result.get('ml_risk_assessment'):
    ml_risk = result['ml_risk_assessment']
    print(f"   ✅ Risk model method: {ml_risk.get('method', 'unknown')}")
    if ml_risk.get('method') == 'ml_fusion':
        print("   ✅ Risk model successfully used ML fusion results!")
    else:
        print("   ❌ Risk model did NOT use ML fusion")
else:
    print("   ❌ No ML risk assessment found")

# Check warnings
print("\n3. Checking for warnings:")
warnings = result.get('warnings', [])
for warning in warnings:
    if 'ML' in warning.get('warning', '') or 'fusion' in warning.get('warning', ''):
        print(f"   ⚠️ {warning['warning']}")

# Show all completed nodes
print("\n4. Completed nodes:")
completed = result.get('completed_nodes', [])
for node in completed:
    print(f"   - {node}")

# Check if all data was properly processed
print("\n5. Data processing summary:")
print(f"   - Total variants: {result.get('variant_count', 0)}")
print(f"   - CADD stats available: {'cadd_stats' in result and result['cadd_stats'] is not None}")
print(f"   - TCGA matches available: {'tcga_matches' in result and result['tcga_matches'] is not None}")
print(f"   - ClinVar annotations available: {'clinvar_annotations' in result and result['clinvar_annotations'] is not None}")
print(f"   - PRS results available: {'prs_results' in result and result['prs_results'] is not None}")
print(f"   - Pathway burden available: {'pathway_burden_results' in result and result['pathway_burden_results'] is not None}")
print(f"   - ML ready variants: {len(result.get('ml_ready_variants', [])) if result.get('ml_ready_variants') else 0}")

# Show risk scores
print("\n6. Final risk scores:")
risk_scores = result.get('risk_scores', {})
for cancer, score in risk_scores.items():
    print(f"   - {cancer}: {score}%")

# Clean up
os.remove('test_ml_fusion_demo.maf')
print("\n✅ Test complete!") 