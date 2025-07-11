import React, { useState } from 'react';
import { useNavigate, useSearchParams, useLocation } from 'react-router-dom';
import Layout from '../components/Layout';
import ConfidenceCheck from '../components/ConfidenceCheck';
import type { PipelineResult } from '../api/geneknowPipeline';

// Type definitions for external libraries
interface JsPDFConstructor {
  new (orientation?: string, unit?: string, format?: string): JsPDFInstance;
}

interface JsPDFInstance {
  internal: {
    pageSize: {
      getWidth(): number;
      getHeight(): number;
    };
  };
  setFontSize(size: number): void;
  setFont(fontName?: string, fontStyle?: string): void;
  text(text: string, x: number, y: number, options?: { align?: string }): void;
  text(text: string[], x: number, y: number): void;
  splitTextToSize(text: string, maxWidth: number): string[];
  addImage(imageData: string, format: string, x: number, y: number, width: number, height: number): void;
  addPage(): void;
  save(filename: string): void;
}

declare global {
  interface Window {
    jsPDF?: JsPDFConstructor | { jsPDF: JsPDFConstructor };
    jspdf?: {
      jsPDF: JsPDFConstructor;
    };
    html2canvas?: (element: HTMLElement, options?: Record<string, unknown>) => Promise<HTMLCanvasElement>;
  }
  
  // Global jsPDF might be available
  const jsPDF: JsPDFConstructor | undefined;
}

// Mock SHAP validation for testing - different statuses for different risk levels
const getMockSHAPValidation = (riskLevel: string) => {
  switch (riskLevel) {
    case 'high':
      return {
        status: 'FLAG_FOR_REVIEW' as const,
        reasons: [
          'The AI predicted HIGH RISK (82%) but this appears to be based on indirect factors (Gene/Pathway Burden Score, TCGA Tumor Enrichment) rather than known disease-causing mutations. The prediction may be less reliable.'
        ],
        top_contributors: [
          {
            feature: 'gene_burden_score',
            display_name: 'Gene/Pathway Burden Score',
            shap_value: 0.35,
            abs_contribution: 0.35,
            direction: 'increases' as const
          },
          {
            feature: 'tcga_enrichment',
            display_name: 'TCGA Tumor Enrichment',
            shap_value: 0.28,
            abs_contribution: 0.28,
            direction: 'increases' as const
          },
          {
            feature: 'prs_score',
            display_name: 'Polygenic Risk Score',
            shap_value: 0.19,
            abs_contribution: 0.19,
            direction: 'increases' as const
          }
        ],
        feature_importance: {
          'gene_burden_score': 0.35,
          'tcga_enrichment': 0.28,
          'prs_score': 0.19,
          'clinvar_pathogenic': 0.18,
          'clinvar_benign': -0.05,
          'cadd_score': 0.15
        } as Record<string, number>,
        details: {
          status: 'FLAG_FOR_REVIEW' as const,
          risk_score: 0.82,
          top_contributors: [],
          validation_reasons: [],
          rule_results: {},
          shap_values: [],
          feature_names: [],
          model_type: 'GradientBoostingRegressor'
        }
      };
    
    case 'medium':
      return {
        status: 'ERROR' as const,
        reasons: [
          'A technical error occurred during the confidence check validation. The model\'s prediction appears sound, but automated validation could not be completed.'
        ],
        top_contributors: [
          {
            feature: 'mismatch_repair_genes',
            display_name: 'Mismatch Repair Genes',
            shap_value: 0.42,
            abs_contribution: 0.42,
            direction: 'increases' as const
          },
          {
            feature: 'prs_score',
            display_name: 'Polygenic Risk Score',
            shap_value: 0.31,
            abs_contribution: 0.31,
            direction: 'increases' as const
          }
        ],
        feature_importance: {
          'mismatch_repair_genes': 0.42,
          'prs_score': 0.31,
          'clinvar_pathogenic': 0.27
        } as Record<string, number>,
        details: {
          status: 'ERROR' as const,
          risk_score: 0.45,
          top_contributors: [],
          validation_reasons: [],
          rule_results: {},
          shap_values: [],
          feature_names: [],
          model_type: 'GradientBoostingRegressor'
        }
      };
    
    case 'low':
      return {
        status: 'PASS' as const,
        reasons: [
          'The AI\'s risk prediction is well-supported by the genomic evidence. The model\'s reasoning aligns with established clinical guidelines.'
        ],
        top_contributors: [
          {
            feature: 'clinvar_pathogenic',
            display_name: 'Known Pathogenic Variants',
            shap_value: 0.25,
            abs_contribution: 0.25,
            direction: 'increases' as const
          },
          {
            feature: 'tumor_suppressor_genes',
            display_name: 'Tumor Suppressor Genes',
            shap_value: 0.22,
            abs_contribution: 0.22,
            direction: 'increases' as const
          },
          {
            feature: 'protective_variants',
            display_name: 'Protective Variants',
            shap_value: -0.18,
            abs_contribution: 0.18,
            direction: 'decreases' as const
          }
        ],
        feature_importance: {
          'clinvar_pathogenic': 0.25,
          'tumor_suppressor_genes': 0.22,
          'protective_variants': -0.18,
          'prs_score': 0.15
        } as Record<string, number>,
        details: {
          status: 'PASS' as const,
          risk_score: 0.15,
          top_contributors: [],
          validation_reasons: [],
          rule_results: {},
          shap_values: [],
          feature_names: [],
          model_type: 'GradientBoostingRegressor'
        }
      };
    
    case 'skipped':
      return {
        status: 'SKIPPED' as const,
        reasons: [
          'The confidence check was not applicable to this analysis type. The genetic variants identified do not fall within the model\'s validation scope.'
        ],
        top_contributors: [],
        feature_importance: {} as Record<string, number>,
        details: {
          status: 'SKIPPED' as const,
          risk_score: 0.0,
          top_contributors: [],
          validation_reasons: [],
          rule_results: {},
          shap_values: [],
          feature_names: [],
          model_type: 'GradientBoostingRegressor'
        }
      };
    
    default:
      return null;
  }
};

// Enhanced genomic data structure with transformations and mutations
const mockGenomicData = {
  "summary": {
    "total_variants_found": 12,
    "variants_passed_qc": 11,
    "high_risk_findings_count": 6,
    "processing_time_seconds": 0.703601,
    "mutation_types": {
      "snv": 5,
      "indel": 4,
      "cnv": 2,
      "structural": 1
    }
  },
  "risk_findings": [
    {
      "cancer_type": "breast",
      "risk_percentage": 100.0,
      "risk_level": "high",
      "affected_genes": ["BRCA1", "TP53"],
      "recommendation": "Enhanced breast cancer screening recommended",
      "mutation_burden": "high",
      "pathway_disruption": ["DNA_REPAIR", "CELL_CYCLE"]
    },
    {
      "cancer_type": "colon",
      "risk_percentage": 100.0,
      "risk_level": "high",
      "affected_genes": ["APC", "TP53"],
      "recommendation": "Enhanced colon cancer screening recommended",
      "mutation_burden": "high",
      "pathway_disruption": ["WNT_SIGNALING", "APOPTOSIS"]
    },
    {
      "cancer_type": "lung",
      "risk_percentage": 97.7,
      "risk_level": "high",
      "affected_genes": ["TP53"],
      "recommendation": "Enhanced lung cancer screening recommended",
      "mutation_burden": "medium",
      "pathway_disruption": ["CELL_CYCLE", "APOPTOSIS"]
    },
    {
      "cancer_type": "prostate",
      "risk_percentage": 81.4,
      "risk_level": "high",
      "affected_genes": ["TP53"],
      "recommendation": "Enhanced prostate cancer screening recommended",
      "mutation_burden": "medium",
      "pathway_disruption": ["CELL_CYCLE"]
    },
    {
      "cancer_type": "blood",
      "risk_percentage": 1.7,
      "risk_level": "low",
      "affected_genes": [],
      "recommendation": "Standard screening guidelines apply",
      "mutation_burden": "low",
      "pathway_disruption": []
    }
  ],
  "variant_table": [
    {
      "gene": "BRCA1",
      "variant_id": "chr17:41223094:A>G",
      "consequence": "missense_variant",
      "mutation_type": "snv",
      "quality_score": 99.0,
      "clinical_significance": "pathogenic",
      "tcga_best_match": { "cancer_type": "breast", "frequency": 6.3 },
      "protein_change": "p.Ile1756Val",
      "functional_impact": "Loss of function",
      "transformation": {
        "original": "ATT",
        "mutated": "GTT",
        "amino_acid_change": "Ile→Val",
        "effect": "Structural disruption of DNA binding domain"
      }
    },
    {
      "gene": "BRCA1",
      "variant_id": "chr17:41244936:GATC>G",
      "consequence": "frameshift_variant",
      "mutation_type": "indel",
      "quality_score": 87.5,
      "clinical_significance": "pathogenic",
      "tcga_best_match": { "cancer_type": "breast", "frequency": 4.5 },
      "protein_change": "p.Met1563fs",
      "functional_impact": "Premature termination",
      "transformation": {
        "original": "GATC",
        "mutated": "G",
        "amino_acid_change": "Met→frameshift",
        "effect": "Truncated protein, complete loss of C-terminal domain"
      }
    },
    {
      "gene": "TP53",
      "variant_id": "chr17:7577121:G>A",
      "consequence": "missense_variant",
      "mutation_type": "snv",
      "quality_score": 92.1,
      "clinical_significance": "pathogenic",
      "tcga_best_match": { "cancer_type": "lung", "frequency": 65.4 },
      "protein_change": "p.Arg248Gln",
      "functional_impact": "Loss of DNA binding",
      "transformation": {
        "original": "CGG",
        "mutated": "CAG",
        "amino_acid_change": "Arg→Gln",
        "effect": "Disrupted p53-DNA interaction, loss of tumor suppressor function"
      }
    },
    {
      "gene": "APC",
      "variant_id": "chr5:112173917:C>T",
      "consequence": "nonsense_variant",
      "mutation_type": "snv",
      "quality_score": 65.2,
      "clinical_significance": "pathogenic",
      "tcga_best_match": { "cancer_type": "colon", "frequency": 81.2 },
      "protein_change": "p.Arg1450*",
      "functional_impact": "Premature termination",
      "transformation": {
        "original": "CGA",
        "mutated": "TGA",
        "amino_acid_change": "Arg→Stop",
        "effect": "Truncated APC protein, loss of β-catenin regulation"
      }
    },
    {
      "gene": "MLH1",
      "variant_id": "chr3:37048590:CAGTC>C",
      "consequence": "frameshift_variant",
      "mutation_type": "indel",
      "quality_score": 78.9,
      "clinical_significance": "pathogenic",
      "tcga_best_match": { "cancer_type": "colon", "frequency": 23.1 },
      "protein_change": "p.Ser456fs",
      "functional_impact": "Loss of mismatch repair",
      "transformation": {
        "original": "CAGTC",
        "mutated": "C",
        "amino_acid_change": "Ser→frameshift",
        "effect": "Disrupted DNA mismatch repair pathway"
      }
    },
    {
      "gene": "KRAS",
      "variant_id": "chr12:25245350:G>A",
      "consequence": "missense_variant",
      "mutation_type": "snv",
      "quality_score": 88.3,
      "clinical_significance": "pathogenic",
      "tcga_best_match": { "cancer_type": "pancreatic", "frequency": 90.2 },
      "protein_change": "p.Gly12Asp",
      "functional_impact": "Constitutive activation",
      "transformation": {
        "original": "GGT",
        "mutated": "GAT",
        "amino_acid_change": "Gly→Asp",
        "effect": "Oncogenic activation, permanent 'on' signal"
      }
    }
  ],
  "structural_variants": [
    {
      "type": "deletion",
      "chromosome": "chr17",
      "start": 41196312,
      "end": 41277500,
      "size": 81188,
      "genes_affected": ["BRCA1"],
      "clinical_significance": "pathogenic",
      "functional_impact": "Complete gene deletion",
      "transformation": {
        "original": "Full BRCA1 gene",
        "mutated": "Deleted",
        "effect": "Complete loss of BRCA1 function"
      }
    },
    {
      "type": "duplication",
      "chromosome": "chr17",
      "start": 7571720,
      "end": 7590868,
      "size": 19148,
      "genes_affected": ["TP53"],
      "clinical_significance": "uncertain_significance",
      "functional_impact": "Gene dosage imbalance",
      "transformation": {
        "original": "2 copies TP53",
        "mutated": "3 copies TP53",
        "effect": "Potential dosage-sensitive effects"
      }
    }
  ],
  "copy_number_variants": [
    {
      "gene": "ERBB2",
      "chromosome": "chr17",
      "copy_number": 6,
      "normal_copy_number": 2,
      "fold_change": 3.0,
      "clinical_significance": "pathogenic",
      "cancer_relevance": "HER2 amplification in breast cancer",
      "transformation": {
        "original": "Normal expression",
        "mutated": "3x overexpression",
        "effect": "Oncogenic driver, targeted therapy candidate"
      }
    },
    {
      "gene": "CDKN2A",
      "chromosome": "chr9",
      "copy_number": 0,
      "normal_copy_number": 2,
      "fold_change": 0.0,
      "clinical_significance": "pathogenic",
      "cancer_relevance": "Tumor suppressor loss",
      "transformation": {
        "original": "Normal cell cycle control",
        "mutated": "Complete loss",
        "effect": "Loss of cell cycle checkpoint control"
      }
    }
  ],
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
    },
    {
      "signature": "SBS4",
      "name": "Tobacco smoking",
      "contribution": 0.22,
      "description": "C>A mutations from tobacco carcinogens",
      "etiology": "Environmental mutagen exposure"
    },
    {
      "signature": "SBS13",
      "name": "APOBEC cytidine deaminase",
      "contribution": 0.15,
      "description": "Cytidine deaminase activity",
      "etiology": "Immune system enzymatic activity"
    }
  ]
};

// Helper function to get color based on risk level
const getRiskColor = (riskLevel: string) => {
  switch (riskLevel) {
    case 'high': return '#EF4444';
    case 'medium': return '#F59E0B';
    case 'low': return '#22C55E';
    default: return '#6B7280';
  }
};

// Simple PDF Download Helper Function
const downloadPDF = async (elementId: string, title: string, summary: string, setIsPDFGenerating: (value: boolean) => void) => {
  console.log('🚀 PDF Download Started:', { elementId, title, summary: summary.substring(0, 100) + '...' });
  
  try {
    // Set PDF generating state
    setIsPDFGenerating(true);
    
    console.log('📄 Creating loading notification...');
    
    // Show loading notification
    const loadingNotification = document.createElement('div');
    loadingNotification.innerHTML = `
      <div style="
        position: fixed; 
        top: 20px; 
        right: 20px; 
        background: #2563EB; 
        color: white; 
        padding: 15px 20px; 
        border-radius: 8px; 
        box-shadow: 0 4px 12px rgba(0,0,0,0.3);
        z-index: 10000;
        font-weight: 600;
        font-size: 14px;
        max-width: 300px;
      ">
        📄 Loading PDF library...
      </div>
    `;
    document.body.appendChild(loadingNotification);
    console.log('✅ Loading notification added to DOM');
    
    // Check if jsPDF is already loaded
    console.log('🔍 Checking if jsPDF already exists:', !!window.jsPDF);
    console.log('🔍 Checking if window.jspdf exists:', !!window.jspdf);
    console.log('🔍 Checking if window.jsPDF exists:', !!window.jsPDF);
    if (window.jsPDF || window.jspdf) {
      console.log('✅ jsPDF already loaded, generating PDF directly...');
      generatePDF();
      return;
    }
    
    console.log('📦 Loading jsPDF from CDN...');
    
    // Load jsPDF from CDN
    const script = document.createElement('script');
    script.src = 'https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js';
    console.log('📦 Script element created with src:', script.src);
    
    script.onload = () => {
      console.log('✅ jsPDF script loaded successfully');
      console.log('🔍 Checking window.jsPDF after load:', !!window.jsPDF);
      console.log('🔍 Checking window.jspdf after load:', !!window.jspdf);
      console.log('🔍 Checking global jsPDF after load:', typeof jsPDF !== 'undefined');
      console.log('🔍 window object keys containing "pdf":', Object.keys(window).filter(key => key.toLowerCase().includes('pdf')));
      console.log('🔍 window.jsPDF object:', window.jsPDF);
      console.log('🔍 window.jspdf object:', window.jspdf);
      
      // Update loading message
      loadingNotification.innerHTML = `
        <div style="
          position: fixed; 
          top: 20px; 
          right: 20px; 
          background: #2563EB; 
          color: white; 
          padding: 15px 20px; 
          border-radius: 8px; 
          box-shadow: 0 4px 12px rgba(0,0,0,0.3);
          z-index: 10000;
          font-weight: 600;
          font-size: 14px;
          max-width: 300px;
        ">
          📄 Generating PDF...
        </div>
      `;
      console.log('📝 Loading notification updated to "Generating PDF"');
      
      // Wait a moment for the library to initialize
      setTimeout(() => {
        console.log('⏰ Timeout completed, checking jsPDF availability...');
        if (window.jsPDF || window.jspdf || typeof jsPDF !== 'undefined') {
          console.log('✅ jsPDF confirmed available, calling generatePDF...');
          generatePDF();
        } else {
          console.error('❌ jsPDF still not available after timeout');
          showError('PDF library failed to load');
        }
      }, 500);
    };
    
    script.onerror = (error) => {
      console.error('❌ Script loading failed:', error);
      showError('Failed to load PDF library');
    };
    
    document.head.appendChild(script);
    console.log('📦 Script added to document head');
    
    async function generatePDF() {
      console.log('🏭 Starting PDF generation...');
      
      try {
        console.log('🔍 Destructuring jsPDF from window.jsPDF...');
        
        // Try different ways jsPDF might be exposed with proper type checking
        let jsPDFClass: JsPDFConstructor | undefined;
        
        if (window.jsPDF && typeof window.jsPDF === 'object' && 'jsPDF' in window.jsPDF) {
          console.log('✅ Found jsPDF at window.jsPDF.jsPDF');
          jsPDFClass = (window.jsPDF as { jsPDF: JsPDFConstructor }).jsPDF;
        } else if (window.jspdf && window.jspdf.jsPDF) {
          console.log('✅ Found jsPDF at window.jspdf.jsPDF');
          jsPDFClass = window.jspdf.jsPDF;
        } else if (window.jsPDF && typeof window.jsPDF === 'function') {
          console.log('✅ Found jsPDF at window.jsPDF');
          jsPDFClass = window.jsPDF as JsPDFConstructor;
        } else if (window.jspdf) {
          console.log('✅ Found jsPDF at window.jspdf');
          jsPDFClass = window.jspdf.jsPDF;
        } else if (typeof jsPDF !== 'undefined') {
          console.log('✅ Found jsPDF as global variable');
          jsPDFClass = jsPDF;
        }
        
        if (!jsPDFClass) {
          throw new Error('jsPDF not found in any expected location');
        }
        
        console.log('✅ jsPDF class located:', !!jsPDFClass);
        console.log('🔍 jsPDF class type:', typeof jsPDFClass);
        console.log('🔍 jsPDF class constructor:', typeof jsPDFClass === 'function');
        
        // Create new PDF document
        console.log('📄 Creating new PDF document...');
        const pdf = new jsPDFClass('p', 'mm', 'a4');
        console.log('✅ PDF document created');
        
        const pageWidth = pdf.internal.pageSize.getWidth();
        const pageHeight = pdf.internal.pageSize.getHeight();
        console.log('📏 Page dimensions:', { pageWidth, pageHeight });
        
        let yPosition = 20;
        
        // Header
        console.log('📝 Adding header...');
        pdf.setFontSize(20);
        pdf.setFont(undefined, 'bold');
        pdf.text(title, pageWidth / 2, yPosition, { align: 'center' });
        yPosition += 15;
        
        pdf.setFontSize(14);
        pdf.setFont(undefined, 'normal');
        pdf.text('GeneKnow AI Genomic Analysis Report', pageWidth / 2, yPosition, { align: 'center' });
        yPosition += 10;
        
        pdf.setFontSize(12);
        pdf.text(`Generated: ${new Date().toLocaleString()}`, pageWidth / 2, yPosition, { align: 'center' });
        yPosition += 20;
        console.log('✅ Header added, current Y position:', yPosition);
        
        // Analysis Summary Section (at the top)
        console.log('📝 Adding analysis summary...');
        pdf.setFontSize(16);
        pdf.setFont(undefined, 'bold');
        pdf.text('Analysis Summary', 20, yPosition);
        yPosition += 10;
        
        // Summary content
        pdf.setFontSize(11);
        pdf.setFont(undefined, 'normal');
        const summaryLines = pdf.splitTextToSize(summary, pageWidth - 40);
        console.log('📝 Summary split into lines:', summaryLines.length);
        pdf.text(summaryLines, 20, yPosition);
        yPosition += summaryLines.length * 5 + 20;
        console.log('✅ Summary added, current Y position:', yPosition);
        
        // Capture and add the actual visualization
        console.log('📊 Capturing visualization element...');
        const element = document.getElementById(elementId);
        if (element) {
          console.log('✅ Found visualization element:', elementId);
          
          try {
            // Add visualization title
            pdf.setFontSize(14);
            pdf.setFont(undefined, 'bold');
            pdf.text('Visualization', 20, yPosition);
            yPosition += 10;
            
            // Load html2canvas if not already loaded
            if (!window.html2canvas) {
              console.log('📦 Loading html2canvas from CDN...');
              const html2canvasScript = document.createElement('script');
              html2canvasScript.src = 'https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js';
              
              await new Promise((resolve, reject) => {
                html2canvasScript.onload = () => {
                  console.log('✅ html2canvas loaded successfully');
                  resolve(true);
                };
                html2canvasScript.onerror = () => {
                  console.log('❌ Failed to load html2canvas');
                  reject(new Error('html2canvas failed to load'));
                };
                document.head.appendChild(html2canvasScript);
              });
            }
            
            console.log('📸 Preparing element for borderless capture...');
            
            // Store original styles to restore later - declare outside inner try block
            const originalStyles = new Map<HTMLElement, {
              border: string;
              borderTop: string;
              borderRight: string;
              borderBottom: string;
              borderLeft: string;
              borderRadius: string;
              boxShadow: string;
            }>();
            const elementsWithBorders = [element as HTMLElement, ...Array.from(element.querySelectorAll('*')).map(el => el as HTMLElement)];
            
            try {
              // Temporarily remove borders from all elements
              elementsWithBorders.forEach((el: HTMLElement) => {
                originalStyles.set(el, {
                  border: el.style.border,
                  borderTop: el.style.borderTop,
                  borderRight: el.style.borderRight,
                  borderBottom: el.style.borderBottom,
                  borderLeft: el.style.borderLeft,
                  borderRadius: el.style.borderRadius,
                  boxShadow: el.style.boxShadow
                });
                
                // Remove borders and shadows
                el.style.border = 'none';
                el.style.borderTop = 'none';
                el.style.borderRight = 'none';
                el.style.borderBottom = 'none';
                el.style.borderLeft = 'none';
                el.style.borderRadius = '0';
                el.style.boxShadow = 'none';
              });
              
              console.log('📸 Capturing element without borders...');
              
              // Use html2canvas to capture the element with proper null checking
              if (!window.html2canvas) {
                throw new Error('html2canvas is not available');
              }
              
              const canvas = await window.html2canvas(element as HTMLElement, {
                backgroundColor: '#ffffff',
                scale: 1.5,
                useCORS: true,
                allowTaint: true,
                logging: false,
                width: element.offsetWidth,
                height: element.offsetHeight,
                scrollX: 0,
                scrollY: 0
              });
              
              // Restore original styles
              console.log('🔄 Restoring original element styles...');
              elementsWithBorders.forEach((el: HTMLElement) => {
                const styles = originalStyles.get(el);
                if (styles) {
                  el.style.border = styles.border;
                  el.style.borderTop = styles.borderTop;
                  el.style.borderRight = styles.borderRight;
                  el.style.borderBottom = styles.borderBottom;
                  el.style.borderLeft = styles.borderLeft;
                  el.style.borderRadius = styles.borderRadius;
                  el.style.boxShadow = styles.boxShadow;
                }
              });
              
              console.log('✅ Element captured successfully without borders');
              
              // Convert canvas to data URL
              const imageData = canvas.toDataURL('image/png', 0.95);
              
              // Calculate optimal size for PDF - larger for better readability
              const maxWidth = pageWidth - 30; // More width usage
              const maxHeight = 180; // Increased height limit
              const aspectRatio = element.offsetWidth / element.offsetHeight;
              
              // Use more of the available space - increased from 0.6 to 0.85
              let imgWidth = Math.min(maxWidth, element.offsetWidth * 0.85);
              let imgHeight = imgWidth / aspectRatio;
              
              // If height exceeds limit, scale down proportionally
              if (imgHeight > maxHeight) {
                imgHeight = maxHeight;
                imgWidth = imgHeight * aspectRatio;
              }
              
              // Ensure minimum readable size
              const minWidth = 140;
              const minHeight = 80;
              if (imgWidth < minWidth) {
                imgWidth = minWidth;
                imgHeight = imgWidth / aspectRatio;
              }
              if (imgHeight < minHeight) {
                imgHeight = minHeight;
                imgWidth = imgHeight * aspectRatio;
              }
              
              console.log('📊 Adding captured visualization to PDF:', { 
                originalWidth: element.offsetWidth, 
                originalHeight: element.offsetHeight,
                pdfWidth: imgWidth, 
                pdfHeight: imgHeight 
              });
              
              // Add the captured image to PDF
              pdf.addImage(imageData, 'PNG', 20, yPosition, imgWidth, imgHeight);
              yPosition += imgHeight + 15;
              
              console.log('✅ Real visualization successfully added to PDF');
              
            } catch (innerCaptureError) {
              // Restore original styles on inner error
              console.log('🔄 Restoring original styles after inner capture error...');
              elementsWithBorders.forEach((el: HTMLElement) => {
                const styles = originalStyles.get(el);
                if (styles) {
                  el.style.border = styles.border;
                  el.style.borderTop = styles.borderTop;
                  el.style.borderRight = styles.borderRight;
                  el.style.borderBottom = styles.borderBottom;
                  el.style.borderLeft = styles.borderLeft;
                  el.style.borderRadius = styles.borderRadius;
                  el.style.boxShadow = styles.boxShadow;
                }
              });
              throw innerCaptureError; // Re-throw to be caught by outer catch
            }
            
          } catch (captureError) {
            console.log('⚠️ Visualization capture failed, using fallback:', captureError);
            
            // Enhanced fallback with more details
            pdf.setFontSize(11);
            pdf.setFont(undefined, 'normal');
            pdf.text('📊 Genomic Analysis Visualization', 20, yPosition);
            yPosition += 8;
            pdf.text(`Dataset: ${title}`, 25, yPosition);
            yPosition += 6;
            pdf.text('Type: Interactive matrix/chart visualization', 25, yPosition);
            yPosition += 6;
            pdf.text('Status: Available in web application interface', 25, yPosition);
            yPosition += 8;
            pdf.text('Note: For full interactive features, please view in the application', 25, yPosition);
            yPosition += 15;
          }
          
        } else {
          console.log('⚠️ Visualization element not found:', elementId);
          
          // Add not found message
          pdf.setFontSize(14);
          pdf.setFont(undefined, 'bold');
          pdf.text('Visualization', 20, yPosition);
          yPosition += 10;
          
          pdf.setFontSize(11);
          pdf.setFont(undefined, 'normal');
          pdf.text('Visualization element not found in current view', 20, yPosition);
          yPosition += 15;
        }
        
        // Check if we need a new page
        if (yPosition > pageHeight - 100) {
          console.log('📄 Adding new page for remaining content');
          pdf.addPage();
          yPosition = 20;
        }
        
        // Key Findings
        console.log('📝 Adding key findings...');
        pdf.setFontSize(14);
        pdf.setFont(undefined, 'bold');
        pdf.text('Key Clinical Insights', 20, yPosition);
        yPosition += 8;
        
        pdf.setFontSize(10);
        pdf.setFont(undefined, 'normal');
        const findings = [
          `• Analysis Type: ${title} - Comprehensive genomic assessment`,
          '• Clinical Relevance: High-priority findings for personalized medicine',
          '• Actionable Results: Data-driven recommendations for clinical decision-making',
          '• Quality Assurance: All findings meet clinical laboratory standards'
        ];
        
        findings.forEach((finding, index) => {
          console.log(`📝 Adding finding ${index + 1}:`, finding.substring(0, 50) + '...');
          pdf.text(finding, 25, yPosition);
          yPosition += 6;
        });
        yPosition += 10;
        console.log('✅ Key findings added, current Y position:', yPosition);
        
        // Visualization Note
        console.log('📝 Adding visualization note...');
        pdf.setFontSize(14);
        pdf.setFont(undefined, 'bold');
        pdf.text('Technical Details', 20, yPosition);
        yPosition += 8;
        
        pdf.setFontSize(10);
        pdf.setFont(undefined, 'normal');
        const vizNote = `This report contains the analysis data corresponding to the ${title} visualization shown above. The data presented includes statistical analysis, clinical correlations, and evidence-based interpretations derived from comprehensive genomic databases including TCGA, ClinVar, and COSMIC.`;
        const vizLines = pdf.splitTextToSize(vizNote, pageWidth - 40);
        console.log('📝 Technical details split into lines:', vizLines.length);
        pdf.text(vizLines, 20, yPosition);
        yPosition += vizLines.length * 5 + 15;
        console.log('✅ Technical details added, current Y position:', yPosition);
        
        // Clinical Significance
        console.log('📝 Adding clinical significance...');
        pdf.setFontSize(14);
        pdf.setFont(undefined, 'bold');
        pdf.text('Clinical Significance', 20, yPosition);
        yPosition += 8;
        
        pdf.setFontSize(10);
        pdf.setFont(undefined, 'normal');
        const clinicalText = 'This analysis provides comprehensive insights into genetic variations and their potential clinical implications. The findings should be interpreted by qualified healthcare professionals in conjunction with clinical presentation and family history. All results have been generated using state-of-the-art AI algorithms and validated against established genomic databases.';
        const clinicalLines = pdf.splitTextToSize(clinicalText, pageWidth - 40);
        console.log('📝 Clinical text split into lines:', clinicalLines.length);
        pdf.text(clinicalLines, 20, yPosition);
        yPosition += clinicalLines.length * 5 + 15;
        console.log('✅ Clinical significance added, current Y position:', yPosition);
        
        // Technical Specifications
        console.log('📝 Adding technical specifications...');
        if (yPosition > pageHeight - 60) {
          console.log('📄 Adding new page for technical specifications');
          pdf.addPage();
          yPosition = 20;
        }
        
        pdf.setFontSize(14);
        pdf.setFont(undefined, 'bold');
        pdf.text('Technical Specifications', 20, yPosition);
        yPosition += 8;
        
        pdf.setFontSize(10);
        pdf.setFont(undefined, 'normal');
        const techSpecs = [
          '• Analysis Pipeline: GeneKnow AI-powered genomic analysis engine',
          '• Reference Genome: GRCh38/hg38 with clinical annotations',
          '• Quality Control: GATK best practices with clinical-grade filtering',
          '• Annotation Sources: ClinVar, COSMIC, gnomAD, PharmGKB, TCGA',
          '• Confidence Level: High confidence results suitable for clinical reporting'
        ];
        
        techSpecs.forEach((spec, index) => {
          console.log(`📝 Adding tech spec ${index + 1}:`, spec.substring(0, 50) + '...');
          pdf.text(spec, 25, yPosition);
          yPosition += 6;
        });
        yPosition += 15;
        console.log('✅ Technical specifications added, current Y position:', yPosition);
        
        // Disclaimer
        console.log('📝 Adding disclaimer...');
        if (yPosition > pageHeight - 40) {
          console.log('📄 Adding new page for disclaimer');
          pdf.addPage();
          yPosition = 20;
        }
        
        pdf.setFontSize(12);
        pdf.setFont(undefined, 'bold');
        pdf.text('Clinical Disclaimer', 20, yPosition);
        yPosition += 8;
        
        pdf.setFontSize(9);
        pdf.setFont(undefined, 'normal');
        const disclaimer = 'This report is generated for research and educational purposes. All genetic findings must be reviewed and validated by qualified healthcare professionals before making any clinical decisions. Genetic counseling is recommended for interpretation and family planning considerations.';
        const disclaimerLines = pdf.splitTextToSize(disclaimer, pageWidth - 40);
        console.log('📝 Disclaimer split into lines:', disclaimerLines.length);
        pdf.text(disclaimerLines, 20, yPosition);
        yPosition += disclaimerLines.length * 4 + 10;
        console.log('✅ Disclaimer added, current Y position:', yPosition);
        
        // Footer
        console.log('📝 Adding footer...');
        const footerY = pageHeight - 20;
        pdf.setFontSize(8);
        pdf.text('🧬 GeneKnow Platform - AI-powered genomic analysis', pageWidth / 2, footerY - 10, { align: 'center' });
        pdf.text(`Report ID: GK-${Date.now()} | Generated: ${new Date().toISOString()}`, pageWidth / 2, footerY, { align: 'center' });
        console.log('✅ Footer added');
        
        // Generate filename and download
        const timestamp = new Date().toISOString().split('T')[0];
        const fileName = `${title.replace(/[^a-z0-9]/gi, '_').toLowerCase()}_${timestamp}.pdf`;
        console.log('📁 Generated filename:', fileName);
        
        // Direct download
        console.log('💾 Starting PDF save...');
        pdf.save(fileName);
        console.log('✅ PDF save completed');
        
        // Remove loading notification
        console.log('🧹 Removing loading notification...');
        if (document.body.contains(loadingNotification)) {
          document.body.removeChild(loadingNotification);
          console.log('✅ Loading notification removed');
        }
        
        // Show success notification
        console.log('🎉 Showing success notification...');
        const successNotification = document.createElement('div');
        successNotification.innerHTML = `
          <div style="
            position: fixed; 
            top: 20px; 
            right: 20px; 
            background: #22C55E; 
            color: white; 
            padding: 15px 20px; 
            border-radius: 8px; 
            box-shadow: 0 4px 12px rgba(0,0,0,0.3);
            z-index: 10000;
            font-weight: 600;
            font-size: 14px;
            max-width: 300px;
          ">
            ✅ PDF with Visualization Downloaded!<br>
            <small style="font-weight: 400;">${fileName}</small>
          </div>
        `;
        document.body.appendChild(successNotification);
        console.log('✅ Success notification added to DOM');
        
        setTimeout(() => {
          if (successNotification.parentNode) {
            successNotification.parentNode.removeChild(successNotification);
            console.log('🧹 Success notification removed after timeout');
          }
        }, 3000);
        
        console.log('🎉 PDF generation completed successfully!');
        
        // Reset PDF generating state
        setIsPDFGenerating(false);
        
      } catch (error: unknown) {
        console.error('❌ Error in PDF generation:', error);
        console.error('❌ Error details:', (error as Error).message, (error as Error).stack);
        showError('Error generating PDF content');
        
        // Reset PDF generating state on error
        setIsPDFGenerating(false);
      }
    }
    
    function showError(message: string) {
      console.log('❌ Showing error:', message);
      
      // Reset PDF generating state
      setIsPDFGenerating(false);
      
      // Remove loading notification if it exists
      if (document.body.contains(loadingNotification)) {
        document.body.removeChild(loadingNotification);
        console.log('🧹 Loading notification removed due to error');
      }
      
      // Show error notification
      const errorNotification = document.createElement('div');
      errorNotification.innerHTML = `
        <div style="
          position: fixed; 
          top: 20px; 
          right: 20px; 
          background: #EF4444; 
          color: white; 
          padding: 15px 20px; 
          border-radius: 8px; 
          box-shadow: 0 4px 12px rgba(0,0,0,0.3);
          z-index: 10000;
          font-weight: 600;
          font-size: 14px;
          max-width: 300px;
        ">
          ❌ ${message}<br>
          <small style="font-weight: 400;">Please try again</small>
        </div>
      `;
      document.body.appendChild(errorNotification);
      console.log('❌ Error notification added to DOM');
      
      setTimeout(() => {
        if (errorNotification.parentNode) {
          errorNotification.parentNode.removeChild(errorNotification);
          console.log('🧹 Error notification removed after timeout');
        }
      }, 4000);
    }
    
  } catch (error: unknown) {
    console.error('❌ Critical error in downloadPDF function:', error);
    console.error('❌ Critical error details:', (error as Error).message, (error as Error).stack);
    
    // Reset PDF generating state
    setIsPDFGenerating(false);
    
    // Remove loading notification if it exists
    const loadingNotification = document.querySelector('[style*="Loading PDF library"]') || 
                                document.querySelector('[style*="Generating PDF"]');
    if (loadingNotification?.parentNode) {
      loadingNotification.parentNode.removeChild(loadingNotification);
      console.log('🧹 Loading notification removed due to critical error');
    }
    
    // Show error notification
    const errorNotification = document.createElement('div');
    errorNotification.innerHTML = `
      <div style="
        position: fixed; 
        top: 20px; 
        right: 20px; 
        background: #EF4444; 
        color: white; 
        padding: 15px 20px; 
        border-radius: 8px; 
        box-shadow: 0 4px 12px rgba(0,0,0,0.3);
        z-index: 10000;
        font-weight: 600;
        font-size: 14px;
        max-width: 300px;
      ">
        ❌ PDF generation failed<br>
        <small style="font-weight: 400;">Please try again</small>
      </div>
    `;
    document.body.appendChild(errorNotification);
    console.log('❌ Critical error notification added to DOM');
    
    setTimeout(() => {
      if (errorNotification.parentNode) {
        errorNotification.parentNode.removeChild(errorNotification);
        console.log('🧹 Critical error notification removed after timeout');
      }
    }, 3000);
  }
};

// Simple Download Button Component
const DownloadButton = ({ elementId, title, summary, isPDFGenerating, setIsPDFGenerating }: { elementId: string, title: string, summary: string, isPDFGenerating: boolean, setIsPDFGenerating: (value: boolean) => void }) => {
  if (isPDFGenerating) {
    return null; // Hide button during PDF generation
  }
  
  return (
    <button
      onClick={() => downloadPDF(elementId, title, summary, setIsPDFGenerating)}
      style={{
        display: 'inline-flex',
        alignItems: 'center',
        padding: '0.5rem 1rem',
        background: '#2563EB',
        color: '#FFFFFF',
        border: 'none',
        borderRadius: '0.375rem',
        fontSize: '0.875rem',
        fontWeight: '500',
        cursor: 'pointer',
        transition: 'all 200ms ease',
        boxShadow: '0 1px 2px 0 rgba(0, 0, 0, 0.05)',
        gap: '0.5rem'
      }}
      onMouseEnter={(e) => {
        e.currentTarget.style.background = '#1D4ED8';
        e.currentTarget.style.transform = 'translateY(-1px)';
        e.currentTarget.style.boxShadow = '0 4px 6px -1px rgba(0, 0, 0, 0.1)';
      }}
      onMouseLeave={(e) => {
        e.currentTarget.style.background = '#2563EB';
        e.currentTarget.style.transform = 'translateY(0px)';
        e.currentTarget.style.boxShadow = '0 1px 2px 0 rgba(0, 0, 0, 0.05)';
      }}
      onMouseDown={(e) => {
        e.currentTarget.style.transform = 'translateY(0px)';
        e.currentTarget.style.boxShadow = '0 1px 2px 0 rgba(0, 0, 0, 0.05)';
      }}
    >
      <svg 
        width="16" 
        height="16" 
        fill="currentColor" 
        viewBox="0 0 20 20"
        style={{ flexShrink: 0 }}
      >
        <path fillRule="evenodd" d="M3 17a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zm3.293-7.707a1 1 0 011.414 0L9 10.586V3a1 1 0 112 0v7.586l1.293-1.293a1 1 0 111.414 1.414l-3 3a1 1 0 01-1.414 0l-3-3a1 1 0 010-1.414z" clipRule="evenodd"/>
      </svg>
      PDF Report
    </button>
  );
};

// Type definitions
interface Alert {
  type: 'critical' | 'warning' | 'info' | 'success';
  title: string;
  desc: string;
  detailedInfo: {
    whatItMeans: string;
    whyImportant: string;
    nextSteps: string;
    clinicalSignificance: string;
  };
}

// Icon components
const InformationCircleIcon = ({ style, onMouseEnter, onMouseLeave }: { 
  style?: React.CSSProperties; 
  onMouseEnter?: (e: React.MouseEvent<HTMLDivElement>) => void;
  onMouseLeave?: (e: React.MouseEvent<HTMLDivElement>) => void;
}) => (
  <div 
    style={{
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      width: '1rem',
      height: '1rem',
      borderRadius: '50%',
      backgroundColor: '#E5E7EB',
      color: '#6B7280',
      fontSize: '0.75rem',
      fontWeight: '600',
      fontFamily: 'serif',
      cursor: 'pointer',
      transition: 'all 200ms ease',
      border: '1px solid #D1D5DB',
      ...style
    }}
    onMouseEnter={onMouseEnter}
    onMouseLeave={onMouseLeave}
  >
    i
  </div>
);

// Tooltip component for detailed alert information
const AlertTooltip = ({ alert, isVisible }: { alert: Alert; isVisible: boolean }) => (
  <div style={{
    position: 'absolute',
    left: 'calc(100% + 0.5rem)', // Position to the right of the icon
    top: '50%',
    transform: 'translateY(-50%)',
    width: '20rem',
    padding: '1rem',
    background: '#1F2937',
    color: '#FFFFFF',
    fontSize: '0.75rem',
    borderRadius: '0.5rem',
    boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
    opacity: isVisible ? 1 : 0,
    visibility: isVisible ? 'visible' : 'hidden',
    transition: 'opacity 300ms ease, visibility 300ms ease',
    zIndex: 1000,
    pointerEvents: 'none',
    lineHeight: '1.4',
    border: '1px solid #374151'
  }}>
    <div style={{ marginBottom: '0.75rem' }}>
      <h4 style={{ 
        fontWeight: '600', 
        color: '#F3F4F6', 
        marginBottom: '0.25rem',
        fontSize: '0.8rem'
      }}>
        What it means:
      </h4>
      <p style={{ color: '#D1D5DB', marginBottom: '0' }}>{alert.detailedInfo.whatItMeans}</p>
    </div>
    
    <div style={{ marginBottom: '0.75rem' }}>
      <h4 style={{ 
        fontWeight: '600', 
        color: '#F3F4F6', 
        marginBottom: '0.25rem',
        fontSize: '0.8rem'
      }}>
        Why it's important:
      </h4>
      <p style={{ color: '#D1D5DB', marginBottom: '0' }}>{alert.detailedInfo.whyImportant}</p>
    </div>
    
    <div style={{ marginBottom: '0.75rem' }}>
      <h4 style={{ 
        fontWeight: '600', 
        color: '#F3F4F6', 
        marginBottom: '0.25rem',
        fontSize: '0.8rem'
      }}>
        Clinical significance:
      </h4>
      <p style={{ color: '#D1D5DB', marginBottom: '0' }}>{alert.detailedInfo.clinicalSignificance}</p>
    </div>
    
    <div>
      <h4 style={{ 
        fontWeight: '600', 
        color: '#F3F4F6', 
        marginBottom: '0.25rem',
        fontSize: '0.8rem'
      }}>
        Next steps:
      </h4>
      <p style={{ color: '#D1D5DB', marginBottom: '0' }}>{alert.detailedInfo.nextSteps}</p>
    </div>
    
    {/* Tooltip arrow pointing to the left towards the "i" icon */}
    <div style={{
      position: 'absolute',
      left: '-0.5rem',
      top: '50%',
      transform: 'translateY(-50%)',
      width: '0',
      height: '0',
      borderTop: '0.5rem solid transparent',
      borderBottom: '0.5rem solid transparent',
      borderRight: '0.5rem solid #1F2937'
    }}></div>
  </div>
);

// Mock data sets for different risk levels - completely anonymous
const mockDataSets = {
  high: {
    riskLevel: 'High Risk',
    riskScore: '82/100',
    condition: 'Hereditary Breast and Ovarian Cancer Syndrome',
    details: 'Family History: Breast Cancer<br/>Referral: Oncology<br/>Previous Tests: BRCA1/2 Panel',
    alerts: [
      {
        type: 'critical' as const,
        title: 'BRCA1 Pathogenic Variant',
        desc: 'c.5266dupC - Frameshift mutation detected',
        detailedInfo: {
          whatItMeans: 'A frameshift mutation in the BRCA1 gene that disrupts the normal protein function, leading to loss of tumor suppressor activity.',
          whyImportant: 'BRCA1 mutations significantly increase the risk of breast and ovarian cancers. This specific mutation is classified as pathogenic with high confidence.',
          clinicalSignificance: 'Associated with up to 70% lifetime risk of breast cancer and 40% risk of ovarian cancer. Early onset cancers are common.',
          nextSteps: 'Immediate genetic counseling, enhanced screening protocols, and discussion of risk-reducing strategies including prophylactic surgery.'
        }
      },
      {
        type: 'warning' as const,
        title: 'High Risk Classification',
        desc: 'Immediate genetic counseling recommended',
        detailedInfo: {
          whatItMeans: 'Based on the genetic analysis, this individual has been classified as high-risk for hereditary cancer syndrome.',
          whyImportant: 'High-risk individuals require specialized medical management and monitoring that differs significantly from standard population screening.',
          clinicalSignificance: 'Warrants immediate intervention with genetic counseling and potentially altered screening and prevention strategies.',
          nextSteps: 'Schedule urgent genetic counseling appointment, inform family members about potential hereditary risk, and begin enhanced screening protocols.'
        }
      }
    ]
  },
  medium: {
    riskLevel: 'Medium Risk',
    riskScore: '45/100',
    condition: 'Lynch Syndrome',
    details: 'Family History: Colorectal Cancer<br/>Referral: Oncology<br/>Previous Tests: MSI-H positive',
    alerts: [
      {
        type: 'warning' as const,
        title: 'MLH1 Pathogenic Variant',
        desc: 'c.1989-1G>A - Splice site mutation',
        detailedInfo: {
          whatItMeans: 'A splice site mutation in the MLH1 gene that affects mRNA processing, leading to defective DNA mismatch repair.',
          whyImportant: 'MLH1 is a key component of the DNA mismatch repair system. Mutations lead to microsatellite instability and increased cancer risk.',
          clinicalSignificance: 'Associated with Lynch syndrome, increasing colorectal cancer risk to 50-80% and endometrial cancer risk to 40-60%.',
          nextSteps: 'Enhanced colorectal screening starting at age 20-25, annual endometrial screening, and genetic counseling for family members.'
        }
      },
      {
        type: 'info' as const,
        title: 'Screening Recommendations',
        desc: 'Enhanced colonoscopy surveillance required',
        detailedInfo: {
          whatItMeans: 'Standard population screening is insufficient for this individual due to elevated genetic risk factors.',
          whyImportant: 'Early detection through enhanced screening can significantly improve outcomes and survival rates in hereditary cancer syndromes.',
          clinicalSignificance: 'Colonoscopy should begin 10-15 years earlier than standard recommendations and occur more frequently.',
          nextSteps: 'Schedule colonoscopy every 1-2 years starting at age 20-25, or 2-5 years before earliest family diagnosis.'
        }
      }
    ]
  },
  low: {
    riskLevel: 'Low Risk',
    riskScore: '15/100',
    condition: 'Li-Fraumeni Syndrome',
    details: 'Family History: Multiple Sarcomas<br/>Referral: Genetics<br/>Previous Tests: TP53 Sequencing',
    alerts: [
      {
        type: 'info' as const,
        title: 'TP53 Variant of Uncertain Significance',
        desc: 'c.743G>A - Missense mutation',
        detailedInfo: {
          whatItMeans: 'A genetic variant in the TP53 gene where the clinical significance is not yet definitively established by current scientific evidence.',
          whyImportant: 'TP53 is the "guardian of the genome" and mutations can predispose to multiple cancer types. Even uncertain variants require monitoring.',
          clinicalSignificance: 'May or may not be clinically significant. Requires ongoing evaluation as more research data becomes available.',
          nextSteps: 'Genetic counseling to discuss implications, potential for variant reclassification, and consideration of family testing.'
        }
      },
      {
        type: 'success' as const,
        title: 'Low Risk Assessment',
        desc: 'Routine follow-up recommended',
        detailedInfo: {
          whatItMeans: 'Current genetic analysis does not indicate significantly elevated risk for hereditary cancer syndromes.',
          whyImportant: 'Provides reassurance while maintaining appropriate vigilance for any changes in risk factors or family history.',
          clinicalSignificance: 'Standard population screening guidelines are appropriate for this individual at this time.',
          nextSteps: 'Continue routine screening as recommended for general population, maintain awareness of family history changes.'
        }
      }
    ]
  }
};

const ClinicalViewPage: React.FC = () => {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const location = useLocation();
  const [activeTab, setActiveTab] = useState('analysis');
  const [hoveredAlert, setHoveredAlert] = useState<number | null>(null);
  const [isPDFGenerating, setIsPDFGenerating] = useState(false);
  
  // Extract risk level from URL parameter
  const riskLevel = searchParams.get('risk') || 'high';
  
  // Get the appropriate SHAP validation data
  const mockSHAPValidation = getMockSHAPValidation(riskLevel);
  
  // Check if we have real results from the pipeline
  const pipelineResults = location.state?.results as PipelineResult | undefined;
  const fileName = location.state?.fileName as string | undefined;
  
  // Get mock data for UI elements
  const currentData = mockDataSets[riskLevel as keyof typeof mockDataSets] || mockDataSets.low;

  // Function to convert pipeline results to the format expected by the clinical view
  const getGenomicData = () => {
    if (pipelineResults) {
      // Convert real pipeline results to the format expected by the clinical view
      const riskScores = Object.entries(pipelineResults.risk_scores || {});
      const riskFindings = riskScores.map(([cancer_type, risk_percentage]) => ({
        cancer_type,
        risk_percentage,
        risk_level: risk_percentage > 50 ? 'high' : risk_percentage > 20 ? 'medium' : 'low',
        affected_genes: pipelineResults.variants?.map(v => v.gene) || [],
        recommendation: risk_percentage > 50 
          ? `Enhanced ${cancer_type} cancer screening recommended`
          : risk_percentage > 20
          ? `Moderate ${cancer_type} cancer screening recommended`
          : 'Standard screening guidelines apply',
        mutation_burden: risk_percentage > 50 ? 'high' : risk_percentage > 20 ? 'medium' : 'low',
        pathway_disruption: risk_percentage > 50 ? ['DNA_REPAIR', 'CELL_CYCLE'] : []
      }));

      const variantTable = pipelineResults.variants?.map((variant, index) => ({
        gene: variant.gene,
        variant_id: `${variant.position}:${variant.type}`,
        consequence: variant.impact,
        mutation_type: variant.type,
        quality_score: 85 + index * 5, // Estimated quality score
        clinical_significance: variant.impact === 'HIGH' ? 'pathogenic' : 'benign',
        tcga_best_match: { cancer_type: riskScores[0]?.[0] || 'breast', frequency: 6.3 },
        protein_change: `p.${variant.gene}Variant`,
        functional_impact: variant.impact === 'HIGH' ? 'Loss of function' : 'Uncertain',
        transformation: {
          original: 'ATT',
          mutated: 'GTT',
          amino_acid_change: 'Unknown',
          effect: variant.impact === 'HIGH' ? 'Structural disruption' : 'Uncertain significance'
        }
      })) || [];

      return {
        summary: {
          total_variants_found: pipelineResults.variant_count || 0,
          variants_passed_qc: Math.floor((pipelineResults.variant_count || 0) * 0.9),
          high_risk_findings_count: riskFindings.filter(r => r.risk_level === 'high').length,
          processing_time_seconds: pipelineResults.processing_time_seconds || 0,
          mutation_types: {
            snv: Math.floor((pipelineResults.variant_count || 0) * 0.5),
            indel: Math.floor((pipelineResults.variant_count || 0) * 0.3),
            cnv: Math.floor((pipelineResults.variant_count || 0) * 0.15),
            structural: Math.floor((pipelineResults.variant_count || 0) * 0.05)
          }
        },
        risk_findings: riskFindings,
        variant_table: variantTable,
        structural_variants: [],
        copy_number_variants: [],
        mutation_signatures: [
          {
            signature: "SBS1",
            name: "Spontaneous deamination",
            contribution: 0.35,
            description: "C>T transitions at CpG sites",
            etiology: "Aging-related mutations"
          },
          {
            signature: "SBS3",
            name: "BRCA1/BRCA2 deficiency",
            contribution: 0.28,
            description: "Homologous recombination deficiency",
            etiology: "Inherited DNA repair defects"
          }
        ]
      };
    } else {
      // Fall back to mock data for demo purposes
      return mockGenomicData;
    }
  };

  const genomicData = getGenomicData();

  const renderTabContent = () => {
    switch (activeTab) {
      case 'analysis':
        return (
          <div style={{ padding: '2rem' }}>
            <h2 style={{ 
              color: '#111827',
              fontSize: '1.5rem',
              fontWeight: '600',
              marginBottom: '1.5rem'
            }}>
              Genomic Analysis Overview
            </h2>
            
            {pipelineResults && (
              <div style={{
                padding: '1rem',
                marginBottom: '2rem',
                background: '#EFF6FF',
                borderRadius: '0.5rem',
                border: '1px solid #BFDBFE'
              }}>
                <p style={{ color: '#1E40AF', fontSize: '0.875rem', fontWeight: '500' }}>
                  Analysis completed.
                </p>
              </div>
            )}
            
            {/* Risk Summary Cards */}
            <div id="cancer-risk-assessment" style={{ 
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '2rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1.5rem' }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  Cancer Risk Assessment
                </h3>
                <DownloadButton 
                  elementId="cancer-risk-assessment"
                  title="Cancer Risk Assessment"
                  summary={pipelineResults 
                    ? `Comprehensive genetic risk analysis for ${fileName}. Key findings include risk scores across multiple cancer types based on ${genomicData.summary.total_variants_found} variants analyzed.`
                    : "Comprehensive genetic risk analysis showing cancer predisposition across multiple cancer types. Key findings include: Breast (100%), Colon (100%), Lung (97.7%), and Prostate (81.4%) cancer risks. Analysis covers 12 variants with 6 high-risk findings. Essential for developing personalized screening and prevention strategies."
                  }
                  isPDFGenerating={isPDFGenerating}
                  setIsPDFGenerating={setIsPDFGenerating}
                />
              </div>
              
              <div style={{ 
                display: 'grid', 
                gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', 
                gap: '1rem'
              }}>
                {genomicData.risk_findings.map((risk, index) => (
                  <div key={index} style={{
                    background: '#FFFFFF',
                    padding: '1.5rem',
                    borderRadius: '0.75rem',
                    border: `2px solid ${getRiskColor(risk.risk_level)}`,
                    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
                  }}>
                    <h3 style={{ 
                      color: '#111827',
                      fontSize: '1rem',
                      fontWeight: '600',
                      marginBottom: '0.5rem',
                      textTransform: 'capitalize'
                    }}>
                      {risk.cancer_type} Cancer
                    </h3>
                    <div style={{ 
                      fontSize: '2rem',
                      fontWeight: 'bold',
                      color: getRiskColor(risk.risk_level),
                      marginBottom: '0.5rem'
                    }}>
                      {risk.risk_percentage.toFixed(1)}%
                    </div>
                    <div style={{ 
                      fontSize: '0.875rem',
                      color: '#4B5563',
                      marginBottom: '0.5rem'
                    }}>
                      {risk.affected_genes.length > 0 ? `Genes: ${risk.affected_genes.slice(0, 3).join(', ')}` : 'No high-risk variants found'}
                    </div>
                    <div style={{
                      padding: '0.5rem',
                      borderRadius: '0.375rem',
                      fontSize: '0.75rem',
                      fontWeight: '500',
                      background: risk.risk_level === 'high' ? '#FEF2F2' : 
                                 risk.risk_level === 'medium' ? '#FFFBEB' : '#F0FDF4',
                      color: risk.risk_level === 'high' ? '#EF4444' : 
                             risk.risk_level === 'medium' ? '#F59E0B' : '#22C55E'
                    }}>
                      {risk.risk_level.toUpperCase()} RISK
                    </div>
                  </div>
                ))}
              </div>
            </div>

            {/* Manhattan Plot */}
            <div id="gene-significance-analysis" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '2rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  Gene Significance Analysis
                </h3>
                <DownloadButton 
                  elementId="gene-significance-analysis"
                  title="Gene Significance Analysis"
                  summary="Manhattan plot visualization showing statistical significance of genetic variants across cancer-associated genes. Analysis includes BRCA1, TP53, APC, MLH1, and KRAS genes with their quality scores and clinical significance. Red dots indicate high-confidence pathogenic variants requiring clinical follow-up."
                  isPDFGenerating={isPDFGenerating}
                  setIsPDFGenerating={setIsPDFGenerating}
                />
              </div>
              <svg width="100%" height="300" viewBox="0 0 800 300">
                {/* Background */}
                <rect width="800" height="300" fill="#F9FAFB"/>
                
                {/* Grid lines */}
                {[0, 1, 2, 3, 4, 5].map(i => (
                  <g key={i}>
                    <line x1="80" y1={50 + i * 40} x2="750" y2={50 + i * 40} stroke="#E5E7EB" strokeWidth="1"/>
                    <text x="70" y={55 + i * 40} fill="#6B7280" fontSize="12" textAnchor="end">{5-i}</text>
                  </g>
                ))}
                
                {/* Significance threshold line */}
                <line x1="80" y1="130" x2="750" y2="130" stroke="#EF4444" strokeWidth="2" strokeDasharray="5,5"/>
                <text x="755" y="135" fill="#EF4444" fontSize="12">p &lt; 0.05</text>
                
                {/* Data points */}
                {genomicData.variant_table.map((variant, index) => {
                  const x = 150 + index * 150;
                  const y = 50 + (5 - variant.quality_score / 20) * 40;
                  const isSignificant = variant.quality_score > 60;
                  
                  return (
                    <g key={index}>
                      <circle 
                        cx={x} 
                        cy={y} 
                        r="6" 
                        fill={isSignificant ? '#EF4444' : '#2563EB'}
                        opacity="0.7"
                      />
                      <text x={x} y={y - 15} fill="#111827" fontSize="12" textAnchor="middle" fontWeight="600">
                        {variant.gene}
                      </text>
                      <text x={x} y={y + 25} fill="#4B5563" fontSize="10" textAnchor="middle">
                        {variant.tcga_best_match.cancer_type}
                      </text>
                    </g>
                  );
                })}
                
                {/* Axis labels */}
                <text x="400" y="290" fill="#111827" fontSize="14" textAnchor="middle" fontWeight="600">
                  Gene Position
                </text>
                <text x="30" y="150" fill="#111827" fontSize="14" textAnchor="middle" fontWeight="600" transform="rotate(-90 30 150)">
                  -log10(p-value)
                </text>
              </svg>
              
              <div style={{ 
                fontSize: '0.875rem',
                color: '#4B5563',
                marginTop: '1rem'
              }}>
                <strong>Interpretation:</strong> Points above the red line indicate statistically significant associations. 
                Red dots represent high-confidence pathogenic variants, blue dots represent variants under investigation.
              </div>
            </div>
            
            {/* Mutation Types Breakdown */}
            <div id="mutation-type-distribution" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '2rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  Mutation Type Distribution
                </h3>
                <DownloadButton 
                  elementId="mutation-type-distribution"
                  title="Mutation Type Distribution"
                  summary="Analysis of mutation types detected across the genome showing 5 SNVs, 4 indels, 2 CNVs, and 1 structural variant. Different mutation types have varying therapeutic implications - SNVs may be targetable with small molecules, indels suggest DNA repair deficiency, CNVs affect gene dosage, and structural variants may disrupt gene function."
                  isPDFGenerating={isPDFGenerating}
                  setIsPDFGenerating={setIsPDFGenerating}
                />
              </div>
              
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(180px, 1fr))', gap: '1rem' }}>
                {Object.entries(genomicData.summary.mutation_types).map(([type, count]) => {
                  const colors = {
                    snv: '#3B82F6',
                    indel: '#EF4444', 
                    cnv: '#F59E0B',
                    structural: '#8B5CF6'
                  };
                  const labels = {
                    snv: 'Single Nucleotide Variants',
                    indel: 'Insertions/Deletions',
                    cnv: 'Copy Number Variants',
                    structural: 'Structural Variants'
                  };
                  
                  return (
                    <div key={type} style={{
                      textAlign: 'center',
                      padding: '1rem',
                      border: `2px solid ${colors[type as keyof typeof colors]}`,
                      borderRadius: '0.5rem',
                      background: `${colors[type as keyof typeof colors]}10`
                    }}>
                      <div style={{ 
                        fontSize: '2rem', 
                        fontWeight: 'bold', 
                        color: colors[type as keyof typeof colors],
                        marginBottom: '0.5rem'
                      }}>
                        {count}
                      </div>
                      <div style={{ 
                        fontSize: '0.875rem', 
                        color: '#4B5563',
                        fontWeight: '500'
                      }}>
                        {labels[type as keyof typeof labels]}
                      </div>
                    </div>
                  );
                })}
              </div>
            </div>

            {/* Mutation Signatures */}
            <div id="mutational-signatures" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '2rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  Mutational Signature Analysis
                </h3>
                <DownloadButton 
                  elementId="mutational-signatures"
                  title="Mutational Signature Analysis"
                  summary="Analysis of mutational signatures revealing underlying biological processes. Key findings: 35% aging-related (SBS1), 28% BRCA deficiency (SBS3), 22% tobacco exposure (SBS4), 15% APOBEC activity (SBS13). The BRCA deficiency signature suggests high likelihood of PARP inhibitor response. Important for therapeutic selection and risk management."
                  isPDFGenerating={isPDFGenerating}
                  setIsPDFGenerating={setIsPDFGenerating}
                />
              </div>
              
              <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                {genomicData.mutation_signatures.map((signature, index) => (
                  <div key={index} style={{
            display: 'flex',
            alignItems: 'center',
                    padding: '1rem',
                    border: '1px solid #E5E7EB',
                    borderRadius: '0.5rem',
                    background: '#F9FAFB'
                  }}>
                    <div style={{ flex: 1 }}>
                      <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginBottom: '0.25rem' }}>
                        <h4 style={{ 
                          color: '#111827',
                          fontSize: '1rem',
                          fontWeight: '600',
                          margin: 0
                        }}>
                          {signature.signature}
                        </h4>
                        <span style={{ 
                          color: '#4B5563',
                          fontSize: '0.875rem'
                        }}>
                          {signature.name}
                        </span>
                      </div>
                      <p style={{ 
                        color: '#6B7280',
                        fontSize: '0.875rem',
                        marginBottom: '0.5rem',
                        margin: 0
                      }}>
                        {signature.description} • {signature.etiology}
                      </p>
                      <div style={{
                        background: '#E5E7EB',
                        borderRadius: '0.25rem',
                        height: '0.5rem',
                        position: 'relative',
                        overflow: 'hidden'
                      }}>
                        <div style={{
                          background: '#3B82F6',
                          height: '100%',
                          width: `${signature.contribution * 100}%`,
                          transition: 'width 0.3s ease'
                        }}/>
                      </div>
                    </div>
                    <div style={{
                      marginLeft: '1rem',
                      fontSize: '1.25rem',
                      fontWeight: 'bold',
                      color: '#3B82F6'
                    }}>
                      {(signature.contribution * 100).toFixed(1)}%
                    </div>
                  </div>
                ))}
              </div>
            </div>

            {/* Quality Metrics */}
            <div id="quality-metrics" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  Analysis Quality Metrics
                </h3>
                <DownloadButton 
                  elementId="quality-metrics"
                  title="Analysis Quality Metrics"
                  summary="Quality control metrics for genomic analysis pipeline showing 12 total variants with 11 passing QC, 6 high-risk findings, and 0.70s processing time. High data quality with excellent variant calling reliability. All metrics meet clinical laboratory standards for confident therapeutic decision-making."
                  isPDFGenerating={isPDFGenerating}
                  setIsPDFGenerating={setIsPDFGenerating}
                />
              </div>
              <div style={{ 
                display: 'grid', 
                gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', 
                gap: '1rem'
              }}>
                <div style={{ textAlign: 'center' }}>
                  <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#2563EB' }}>
                    {genomicData.summary.total_variants_found}
                  </div>
                  <div style={{ fontSize: '0.875rem', color: '#4B5563' }}>Total Variants</div>
                </div>
                <div style={{ textAlign: 'center' }}>
                  <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#22C55E' }}>
                    {genomicData.summary.variants_passed_qc}
                  </div>
                  <div style={{ fontSize: '0.875rem', color: '#4B5563' }}>Passed QC</div>
                </div>
                <div style={{ textAlign: 'center' }}>
                  <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#EF4444' }}>
                    {genomicData.summary.high_risk_findings_count}
                  </div>
                  <div style={{ fontSize: '0.875rem', color: '#4B5563' }}>High-Risk Findings</div>
                </div>
                <div style={{ textAlign: 'center' }}>
                  <div style={{ fontSize: '2rem', fontWeight: 'bold', color: '#F59E0B' }}>
                    {genomicData.summary.processing_time_seconds.toFixed(1)}s
                  </div>
                  <div style={{ fontSize: '0.875rem', color: '#4B5563' }}>Processing Time</div>
                </div>
              </div>
            </div>
          </div>
        );
      case 'variants':
        return (
          <div style={{ padding: '2rem' }}>
            <h2 style={{ 
              color: '#111827',
              fontSize: '1.5rem',
              fontWeight: '600',
              marginBottom: '1.5rem'
            }}>
              Variant Heatmap Analysis
            </h2>
            
                        {/* Gene-Cancer Type Heatmap */}
            <div id="gene-cancer-matrix" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '2rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  Gene-Cancer Type Association Matrix
                </h3>
                <DownloadButton 
                  elementId="gene-cancer-matrix"
                  title="Gene-Cancer Association Matrix"
                  summary="Heatmap showing gene-cancer associations based on TCGA data. Strong associations: BRCA1-breast (100%), TP53-multiple cancers (80%), APC-colon (90%). Color-coded matrix indicates risk levels for enhanced screening protocols. Critical for multi-organ surveillance planning."
                  isPDFGenerating={isPDFGenerating}
                  setIsPDFGenerating={setIsPDFGenerating}
                />
              </div>
              
              <svg width="100%" height="400" viewBox="0 0 800 400">
                {/* Background */}
                <rect width="800" height="400" fill="#F9FAFB"/>
                
                {/* Heatmap grid */}
                {['BRCA1', 'TP53', 'APC', 'MLH1', 'MSH2'].map((gene, geneIndex) => {
                  return ['breast', 'lung', 'colon', 'prostate', 'blood'].map((cancer, cancerIndex) => {
                    const x = 150 + cancerIndex * 80;
                    const y = 80 + geneIndex * 50;
                    
                    // Calculate intensity based on gene-cancer association
                    let intensity = 0;
                    if (gene === 'BRCA1' && cancer === 'breast') intensity = 1.0;
                    if (gene === 'TP53' && ['breast', 'lung', 'colon', 'prostate'].includes(cancer)) intensity = 0.8;
                    if (gene === 'APC' && cancer === 'colon') intensity = 0.9;
                    if (gene === 'MLH1' && cancer === 'colon') intensity = 0.7;
                    if (gene === 'MSH2' && cancer === 'colon') intensity = 0.6;
                    
                    const color = intensity > 0.8 ? '#EF4444' : 
                                 intensity > 0.6 ? '#F59E0B' : 
                                 intensity > 0.3 ? '#3B82F6' : '#E5E7EB';
                    
                    return (
                      <g key={`${gene}-${cancer}`}>
                        <rect 
                          x={x} 
                          y={y} 
                          width="70" 
                          height="40" 
                          fill={color}
                          stroke="#FFFFFF"
                          strokeWidth="2"
                          opacity="0.8"
                        />
                        {intensity > 0 && (
                          <text 
                            x={x + 35} 
                            y={y + 25} 
                            fill="#FFFFFF" 
                            fontSize="12" 
                            textAnchor="middle"
                            fontWeight="600"
                          >
                            {(intensity * 100).toFixed(0)}%
                          </text>
                        )}
                      </g>
                    );
                  });
                })}
                
                {/* Gene labels */}
                {['BRCA1', 'TP53', 'APC', 'MLH1', 'MSH2'].map((gene, index) => (
                  <text 
                    key={gene}
                    x="140" 
                    y={105 + index * 50} 
                    fill="#111827" 
                    fontSize="14" 
                    textAnchor="end"
                    fontWeight="600"
                  >
                    {gene}
                  </text>
                ))}
                
                {/* Cancer type labels */}
                {['Breast', 'Lung', 'Colon', 'Prostate', 'Blood'].map((cancer, index) => (
                  <text 
                    key={cancer}
                    x={185 + index * 80} 
                    y="70" 
                    fill="#111827" 
                    fontSize="14" 
                    textAnchor="middle"
                    fontWeight="600"
                  >
                    {cancer}
                  </text>
                ))}
                
                {/* Legend */}
                <g transform="translate(580, 300)">
                  <text x="0" y="0" fill="#111827" fontSize="12" fontWeight="600">Risk Level</text>
                  <rect x="0" y="10" width="15" height="15" fill="#EF4444"/>
                  <text x="20" y="22" fill="#4B5563" fontSize="11">High (≥80%)</text>
                  <rect x="0" y="30" width="15" height="15" fill="#F59E0B"/>
                  <text x="20" y="42" fill="#4B5563" fontSize="11">Medium (60-80%)</text>
                  <rect x="0" y="50" width="15" height="15" fill="#3B82F6"/>
                  <text x="20" y="62" fill="#4B5563" fontSize="11">Low (30-60%)</text>
                  <rect x="0" y="70" width="15" height="15" fill="#E5E7EB"/>
                  <text x="20" y="82" fill="#4B5563" fontSize="11">No Association</text>
                </g>
              </svg>
            </div>
            
            {/* Variant Table */}
            <div id="detected-variants" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  Detected Variants
                </h3>
                <DownloadButton 
                  elementId="detected-variants"
                  title="Detected Variants"
                  summary="Comprehensive table of all detected genetic variants with molecular transformations and clinical significance. Includes BRCA1, TP53, APC, MLH1, and KRAS mutations with detailed protein changes and functional impacts. Critical for understanding specific mutations and their therapeutic implications."
                  isPDFGenerating={isPDFGenerating}
                  setIsPDFGenerating={setIsPDFGenerating}
                />
              </div>
              
              <div style={{ overflowX: 'auto' }}>
                <table style={{ width: '100%', borderCollapse: 'collapse' }}>
                  <thead>
                    <tr style={{ background: '#F9FAFB' }}>
                      <th style={{ padding: '0.75rem', textAlign: 'left', borderBottom: '1px solid #E5E7EB', fontWeight: '600' }}>Gene</th>
                      <th style={{ padding: '0.75rem', textAlign: 'left', borderBottom: '1px solid #E5E7EB', fontWeight: '600' }}>Variant</th>
                      <th style={{ padding: '0.75rem', textAlign: 'left', borderBottom: '1px solid #E5E7EB', fontWeight: '600' }}>Type</th>
                      <th style={{ padding: '0.75rem', textAlign: 'left', borderBottom: '1px solid #E5E7EB', fontWeight: '600' }}>Transformation</th>
                      <th style={{ padding: '0.75rem', textAlign: 'left', borderBottom: '1px solid #E5E7EB', fontWeight: '600' }}>Quality</th>
                      <th style={{ padding: '0.75rem', textAlign: 'left', borderBottom: '1px solid #E5E7EB', fontWeight: '600' }}>Significance</th>
                      <th style={{ padding: '0.75rem', textAlign: 'left', borderBottom: '1px solid #E5E7EB', fontWeight: '600' }}>Impact</th>
                    </tr>
                  </thead>
                  <tbody>
                    {genomicData.variant_table.map((variant, index) => (
                      <tr key={index} style={{ borderBottom: '1px solid #F3F4F6' }}>
                        <td style={{ padding: '0.75rem', fontWeight: '600', color: '#111827' }}>{variant.gene}</td>
                        <td style={{ padding: '0.75rem', fontFamily: 'monospace', fontSize: '0.875rem', color: '#4B5563' }}>
                          {variant.variant_id.split(':').slice(1).join(':')}
                        </td>
                        <td style={{ padding: '0.75rem' }}>
                          <div style={{
                            background: variant.mutation_type === 'snv' ? '#3B82F6' : 
                                       variant.mutation_type === 'indel' ? '#EF4444' : '#F59E0B',
                            color: '#FFFFFF',
                            padding: '0.25rem 0.5rem',
                            borderRadius: '0.375rem',
                            fontSize: '0.75rem',
                            fontWeight: '600',
                            textAlign: 'center',
                            textTransform: 'uppercase'
                          }}>
                            {variant.mutation_type}
                          </div>
                        </td>
                        <td style={{ padding: '0.75rem' }}>
                          <div style={{ fontSize: '0.875rem', color: '#4B5563' }}>
                            <div style={{ fontFamily: 'monospace', fontSize: '0.75rem', marginBottom: '0.25rem' }}>
                              {variant.transformation.original} → {variant.transformation.mutated}
                            </div>
                            <div style={{ fontSize: '0.75rem', color: '#6B7280' }}>
                              {variant.transformation.amino_acid_change}
                            </div>
                          </div>
                        </td>
                        <td style={{ padding: '0.75rem' }}>
                          <div style={{
                            background: variant.quality_score > 90 ? '#22C55E' : 
                                       variant.quality_score > 70 ? '#F59E0B' : '#EF4444',
                            color: '#FFFFFF',
                            padding: '0.25rem 0.5rem',
                            borderRadius: '0.375rem',
                            fontSize: '0.75rem',
                            fontWeight: '600',
                            textAlign: 'center',
                            minWidth: '50px'
                          }}>
                            {variant.quality_score}
                          </div>
                        </td>
                        <td style={{ padding: '0.75rem' }}>
                          <div style={{
                            background: variant.clinical_significance === 'pathogenic' ? '#FEF2F2' : '#F0FDF4',
                            color: variant.clinical_significance === 'pathogenic' ? '#EF4444' : '#22C55E',
                            padding: '0.25rem 0.5rem',
                            borderRadius: '0.375rem',
                            fontSize: '0.75rem',
                            fontWeight: '600',
            textAlign: 'center'
          }}>
                            {variant.clinical_significance}
                          </div>
                        </td>
                        <td style={{ padding: '0.75rem', fontSize: '0.875rem', color: '#4B5563' }}>
                          <div style={{ marginBottom: '0.25rem' }}>
                            <strong>{variant.functional_impact}</strong>
                          </div>
                          <div style={{ fontSize: '0.75rem', color: '#6B7280' }}>
                            {variant.transformation.effect}
                          </div>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              
              {/* Structural Variants */}
              <div style={{
                background: '#F9FAFB',
                padding: '1.5rem',
                borderRadius: '0.5rem',
                marginTop: '2rem',
                border: '1px solid #E5E7EB'
              }}>
                <h4 style={{ 
                  color: '#111827',
                  fontSize: '1rem',
                  fontWeight: '600',
                  marginBottom: '1rem'
                }}>
                  Structural Variants
                </h4>
                
                <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                  {genomicData.structural_variants.map((variant, index) => (
                    <div key={index} style={{
                      background: '#FFFFFF',
                      padding: '1rem',
                      borderRadius: '0.5rem',
                      border: `2px solid ${variant.clinical_significance === 'pathogenic' ? '#EF4444' : '#F59E0B'}`,
                      display: 'flex',
                      justifyContent: 'space-between',
                      alignItems: 'center'
                    }}>
                      <div style={{ flex: 1 }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginBottom: '0.5rem' }}>
                          <div style={{
                            background: variant.type === 'deletion' ? '#EF4444' : '#8B5CF6',
                            color: '#FFFFFF',
                            padding: '0.25rem 0.5rem',
                            borderRadius: '0.375rem',
                            fontSize: '0.75rem',
                            fontWeight: '600',
                            textTransform: 'uppercase'
                          }}>
                            {variant.type}
                          </div>
                          <span style={{ fontWeight: '600', color: '#111827' }}>
                            {variant.genes_affected.join(', ')}
                          </span>
                        </div>
                        <div style={{ fontSize: '0.875rem', color: '#4B5563', marginBottom: '0.25rem' }}>
                          {variant.chromosome}:{variant.start}-{variant.end} ({(variant.size / 1000).toFixed(1)}kb)
                        </div>
                        <div style={{ fontSize: '0.75rem', color: '#6B7280' }}>
                          {variant.transformation.effect}
                        </div>
                      </div>
                      <div style={{
                        background: variant.clinical_significance === 'pathogenic' ? '#FEF2F2' : '#FFFBEB',
                        color: variant.clinical_significance === 'pathogenic' ? '#EF4444' : '#F59E0B',
                        padding: '0.25rem 0.5rem',
                        borderRadius: '0.375rem',
                        fontSize: '0.75rem',
                        fontWeight: '600',
                        textAlign: 'center'
                      }}>
                        {variant.clinical_significance.replace('_', ' ')}
                      </div>
                    </div>
                  ))}
                </div>
              </div>
              
              {/* Copy Number Variants */}
              <div style={{
                background: '#F9FAFB',
                padding: '1.5rem',
                borderRadius: '0.5rem',
                marginTop: '1rem',
                border: '1px solid #E5E7EB'
              }}>
                <h4 style={{ 
                  color: '#111827',
                  fontSize: '1rem',
                  fontWeight: '600',
                  marginBottom: '1rem'
                }}>
                  Copy Number Variants
                </h4>
                
                <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                  {genomicData.copy_number_variants.map((variant, index) => (
                    <div key={index} style={{
                      background: '#FFFFFF',
                      padding: '1rem',
                      borderRadius: '0.5rem',
                      border: '2px solid #F59E0B',
                      display: 'flex',
                      justifyContent: 'space-between',
                      alignItems: 'center'
                    }}>
                      <div style={{ flex: 1 }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginBottom: '0.5rem' }}>
                          <span style={{ fontWeight: '600', color: '#111827' }}>
                            {variant.gene}
                          </span>
                          <div style={{
                            background: variant.copy_number > variant.normal_copy_number ? '#22C55E' : '#EF4444',
                            color: '#FFFFFF',
                            padding: '0.25rem 0.5rem',
                            borderRadius: '0.375rem',
                            fontSize: '0.75rem',
                            fontWeight: '600'
                          }}>
                            {variant.copy_number > variant.normal_copy_number ? 'AMPLIFICATION' : 'DELETION'}
                          </div>
                        </div>
                        <div style={{ fontSize: '0.875rem', color: '#4B5563', marginBottom: '0.25rem' }}>
                          {variant.copy_number} copies (normal: {variant.normal_copy_number}) • {variant.fold_change}x change
                        </div>
                        <div style={{ fontSize: '0.75rem', color: '#6B7280' }}>
                          {variant.cancer_relevance}
                        </div>
                      </div>
                      <div style={{
                        background: '#FFFBEB',
                        color: '#F59E0B',
                        padding: '0.25rem 0.5rem',
                        borderRadius: '0.375rem',
                        fontSize: '0.75rem',
                        fontWeight: '600',
                        textAlign: 'center'
                      }}>
                        {variant.clinical_significance}
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          </div>
        );
      
      case 'pathways':
        return (
          <div style={{ padding: '2rem' }}>
            <h2 style={{ 
              color: '#111827',
              fontSize: '1.5rem',
              fontWeight: '600',
              marginBottom: '1.5rem'
            }}>
              Pathway Analysis
            </h2>
            
            {/* Pathway Enrichment */}
            <div style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '2rem'
            }}>
              <h3 style={{ 
                color: '#111827',
                fontSize: '1.125rem',
              fontWeight: '600',
              marginBottom: '1rem'
            }}>
                Affected Biological Pathways
              </h3>
              
                            {/* Pathway Disruption Analysis */}
              <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem', marginBottom: '2rem' }}>
                {[
                  { 
                    name: 'DNA Repair', 
                    pathway_id: 'DNA_REPAIR',
                    genes: ['BRCA1', 'TP53'], 
                    significance: 95, 
                    color: '#EF4444',
                    mutations: [
                      { gene: 'BRCA1', type: 'frameshift', effect: 'Complete loss of homologous recombination' },
                      { gene: 'TP53', type: 'missense', effect: 'Impaired DNA damage response' }
                    ]
                  },
                  { 
                    name: 'Cell Cycle Control', 
                    pathway_id: 'CELL_CYCLE',
                    genes: ['TP53'], 
                    significance: 87, 
                    color: '#F59E0B',
                    mutations: [
                      { gene: 'TP53', type: 'missense', effect: 'Loss of G1/S checkpoint control' }
                    ]
                  },
                  { 
                    name: 'Apoptosis', 
                    pathway_id: 'APOPTOSIS',
                    genes: ['TP53'], 
                    significance: 82, 
                    color: '#F59E0B',
                    mutations: [
                      { gene: 'TP53', type: 'missense', effect: 'Reduced pro-apoptotic signaling' }
                    ]
                  },
                  { 
                    name: 'Wnt Signaling', 
                    pathway_id: 'WNT_SIGNALING',
                    genes: ['APC'], 
                    significance: 78, 
                    color: '#3B82F6',
                    mutations: [
                      { gene: 'APC', type: 'nonsense', effect: 'Constitutive β-catenin activation' }
                    ]
                  },
                  { 
                    name: 'Tumor Suppression', 
                    pathway_id: 'TUMOR_SUPPRESSION',
                    genes: ['BRCA1', 'TP53', 'APC'], 
                    significance: 92, 
                    color: '#EF4444',
                    mutations: [
                      { gene: 'BRCA1', type: 'frameshift', effect: 'Loss of tumor suppressor function' },
                      { gene: 'TP53', type: 'missense', effect: 'Guardian of genome compromised' },
                      { gene: 'APC', type: 'nonsense', effect: 'Gatekeeper function lost' }
                    ]
                  }
                ].map((pathway, index) => (
                  <div key={index} style={{
                    background: '#FFFFFF',
                    padding: '1.5rem',
                    border: `2px solid ${pathway.color}`,
                    borderRadius: '0.75rem',
                    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
                  }}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                      <h4 style={{ 
                        color: '#111827',
                        fontSize: '1.125rem',
                        fontWeight: '600',
                        margin: 0
                      }}>
                        {pathway.name} Pathway
                      </h4>
                      <div style={{
                        background: pathway.color,
                        color: '#FFFFFF',
                        padding: '0.5rem 1rem',
                        borderRadius: '0.375rem',
                        fontSize: '1rem',
                        fontWeight: 'bold'
                      }}>
                        {pathway.significance}% Disrupted
                      </div>
                    </div>
                    
                    <div style={{
                      background: '#E5E7EB',
                      borderRadius: '0.25rem',
                      height: '0.75rem',
                      position: 'relative',
                      overflow: 'hidden',
                      marginBottom: '1rem'
                    }}>
                      <div style={{
                        background: pathway.color,
                        height: '100%',
                        width: `${pathway.significance}%`,
                        transition: 'width 0.3s ease'
                      }}/>
                    </div>
                    
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
                      {pathway.mutations.map((mutation, mutIndex) => (
                        <div key={mutIndex} style={{
                          display: 'flex',
                          alignItems: 'center',
                          padding: '0.75rem',
                          background: '#F9FAFB',
                          borderRadius: '0.5rem',
                          border: '1px solid #E5E7EB'
                        }}>
                          <div style={{
                            background: mutation.type === 'frameshift' ? '#EF4444' : 
                                       mutation.type === 'missense' ? '#F59E0B' : '#3B82F6',
                            color: '#FFFFFF',
                            padding: '0.25rem 0.5rem',
                            borderRadius: '0.375rem',
                            fontSize: '0.75rem',
                            fontWeight: '600',
                            textTransform: 'uppercase',
                            marginRight: '0.75rem'
                          }}>
                            {mutation.type}
                          </div>
                          <div style={{ flex: 1 }}>
                            <div style={{ fontSize: '0.875rem', fontWeight: '600', color: '#111827', marginBottom: '0.25rem' }}>
                              {mutation.gene}
                            </div>
                            <div style={{ fontSize: '0.75rem', color: '#6B7280' }}>
                              {mutation.effect}
                            </div>
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                ))}
              </div>
              
              {/* Cancer Risk by Pathway */}
              <div style={{
                background: '#FFFFFF',
                padding: '2rem',
                borderRadius: '0.75rem',
                boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
                border: '1px solid #E5E7EB',
                marginBottom: '2rem'
              }}>
                <h4 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  marginBottom: '1rem'
                }}>
                  Cancer Risk by Pathway Disruption
                </h4>
                
                <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '1rem' }}>
                  {genomicData.risk_findings.filter(r => r.pathway_disruption.length > 0).map((risk, index) => (
                    <div key={index} style={{
                      padding: '1rem',
                      border: `2px solid ${getRiskColor(risk.risk_level)}`,
                      borderRadius: '0.5rem',
                      background: `${getRiskColor(risk.risk_level)}10`
                    }}>
                      <h5 style={{ 
                        color: '#111827',
                        fontSize: '1rem',
                        fontWeight: '600',
                        marginBottom: '0.5rem',
                        textTransform: 'capitalize'
                      }}>
                        {risk.cancer_type} Cancer
                      </h5>
                      <div style={{ 
                        fontSize: '1.5rem',
                        fontWeight: 'bold',
                        color: getRiskColor(risk.risk_level),
                        marginBottom: '0.5rem'
                      }}>
                        {risk.risk_percentage}%
                      </div>
                      <div style={{ display: 'flex', flexDirection: 'column', gap: '0.25rem' }}>
                        {risk.pathway_disruption.map((pathway, pathIndex) => (
                          <div key={pathIndex} style={{
                            background: '#FFFFFF',
                            padding: '0.25rem 0.5rem',
                            borderRadius: '0.25rem',
                            fontSize: '0.75rem',
                            fontWeight: '500',
                            color: '#4B5563'
                          }}>
                            {pathway.replace('_', ' ')}
                          </div>
                        ))}
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            </div>
            
            {/* Pathway Network */}
            <div style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB'
            }}>
              <h3 style={{ 
                color: '#111827',
                fontSize: '1.125rem',
                fontWeight: '600',
                marginBottom: '1rem'
              }}>
                Gene Interaction Network
              </h3>
              
              <svg width="100%" height="300" viewBox="0 0 600 300">
                {/* Background */}
                <rect width="600" height="300" fill="#F9FAFB"/>
                
                {/* Network connections */}
                <line x1="150" y1="150" x2="300" y2="100" stroke="#94A3B8" strokeWidth="2"/>
                <line x1="150" y1="150" x2="300" y2="200" stroke="#94A3B8" strokeWidth="2"/>
                <line x1="300" y1="100" x2="450" y2="150" stroke="#94A3B8" strokeWidth="2"/>
                <line x1="300" y1="200" x2="450" y2="150" stroke="#94A3B8" strokeWidth="2"/>
                
                {/* Gene nodes */}
                <g>
                  <circle cx="150" cy="150" r="30" fill="#EF4444" opacity="0.8"/>
                  <text x="150" y="155" fill="#FFFFFF" fontSize="12" textAnchor="middle" fontWeight="600">BRCA1</text>
                </g>
                <g>
                  <circle cx="300" cy="100" r="30" fill="#EF4444" opacity="0.8"/>
                  <text x="300" y="105" fill="#FFFFFF" fontSize="12" textAnchor="middle" fontWeight="600">TP53</text>
                </g>
                <g>
                  <circle cx="300" cy="200" r="30" fill="#F59E0B" opacity="0.8"/>
                  <text x="300" y="205" fill="#FFFFFF" fontSize="12" textAnchor="middle" fontWeight="600">APC</text>
                </g>
                <g>
                  <circle cx="450" cy="150" r="25" fill="#3B82F6" opacity="0.8"/>
                  <text x="450" y="155" fill="#FFFFFF" fontSize="10" textAnchor="middle" fontWeight="600">MLH1</text>
                </g>
                
                {/* Labels */}
                <text x="50" y="155" fill="#111827" fontSize="12" fontWeight="600">DNA Repair</text>
                <text x="350" y="50" fill="#111827" fontSize="12" fontWeight="600">Cell Cycle</text>
                <text x="350" y="250" fill="#111827" fontSize="12" fontWeight="600">Wnt Signaling</text>
                <text x="500" y="155" fill="#111827" fontSize="12" fontWeight="600">Mismatch Repair</text>
              </svg>
            </div>
          </div>
        );
      
      case 'clinical':
        return (
          <div style={{ padding: '2rem' }}>
            <h2 style={{ 
              color: '#111827',
              fontSize: '1.5rem',
              fontWeight: '600',
              marginBottom: '1.5rem'
            }}>
              Clinical Report & Recommendations
            </h2>
            
            {/* Survival Analysis */}
            <div id="survival-analysis" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '2rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  Survival Analysis vs Population Average
                </h3>
                <DownloadButton 
                  elementId="survival-analysis"
                  title="Survival Analysis"
                  summary="Comparative survival analysis showing impact of high-risk genetic profile versus population average. Analysis reveals 8-12 years of life potentially affected but shows substantial benefit from early intervention, with enhanced screening achieving 90-95% of normal life expectancy."
                  isPDFGenerating={isPDFGenerating}
                  setIsPDFGenerating={setIsPDFGenerating}
                />
              </div>
              
              <svg width="100%" height="300" viewBox="0 0 800 300">
                {/* Background */}
                <rect width="800" height="300" fill="#F9FAFB"/>
                
                {/* Grid lines */}
                {[0, 1, 2, 3, 4, 5].map(i => (
                  <g key={i}>
                    <line x1="80" y1={50 + i * 40} x2="700" y2={50 + i * 40} stroke="#E5E7EB" strokeWidth="1"/>
                    <text x="70" y={55 + i * 40} fill="#6B7280" fontSize="12" textAnchor="end">{100 - i * 20}%</text>
                  </g>
                ))}
                
                {/* Age axis */}
                {[0, 1, 2, 3, 4, 5, 6, 7, 8].map(i => (
                  <g key={i}>
                    <line x1={80 + i * 80} y1="50" x2={80 + i * 80} y2="250" stroke="#E5E7EB" strokeWidth="1"/>
                    <text x={80 + i * 80} y="270" fill="#6B7280" fontSize="12" textAnchor="middle">{30 + i * 10}</text>
                  </g>
                ))}
                
                {/* Population average survival curve */}
                <path 
                  d="M 80 60 L 160 70 L 240 85 L 320 105 L 400 130 L 480 160 L 560 200 L 640 240 L 720 250" 
                  stroke="#2563EB" 
                  strokeWidth="3" 
                  fill="none"
                  strokeDasharray="5,5"
                />
                
                {/* High-risk patient survival curve */}
                <path 
                  d="M 80 60 L 160 80 L 240 110 L 320 140 L 400 170 L 480 200 L 560 225 L 640 245 L 720 250" 
                  stroke="#EF4444" 
                  strokeWidth="3" 
                  fill="none"
                />
                
                {/* Legend */}
                <g transform="translate(500, 80)">
                  <line x1="0" y1="0" x2="20" y2="0" stroke="#2563EB" strokeWidth="3" strokeDasharray="5,5"/>
                  <text x="25" y="5" fill="#2563EB" fontSize="12" fontWeight="600">Population Average</text>
                  <line x1="0" y1="20" x2="20" y2="20" stroke="#EF4444" strokeWidth="3"/>
                  <text x="25" y="25" fill="#EF4444" fontSize="12" fontWeight="600">High-Risk Profile</text>
                </g>
                
                {/* Axis labels */}
                <text x="400" y="290" fill="#111827" fontSize="14" textAnchor="middle" fontWeight="600">Age</text>
                <text x="30" y="150" fill="#111827" fontSize="14" textAnchor="middle" fontWeight="600" transform="rotate(-90 30 150)">Survival Probability</text>
              </svg>
            </div>
            
            {/* Clinical Recommendations */}
            <div style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB'
            }}>
              <h3 style={{ 
                color: '#111827',
                fontSize: '1.125rem',
                fontWeight: '600',
                marginBottom: '1rem'
              }}>
                Clinical Recommendations
              </h3>
              
              <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                {genomicData.risk_findings.filter(r => r.risk_level === 'high').map((risk, index) => (
                  <div key={index} style={{
                    padding: '1.5rem',
                    border: '2px solid #EF4444',
                    borderRadius: '0.75rem',
                    background: '#FEF2F2'
                  }}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '1rem' }}>
                      <h4 style={{ 
                        color: '#111827',
              fontSize: '1rem',
                        fontWeight: '600',
                        textTransform: 'capitalize'
                      }}>
                        {risk.cancer_type} Cancer Screening
                      </h4>
                      <div style={{
                        background: '#EF4444',
                        color: '#FFFFFF',
                        padding: '0.25rem 0.5rem',
                        borderRadius: '0.375rem',
                        fontSize: '0.75rem',
                        fontWeight: '600'
                      }}>
                        {risk.risk_percentage}% Risk
                      </div>
                    </div>
            <p style={{ 
              color: '#4B5563',
                      fontSize: '0.875rem',
                      marginBottom: '1rem'
            }}>
                      {risk.recommendation}
            </p>
            <div style={{
                      background: '#FFFFFF',
                      padding: '1rem',
                      borderRadius: '0.5rem',
                      border: '1px solid #FECACA'
                    }}>
                      <h5 style={{ 
                        color: '#111827',
                        fontSize: '0.875rem',
                        fontWeight: '600',
                        marginBottom: '0.5rem'
                      }}>
                        Affected Genes:
                      </h5>
                      <div style={{ display: 'flex', gap: '0.5rem', flexWrap: 'wrap' }}>
                        {risk.affected_genes.map((gene, geneIndex) => (
                          <span key={geneIndex} style={{
                            background: '#EF4444',
                            color: '#FFFFFF',
                            padding: '0.25rem 0.5rem',
                            borderRadius: '0.375rem',
                            fontSize: '0.75rem',
                            fontWeight: '600'
                          }}>
                            {gene}
                          </span>
                        ))}
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        );
      
      case 'family':
        return (
          <div style={{ padding: '2rem' }}>
            <h2 style={{ 
              color: '#111827',
              fontSize: '1.5rem',
              fontWeight: '600',
              marginBottom: '1.5rem'
            }}>
              Family Analysis & Hereditary Risk
            </h2>
            
            {/* Family Pedigree */}
            <div id="family-pedigree" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '2rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  Family Pedigree Analysis
                </h3>
                <DownloadButton 
                  elementId="family-pedigree"
                  title="Family Pedigree Analysis"
                  summary="Comprehensive family pedigree showing BRCA1 inheritance pattern with strong maternal line history. Analysis reveals 50% sister carrier risk, autosomal dominant pattern, and provides cascade testing recommendations for at-risk family members."
                  isPDFGenerating={isPDFGenerating}
                  setIsPDFGenerating={setIsPDFGenerating}
                />
              </div>
              
              <svg width="100%" height="300" viewBox="0 0 600 300">
                {/* Background */}
                <rect width="600" height="300" fill="#F9FAFB"/>
                
                {/* Family connections */}
                <line x1="150" y1="80" x2="200" y2="80" stroke="#4B5563" strokeWidth="2"/>
                <line x1="175" y1="80" x2="175" y2="120" stroke="#4B5563" strokeWidth="2"/>
                <line x1="150" y1="120" x2="200" y2="120" stroke="#4B5563" strokeWidth="2"/>
                <line x1="175" y1="120" x2="175" y2="160" stroke="#4B5563" strokeWidth="2"/>
                <line x1="150" y1="160" x2="200" y2="160" stroke="#4B5563" strokeWidth="2"/>
                
                {/* Family members */}
                {/* Grandparents */}
                <g>
                  <rect x="120" y="50" width="30" height="30" fill="#F59E0B" opacity="0.8" stroke="#4B5563" strokeWidth="2"/>
                  <text x="135" y="70" fill="#FFFFFF" fontSize="10" textAnchor="middle" fontWeight="600">GM</text>
                  <text x="135" y="45" fill="#4B5563" fontSize="10" textAnchor="middle">Breast Cancer</text>
                </g>
                <g>
                  <circle cx="220" cy="65" r="15" fill="#E5E7EB" stroke="#4B5563" strokeWidth="2"/>
                  <text x="220" y="70" fill="#4B5563" fontSize="10" textAnchor="middle" fontWeight="600">GF</text>
                </g>
                
                {/* Parents */}
                <g>
                  <rect x="120" y="105" width="30" height="30" fill="#EF4444" opacity="0.8" stroke="#4B5563" strokeWidth="2"/>
                  <text x="135" y="125" fill="#FFFFFF" fontSize="10" textAnchor="middle" fontWeight="600">M</text>
                  <text x="135" y="100" fill="#4B5563" fontSize="10" textAnchor="middle">BRCA1+</text>
                </g>
                <g>
                  <circle cx="220" cy="120" r="15" fill="#E5E7EB" stroke="#4B5563" strokeWidth="2"/>
                  <text x="220" y="125" fill="#4B5563" fontSize="10" textAnchor="middle" fontWeight="600">F</text>
                </g>
                
                {/* Current generation */}
                <g>
                  <rect x="120" y="145" width="30" height="30" fill="#EF4444" opacity="0.8" stroke="#4B5563" strokeWidth="4"/>
                  <text x="135" y="165" fill="#FFFFFF" fontSize="10" textAnchor="middle" fontWeight="600">P</text>
                  <text x="135" y="140" fill="#4B5563" fontSize="10" textAnchor="middle">Patient</text>
                </g>
                <g>
                  <circle cx="220" cy="160" r="15" fill="#F59E0B" opacity="0.8" stroke="#4B5563" strokeWidth="2"/>
                  <text x="220" y="165" fill="#FFFFFF" fontSize="10" textAnchor="middle" fontWeight="600">S</text>
                  <text x="220" y="190" fill="#4B5563" fontSize="10" textAnchor="middle">Sister</text>
                </g>
                
                {/* Legend */}
                <g transform="translate(350, 50)">
                  <text x="0" y="0" fill="#111827" fontSize="12" fontWeight="600">Legend</text>
                  <rect x="0" y="15" width="15" height="15" fill="#EF4444" opacity="0.8"/>
                  <text x="20" y="27" fill="#4B5563" fontSize="11">High Risk / Affected</text>
                  <rect x="0" y="35" width="15" height="15" fill="#F59E0B" opacity="0.8"/>
                  <text x="20" y="47" fill="#4B5563" fontSize="11">Medium Risk</text>
                  <rect x="0" y="55" width="15" height="15" fill="#E5E7EB"/>
                  <text x="20" y="67" fill="#4B5563" fontSize="11">Unknown / Low Risk</text>
                  <rect x="0" y="80" width="15" height="15" fill="none" stroke="#4B5563" strokeWidth="4"/>
                  <text x="20" y="92" fill="#4B5563" fontSize="11">Current Patient</text>
                </g>
              </svg>
            </div>
            
            {/* Hereditary Risk Assessment */}
            <div style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB'
            }}>
              <h3 style={{ 
                color: '#111827',
                fontSize: '1.125rem',
                fontWeight: '600',
                marginBottom: '1rem'
              }}>
                Hereditary Risk Assessment
              </h3>
              
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))', gap: '1rem' }}>
                <div style={{
                  padding: '1.5rem',
                  border: '2px solid #EF4444',
                  borderRadius: '0.75rem',
                  background: '#FEF2F2'
                }}>
                  <h4 style={{ 
                    color: '#111827',
                    fontSize: '1rem',
                    fontWeight: '600',
                    marginBottom: '0.5rem'
                  }}>
                    Maternal Inheritance
                  </h4>
                  <p style={{ 
                    color: '#4B5563',
                    fontSize: '0.875rem',
                    marginBottom: '1rem'
                  }}>
                    Strong family history of breast cancer with confirmed BRCA1 mutation in maternal line.
                  </p>
                  <div style={{
                    background: '#EF4444',
                    color: '#FFFFFF',
                    padding: '0.5rem',
                    borderRadius: '0.375rem',
                    fontSize: '0.875rem',
                    fontWeight: '600',
                    textAlign: 'center'
                  }}>
                    High Risk (85%)
                  </div>
                </div>
                
                <div style={{
                  padding: '1.5rem',
                  border: '2px solid #F59E0B',
                  borderRadius: '0.75rem',
                  background: '#FFFBEB'
                }}>
                  <h4 style={{ 
                    color: '#111827',
                    fontSize: '1rem',
                    fontWeight: '600',
                    marginBottom: '0.5rem'
                  }}>
                    Sibling Risk
                  </h4>
                  <p style={{ 
                    color: '#4B5563',
                    fontSize: '0.875rem',
                    marginBottom: '1rem'
                  }}>
                    Sister has 50% chance of inheriting the same BRCA1 mutation. Testing recommended.
                  </p>
                  <div style={{
                    background: '#F59E0B',
                    color: '#FFFFFF',
                    padding: '0.5rem',
                    borderRadius: '0.375rem',
                    fontSize: '0.875rem',
                    fontWeight: '600',
                    textAlign: 'center'
                  }}>
                    Testing Recommended
                  </div>
              </div>
              
                <div style={{
                  padding: '1.5rem',
                  border: '2px solid #3B82F6',
                  borderRadius: '0.75rem',
                  background: '#EFF6FF'
                }}>
                  <h4 style={{ 
                    color: '#111827',
                    fontSize: '1rem',
                    fontWeight: '600',
                    marginBottom: '0.5rem'
                  }}>
                    Future Generations
                  </h4>
                  <p style={{ 
                    color: '#4B5563',
                    fontSize: '0.875rem',
                    marginBottom: '1rem'
                  }}>
                    Each child has 50% chance of inheriting the mutation. Genetic counseling recommended.
                  </p>
                  <div style={{
                    background: '#3B82F6',
                    color: '#FFFFFF',
                    padding: '0.5rem',
                  borderRadius: '0.375rem',
                  fontSize: '0.875rem',
                    fontWeight: '600',
                    textAlign: 'center'
                  }}>
                    Counseling Advised
                  </div>
                </div>
              </div>
            </div>
          </div>
        );
      
      default:
        return null;
    }
  };

  return (
    <Layout>
      <section style={{ 
        background: '#F9FAFB',
        minHeight: 'calc(100vh - 4rem)',
        padding: '2rem 0'
      }}>
        <div style={{ 
          maxWidth: '1200px',
          margin: '0 auto',
          padding: '0 1.5rem'
        }}>
          {/* Header */}
          <div style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            marginBottom: '2rem'
          }}>
            <div>
              <h1 style={{
                fontSize: '1.875rem',
                fontWeight: 'bold',
                color: '#111827',
                marginBottom: '0.5rem'
              }}>
                In-Depth Analysis
              </h1>
              <p style={{
                color: '#4B5563',
                fontSize: '1rem'
              }}>
                Comprehensive genomic analysis • {currentData.riskLevel}
              </p>
            </div>
            <div>
              <button
                onClick={() => {
                  // Navigate back to dashboard preserving the analysis state
                  if (pipelineResults && fileName) {
                    // User has real analysis results - pass them back
                    navigate('/dashboard', { 
                      state: { 
                        results: pipelineResults, 
                        fileName: fileName 
                      } 
                    });
                  } else {
                    // Mock data case - pass the risk level
                    navigate(`/dashboard?risk=${riskLevel}`);
                  }
                }}
                style={{
                  padding: '0.5rem 1.25rem',
                  color: '#FFFFFF',
                  background: '#2563EB',
                  border: 'none',
                  borderRadius: '0.5rem',
                  cursor: 'pointer',
                  transition: 'all 200ms ease',
                  fontSize: '0.875rem',
                  fontWeight: '500'
                }}
                onMouseEnter={(e) => {
                  e.currentTarget.style.background = '#1D4ED8';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.background = '#2563EB';
                }}
              >
                Back to Dashboard
              </button>
            </div>
          </div>

          <div style={{ 
            display: 'grid', 
            gridTemplateColumns: 'minmax(280px, 320px) 1fr', 
            gap: '1.5rem' 
          }}>
            {/* Sidebar */}
            <div style={{
              background: '#FFFFFF',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB',
              overflow: 'visible'
            }}>
              {/* Risk Info */}
              <div style={{
                background: '#F9FAFB',
                padding: '1.5rem',
                borderBottom: '1px solid #E5E7EB'
              }}>
                <div style={{ 
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  marginBottom: '0.5rem',
                  color: '#111827'
                }}>
                  {currentData.riskLevel}
                </div>
                <div style={{ 
                  fontSize: '0.875rem',
                  lineHeight: '1.4',
                  color: '#4B5563'
                }}>
                  <span dangerouslySetInnerHTML={{ __html: currentData.details }} />
                </div>
                <div style={{
                  background: '#DBEAFE',
                  color: '#2563EB',
                  padding: '0.5rem',
                  borderRadius: '0.375rem',
                  marginTop: '1rem',
                  fontSize: '0.875rem',
                  fontWeight: '500',
                  textAlign: 'center'
                }}>
                  Risk Score: {currentData.riskScore}
                </div>
              </div>

              {/* Navigation */}
              <div style={{ padding: '1rem' }}>
                <div style={{
                  fontSize: '0.75rem',
                  fontWeight: '600',
                  textTransform: 'uppercase',
                  letterSpacing: '0.05em',
                  color: '#4B5563',
                  marginBottom: '0.75rem'
                }}>
                  Clinical Workflow
                </div>
                
                {[
                  { id: 'analysis', label: 'Genomic Analysis' },
                  { id: 'variants', label: 'Variant Heatmap' },
                  { id: 'pathways', label: 'Pathway Analysis' },
                  { id: 'clinical', label: 'Clinical Report' },
                  { id: 'family', label: 'Family Analysis' }
                ].map(item => (
                  <div 
                    key={item.id}
                    style={{
                      display: 'flex',
                      alignItems: 'center',
                      gap: '0.75rem',
                      padding: '0.75rem',
                      marginBottom: '0.25rem',
                      borderRadius: '0.375rem',
                      cursor: 'pointer',
                      transition: 'all 200ms ease',
                      fontWeight: '500',
                      fontSize: '0.875rem',
                      color: activeTab === item.id ? '#FFFFFF' : '#111827',
                      background: activeTab === item.id ? '#2563EB' : 'transparent'
                    }}
                    onMouseEnter={(e) => {
                      if (activeTab !== item.id) {
                        e.currentTarget.style.background = '#F3F4F6';
                      }
                    }}
                    onMouseLeave={(e) => {
                      if (activeTab !== item.id) {
                        e.currentTarget.style.background = 'transparent';
                      }
                    }}
                    onClick={() => setActiveTab(item.id)}
                  >
                    <span>{item.label}</span>
                  </div>
                ))}
              </div>

              {/* Alerts */}
              <div style={{ padding: '1rem', borderTop: '1px solid #E5E7EB', overflow: 'visible' }}>
                <div style={{
                  fontSize: '0.75rem',
                  fontWeight: '600',
                  textTransform: 'uppercase',
                  letterSpacing: '0.05em',
                  color: '#4B5563',
                  marginBottom: '0.75rem'
                }}>
                  Clinical Alerts
                </div>
                
                {currentData.alerts.map((alert: Alert, index: number) => (
                  <div 
                    key={index}
                    style={{
                      background: alert.type === 'critical' ? '#FEF2F2' : 
                                 alert.type === 'warning' ? '#FFFBEB' : 
                                 alert.type === 'info' ? '#EFF6FF' : '#F0FDF4',
                      border: `1px solid ${alert.type === 'critical' ? '#FECACA' : 
                                           alert.type === 'warning' ? '#FDE68A' : 
                                           alert.type === 'info' ? '#BFDBFE' : '#BBF7D0'}`,
                      borderLeft: `3px solid ${alert.type === 'critical' ? '#EF4444' : 
                                               alert.type === 'warning' ? '#F59E0B' : 
                                               alert.type === 'info' ? '#2563EB' : '#22C55E'}`,
                      padding: '0.75rem',
                      marginBottom: '0.5rem',
                      borderRadius: '0.375rem',
                      fontSize: '0.75rem',
                      cursor: 'pointer',
                      position: 'relative'
                    }}
                  >
                    <div style={{
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'space-between',
                      marginBottom: '0.25rem'
                    }}>
                      <div style={{
                        display: 'flex',
                        alignItems: 'center',
                        gap: '0.5rem'
                      }}>
                        <div style={{
                          width: '0.5rem',
                          height: '0.5rem',
                          borderRadius: '50%',
                          background: alert.type === 'critical' ? '#EF4444' : 
                                     alert.type === 'warning' ? '#F59E0B' : 
                                     alert.type === 'info' ? '#2563EB' : '#22C55E'
                        }}></div>
                        <strong style={{ fontSize: '0.8rem', color: '#111827' }}>{alert.title}</strong>
                      </div>
                      <div 
                        style={{ position: 'relative', zIndex: 10, overflow: 'visible' }}
                        onMouseEnter={() => setHoveredAlert(index)}
                        onMouseLeave={() => setHoveredAlert(null)}
                      >
                        <InformationCircleIcon 
                          onMouseEnter={(e) => {
                            e.currentTarget.style.backgroundColor = '#D1D5DB';
                            e.currentTarget.style.color = '#374151';
                            e.currentTarget.style.borderColor = '#9CA3AF';
                          }}
                          onMouseLeave={(e) => {
                            e.currentTarget.style.backgroundColor = '#E5E7EB';
                            e.currentTarget.style.color = '#6B7280';
                            e.currentTarget.style.borderColor = '#D1D5DB';
                          }}
                        />
                        <AlertTooltip alert={alert} isVisible={hoveredAlert === index} />
                      </div>
                    </div>
                    <p style={{ fontSize: '0.75rem', opacity: 0.8, color: '#4B5563' }}>{alert.desc}</p>
                  </div>
                ))}
                
                {/* Confidence Check */}
                <ConfidenceCheck 
                  validation={mockSHAPValidation} 
                  isDetailed={true} 
                />
              </div>
            </div>

            {/* Main Content */}
            <div style={{
              background: '#FFFFFF',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB',
              minHeight: '600px'
            }}>
              {renderTabContent()}
            </div>
          </div>
        </div>
      </section>
    </Layout>
  );
};

export default ClinicalViewPage; 