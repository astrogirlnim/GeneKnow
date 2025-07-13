import React, { useState, useRef, useEffect, useMemo, useCallback } from 'react';
import { useNavigate, useSearchParams, useLocation } from 'react-router-dom';
import Layout from '../components/Layout';
import ConfidenceCheck from '../components/ConfidenceCheck';
import DetailedDisclaimer from '../components/DetailedDisclaimer';
import type { PipelineResult } from '../api/geneknowPipeline';
import html2canvas from 'html2canvas';
import jsPDF from 'jspdf';
import { invoke } from '@tauri-apps/api/core';



// Type definitions for genomic data structures
// (Structural and Copy Number variant interfaces removed as they're not currently used)

// Type definitions for subtab content
interface SubtabContent {
  id: string;
  title: string;
  elementIds: string[];
}

interface VariantTableRow {
  gene: string;
  variant: string;
  type: string;
  quality: string;
  significance: string;
  mainImpact: string;
  subImpact: string;
}

// Mock SHAP validation function removed - only real pipeline data is used

// Mock genomic data removed - only real pipeline data is used

// Helper function to get color based on risk level
const getRiskColor = (riskLevel: string) => {
  switch (riskLevel) {
    case 'high': return '#EF4444';
    case 'medium': return '#F59E0B';
    case 'low': return '#22C55E';
    default: return '#6B7280';
  }
};

// Enhanced PDF Download Function - Consolidated per subtab
const downloadSubtabPDF = async (subtabContent: SubtabContent, setIsPDFGenerating: (value: boolean) => void) => {
  console.log('🚀 PDF Download Started for subtab:', subtabContent.title);
  
  try {
    setIsPDFGenerating(true);
    
    // Show loading notification
    const showNotification = (message: string, bgColor: string = '#2563EB') => {
      // Remove any existing notifications
      const existingNotifications = document.querySelectorAll('[data-pdf-notification]');
      existingNotifications.forEach(n => n.remove());
      
      const notification = document.createElement('div');
      notification.setAttribute('data-pdf-notification', 'true');
      notification.innerHTML = `
        <div style="
          position: fixed; 
          top: 20px; 
          right: 20px; 
          background: ${bgColor}; 
          color: white; 
          padding: 15px 20px; 
          border-radius: 8px; 
          box-shadow: 0 4px 12px rgba(0,0,0,0.3);
          z-index: 10000;
          font-weight: 600;
          font-size: 14px;
          max-width: 300px;
        ">
          ${message}
        </div>
      `;
      document.body.appendChild(notification);
      return notification;
    };
    
    showNotification('📄 Generating PDF...');
    
    // Create PDF document
    const pdf = new jsPDF('p', 'mm', 'a4');
        const pageWidth = pdf.internal.pageSize.getWidth();
        const pageHeight = pdf.internal.pageSize.getHeight();
        let yPosition = 20;
        
    // Add header
        pdf.setFontSize(20);
    pdf.setFont('helvetica', 'bold');
    pdf.text(subtabContent.title, pageWidth / 2, yPosition, { align: 'center' });
        yPosition += 15;
        
        pdf.setFontSize(14);
    pdf.setFont('helvetica', 'normal');
        pdf.text('GeneKnow AI Genomic Analysis Report', pageWidth / 2, yPosition, { align: 'center' });
        yPosition += 10;
        
        pdf.setFontSize(12);
        pdf.text(`Generated: ${new Date().toLocaleString()}`, pageWidth / 2, yPosition, { align: 'center' });
        yPosition += 20;
    
    // Capture and add content from each element
    for (const elementId of subtabContent.elementIds) {
        const element = document.getElementById(elementId);
      if (!element) {
        console.warn(`Element ${elementId} not found, skipping`);
        continue;
      }
      
      try {
        // Add section title
        pdf.setFontSize(16);
        pdf.setFont('helvetica', 'bold');
        const sectionTitle = element.querySelector('h3')?.textContent || elementId.replace(/-/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
        
        // Check if we need a new page (more space needed for gene significance analysis)
        const spaceNeeded = elementId.includes('gene-significance-analysis') ? 150 : 60;
        if (yPosition > pageHeight - spaceNeeded) {
          pdf.addPage();
          yPosition = 20;
        }
        
        pdf.text(sectionTitle, 20, yPosition);
        yPosition += 12;
        
        // Special handling for detected-variants - create a text-based table
        if (elementId === 'detected-variants') {
          console.log(`📊 Creating text-based table for detected variants`);
          
          // Extract variant data from the table in the DOM
          const table = element.querySelector('table');
          const variantData: VariantTableRow[] = [];
          
          if (table) {
            const rows = table.querySelectorAll('tbody tr');
            rows.forEach((row) => {
              const cells = row.querySelectorAll('td');
              if (cells.length >= 7) {
                // Extract quality score from the colored div
                const qualityCell = cells[4];
                const qualityDiv = qualityCell.querySelector('div');
                const qualityScore = qualityDiv?.textContent?.trim() || 'N/A';
                
                // Extract clinical significance from the colored div
                const significanceCell = cells[5];
                const significanceDiv = significanceCell.querySelector('div');
                const significance = significanceDiv?.textContent?.trim() || 'N/A';
                
                // Extract impact with both main text and sub-text
                const impactCell = cells[6];
                const mainImpact = impactCell.querySelector('div:first-child strong')?.textContent?.trim() || '';
                const subImpact = impactCell.querySelector('div:last-child')?.textContent?.trim() || '';
                
                variantData.push({
                  gene: cells[0]?.textContent?.trim() || 'N/A',
                  variant: cells[1]?.textContent?.trim() || 'N/A',
                  type: cells[2]?.textContent?.trim() || 'N/A',
                  quality: qualityScore,
                  significance: significance,
                  mainImpact: mainImpact,
                  subImpact: subImpact
                });
              }
            });
          }
          
          if (variantData.length === 0) {
            // Add "no variants found" message
            pdf.setFontSize(12);
            pdf.setFont('helvetica', 'normal');
            pdf.text('No variants detected in the analysis.', 20, yPosition);
            yPosition += 15;
          } else {
            // Simplified column structure - removed problematic DNA/Protein Change columns
            const headers = ['Gene', 'Variant', 'Type', 'Quality', 'Significance', 'Impact'];
            const colWidths = [22, 35, 18, 18, 25, 52]; // Total: 170mm
            const startX = 20;
            const baseRowHeight = 8;
            
            // Helper function to calculate required height for multi-line text
            const calculateCellHeight = (text: string, width: number, fontSize: number = 8): number => {
              pdf.setFontSize(fontSize);
              const lines = pdf.splitTextToSize(text, width - 4);
              return Math.max(baseRowHeight, lines.length * 3 + 4);
            };
            
            // Helper function to draw multi-line text in a cell
            const drawCellText = (text: string, x: number, y: number, width: number, height: number, fontSize: number = 8) => {
              pdf.setFontSize(fontSize);
              const lines = pdf.splitTextToSize(text, width - 4);
              let lineY = y + 3;
              
              lines.forEach((line: string) => {
                if (lineY < y + height - 1) { // Ensure we don't overflow the cell
                  pdf.text(line, x + 2, lineY);
                  lineY += 3;
                }
              });
            };
            
            // Draw header row
            let currentX = startX;
            pdf.setFontSize(9);
            pdf.setFont('helvetica', 'bold');
            pdf.setFillColor(249, 250, 251); // Light gray background
            pdf.rect(startX, yPosition, colWidths.reduce((sum, width) => sum + width, 0), 10, 'F');
            
            for (let i = 0; i < headers.length; i++) {
              drawCellText(headers[i], currentX, yPosition, colWidths[i], 10, 9);
              currentX += colWidths[i];
            }
            yPosition += 10;
            
            // Draw data rows with proper text wrapping
            pdf.setFont('helvetica', 'normal');
            
            variantData.slice(0, 20).forEach((variant: VariantTableRow, index: number) => {
              // Check if we need a new page
              if (yPosition > pageHeight - 60) {
                pdf.addPage();
                yPosition = 20;
                
                // Re-draw headers on new page
                currentX = startX;
                pdf.setFontSize(9);
                pdf.setFont('helvetica', 'bold');
                pdf.setFillColor(249, 250, 251);
                pdf.rect(startX, yPosition, colWidths.reduce((sum, width) => sum + width, 0), 10, 'F');
                
                for (let i = 0; i < headers.length; i++) {
                  drawCellText(headers[i], currentX, yPosition, colWidths[i], 10, 9);
                  currentX += colWidths[i];
                }
                yPosition += 10;
                pdf.setFont('helvetica', 'normal');
              }
              
              // Prepare row data with clean formatting (removed problematic DNA/Protein columns)
              const cellData = [
                variant.gene,
                variant.variant.length > 20 ? variant.variant.substring(0, 20) + '...' : variant.variant,
                variant.type,
                variant.quality,
                variant.significance.replace('_', ' '),
                variant.mainImpact + (variant.subImpact ? '\n' + variant.subImpact : '')
              ];
              
              // Calculate row height based on the tallest cell
              let maxRowHeight = baseRowHeight;
              cellData.forEach((text, i) => {
                const cellHeight = calculateCellHeight(text.toString(), colWidths[i]);
                maxRowHeight = Math.max(maxRowHeight, cellHeight);
              });
              
              // Ensure minimum readable height
              maxRowHeight = Math.max(maxRowHeight, 12);
              
              // Draw row background (alternate colors)
              if (index % 2 === 1) {
                pdf.setFillColor(249, 250, 251);
                pdf.rect(startX, yPosition, colWidths.reduce((sum, width) => sum + width, 0), maxRowHeight, 'F');
              }
              
              // Draw cell borders and content
              currentX = startX;
              cellData.forEach((text, i) => {
                // Draw cell border
                pdf.setDrawColor(230, 230, 230);
                pdf.rect(currentX, yPosition, colWidths[i], maxRowHeight);
                
                // Draw cell content
                pdf.setTextColor(0, 0, 0);
                drawCellText(text.toString(), currentX, yPosition, colWidths[i], maxRowHeight, 8);
                
                currentX += colWidths[i];
              });
              
              yPosition += maxRowHeight;
            });
            
            // Add summary if we truncated
            if (variantData.length > 20) {
              yPosition += 10;
              pdf.setFontSize(10);
              pdf.setFont('helvetica', 'italic');
              pdf.setTextColor(100, 100, 100);
              pdf.text(`Note: Showing first 20 of ${variantData.length} total variants.`, 20, yPosition);
              pdf.text('See full analysis in the application for complete variant list.', 20, yPosition + 5);
              yPosition += 20;
            }
          }
          
          console.log(`✅ Generated text-based variants table with improved formatting`);
        } else {
          // Original image capture for other elements
          const canvas = await html2canvas(element, {
            logging: false,
            // Only apply specific settings to gene significance analysis
            ...(elementId.includes('gene-significance-analysis') ? {
              backgroundColor: '#FFFFFF',
              width: element.scrollWidth,
              height: element.scrollHeight
            } : {
              // For other elements, use transparent background and let browser handle it naturally
              backgroundColor: 'transparent'
            })
          });
          
          const imageData = canvas.toDataURL('image/png', 0.9);
          
          // Calculate optimal size for PDF with better scaling for different content types
          const maxWidth = pageWidth - 40;
          const maxHeight = pageHeight - yPosition - 40;
          const aspectRatio = canvas.width / canvas.height;
          
          let imgWidth: number;
          let imgHeight: number;
        
          // Special handling for gene significance analysis - use full available width
          if (elementId.includes('gene-significance-analysis')) {
            console.log(`📊 Auto-sizing gene significance analysis to full PDF width`);
            console.log(`📐 Canvas dimensions: ${canvas.width}x${canvas.height}, aspect ratio: ${aspectRatio.toFixed(2)}`);
            console.log(`📄 Available space: ${maxWidth}x${maxHeight}`);
            
            // For gene significance analysis, prioritize readability over fitting perfectly
            // Use at least 80% of page width, but allow it to be larger if needed
            const minDesiredWidth = maxWidth * 0.8;
            const maxDesiredWidth = maxWidth;
            
            // Calculate height based on desired width
            const desiredHeight = minDesiredWidth / aspectRatio;
            
            // If the minimum desired size fits, use it
            if (desiredHeight <= maxHeight) {
              imgWidth = minDesiredWidth;
              imgHeight = desiredHeight;
              console.log(`✅ Using minimum desired size: ${imgWidth.toFixed(1)}x${imgHeight.toFixed(1)}`);
            } else {
              // If minimum size is too tall, try to find the largest size that fits
              imgHeight = maxHeight * 0.9; // Use 90% of available height to leave some margin
              imgWidth = imgHeight * aspectRatio;
              
              // Ensure we don't exceed page width
              if (imgWidth > maxDesiredWidth) {
                imgWidth = maxDesiredWidth;
                imgHeight = imgWidth / aspectRatio;
              }
              
              console.log(`⚠️ Using height-constrained size: ${imgWidth.toFixed(1)}x${imgHeight.toFixed(1)}`);
            }
            
            // Ensure minimum readable size - gene significance should be at least 120mm wide
            const minReadableWidth = 120;
            if (imgWidth < minReadableWidth) {
              console.log(`📏 Enforcing minimum readable width of ${minReadableWidth}mm`);
              imgWidth = Math.min(minReadableWidth, maxWidth);
              imgHeight = imgWidth / aspectRatio;
              
              // If this makes it too tall, we'll let it exceed the height limit for readability
              if (imgHeight > maxHeight) {
                console.log(`📏 Gene significance exceeds height limit for readability - allowing overflow`);
              }
            }
            
            console.log(`✅ Gene significance final size: ${imgWidth.toFixed(1)}x${imgHeight.toFixed(1)}`);
            
            // For gene significance analysis, if it's large, consider adding a new page for better presentation
            if (imgHeight > 100 && (yPosition + imgHeight) > (pageHeight - 40)) {
              console.log(`📄 Adding new page for better gene significance presentation`);
              pdf.addPage();
              yPosition = 20;
              // Re-add section title on new page
              pdf.setFontSize(16);
              pdf.setFont('helvetica', 'bold');
              pdf.text(sectionTitle, 20, yPosition);
              yPosition += 12;
            }
          } else {
            // Use scaling factors for other content types
            let scaleFactor = 0.3; // Default scale factor (increased from 0.2)
            
            // Special handling for other chart/visualization elements
            if (elementId.includes('significance') || elementId.includes('matrix') || elementId.includes('chart')) {
              scaleFactor = 0.5; // Higher scale for charts
            }
            
            // For very wide elements (like tables), use adaptive scaling
            if (canvas.width > 1000) {
              scaleFactor = Math.min(0.6, maxWidth / canvas.width);
            }
            
            imgWidth = Math.min(maxWidth, canvas.width * scaleFactor);
            imgHeight = imgWidth / aspectRatio;
            
            // Ensure minimum size for important visualizations
            const minWidth = elementId.includes('significance') ? 120 : 80;
            if (imgWidth < minWidth) {
              imgWidth = Math.min(minWidth, maxWidth);
              imgHeight = imgWidth / aspectRatio;
            }
            
            // If height exceeds available space, scale down
            if (imgHeight > maxHeight) {
              imgHeight = maxHeight;
              imgWidth = imgHeight * aspectRatio;
            }
          }
                
          // Special handling for gene significance analysis - use full available width
          if (elementId.includes('gene-significance-analysis')) {
            console.log(`📊 Auto-sizing gene significance analysis to full PDF width`);
            console.log(`📐 Canvas dimensions: ${canvas.width}x${canvas.height}, aspect ratio: ${aspectRatio.toFixed(2)}`);
            console.log(`📄 Available space: ${maxWidth}x${maxHeight}`);
            
            // For gene significance analysis, prioritize readability over fitting perfectly
            // Use at least 80% of page width, but allow it to be larger if needed
            const minDesiredWidth = maxWidth * 0.8;
            const maxDesiredWidth = maxWidth;
            
            // Calculate height based on desired width
            const desiredHeight = minDesiredWidth / aspectRatio;
            
            // If the minimum desired size fits, use it
            if (desiredHeight <= maxHeight) {
              imgWidth = minDesiredWidth;
              imgHeight = desiredHeight;
              console.log(`✅ Using minimum desired size: ${imgWidth.toFixed(1)}x${imgHeight.toFixed(1)}`);
            } else {
              // If minimum size is too tall, try to find the largest size that fits
              imgHeight = maxHeight * 0.9; // Use 90% of available height to leave some margin
              imgWidth = imgHeight * aspectRatio;
              
              // Ensure we don't exceed page width
              if (imgWidth > maxDesiredWidth) {
                imgWidth = maxDesiredWidth;
                imgHeight = imgWidth / aspectRatio;
              }
              
              console.log(`⚠️ Using height-constrained size: ${imgWidth.toFixed(1)}x${imgHeight.toFixed(1)}`);
            }
            
            // Ensure minimum readable size - gene significance should be at least 120mm wide
            const minReadableWidth = 120;
            if (imgWidth < minReadableWidth) {
              console.log(`📏 Enforcing minimum readable width of ${minReadableWidth}mm`);
              imgWidth = Math.min(minReadableWidth, maxWidth);
              imgHeight = imgWidth / aspectRatio;
              
              // If this makes it too tall, we'll let it exceed the height limit for readability
              if (imgHeight > maxHeight) {
                console.log(`📏 Gene significance exceeds height limit for readability - allowing overflow`);
              }
            }
            
            console.log(`✅ Gene significance final size: ${imgWidth.toFixed(1)}x${imgHeight.toFixed(1)}`);
            
            // For gene significance analysis, if it's large, consider adding a new page for better presentation
            if (imgHeight > 100 && (yPosition + imgHeight) > (pageHeight - 40)) {
              console.log(`📄 Adding new page for better gene significance presentation`);
              pdf.addPage();
              yPosition = 20;
              // Re-add section title on new page
              pdf.setFontSize(16);
              pdf.setFont('helvetica', 'bold');
              pdf.text(sectionTitle, 20, yPosition);
              yPosition += 12;
            }
          } else {
            // Use scaling factors for other content types
            let scaleFactor = 0.3; // Default scale factor (increased from 0.2)
            
            // Special handling for other chart/visualization elements
            if (elementId.includes('significance') || elementId.includes('matrix') || elementId.includes('chart')) {
              scaleFactor = 0.5; // Higher scale for charts
            }
            
            // For very wide elements (like tables), use adaptive scaling
            if (canvas.width > 1000) {
              scaleFactor = Math.min(0.6, maxWidth / canvas.width);
            }
            
            imgWidth = Math.min(maxWidth, canvas.width * scaleFactor);
            imgHeight = imgWidth / aspectRatio;
            
            // Ensure minimum size for important visualizations
            const minWidth = elementId.includes('significance') ? 120 : 80;
            if (imgWidth < minWidth) {
              imgWidth = Math.min(minWidth, maxWidth);
              imgHeight = imgWidth / aspectRatio;
            }
            
            // If height exceeds available space, scale down
            if (imgHeight > maxHeight) {
              imgHeight = maxHeight;
              imgWidth = imgHeight * aspectRatio;
            }
          }
                
          // Add the image to PDF
          pdf.addImage(imageData, 'PNG', 20, yPosition, imgWidth, imgHeight);
          yPosition += imgHeight + 15;
        }
              
        console.log(`✅ Added section: ${sectionTitle}`);
        
      } catch (error) {
        console.warn(`Failed to capture element ${elementId}:`, error);
        
        // Add fallback text
            pdf.setFontSize(11);
        pdf.setFont('helvetica', 'normal');
        pdf.text(`[Visualization from ${elementId} - see application for full details]`, 20, yPosition);
          yPosition += 10;
      }
    }
    
    // Add footer
    const footerY = pageHeight - 15;
    pdf.setFontSize(8);
    pdf.setFont('helvetica', 'normal');
    pdf.text('GeneKnow Platform - AI-powered genomic analysis', pageWidth / 2, footerY - 5, { align: 'center' });
    
    // Format the timestamp nicely
    const now = new Date();
    const formattedDate = now.toLocaleDateString('en-US', { 
      year: 'numeric', 
      month: 'long', 
      day: 'numeric',
      hour: '2-digit',
      minute: '2-digit'
    });
    pdf.text(`Report ID: GK-${Date.now()} | Generated: ${formattedDate}`, pageWidth / 2, footerY, { align: 'center' });
    
    // Generate filename
    const timestamp = new Date().toISOString().split('T')[0];
    const fileName = `geneknow_${subtabContent.id}_${timestamp}.pdf`;
    
    // Get PDF as blob and convert to bytes
    const pdfBlob = pdf.output('blob');
    const pdfBytes = await pdfBlob.arrayBuffer();
    const pdfByteArray = Array.from(new Uint8Array(pdfBytes));
    
    // Use Tauri save_file_dialog to save the PDF
    try {
      const savedPath = await invoke<string>('save_file_dialog', {
        filename: fileName,
        fileContent: pdfByteArray
      });
      
      console.log('✅ PDF saved successfully to:', savedPath);
        
        // Show success notification
      showNotification(`✅ PDF Downloaded!<br><small>Saved successfully</small>`, '#22C55E');
        setTimeout(() => {
        const notifications = document.querySelectorAll('[data-pdf-notification]');
        notifications.forEach(n => n.remove());
      }, 4000);
      
    } catch (saveError) {
      console.error('❌ Error saving PDF:', saveError);
      throw new Error(`Failed to save PDF: ${saveError}`);
    }
        
        console.log('🎉 PDF generation completed successfully!');
        
  } catch (error) {
    console.error('❌ Error generating PDF:', error);
      
      // Show error notification
    const showErrorNotification = (message: string) => {
      const notifications = document.querySelectorAll('[data-pdf-notification]');
      notifications.forEach(n => n.remove());
      
      const notification = document.createElement('div');
      notification.setAttribute('data-pdf-notification', 'true');
      notification.innerHTML = `
        <div style="
          position: fixed; 
          top: 20px; 
          right: 20px; 
          background: #EF4444; 
          color: white; 
          padding: 15px 20px; 
          border-radius: 8px; 
          box-shadow: 0 4px 12px rgba(0, 0, 0, 0.3);
          z-index: 10000;
          font-weight: 600;
          font-size: 14px;
          max-width: 300px;
        ">
          ❌ ${message}<br>
          <small style="font-weight: 400;">Please try again</small>
        </div>
      `;
      document.body.appendChild(notification);
      
      setTimeout(() => {
        if (notification.parentNode) {
          notification.parentNode.removeChild(notification);
        }
      }, 4000);
    };
    
    showErrorNotification('PDF generation failed');
  } finally {
    setIsPDFGenerating(false);
  }
};

// Consolidated Download Button Component - One per subtab
const SubtabDownloadButton = ({ subtabContent, isPDFGenerating, setIsPDFGenerating }: { 
  subtabContent: SubtabContent; 
  isPDFGenerating: boolean; 
  setIsPDFGenerating: (value: boolean) => void;
}) => {
  if (isPDFGenerating) {
    return (
      <div style={{
        display: 'inline-flex',
        alignItems: 'center',
        padding: '0.5rem 1rem',
        background: '#9CA3AF',
        color: '#FFFFFF',
        border: 'none',
        borderRadius: '0.375rem',
        fontSize: '0.875rem',
        fontWeight: '500',
        gap: '0.5rem'
      }}>
        <div style={{
          width: '16px',
          height: '16px',
          border: '2px solid #ffffff',
          borderTop: '2px solid transparent',
          borderRadius: '50%',
          animation: 'spin 1s linear infinite'
        }}></div>
        Generating PDF...
      </div>
    );
  }
  
  return (
    <button
              onClick={() => downloadSubtabPDF(subtabContent, setIsPDFGenerating)}
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
      PDF Download
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

interface QualityMetrics {
  quality?: number;
  depth?: number;
  allele_freq?: number;
}

interface PathwayResult {
  burden_score?: number;
  contributing_genes?: string[];
  damaging_genes?: string[];
  gene_damaging_counts?: Record<string, number>;
}

interface VariantTransformation {
  original?: string;
  mutated?: string;
  amino_acid_change?: string;
  effect?: string;
}

interface ExtendedVariant {
  transformation?: VariantTransformation;
  functional_impact?: string;
  [key: string]: unknown;
}

interface ClinicalRecommendation {
  title: string;
  description: string;
  priority: string;
  category: string;
  cancer_type: string;
  risk_level: string;
  risk_percentage: number;
  recommendation: string;
  screening_protocol?: {
    test: string;
    frequency: string;
    start_age: string;
  };
  prevention_options?: string[];
}

interface PathwayMutation {
  gene: string;
  type: string;
  effect: string;
}

interface DisruptedPathway {
  name: string;
  pathway_id: string;
  significance: number;
  affected_genes: string[];
  mutations: PathwayMutation[];
  description?: string;
  genes_affected_ratio?: string;
}

interface PathwayAnalysisSummary {
  total_pathways_disrupted: number;
  highly_disrupted_pathways: number;
  total_genes_affected: number;
  pathway_interaction_count: number;
  overall_burden_score: number;
  high_burden_pathways: string[];
}

interface PathwayAnalysisData {
  disrupted_pathways: DisruptedPathway[];
  cancer_pathway_associations: Record<string, string[]>;
  summary?: PathwayAnalysisSummary;
}

interface PathwayBurdenResult {
  burden_score: number;
  damaging_genes?: string[];
  genes_with_damaging?: number;
  genes_in_pathway?: number;
  description?: string;
}

interface PipelineWarning {
  message?: string;
  [key: string]: unknown;
}

type WarningType = string | PipelineWarning;

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

// Enhanced tooltip component with smart positioning
const SmartTooltip = ({ content, isVisible, triggerRef }: { 
  content: string; 
  isVisible: boolean; 
  triggerRef?: React.RefObject<HTMLDivElement> | null;
}) => {
  const [position, setPosition] = useState<'right' | 'above' | 'below'>('right');
  const tooltipRef = useRef<HTMLDivElement>(null);
  
  useEffect(() => {
    if (isVisible && triggerRef?.current && tooltipRef.current) {
      const triggerRect = triggerRef.current.getBoundingClientRect();
      // const tooltipRect = tooltipRef.current.getBoundingClientRect(); // Removed: unused variable
      const viewportWidth = window.innerWidth;
      const viewportHeight = window.innerHeight;
      
      // Check if there's enough space to the right
      const spaceToRight = viewportWidth - triggerRect.right;
      const tooltipWidth = 288; // 18rem = 288px
      
      // Check if there's enough space above/below
      const spaceAbove = triggerRect.top;
      const spaceBelow = viewportHeight - triggerRect.bottom;
      const tooltipHeight = 100; // Approximate height
      
      // Determine best position
      if (spaceToRight >= tooltipWidth + 16) {
        // Enough space to the right
        setPosition('right');
      } else if (spaceAbove >= tooltipHeight + 16) {
        // Not enough space to right, but enough space above
        setPosition('above');
      } else if (spaceBelow >= tooltipHeight + 16) {
        // Not enough space to right or above, try below
        setPosition('below');
      } else {
        // Default to right if no good position found
        setPosition('right');
      }
    }
  }, [isVisible, triggerRef]);
  
  const getTooltipStyles = () => {
    const baseStyles = {
      position: 'absolute' as const,
      width: '18rem',
      padding: '0.75rem',
      background: '#1F2937',
      color: '#FFFFFF',
      fontSize: '0.75rem',
      borderRadius: '0.5rem',
      boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
      opacity: isVisible ? 1 : 0,
      visibility: isVisible ? 'visible' as const : 'hidden' as const,
      transition: 'opacity 300ms ease, visibility 300ms ease',
      zIndex: 1000,
      pointerEvents: 'none' as const,
      lineHeight: '1.4',
      border: '1px solid #374151'
    };
    
    switch (position) {
      case 'right':
        return {
          ...baseStyles,
          left: 'calc(100% + 0.5rem)',
          top: '50%',
          transform: 'translateY(-50%)',
        };
      case 'above':
        return {
          ...baseStyles,
          left: '50%',
          bottom: 'calc(100% + 0.5rem)',
          transform: 'translateX(-50%)',
        };
      case 'below':
        return {
          ...baseStyles,
          left: '50%',
          top: 'calc(100% + 0.5rem)',
          transform: 'translateX(-50%)',
        };
      default:
        return baseStyles;
    }
  };
  
  const getArrowStyles = () => {
    const baseArrowStyles = {
      position: 'absolute' as const,
      width: '0',
      height: '0',
    };
    
    switch (position) {
      case 'right':
        return {
          ...baseArrowStyles,
          left: '-0.5rem',
          top: '50%',
          transform: 'translateY(-50%)',
          borderTop: '0.5rem solid transparent',
          borderBottom: '0.5rem solid transparent',
          borderRight: '0.5rem solid #1F2937'
        };
      case 'above':
        return {
          ...baseArrowStyles,
          left: '50%',
          bottom: '-0.5rem',
          transform: 'translateX(-50%)',
          borderLeft: '0.5rem solid transparent',
          borderRight: '0.5rem solid transparent',
          borderTop: '0.5rem solid #1F2937'
        };
      case 'below':
        return {
          ...baseArrowStyles,
          left: '50%',
          top: '-0.5rem',
          transform: 'translateX(-50%)',
          borderLeft: '0.5rem solid transparent',
          borderRight: '0.5rem solid transparent',
          borderBottom: '0.5rem solid #1F2937'
        };
      default:
        return baseArrowStyles;
    }
  };
  
  return (
    <div ref={tooltipRef} style={getTooltipStyles()}>
      <p style={{ color: '#D1D5DB', marginBottom: '0' }}>{content}</p>
      <div style={getArrowStyles()}></div>
    </div>
  );
};

// Legacy SectionTooltip component removed - was unused

// Function to generate real clinical alerts from pipeline data
const generateClinicalAlerts = (pipelineResults: PipelineResult | undefined): Alert[] => {
  console.log('🔍 ALERT DEBUG: generateClinicalAlerts called with:', !!pipelineResults);
  if (!pipelineResults) {
    console.log('🔍 ALERT DEBUG: No pipeline results, returning empty array');
    return [];
  }

  const alerts: Alert[] = [];
  const structuredJson = pipelineResults.structured_json;
  console.log('🔍 ALERT DEBUG: structuredJson available:', !!structuredJson);
  console.log('🔍 ALERT DEBUG: structuredJson keys:', structuredJson ? Object.keys(structuredJson) : 'none');

  // 1. Check for pathogenic and likely pathogenic variants
  const variantDetails = structuredJson?.variant_details || [];
  const pathogenicVariants = variantDetails.filter(variant => 
    variant.clinical_significance === 'pathogenic' || 
    variant.clinical_significance === 'likely_pathogenic'
  );

  // Generate alerts for pathogenic variants
  pathogenicVariants.forEach(variant => {
    const isPathogenic = variant.clinical_significance === 'pathogenic';
    const gene = variant.gene || 'Unknown';
    const variantId = variant.variant_id || variant.variant || 'Unknown';
    
    // Get protein change or use variant ID
    const proteinChange = variant.protein_change || variant.hgvs_p || variantId;
    
    alerts.push({
      type: isPathogenic ? 'critical' : 'warning',
      title: `${gene} ${isPathogenic ? 'Pathogenic' : 'Likely Pathogenic'} Variant`,
      desc: `${proteinChange} - ${variant.consequence || 'mutation'} detected`,
      detailedInfo: {
        whatItMeans: `A ${isPathogenic ? 'pathogenic' : 'likely pathogenic'} variant in the ${gene} gene has been identified. This variant is ${isPathogenic ? 'known to cause' : 'likely to cause'} disease or increase cancer risk.`,
        whyImportant: `${gene} variants can significantly impact cancer risk. ${isPathogenic ? 'Pathogenic' : 'Likely pathogenic'} variants have strong evidence for clinical significance and require immediate attention.`,
        clinicalSignificance: variant.functional_impact || `${isPathogenic ? 'High' : 'Moderate'} research significance. This variant requires medical evaluation and may affect treatment decisions.`,
        nextSteps: `Genetic counseling is recommended. Consider enhanced screening protocols and discuss with oncology team. Family members should be informed about potential hereditary risk.`
      }
    });
  });

  // 2. Check for variants of uncertain significance (VUS) in important genes
  const importantGenes = ['TP53', 'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'APC'];
  const vusInImportantGenes = variantDetails.filter(variant => 
    variant.clinical_significance === 'uncertain_significance' && 
    importantGenes.includes(variant.gene)
  );

  vusInImportantGenes.forEach(variant => {
    const gene = variant.gene || 'Unknown';
    const proteinChange = variant.protein_change || variant.hgvs_p || variant.variant_id || 'Unknown';
    
    alerts.push({
      type: 'info',
      title: `${gene} Variant of Uncertain Significance`,
      desc: `${proteinChange} - Research significance unclear`,
      detailedInfo: {
        whatItMeans: `A genetic variant in the ${gene} gene where the research significance is not yet definitively established by current scientific evidence.`,
        whyImportant: `${gene} is an important cancer-related gene. Even uncertain variants require monitoring as new research may reclassify them as pathogenic or benign.`,
        clinicalSignificance: 'May or may not be medically significant. Requires ongoing evaluation as more research data becomes available.',
        nextSteps: 'Genetic counseling to discuss implications, potential for variant reclassification, and consideration of family testing.'
      }
    });
  });

  // 3. Check SHAP validation results for model reliability alerts
  const shapValidation = structuredJson?.shap_validation;
  if (shapValidation) {
    if (shapValidation.status === 'FLAG_FOR_REVIEW') {
      alerts.push({
        type: 'warning',
        title: 'AI Model Review Required',
        desc: 'Machine learning prediction requires manual verification',
        detailedInfo: {
          whatItMeans: 'The AI model has flagged this analysis for review due to potential inconsistencies in the prediction logic or unusual patterns in the data.',
          whyImportant: 'This flag indicates the model prediction may be less reliable than usual. Manual review by a healthcare provider is recommended to ensure accuracy.',
          clinicalSignificance: `Review reasons: ${shapValidation.reasons?.join('; ') || 'Model validation concerns detected'}`,
          nextSteps: 'Have a healthcare provider review the analysis results. Consider additional testing or consultation with a genetic counselor.'
        }
      });
    } else if (shapValidation.status === 'ERROR') {
      alerts.push({
        type: 'critical',
        title: 'Analysis Validation Error',
        desc: 'Technical issue with risk prediction model',
        detailedInfo: {
          whatItMeans: 'The analysis encountered a technical error during the validation process. This may affect the reliability of the risk predictions.',
          whyImportant: 'Validation errors can indicate problems with the analysis that may impact the accuracy of the results.',
          clinicalSignificance: 'Results should be interpreted with caution due to validation errors.',
          nextSteps: 'Consult with a healthcare provider for proper interpretation. Consider re-running the analysis or seeking additional genetic testing.'
        }
      });
    }
  }

  // 4. Check for elevated cancer risk (>20% for clinical significance)
  const riskScores = pipelineResults.risk_scores || {};
  console.log('🔍 ALERT DEBUG: All risk scores:', riskScores);
  
  const elevatedRiskCancers = Object.entries(riskScores).filter(([, score]) => score >= 20);
  console.log('🔍 ALERT DEBUG: Elevated risk cancers (>=20%):', elevatedRiskCancers);
  
  if (elevatedRiskCancers.length > 0) {
    const topRisk = elevatedRiskCancers.reduce((prev, [cancer, score]) => 
      score > prev.score ? { cancer, score } : prev,
      { cancer: '', score: 0 }
    );

    console.log('🔍 ALERT DEBUG: Adding elevated risk alert for:', topRisk);
    
    alerts.push({
      type: topRisk.score >= 30 ? 'critical' : 'warning',
      title: `Elevated ${topRisk.cancer} Cancer Risk`,
      desc: `${topRisk.score.toFixed(1)}% risk detected - Enhanced screening recommended`,
      detailedInfo: {
        whatItMeans: `Based on the genetic analysis, you have been identified as having elevated risk for ${topRisk.cancer} cancer (${topRisk.score.toFixed(1)}%).`,
        whyImportant: `This elevated risk is significantly higher than the general population average. Early detection and prevention strategies can significantly improve outcomes.`,
        clinicalSignificance: `This risk level warrants ${topRisk.score >= 30 ? 'immediate medical attention' : 'enhanced screening protocols'} and specialized medical management.`,
        nextSteps: `Schedule consultation with oncology team. Begin enhanced screening protocols. Consider genetic counseling to discuss risk management strategies and family implications.`
      }
    });
  }

  // 5. Check for quality or pipeline warnings
  const warnings = structuredJson?.warnings || [];
  if (Array.isArray(warnings) && warnings.length > 0) {
    // Only show the most relevant warnings (limit to avoid cluttering)
    const relevantWarnings = warnings.slice(0, 2);
    
    relevantWarnings.forEach((warning: WarningType) => {
      alerts.push({
        type: 'info',
        title: 'Analysis Quality Notice',
        desc: typeof warning === 'string' ? warning : (warning?.message || 'Quality check flagged'),
        detailedInfo: {
          whatItMeans: 'A quality check during the analysis process has identified a potential issue or limitation in the data.',
          whyImportant: 'Quality issues can affect the reliability of the analysis results. Being aware of these limitations helps with proper interpretation.',
          clinicalSignificance: 'Results should be interpreted considering these quality factors.',
          nextSteps: 'Discuss these quality considerations with your healthcare provider. Additional testing may be recommended if needed.'
        }
      });
    });
  }

  // 6. Default success alert for low-risk cases with no significant findings
  if (alerts.length === 0) {
    const overallRiskScore = Math.max(...Object.values(riskScores));
    console.log('🔍 ALERT DEBUG: No alerts generated, adding default success alert. Overall risk score:', overallRiskScore);
    
    alerts.push({
      type: 'success',
      title: 'Low Risk Assessment',
      desc: 'No significant genetic risk factors identified',
      detailedInfo: {
        whatItMeans: 'Current genetic analysis does not indicate significantly elevated risk for hereditary cancer syndromes based on the variants analyzed.',
        whyImportant: 'This provides reassurance while maintaining appropriate vigilance for any changes in risk factors or family history.',
        clinicalSignificance: 'Standard population screening guidelines are appropriate at this time.',
        nextSteps: 'Continue routine screening as recommended for general population. Maintain awareness of family history changes and discuss with healthcare provider.'
      }
    });
  }

  console.log('🔍 ALERT DEBUG: Final alerts being returned:', alerts);
  console.log('🔍 ALERT DEBUG: Total alert count:', alerts.length);
  return alerts;
};

// Mock data sets for different risk levels - completely anonymous
const mockDataSets = {
  high: {
    riskLevel: 'High Risk',
    riskScore: '82/100',
    condition: 'Hereditary Breast and Ovarian Cancer Syndrome',
    details: 'Analysis Type: Genomic Variant Analysis<br/>Method: Machine Learning Risk Assessment<br/>Data Source: Uploaded VCF File',
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
    details: 'Analysis Type: Genomic Variant Analysis<br/>Method: Machine Learning Risk Assessment<br/>Data Source: Uploaded VCF File',
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
    details: 'Analysis Type: Genomic Variant Analysis<br/>Method: Machine Learning Risk Assessment<br/>Data Source: Uploaded VCF File',
    alerts: [
      {
        type: 'info' as const,
        title: 'TP53 Variant of Uncertain Significance',
        desc: 'c.743G>A - Missense mutation',
        detailedInfo: {
          whatItMeans: 'A genetic variant in the TP53 gene where the research significance is not yet definitively established by current scientific evidence.',
          whyImportant: 'TP53 is the "guardian of the genome" and mutations can predispose to multiple cancer types. Even uncertain variants require monitoring.',
          clinicalSignificance: 'May or may not be medically significant. Requires ongoing evaluation as more research data becomes available.',
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
    visibility: isVisible ? 'visible' as const : 'hidden' as const,
    transition: 'opacity 300ms ease, visibility 300ms ease',
    zIndex: 1000,
    pointerEvents: 'none' as const,
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
        Research significance:
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

const ClinicalViewPage: React.FC = () => {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const location = useLocation();
  const [activeTab, setActiveTab] = useState('analysis');
  const [hoveredAlert, setHoveredAlert] = useState<number | null>(null);
  const [isPDFGenerating, setIsPDFGenerating] = useState(false);
  const [hoveredTooltip, setHoveredTooltip] = useState<string | null>(null);
  
  // Add pagination state for variants table
  const [currentPage, setCurrentPage] = useState(1);
  const [variantsPerPage] = useState(20); // Show 20 variants per page
  
  // Check if we have real results from the pipeline
  const pipelineResults = location.state?.results as PipelineResult | undefined;
  const fileName = location.state?.fileName as string | undefined;
  
  // Extract risk level from URL parameter - use same default as dashboard for consistency
  const riskLevel = searchParams.get('risk') || 'low';
  
  // Use real SHAP validation from pipeline results only
  const getConfidenceCheckValidation = () => {
    // Use real SHAP validation if available
    if (pipelineResults?.structured_json?.shap_validation) {
      return pipelineResults.structured_json.shap_validation;
    }
    
    // If no pipeline results or no SHAP validation, return SKIPPED status
    return {
      status: 'SKIPPED' as const,
      reasons: ['No pipeline results available'],
      top_contributors: [],
      feature_importance: {},
      details: {
        status: 'SKIPPED' as const,
        risk_score: 0.0,
        top_contributors: [],
        validation_reasons: ['No pipeline results available'],
        rule_results: {},
        shap_values: [],
        feature_names: [],
        model_type: 'None'
      }
    };
  };
  
  const confidenceCheckValidation = getConfidenceCheckValidation();
  
  // Get mock data for UI elements (only used when no real pipeline results)
  const currentData = mockDataSets[riskLevel as keyof typeof mockDataSets] || mockDataSets.low;

  // Helper function to determine file type from filename
  const getFileTypeDescription = (filename?: string): string => {
    if (!filename) return 'Unknown File Type';
    
    const lowercaseFile = filename.toLowerCase();
    
    if (lowercaseFile.endsWith('.maf') || lowercaseFile.endsWith('.maf.gz')) {
      return 'Uploaded MAF File';
    } else if (lowercaseFile.endsWith('.vcf') || lowercaseFile.endsWith('.vcf.gz')) {
      return 'Uploaded VCF File';
    } else if (lowercaseFile.endsWith('.fastq') || lowercaseFile.endsWith('.fq') || 
               lowercaseFile.endsWith('.fastq.gz') || lowercaseFile.endsWith('.fq.gz')) {
      return 'Uploaded FASTQ File';
    } else if (lowercaseFile.endsWith('.bam')) {
      return 'Uploaded BAM File';
    } else if (lowercaseFile.endsWith('.sam')) {
      return 'Uploaded SAM File';
    } else {
      return 'Uploaded Genomic File';
    }
  };

  // Get real data for sidebar when available
  const getSidebarData = () => {
    if (pipelineResults) {
      console.log('🔍 SIDEBAR DEBUG: pipelineResults available, generating real alerts');
      console.log('🔍 Risk scores:', pipelineResults.risk_scores);
      console.log('🔍 Structured JSON:', pipelineResults.structured_json);
      
      // Calculate real risk level and score from pipeline results
      const riskScores = Object.entries(pipelineResults.risk_scores || {});
      const highestRisk = riskScores.reduce((prev, [cancer, score]) => 
        score > prev.score ? { cancer, score } : prev,
        { cancer: '', score: 0 }
      );
      
      const overallRiskScore = Math.round(highestRisk.score);
      const riskCategory = overallRiskScore >= 30 ? 'High Risk' : 
                          overallRiskScore >= 15 ? 'Moderate Risk' : 'Low Risk';
      
      // Generate real clinical alerts from pipeline data
      const realAlerts = generateClinicalAlerts(pipelineResults);
      console.log('🔍 SIDEBAR DEBUG: Generated alerts:', realAlerts);
      console.log('🔍 SIDEBAR DEBUG: Alert count:', realAlerts.length);
      
      return {
        riskLevel: riskCategory,
        riskScore: `${overallRiskScore}/100`,
        condition: `${highestRisk.cancer} Cancer Risk Assessment`,
        details: `Analysis Type: Genomic Variant Analysis<br/>Method: Machine Learning Risk Assessment<br/>Data Source: ${getFileTypeDescription(fileName)}`,
        alerts: realAlerts
      };
    }
    
    console.log('🔍 SIDEBAR DEBUG: No pipelineResults, falling back to mock data');
    // Fallback to mock data if no pipeline results (but update data source to be accurate)
    return {
      ...currentData,
      details: `Analysis Type: Genomic Variant Analysis<br/>Method: Machine Learning Risk Assessment<br/>Data Source: ${getFileTypeDescription(fileName)}`
    };
  };

  const sidebarData = getSidebarData();

  // Define subtab content for consolidated PDF downloads
  const subtabContents: Record<string, SubtabContent> = {
    analysis: {
      id: 'analysis',
      title: 'Genomic Analysis',
      elementIds: ['cancer-risk-assessment', 'gene-significance-analysis', 'mutation-type-distribution', 'mutational-signatures', 'quality-metrics']
    },
    variants: {
      id: 'variants',
      title: 'Variant Heatmap',
      elementIds: ['gene-cancer-matrix', 'detected-variants']
    },
    pathways: {
      id: 'pathways',
      title: 'Pathway Analysis',
      elementIds: ['pathway-disruption-analysis', 'cancer-risk-by-pathways']
    },
    clinical: {
      id: 'clinical',
      title: 'In-Depth Report',
      elementIds: ['survival-analysis', 'clinical-recommendations']
    }
  };

  // Function to convert pipeline results to the format expected by the clinical view
  const getGenomicData = () => {
    if (pipelineResults) {
      console.log('🔍 Using real pipeline results:', pipelineResults);
      
      // Extract risk scores and findings
      const riskScores = Object.entries(pipelineResults.risk_scores || {});
      const pathwayBurdenResults = pipelineResults.pathway_burden_results || {};
      
      // Build cancer-gene associations from pathway burden data
      const cancerGeneMapping: Record<string, string[]> = {
        'breast': [],
        'lung': [],
        'colon': [],
        'prostate': [],
        'blood': []
      };
      
      // Map pathways to cancer types and collect genes
      const pathwayCancerMapping: Record<string, string[]> = {
        'dna_repair': ['breast', 'colon'],
        'tumor_suppressors': ['breast', 'lung', 'colon', 'prostate'],
        'oncogenes': ['lung', 'colon', 'breast'],
        'ras_mapk': ['lung', 'colon', 'prostate'],
        'cell_cycle': ['breast', 'lung', 'prostate'],
        'apoptosis': ['breast', 'lung', 'colon'],
        'chromatin_remodeling': ['blood', 'lung'],
        'mismatch_repair': ['colon'],
        'wnt_signaling': ['colon'],
        'pi3k_akt': ['breast', 'prostate']
      };
      
      // Populate cancer-gene associations from pathway burden results
      Object.entries(pathwayBurdenResults).forEach(([pathwayName, pathwayResult]) => {
        const associatedCancers = pathwayCancerMapping[pathwayName] || [];
        const burdenScore = pathwayResult.burden_score || 0;
        
        // Only include genes from high-burden pathways
        if (burdenScore > 0.5 && pathwayResult.damaging_genes) {
          pathwayResult.damaging_genes.forEach(gene => {
            associatedCancers.forEach(cancer => {
              if (cancerGeneMapping[cancer] && !cancerGeneMapping[cancer].includes(gene)) {
                cancerGeneMapping[cancer].push(gene);
              }
            });
          });
        }
      });
      
      // Also add genes from multi-pathway associations
      const multiPathwayGenes = pipelineResults.pathway_burden_summary?.multi_pathway_genes || {};
      Object.entries(multiPathwayGenes).forEach(([gene, pathways]) => {
        if (pathways.length > 1) {
          // Gene appears in multiple pathways, add to all relevant cancer types
          pathways.forEach(pathwayName => {
            const associatedCancers = pathwayCancerMapping[pathwayName] || [];
            associatedCancers.forEach(cancer => {
              if (cancerGeneMapping[cancer] && !cancerGeneMapping[cancer].includes(gene)) {
                cancerGeneMapping[cancer].push(gene);
              }
            });
          });
        }
      });
      
      console.log('🧬 Built cancer-gene associations:', cancerGeneMapping);
      
      const riskFindings = riskScores.map(([cancer_type, risk_percentage]) => ({
        cancer_type,
        risk_percentage,
        risk_level: risk_percentage > 50 ? 'high' : risk_percentage > 20 ? 'medium' : 'low',
        affected_genes: cancerGeneMapping[cancer_type] || [], // Use pathway-derived genes
        recommendation: risk_percentage > 50 
          ? `Enhanced ${cancer_type} cancer screening recommended`
          : risk_percentage > 20
          ? `Moderate ${cancer_type} cancer screening recommended`
          : 'Standard screening guidelines apply',
        mutation_burden: risk_percentage > 50 ? 'high' : risk_percentage > 20 ? 'medium' : 'low',
        pathway_disruption: [] // Backend doesn't provide this yet
      }));

      // Extract variant details from structured_json if available
      const structuredJson = pipelineResults.structured_json || {};
      const variantDetails = Array.isArray(structuredJson.variant_details) ? structuredJson.variant_details : [];
      
      // Transform variants to expected format
      const variantTable = variantDetails.map((variant, index) => {
        // Log missing data
        if (!variant.transformation) {
          console.warn(`⚠️ Missing transformation data for variant ${variant.gene}`);
        }
        if (!variant.protein_change) {
          console.warn(`⚠️ Missing protein change for variant ${variant.gene}`);
        }
        
        // Calculate more detailed quality score
        const qualityMetrics = variant.quality_metrics || {} as QualityMetrics;
        const baseQuality = qualityMetrics.quality || 85;
        const depth = qualityMetrics.depth || 0;
        const alleleFreq = qualityMetrics.allele_freq || 0;
        
        // Adjust quality score based on multiple factors
        let adjustedQuality = baseQuality;
        
        // Penalize low depth
        if (depth < 10) {
          adjustedQuality = Math.max(adjustedQuality - 20, 0);
        } else if (depth < 20) {
          adjustedQuality = Math.max(adjustedQuality - 10, 0);
        }
        
        // Penalize extreme allele frequencies
        if (alleleFreq < 0.1 || alleleFreq > 0.9) {
          adjustedQuality = Math.max(adjustedQuality - 15, 0);
        }
        
        // Boost quality for high-confidence clinical significance
        if (variant.clinical_significance === 'pathogenic') {
          adjustedQuality = Math.min(adjustedQuality + 5, 100);
        } else if (variant.clinical_significance === 'likely_pathogenic') {
          adjustedQuality = Math.min(adjustedQuality + 3, 100);
        }
        
        // Add quality details for tooltip
        const qualityDetails = {
          base_quality: baseQuality,
          depth: depth,
          allele_frequency: alleleFreq,
          adjusted_quality: Math.round(adjustedQuality),
          quality_factors: [] as string[]
        };
        
        if (depth < 10) qualityDetails.quality_factors.push('Low read depth');
        if (alleleFreq < 0.1) qualityDetails.quality_factors.push('Low allele frequency');
        if (alleleFreq > 0.9) qualityDetails.quality_factors.push('High allele frequency');
        if (variant.clinical_significance === 'pathogenic') qualityDetails.quality_factors.push('Pathogenic variant');

        // Better handling of TCGA matches
        const tcgaMatch = variant.tcga_matches ? 
          Object.entries(variant.tcga_matches)[0]?.[1] || null : null;
        
        return {
          gene: variant.gene || 'Unknown',
          variant_id: variant.variant || variant.variant_id || `Unknown_${index}`,
          consequence: variant.consequence || 'unknown',
          mutation_type: variant.mutation_type || 'snv', // Default to SNV if not specified
          quality_score: Math.round(adjustedQuality),
          quality_details: qualityDetails,
          clinical_significance: variant.clinical_significance || 'uncertain_significance',
          tcga_best_match: tcgaMatch || { cancer_type: 'not_cancer_associated', frequency: 0 },
          protein_change: variant.protein_change || variant.hgvs_p || `p.Unknown${index}`,
          functional_impact: variant.functional_impact || 'Unknown',
          // Use real transformation data if available, otherwise provide placeholder
          transformation: variant.transformation || {
            original: 'N/A',
            mutated: 'N/A',
            amino_acid_change: 'Unknown',
            effect: 'Effect not determined'
          }
        };
      });

      // Log what backend data we're missing
      console.log('📊 Backend data analysis:');
      console.log('✅ Available:', {
        riskScores: !!pipelineResults.risk_scores,
        variantCount: !!pipelineResults.variant_count,
        processingTime: !!pipelineResults.processing_time_seconds,
        variantDetails: variantDetails.length > 0,
        tcgaData: !!pipelineResults.tcga_matches,
        caddStats: !!pipelineResults.cadd_stats
      });
      
      console.log('❌ Missing:', {
        mutationTypes: !structuredJson.summary?.mutation_types,
        mutationSignatures: !structuredJson.mutation_signatures,
        structuralVariants: !structuredJson.structural_variants,
        copyNumberVariants: !structuredJson.copy_number_variants,
        pathwayAnalysis: !structuredJson.pathway_analysis,
        survivalAnalysis: !structuredJson.survival_analysis
      });

      // Build the genomic data object using real data where available
      return {
        summary: {
          total_variants_found: pipelineResults.variant_count || structuredJson.summary?.total_variants_found || 0,
          variants_passed_qc: structuredJson.summary?.variants_passed_qc || Math.floor((pipelineResults.variant_count || 0) * 0.9),
          high_risk_findings_count: structuredJson.summary?.high_risk_findings || riskFindings.filter(r => r.risk_level === 'high').length,
          processing_time_seconds: pipelineResults.processing_time_seconds || 0,
          // Use real mutation types if available, otherwise estimate
          mutation_types: structuredJson.summary?.mutation_types || {
            snv: Math.floor((pipelineResults.variant_count || 0) * 0.5),
            indel: Math.floor((pipelineResults.variant_count || 0) * 0.3),
            cnv: Math.floor((pipelineResults.variant_count || 0) * 0.15),
            structural: Math.floor((pipelineResults.variant_count || 0) * 0.05)
          }
        },
        risk_findings: riskFindings,
        variant_table: variantTable,
        // Use real data if available, otherwise empty arrays (not mock data)
        structural_variants: structuredJson.structural_variants || [],
        copy_number_variants: structuredJson.copy_number_variants || [],
        // Use real mutation signatures if available, otherwise show message
        mutation_signatures: structuredJson.mutation_signatures || [
          {
            signature: "Analysis Pending",
            name: "Mutation signature analysis not yet available",
            contribution: 0,
            description: "This feature requires additional backend implementation",
            etiology: "Coming soon"
          }
        ]
      };
    } else {
      // No pipeline results - return empty state
      return {
        summary: {
          total_variants_found: 0,
          variants_passed_qc: 0,
          high_risk_findings_count: 0,
          processing_time_seconds: 0,
          mutation_types: {
            snv: 0,
            indel: 0,
            cnv: 0,
            structural: 0
          }
        },
        risk_findings: [],
        variant_table: [],
        structural_variants: [],
        copy_number_variants: [],
        mutation_signatures: []
      };
    }
  };

  // Memoize genomic data to prevent recalculation
  const genomicData = useMemo(() => getGenomicData(), [pipelineResults, fileName]);

  // Memoize gene associations calculation to prevent lag
  const geneAssociations = useMemo(() => {
    console.log('🔄 Calculating gene associations (memoized)');
    
    if (!pipelineResults?.pathway_burden_results) {
      return { associations: {}, genes: [], cancerTypes: ['breast', 'lung', 'colon', 'prostate', 'blood'] };
    }

    const pathwayBurdenResults = pipelineResults.pathway_burden_results;
    const allGenes = new Set<string>();
    
    // Collect genes from pathway burden analysis
    Object.values(pathwayBurdenResults).forEach((pathwayResult: PathwayResult) => {
      if (pathwayResult.contributing_genes) {
        pathwayResult.contributing_genes.forEach((gene: string) => allGenes.add(gene));
      }
      if (pathwayResult.damaging_genes) {
        pathwayResult.damaging_genes.forEach((gene: string) => allGenes.add(gene));
      }
    });

    // Limit to top 8 genes for performance
    const limitedGenes = Array.from(allGenes).slice(0, 8);
    const cancerTypes = ['breast', 'lung', 'colon', 'prostate', 'blood'];
    
    const associations: Record<string, Record<string, number>> = {};
    
    // Initialize with zeros
    limitedGenes.forEach(gene => {
      associations[gene] = {};
      cancerTypes.forEach(cancer => {
        associations[gene][cancer] = 0;
      });
    });

    // Build associations from pathway burden data
    const pathwayCancerMapping: Record<string, string[]> = {
      'dna_repair': ['breast', 'colon'],
      'tumor_suppressors': ['breast', 'lung', 'colon', 'prostate'],
      'oncogenes': ['lung', 'colon', 'breast'],
      'ras_mapk': ['lung', 'colon', 'prostate'],
      'cell_cycle': ['breast', 'lung', 'prostate'],
      'apoptosis': ['breast', 'lung', 'colon'],
      'chromatin_remodeling': ['blood', 'lung'],
      'mismatch_repair': ['colon'],
      'wnt_signaling': ['colon'],
      'pi3k_akt': ['breast', 'prostate']
    };

    Object.entries(pathwayBurdenResults).forEach(([pathwayName, pathwayResult]: [string, PathwayResult]) => {
      const associatedCancers = pathwayCancerMapping[pathwayName] || [];
      const burdenScore = pathwayResult.burden_score || 0;
      
      if (pathwayResult.damaging_genes) {
        pathwayResult.damaging_genes.forEach((gene: string) => {
          if (associations[gene]) {
            associatedCancers.forEach(cancer => {
              const geneVariantCount = pathwayResult.gene_damaging_counts?.[gene] || 1;
              const associationScore = burdenScore * 100 * Math.min(geneVariantCount, 2);
              associations[gene][cancer] = Math.max(
                associations[gene][cancer] || 0,
                Math.min(associationScore, 100)
              );
            });
          }
        });
      }
    });

    console.log('✅ Gene associations calculated:', Object.keys(associations).length, 'genes');
    return { associations, genes: limitedGenes, cancerTypes };
  }, [pipelineResults?.pathway_burden_results]);

  // Memoize paginated variants
  const paginatedVariants = useMemo(() => {
    const startIndex = (currentPage - 1) * variantsPerPage;
    const endIndex = startIndex + variantsPerPage;
    return genomicData.variant_table.slice(startIndex, endIndex);
  }, [genomicData.variant_table, currentPage, variantsPerPage]);

  const totalPages = Math.ceil(genomicData.variant_table.length / variantsPerPage);

  // Pagination component
  const PaginationControls = useCallback(() => (
    <div style={{
      display: 'flex',
      justifyContent: 'center',
      alignItems: 'center',
      gap: '1rem',
      padding: '1rem',
      borderTop: '1px solid #E5E7EB'
    }}>
      <button
        onClick={() => setCurrentPage(prev => Math.max(prev - 1, 1))}
        disabled={currentPage === 1}
        style={{
          padding: '0.5rem 1rem',
          background: currentPage === 1 ? '#F3F4F6' : '#2563EB',
          color: currentPage === 1 ? '#9CA3AF' : '#FFFFFF',
          border: 'none',
          borderRadius: '0.375rem',
          cursor: currentPage === 1 ? 'not-allowed' : 'pointer',
          fontSize: '0.875rem'
        }}
      >
        Previous
      </button>
      
      <span style={{ fontSize: '0.875rem', color: '#4B5563' }}>
        Page {currentPage} of {totalPages} ({genomicData.variant_table.length} total variants)
      </span>
      
      <button
        onClick={() => setCurrentPage(prev => Math.min(prev + 1, totalPages))}
        disabled={currentPage === totalPages}
        style={{
          padding: '0.5rem 1rem',
          background: currentPage === totalPages ? '#F3F4F6' : '#2563EB',
          color: currentPage === totalPages ? '#9CA3AF' : '#FFFFFF',
          border: 'none',
          borderRadius: '0.375rem',
          cursor: currentPage === totalPages ? 'not-allowed' : 'pointer',
          fontSize: '0.875rem'
        }}
      >
        Next
      </button>
    </div>
  ), [currentPage, totalPages, genomicData.variant_table.length]);

  const renderTabContent = () => {
    switch (activeTab) {
      case 'analysis':
        return (
          <div style={{ padding: '1rem' }}>
            <div style={{ 
              display: 'flex',
              justifyContent: 'space-between',
              alignItems: 'center',
              marginBottom: '1.5rem'
            }}>
            <h2 style={{ 
              color: '#111827',
              fontSize: '1.25rem',
              fontWeight: '600',
              margin: 0
            }}>
              Genomic Analysis Overview
            </h2>
              <SubtabDownloadButton 
                subtabContent={subtabContents.analysis}
                isPDFGenerating={isPDFGenerating}
                setIsPDFGenerating={setIsPDFGenerating}
              />
            </div>
            
            {pipelineResults && (
              <div style={{
                padding: '0.75rem',
                marginBottom: '1rem',
                background: '#EFF6FF',
                borderRadius: '0.5rem',
                border: '1px solid #BFDBFE'
              }}>
                <p style={{ color: '#1E40AF', fontSize: '0.8rem', fontWeight: '500' }}>
                  Analysis completed.
                </p>
              </div>
            )}
            
            {/* Risk Summary Cards */}
            <div id="cancer-risk-assessment" style={{ 
              background: '#FFFFFF',
              padding: '1rem',
              borderRadius: '0.5rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '1rem'
            }}>
                              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <div style={{ 
                  display: 'flex',
                  alignItems: 'center',
                  gap: '0.5rem'
                }}>
                  <h3 style={{ 
                    color: '#111827',
                    fontSize: '1.125rem',
                    fontWeight: '600',
                    margin: 0
                  }}>
                    Cancer Risk Assessment
                  </h3>
                  <div 
                    style={{ position: 'relative', display: 'inline-flex' }}
                    onMouseEnter={() => setHoveredTooltip('risk-assessment')}
                    onMouseLeave={() => setHoveredTooltip(null)}
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
                    <SmartTooltip content="Shows calculated cancer risk percentages based on genetic variants and pathway analysis. Each cancer type is assessed individually using machine learning models trained on population data and clinical outcomes." isVisible={hoveredTooltip === 'risk-assessment'} triggerRef={null} />
                  </div>
                </div>

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
              padding: '1rem',
              borderRadius: '0.5rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '1rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <div style={{ 
                  display: 'flex',
                  alignItems: 'center',
                  gap: '0.5rem'
                }}>
                  <h3 style={{ 
                    color: '#111827',
                    fontSize: '1.125rem',
                    fontWeight: '600',
                    margin: 0
                  }}>
                    Gene Significance Analysis
                  </h3>
                  <div 
                    style={{ position: 'relative', display: 'inline-flex' }}
                    onMouseEnter={() => setHoveredTooltip('gene-significance')}
                    onMouseLeave={() => setHoveredTooltip(null)}
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
                    <SmartTooltip content="Manhattan plot showing statistical significance of genetic variants. Higher points indicate stronger associations with cancer risk. Red dots represent pathogenic variants, blue dots show variants under investigation." isVisible={hoveredTooltip === 'gene-significance'} triggerRef={null} />
                  </div>
                </div>

              </div>
              <svg width="100%" height="200" viewBox="0 0 800 200">
                {/* Background */}
                <rect width="800" height="200" fill="#F9FAFB"/>
                
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
                        {(() => {
                          const cancerType = (variant.tcga_best_match as { cancer_type?: string })?.cancer_type;
                          if (cancerType === 'not_cancer_associated') {
                            return 'No cancer association';
                          }
                          return cancerType || 'Unknown';
                        })()}
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
                <strong>Interpretation:</strong> This plot shows gene variants found in your analysis. 
                Red dots represent high-confidence pathogenic variants, blue dots represent variants under investigation. 
                "No cancer association" indicates genes not found in cancer databases (this is normal for most genes).
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
                <div style={{ 
                  display: 'flex',
                  alignItems: 'center',
                  gap: '0.5rem'
                }}>
                  <h3 style={{ 
                    color: '#111827',
                    fontSize: '1.125rem',
                    fontWeight: '600',
                    margin: 0
                  }}>
                    Mutation Type Distribution
                  </h3>
                  <div 
                    style={{ position: 'relative', display: 'inline-flex' }}
                    onMouseEnter={() => setHoveredTooltip('mutation-types')}
                    onMouseLeave={() => setHoveredTooltip(null)}
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
                    <SmartTooltip content="Breakdown of mutation types found in the genetic analysis. SNVs (single nucleotide variants) are point mutations, INDELs are insertions/deletions, CNVs are copy number variations, and structural variants are large chromosomal rearrangements." isVisible={hoveredTooltip === 'mutation-types'} triggerRef={null} />
                  </div>
                </div>

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
                <div style={{ 
                  display: 'flex',
                  alignItems: 'center',
                  gap: '0.5rem'
                }}>
                  <h3 style={{ 
                    color: '#111827',
                    fontSize: '1.125rem',
                    fontWeight: '600',
                    margin: 0
                  }}>
                    Mutational Signature Analysis
                  </h3>
                  <div 
                    style={{ position: 'relative', display: 'inline-flex' }}
                    onMouseEnter={() => setHoveredTooltip('mutational-signatures')}
                    onMouseLeave={() => setHoveredTooltip(null)}
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
                    <SmartTooltip content="Mutational signatures reveal the underlying biological processes that caused DNA damage. Each signature represents a specific pattern of mutations caused by factors like aging, smoking, UV exposure, DNA repair defects, or chemotherapy. The percentage shows each signature's contribution to the overall mutational profile." isVisible={hoveredTooltip === 'mutational-signatures'} triggerRef={null} />
                  </div>
                </div>
              </div>
              
              <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                {genomicData.mutation_signatures.length === 0 ? (
                  <div style={{
                    textAlign: 'center',
                    padding: '3rem',
                    color: '#6B7280'
                  }}>
                    <h4 style={{ 
                      color: '#111827',
                      fontSize: '1.125rem',
                      fontWeight: '600',
                      marginBottom: '0.5rem'
                    }}>
                      No Mutational Signatures Found
                    </h4>
                    <p style={{ 
                      color: '#6B7280',
                      fontSize: '0.875rem',
                      maxWidth: '400px',
                      margin: '0 auto'
                    }}>
                      The analysis did not detect any significant mutational signatures in the provided genomic data. This could indicate insufficient data or variants not matching known signature patterns.
                    </p>
                  </div>
                ) : (
                  genomicData.mutation_signatures.map((signature, index) => (
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
                  ))
                )}
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
                <div style={{ 
                  display: 'flex',
                  alignItems: 'center',
                  gap: '0.5rem'
                }}>
                  <h3 style={{ 
                    color: '#111827',
                    fontSize: '1.125rem',
                    fontWeight: '600',
                    margin: 0
                  }}>
                    Analysis Quality Metrics
                  </h3>
                  <div 
                    style={{ position: 'relative', display: 'inline-flex' }}
                    onMouseEnter={() => setHoveredTooltip('quality-metrics')}
                    onMouseLeave={() => setHoveredTooltip(null)}
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
                    <SmartTooltip content="Quality control metrics for the genomic analysis. Shows total variants found, those passing quality filters, high-risk findings requiring clinical attention, and processing time for the analysis pipeline." isVisible={hoveredTooltip === 'quality-metrics'} triggerRef={null} />
                  </div>
                </div>
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
          <div style={{ 
            padding: '2rem',
            maxWidth: '100%',
            overflow: 'hidden'
          }}>
            <div style={{ 
              display: 'flex',
              justifyContent: 'space-between',
              alignItems: 'center',
              marginBottom: '1.5rem'
            }}>
              <div style={{
                display: 'flex',
                alignItems: 'center',
                gap: '0.5rem'
            }}>
              <h2 style={{ 
                color: '#111827',
                fontSize: '1.5rem',
                fontWeight: '600',
                margin: 0
              }}>
                Variant Heatmap Analysis
              </h2>
              <div 
                style={{ position: 'relative', display: 'inline-flex' }}
                onMouseEnter={() => setHoveredTooltip('variant-heatmap')}
                onMouseLeave={() => setHoveredTooltip(null)}
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
                <SmartTooltip content="Interactive heatmap showing gene-cancer associations based on pathway burden analysis. Darker colors indicate stronger associations between specific genes and cancer types based on genetic variant patterns." isVisible={hoveredTooltip === 'variant-heatmap'} triggerRef={null} />
              </div>
              </div>
              <SubtabDownloadButton 
                subtabContent={subtabContents.variants}
                isPDFGenerating={isPDFGenerating}
                setIsPDFGenerating={setIsPDFGenerating}
              />
            </div>
            
            {/* Performance-optimized Gene-Cancer Type Heatmap */}
            <div id="gene-cancer-matrix" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB',
              marginBottom: '2rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <div style={{ 
                  display: 'flex',
                  alignItems: 'center',
                  gap: '0.5rem'
                }}>
                  <h3 style={{ 
                    color: '#111827',
                    fontSize: '1.125rem',
                    fontWeight: '600',
                    margin: 0
                  }}>
                    Gene-Cancer Type Association Matrix
                  </h3>
                  <div 
                    style={{ position: 'relative', display: 'inline-flex' }}
                    onMouseEnter={() => setHoveredTooltip('gene-cancer-matrix')}
                    onMouseLeave={() => setHoveredTooltip(null)}
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
                    <SmartTooltip content="Color-coded heatmap showing the strength of association between specific genes and cancer types. 
                        Red indicates high association (≥80%), orange medium (60-80%), blue low (30-60%), light blue very low (1-30%), 
                        and gray no association. Based on pathway burden analysis and research literature. 
                        Use this to understand which genes are most relevant for specific cancer screening protocols." isVisible={hoveredTooltip === 'gene-cancer-matrix'} triggerRef={null} />
                  </div>
                </div>
              </div>
              
              {/* Optimized heatmap rendering */}
              {(() => {
                if (!geneAssociations.genes || geneAssociations.genes.length === 0) {
                  return (
                    <div style={{ 
                      padding: '2rem', 
                      textAlign: 'center',
                      color: '#6B7280',
                      fontSize: '0.875rem'
                    }}>
                      No variant data available for gene-cancer association matrix
                    </div>
                  );
                }

                const { associations = {}, genes = [], cancerTypes = ['breast', 'lung', 'colon', 'prostate', 'blood'] } = geneAssociations;
                
                // Calculate summary statistics
                const totalAssociations = genes.reduce((sum, gene) => {
                  return sum + cancerTypes.reduce((geneSum, cancer) => {
                    return geneSum + (associations[gene]?.[cancer] > 0 ? 1 : 0);
                  }, 0);
                }, 0);

                const highRiskAssociations = genes.reduce((sum, gene) => {
                  return sum + cancerTypes.reduce((geneSum, cancer) => {
                    return geneSum + (associations[gene]?.[cancer] > 80 ? 1 : 0);
                  }, 0);
                }, 0);

                const pathwayBurdenSummary = pipelineResults?.pathway_burden_summary;
                
                return (
                  <div>
                    {/* Summary Statistics */}
                    <div style={{ 
                      marginBottom: '1.5rem',
                      padding: '1rem',
                      background: '#F9FAFB',
                      borderRadius: '0.5rem',
                      border: '1px solid #E5E7EB'
                    }}>
                      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '1rem' }}>
                        <div style={{ textAlign: 'center' }}>
                          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#3B82F6' }}>
                            {genes.length}
                          </div>
                          <div style={{ fontSize: '0.875rem', color: '#6B7280' }}>
                            Top Genes Analyzed
                          </div>
                        </div>
                        <div style={{ textAlign: 'center' }}>
                          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#10B981' }}>
                            {totalAssociations}
                          </div>
                          <div style={{ fontSize: '0.875rem', color: '#6B7280' }}>
                            Total Associations
                          </div>
                        </div>
                        <div style={{ textAlign: 'center' }}>
                          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#EF4444' }}>
                            {highRiskAssociations}
                          </div>
                          <div style={{ fontSize: '0.875rem', color: '#6B7280' }}>
                            High Risk (&gt;80%)
                          </div>
                        </div>
                        {pathwayBurdenSummary && (
                          <div style={{ textAlign: 'center' }}>
                            <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#8B5CF6' }}>
                              {pathwayBurdenSummary.high_burden_pathways?.length || 0}
                            </div>
                            <div style={{ fontSize: '0.875rem', color: '#6B7280' }}>
                              High-Burden Pathways
                            </div>
                          </div>
                        )}
                      </div>
                      
                      {genes.length >= 8 && (
                        <div style={{ 
                          marginTop: '0.75rem', 
                          padding: '0.5rem',
                          background: '#EFF6FF',
                          borderRadius: '0.375rem',
                          border: '1px solid #BFDBFE'
                        }}>
                          <div style={{ fontSize: '0.75rem', color: '#1E40AF', fontWeight: '500' }}>
                            ℹ️ Showing top {genes.length} genes for optimal performance. {genomicData.variant_table.length > genes.length ? `${genomicData.variant_table.length - genes.length} additional genes available in variant table below.` : ''}
                          </div>
                        </div>
                      )}
                    </div>
                    
                    {/* Optimized SVG Heatmap */}
                    <div style={{
                      width: '100%',
                      background: '#F9FAFB',
                      borderRadius: '0.5rem',
                      padding: '1rem'
                    }}>
                      <svg 
                        width="100%" 
                        height="300" 
                        viewBox="0 0 700 250" 
                        style={{ 
                          width: '100%',
                          height: 'auto'
                        }}
                        preserveAspectRatio="xMidYMid meet"
                      >
                        {/* Background */}
                        <rect width="700" height="250" fill="#F9FAFB"/>
                        
                        {/* Optimized heatmap grid - only render visible cells */}
                        {genes.map((gene, geneIndex) => {
                          return cancerTypes.map((cancer, cancerIndex) => {
                            const cellWidth = 80;
                            const cellHeight = 35;
                            const startX = 120;
                            const startY = 50;
                            
                            const x = startX + cancerIndex * (cellWidth + 5);
                            const y = startY + geneIndex * (cellHeight + 5);
                            
                            const intensity = (associations[gene]?.[cancer] || 0) / 100;
                            
                            const color = intensity > 0.8 ? '#EF4444' : 
                                         intensity > 0.6 ? '#F59E0B' : 
                                         intensity > 0.3 ? '#3B82F6' : 
                                         intensity > 0 ? '#93C5FD' : '#E5E7EB';
                            
                            return (
                              <g key={`${gene}-${cancer}`}>
                                <rect 
                                  x={x} 
                                  y={y} 
                                  width={cellWidth} 
                                  height={cellHeight} 
                                  fill={color}
                                  stroke="#FFFFFF"
                                  strokeWidth="2"
                                  opacity="0.8"
                                  rx="4"
                                />
                                {intensity > 0 && (
                                  <text 
                                    x={x + cellWidth/2} 
                                    y={y + cellHeight/2} 
                                    fill="#FFFFFF" 
                                    fontSize="12" 
                                    textAnchor="middle"
                                    fontWeight="600"
                                    dominantBaseline="middle"
                                  >
                                    {(intensity * 100).toFixed(0)}%
                                  </text>
                                )}
                              </g>
                            );
                          });
                        })}
                      
                        {/* Gene labels */}
                        {genes.map((gene, index) => (
                          <text 
                            key={gene}
                            x="110" 
                            y={50 + 17.5 + index * 40} 
                            fill="#111827" 
                            fontSize="14" 
                            textAnchor="end"
                            fontWeight="600"
                            dominantBaseline="middle"
                          >
                            {gene}
                          </text>
                        ))}
                        
                        {/* Cancer type labels */}
                        {cancerTypes.map((cancer, index) => (
                          <text 
                            key={cancer}
                            x={120 + 40 + index * 85} 
                            y="35" 
                            fill="#111827" 
                            fontSize="14" 
                            textAnchor="middle"
                            fontWeight="600"
                            dominantBaseline="middle"
                          >
                            {cancer.charAt(0).toUpperCase() + cancer.slice(1)}
                          </text>
                        ))}
                      
                        {/* Legend */}
                        <g transform="translate(550, 60)">
                          <text x="0" y="0" fill="#111827" fontSize="12" fontWeight="600">Risk Level</text>
                          <rect x="0" y="12" width="12" height="12" fill="#EF4444" rx="2"/>
                          <text x="16" y="23" fill="#4B5563" fontSize="10" dominantBaseline="middle">High (≥80%)</text>
                          <rect x="0" y="30" width="12" height="12" fill="#F59E0B" rx="2"/>
                          <text x="16" y="41" fill="#4B5563" fontSize="10" dominantBaseline="middle">Medium (60-80%)</text>
                          <rect x="0" y="48" width="12" height="12" fill="#3B82F6" rx="2"/>
                          <text x="16" y="59" fill="#4B5563" fontSize="10" dominantBaseline="middle">Low (30-60%)</text>
                          <rect x="0" y="66" width="12" height="12" fill="#93C5FD" rx="2"/>
                          <text x="16" y="77" fill="#4B5563" fontSize="10" dominantBaseline="middle">Very Low (1-30%)</text>
                          <rect x="0" y="84" width="12" height="12" fill="#E5E7EB" rx="2"/>
                          <text x="16" y="95" fill="#4B5563" fontSize="10" dominantBaseline="middle">No Association</text>
                        </g>
                      </svg>
                    </div>
                  </div>
                );
              })()}
            </div>
            
            {/* Performance-optimized Variant Table with Pagination */}
            <div id="detected-variants" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <div style={{ 
                  display: 'flex',
                  alignItems: 'center',
                  gap: '0.5rem'
                }}>
                  <h3 style={{ 
                    color: '#111827',
                    fontSize: '1.125rem',
                    fontWeight: '600',
                    margin: 0
                  }}>
                    Detected Variants
                  </h3>
                  <div 
                    style={{ position: 'relative', display: 'inline-flex' }}
                    onMouseEnter={() => setHoveredTooltip('detected-variants')}
                    onMouseLeave={() => setHoveredTooltip(null)}
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
                    <SmartTooltip content="Comprehensive table of all genetic variants found in the analysis. Includes gene names, genomic positions, mutation types, protein changes, quality scores, research significance, and functional impact assessments. Paginated for optimal performance." isVisible={hoveredTooltip === 'detected-variants'} triggerRef={null} />
                  </div>
                </div>
              </div>
              
              {/* Table stats */}
              <div style={{ 
                display: 'flex', 
                justifyContent: 'space-between', 
                alignItems: 'center',
                padding: '0.5rem 0',
                fontSize: '0.875rem',
                color: '#6B7280',
                borderBottom: '1px solid #E5E7EB',
                marginBottom: '0.5rem'
              }}>
                <span>
                  Showing page {currentPage} of {totalPages} ({genomicData.variant_table.length} total variants)
                </span>
                <span style={{ 
                  display: 'flex', 
                  alignItems: 'center', 
                  gap: '0.25rem',
                  fontStyle: 'italic'
                }}>
                  {variantsPerPage} variants per page
                </span>
              </div>
              
              <div style={{ 
                overflowX: 'visible', 
                border: '1px solid #E5E7EB',
                borderRadius: '0.5rem'
              }}>
                <table style={{ 
                  width: '100%', 
                  borderCollapse: 'collapse', 
                  tableLayout: 'fixed',
                  minWidth: '800px'
                }}>
                  <thead style={{ 
                    position: 'sticky', 
                    top: 0, 
                    zIndex: 1,
                    background: '#F9FAFB' // Ensure header background stays visible when scrolling
                  }}>
                    <tr style={{ background: '#F9FAFB', borderBottom: '2px solid #E5E7EB' }}>
                      <th style={{ 
                        padding: '1rem', 
                        textAlign: 'left', 
                        borderBottom: '1px solid #E5E7EB', 
                        fontWeight: '600',
                        color: '#111827',
                        fontSize: '0.875rem',
                        width: '10%'
                      }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                          Gene
                          <div 
                            style={{ position: 'relative', display: 'inline-flex' }}
                            onMouseEnter={() => setHoveredTooltip('header-gene')}
                            onMouseLeave={() => setHoveredTooltip(null)}
                          >
                            <InformationCircleIcon 
                              style={{ width: '14px', height: '14px' }}
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
                            <SmartTooltip content="Gene symbol where the variant is located. These are typically cancer-associated genes like BRCA1, TP53, KRAS, etc. that are important for cancer risk assessment." isVisible={hoveredTooltip === 'header-gene'} triggerRef={null} />
                          </div>
                        </div>
                      </th>
                      <th style={{ 
                        padding: '1rem', 
                        textAlign: 'left', 
                        borderBottom: '1px solid #E5E7EB', 
                        fontWeight: '600',
                        color: '#111827',
                        fontSize: '0.875rem',
                        width: '15%'
                      }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                          Variant
                          <div 
                            style={{ position: 'relative', display: 'inline-flex' }}
                            onMouseEnter={() => setHoveredTooltip('header-variant')}
                            onMouseLeave={() => setHoveredTooltip(null)}
                          >
                            <InformationCircleIcon 
                              style={{ width: '14px', height: '14px' }}
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
                            <div style={{
                              position: 'absolute',
                              left: '0',
                              top: 'calc(100% + 0.5rem)',
                              width: '18rem',
                              padding: '0.75rem',
                              background: '#1F2937',
                              color: '#FFFFFF',
                              fontSize: '0.75rem',
                              borderRadius: '0.5rem',
                              boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
                              zIndex: 1000,
                              visibility: hoveredTooltip === 'header-variant' ? 'visible' as const : 'hidden' as const,
                              opacity: hoveredTooltip === 'header-variant' ? 1 : 0,
                              transition: 'opacity 0.2s, visibility 0.2s',
                              pointerEvents: 'none',
                              border: '1px solid #374151'
                            }}>
                              <p style={{ color: '#D1D5DB', marginBottom: '0' }}>Genomic coordinates and nucleotide change. Format: position:reference→alternate (e.g., 17:41223094:A→G). This uniquely identifies the DNA change.</p>
                              <div style={{
                                position: 'absolute',
                                left: '0.75rem',
                                bottom: '100%',
                                width: '0',
                                height: '0',
                                borderLeft: '0.5rem solid transparent',
                                borderRight: '0.5rem solid transparent',
                                borderBottom: '0.5rem solid #1F2937'
                              }}></div>
                            </div>
                          </div>
                        </div>
                      </th>
                      <th style={{ 
                        padding: '1rem', 
                        textAlign: 'left', 
                        borderBottom: '1px solid #E5E7EB', 
                        fontWeight: '600',
                        color: '#111827',
                        fontSize: '0.875rem',
                        width: '12%'
                      }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                          Type
                          <div 
                            style={{ position: 'relative', display: 'inline-flex' }}
                            onMouseEnter={() => setHoveredTooltip('header-type')}
                            onMouseLeave={() => setHoveredTooltip(null)}
                          >
                            <InformationCircleIcon 
                              style={{ width: '14px', height: '14px' }}
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
                            <div style={{
                              position: 'absolute',
                              left: '50%',
                              top: 'calc(100% + 0.5rem)',
                              transform: 'translateX(-50%)',
                              width: '18rem',
                              padding: '0.75rem',
                              background: '#1F2937',
                              color: '#FFFFFF',
                              fontSize: '0.75rem',
                              borderRadius: '0.5rem',
                              boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
                              zIndex: 1000,
                              visibility: hoveredTooltip === 'header-type' ? 'visible' as const : 'hidden' as const,
                              opacity: hoveredTooltip === 'header-type' ? 1 : 0,
                              transition: 'opacity 0.2s, visibility 0.2s',
                              pointerEvents: 'none',
                              border: '1px solid #374151'
                            }}>
                              <p style={{ color: '#D1D5DB', marginBottom: '0' }}>Mutation type: SNV (Single Nucleotide Variant - blue), INDEL (Insertion/Deletion - red), or CNV (Copy Number Variant - yellow). Different types have different research implications.</p>
                              <div style={{
                                position: 'absolute',
                                bottom: '100%',
                                left: '50%',
                                transform: 'translateX(-50%)',
                                width: 0,
                                height: 0,
                                borderLeft: '6px solid transparent',
                                borderRight: '6px solid transparent',
                                borderBottom: '6px solid #1F2937'
                              }} />
                            </div>
                          </div>
                        </div>
                      </th>
                      <th style={{ 
                        padding: '0.75rem', 
                        textAlign: 'left', 
                        borderBottom: '1px solid #E5E7EB', 
                        fontWeight: '600',
                        color: '#111827',
                        fontSize: '0.875rem',
                        width: '14%'
                      }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                          Transformation
                          <div 
                            style={{ position: 'relative', display: 'inline-flex' }}
                            onMouseEnter={() => setHoveredTooltip('header-transformation')}
                            onMouseLeave={() => setHoveredTooltip(null)}
                          >
                            <InformationCircleIcon 
                              style={{ width: '14px', height: '14px' }}
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
                            <div style={{
                              position: 'absolute',
                              left: '50%',
                              top: 'calc(100% + 0.5rem)',
                              transform: 'translateX(-50%)',
                              width: '18rem',
                              padding: '0.75rem',
                              background: '#1F2937',
                              color: '#FFFFFF',
                              fontSize: '0.75rem',
                              borderRadius: '0.5rem',
                              boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
                              zIndex: 1000,
                              visibility: hoveredTooltip === 'header-transformation' ? 'visible' as const : 'hidden' as const,
                              opacity: hoveredTooltip === 'header-transformation' ? 1 : 0,
                              transition: 'opacity 0.2s, visibility 0.2s',
                              pointerEvents: 'none',
                              border: '1px solid #374151'
                            }}>
                              <p style={{ color: '#D1D5DB', marginBottom: '0' }}>Protein-level change caused by the DNA variant. Shows original codon → mutated codon and resulting amino acid change (e.g., Arg→His). Critical for understanding functional impact.</p>
                              <div style={{
                                position: 'absolute',
                                bottom: '100%',
                                left: '50%',
                                transform: 'translateX(-50%)',
                                width: 0,
                                height: 0,
                                borderLeft: '6px solid transparent',
                                borderRight: '6px solid transparent',
                                borderBottom: '6px solid #1F2937'
                              }} />
                            </div>
                          </div>
                        </div>
                      </th>
                      <th style={{ 
                        padding: '0.75rem', 
                        textAlign: 'left', 
                        borderBottom: '1px solid #E5E7EB', 
                        fontWeight: '600',
                        color: '#111827',
                        fontSize: '0.875rem',
                        width: '10%'
                      }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                          Quality
                          <div 
                            style={{ position: 'relative', display: 'inline-flex' }}
                            onMouseEnter={() => setHoveredTooltip('header-quality')}
                            onMouseLeave={() => setHoveredTooltip(null)}
                          >
                            <InformationCircleIcon 
                              style={{ width: '14px', height: '14px' }}
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
                            <div style={{
                              position: 'absolute',
                              left: '50%',
                              top: 'calc(100% + 0.5rem)',
                              transform: 'translateX(-50%)',
                              width: '18rem',
                              padding: '0.75rem',
                              background: '#1F2937',
                              color: '#FFFFFF',
                              fontSize: '0.75rem',
                              borderRadius: '0.5rem',
                              boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
                              zIndex: 1000,
                              visibility: hoveredTooltip === 'header-quality' ? 'visible' as const : 'hidden' as const,
                              opacity: hoveredTooltip === 'header-quality' ? 1 : 0,
                              transition: 'opacity 0.2s, visibility 0.2s',
                              pointerEvents: 'none',
                              border: '1px solid #374151'
                            }}>
                              <p style={{ color: '#D1D5DB', marginBottom: '0' }}>Variant call quality score (0-100). Green (&gt;90) = high confidence, Yellow (70-90) = moderate confidence, Red (&lt;70) = low confidence. Score of 100 indicates variants that passed all quality filters.</p>
                              <div style={{
                                position: 'absolute',
                                bottom: '100%',
                                left: '50%',
                                transform: 'translateX(-50%)',
                                width: 0,
                                height: 0,
                                borderLeft: '6px solid transparent',
                                borderRight: '6px solid transparent',
                                borderBottom: '6px solid #1F2937'
                              }} />
                            </div>
                          </div>
                        </div>
                      </th>
                      <th style={{ 
                        padding: '0.75rem', 
                        textAlign: 'left', 
                        borderBottom: '1px solid #E5E7EB', 
                        fontWeight: '600',
                        color: '#111827',
                        fontSize: '0.875rem',
                        width: '15%'
                      }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                          Significance
                          <div 
                            style={{ position: 'relative', display: 'inline-flex' }}
                            onMouseEnter={() => setHoveredTooltip('header-significance')}
                            onMouseLeave={() => setHoveredTooltip(null)}
                          >
                            <InformationCircleIcon 
                              style={{ width: '14px', height: '14px' }}
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
                            <div style={{
                              position: 'absolute',
                              left: '50%',
                              top: 'calc(100% + 0.5rem)',
                              transform: 'translateX(-50%)',
                              width: '18rem',
                              padding: '0.75rem',
                              background: '#1F2937',
                              color: '#FFFFFF',
                              fontSize: '0.75rem',
                              borderRadius: '0.5rem',
                              boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
                              zIndex: 1000,
                              visibility: hoveredTooltip === 'header-significance' ? 'visible' as const : 'hidden' as const,
                              opacity: hoveredTooltip === 'header-significance' ? 1 : 0,
                              transition: 'opacity 0.2s, visibility 0.2s',
                              pointerEvents: 'none',
                              border: '1px solid #374151'
                            }}>
                              <p style={{ color: '#D1D5DB', marginBottom: '0' }}>Research significance based on ClinVar and population data. Pathogenic (red) = disease-causing, Likely pathogenic (orange) = probably disease-causing, Uncertain significance (green) = unknown research impact.</p>
                              <div style={{
                                position: 'absolute',
                                bottom: '100%',
                                left: '50%',
                                transform: 'translateX(-50%)',
                                width: 0,
                                height: 0,
                                borderLeft: '6px solid transparent',
                                borderRight: '6px solid transparent',
                                borderBottom: '6px solid #1F2937'
                              }} />
                            </div>
                          </div>
                        </div>
                      </th>
                      <th style={{ 
                        padding: '0.75rem', 
                        textAlign: 'left', 
                        borderBottom: '1px solid #E5E7EB', 
                        fontWeight: '600',
                        color: '#111827',
                        fontSize: '0.875rem',
                        width: '24%'
                      }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                          Impact
                          <div 
                            style={{ position: 'relative', display: 'inline-flex' }}
                            onMouseEnter={() => setHoveredTooltip('header-impact')}
                            onMouseLeave={() => setHoveredTooltip(null)}
                          >
                            <InformationCircleIcon 
                              style={{ width: '14px', height: '14px' }}
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
                            <div style={{
                              position: 'absolute',
                              left: '50%',
                              top: 'calc(100% + 0.5rem)',
                              transform: 'translateX(-50%)',
                              width: '18rem',
                              padding: '0.75rem',
                              background: '#1F2937',
                              color: '#FFFFFF',
                              fontSize: '0.75rem',
                              borderRadius: '0.5rem',
                              boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
                              zIndex: 1000,
                              visibility: hoveredTooltip === 'header-impact' ? 'visible' as const : 'hidden' as const,
                              opacity: hoveredTooltip === 'header-impact' ? 1 : 0,
                              transition: 'opacity 0.2s, visibility 0.2s',
                              pointerEvents: 'none',
                              border: '1px solid #374151'
                            }}>
                              <p style={{ color: '#D1D5DB', marginBottom: '0' }}>Predicted functional impact of the variant on protein function. Includes details about amino acid property changes, charge effects, and potential disruption to protein structure or function.</p>
                              <div style={{
                                position: 'absolute',
                                bottom: '100%',
                                left: '50%',
                                transform: 'translateX(-50%)',
                                width: 0,
                                height: 0,
                                borderLeft: '6px solid transparent',
                                borderRight: '6px solid transparent',
                                borderBottom: '6px solid #1F2937'
                              }} />
                            </div>
                          </div>
                        </div>
                      </th>
                    </tr>
                  </thead>
                  <tbody>
                    {/* Render only paginated variants for performance */}
                    {paginatedVariants.map((variant, index) => (
                      <tr key={`${currentPage}-${index}`} style={{ 
                        borderBottom: '1px solid #F3F4F6',
                        minHeight: '60px' // Ensure minimum row height
                      }}>
                        <td style={{ 
                          padding: '1rem', 
                          fontWeight: '600', 
                          color: '#111827',
                          wordWrap: 'break-word',
                          overflowWrap: 'break-word',
                          width: '10%',
                          lineHeight: '1.5',
                          verticalAlign: 'middle'
                        }}>
                          {variant.gene}
                        </td>
                        <td style={{ 
                          padding: '1rem', 
                          fontFamily: 'monospace', 
                          fontSize: '0.875rem', 
                          color: '#4B5563',
                          wordWrap: 'break-word',
                          overflowWrap: 'break-word',
                          width: '15%',
                          lineHeight: '1.5',
                          verticalAlign: 'middle'
                        }}>
                          {variant.variant_id.split(':').slice(1).join(':')}
                        </td>
                        <td style={{ 
                          padding: '1rem',
                          width: '12%',
                          verticalAlign: 'middle'
                        }}>
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
                        <td style={{ 
                          padding: '1rem',
                          width: '14%',
                          verticalAlign: 'middle'
                        }}>
                          <div style={{ 
                            fontSize: '0.875rem', 
                            color: '#4B5563',
                            wordWrap: 'break-word',
                            overflowWrap: 'break-word'
                          }}>
                            <div style={{ 
                              fontFamily: 'monospace', 
                              fontSize: '0.75rem', 
                              marginBottom: '0.25rem',
                              wordWrap: 'break-word',
                              overflowWrap: 'break-word'
                            }}>
                              {variant.transformation.original} → {variant.transformation.mutated}
                            </div>
                            <div style={{ 
                              fontSize: '0.75rem', 
                              color: '#6B7280',
                              wordWrap: 'break-word',
                              overflowWrap: 'break-word'
                            }}>
                              {variant.transformation.amino_acid_change}
                            </div>
                          </div>
                        </td>
                        <td style={{ 
                          padding: '1rem',
                          width: '10%',
                          verticalAlign: 'middle'
                        }}>
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
                        <td style={{ 
                          padding: '1rem',
                          width: '15%',
                          verticalAlign: 'middle'
                        }}>
                          <div style={{
                            background: variant.clinical_significance === 'pathogenic' ? '#FEF2F2' : '#F0FDF4',
                            color: variant.clinical_significance === 'pathogenic' ? '#EF4444' : '#22C55E',
                            padding: '0.25rem 0.5rem',
                            borderRadius: '0.375rem',
                            fontSize: '0.75rem',
                            fontWeight: '600',
                            textAlign: 'center',
                            wordWrap: 'break-word',
                            overflowWrap: 'break-word'
                          }}>
                            {variant.clinical_significance}
                          </div>
                        </td>
                        <td style={{ 
                          padding: '1rem', 
                          fontSize: '0.875rem', 
                          color: '#4B5563',
                          width: '24%',
                          lineHeight: '1.5',
                          verticalAlign: 'middle'
                        }}>
                          <div style={{ 
                            marginBottom: '0.25rem',
                            wordWrap: 'break-word',
                            overflowWrap: 'break-word'
                          }}>
                            <strong>{variant.functional_impact}</strong>
                          </div>
                          <div style={{ 
                            fontSize: '0.75rem', 
                            color: '#6B7280',
                            wordWrap: 'break-word',
                            overflowWrap: 'break-word'
                          }}>
                            {(variant as ExtendedVariant).transformation?.effect || variant.functional_impact || 'Impact not determined'}
                          </div>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
                </div>              
              {/* Pagination Controls */}
              {totalPages > 1 && <PaginationControls />}
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
                {genomicData.structural_variants.length === 0 ? (
                  <div style={{
                    textAlign: 'center',
                    padding: '2rem',
                    color: '#6B7280'
                  }}>
                    <h4 style={{ 
                      color: '#111827',
                      fontSize: '1rem',
                      fontWeight: '600',
                      marginBottom: '0.5rem'
                    }}>
                      No Structural Variants Detected
                    </h4>
                    <p style={{ 
                      color: '#6B7280',
                      fontSize: '0.875rem',
                      maxWidth: '350px',
                      margin: '0 auto'
                    }}>
                      No large-scale genomic rearrangements were found in the analysis.
                    </p>
                  </div>
                ) : (
                  genomicData.structural_variants.map((variant, index) => (
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
                          {(variant as ExtendedVariant).transformation?.effect || variant.functional_impact || 'Impact not determined'}
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
                  ))
                )}
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
                {genomicData.copy_number_variants.length === 0 ? (
                  <div style={{
                    textAlign: 'center',
                    padding: '2rem',
                    color: '#6B7280'
                  }}>
                    <h4 style={{ 
                      color: '#111827',
                      fontSize: '1rem',
                      fontWeight: '600',
                      marginBottom: '0.5rem'
                    }}>
                      No Copy Number Variants Detected
                    </h4>
                    <p style={{ 
                      color: '#6B7280',
                      fontSize: '0.875rem',
                      maxWidth: '350px',
                      margin: '0 auto'
                    }}>
                      No significant gene amplifications or deletions were detected in the analysis.
                    </p>
                  </div>
                ) : (
                  genomicData.copy_number_variants.map((variant, index) => (
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
                  ))
                )}
              </div>
            </div>
          </div>
        );
      
      case 'pathways':
        return (
          <div style={{ padding: '2rem' }}>
            <div style={{ 
              display: 'flex',
              justifyContent: 'space-between',
              alignItems: 'center',
              marginBottom: '1.5rem'
            }}>
              <div style={{ 
                display: 'flex',
                alignItems: 'center',
                gap: '0.5rem'
            }}>
              <h2 style={{ 
                color: '#111827',
                fontSize: '1.5rem',
                fontWeight: '600',
                margin: 0
              }}>
                Pathway Analysis
              </h2>
              <div 
                style={{ position: 'relative', display: 'inline-flex' }}
                onMouseEnter={() => setHoveredTooltip('pathway-analysis')}
                onMouseLeave={() => setHoveredTooltip(null)}
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
                <SmartTooltip content="Comprehensive analysis of biological pathways affected by genetic variants. Shows which cellular processes are disrupted and their impact on cancer risk. Pathway burden scores indicate the severity of disruption, with higher percentages representing greater dysfunction." isVisible={hoveredTooltip === 'pathway-analysis'} triggerRef={null} />
              </div>
              </div>
              <SubtabDownloadButton 
                subtabContent={subtabContents.pathways}
                isPDFGenerating={isPDFGenerating}
                setIsPDFGenerating={setIsPDFGenerating}
              />
            </div>
            
            {/* Check if pathway analysis data is available */}
            {(() => {
              // First check structured_json for pathway analysis
              let pathwayData: PathwayAnalysisData | null = pipelineResults?.structured_json?.pathway_analysis || null;
              
              // If not found in structured_json, try to construct from pathway_burden_results
              if (!pathwayData && pipelineResults?.pathway_burden_results) {
                console.log('🔍 Constructing pathway data from pathway_burden_results:', pipelineResults.pathway_burden_results);
                
                const pathwayBurdenResults = pipelineResults.pathway_burden_results;
                const pathwayBurdenSummary = pipelineResults.pathway_burden_summary;
                
                // Transform pathway burden results into the format expected by frontend
                const disrupted_pathways: DisruptedPathway[] = [];
                const cancer_pathway_associations: Record<string, string[]> = {};
                
                // Convert pathway burden results to disrupted pathways format
                for (const [pathway_name, burden_result] of Object.entries(pathwayBurdenResults)) {
                  const typedBurdenResult = burden_result as PathwayBurdenResult;
                  const burden_score = typedBurdenResult.burden_score || 0;
                  if (burden_score > 0.1) {
                    // Create mutations list from damaging genes
                    const mutations: PathwayMutation[] = [];
                    if (typedBurdenResult.damaging_genes) {
                      typedBurdenResult.damaging_genes.forEach((gene: string) => {
                        mutations.push({
                          gene: gene,
                          type: "missense",
                          effect: `Damaging variant in ${gene}`
                        });
                      });
                    }
                    
                    disrupted_pathways.push({
                      name: pathway_name.replace("_", " ").replace(/\b\w/g, l => l.toUpperCase()),
                      pathway_id: pathway_name,
                      significance: Math.round(burden_score * 100 * 10) / 10, // Convert to percentage
                      affected_genes: typedBurdenResult.damaging_genes || [],
                      mutations: mutations,
                      description: typedBurdenResult.description || `${pathway_name} pathway`,
                      genes_affected_ratio: `${typedBurdenResult.genes_with_damaging || 0}/${typedBurdenResult.genes_in_pathway || 0}`
                    });
                  }
                }
                
                // Create cancer pathway associations based on high burden pathways
                const high_burden_pathways = pathwayBurdenSummary?.high_burden_pathways || [];
                if (high_burden_pathways.length > 0) {
                  // Map pathways to cancer types based on common associations
                  const pathway_cancer_mapping: Record<string, string[]> = {
                    "oncogenes": ["lung", "colon", "breast"],
                    "tumor_suppressors": ["breast", "lung", "colon", "prostate"],
                    "dna_repair": ["breast", "colon"],
                    "chromatin_remodeling": ["blood", "lung"],
                    "ras_mapk": ["lung", "colon", "prostate"],
                    "cell_cycle": ["breast", "lung", "prostate"],
                    "apoptosis": ["breast", "lung", "colon"],
                    "mismatch_repair": ["colon"],
                    "wnt_signaling": ["colon"],
                    "pi3k_akt": ["breast", "prostate"]
                  };
                  
                  high_burden_pathways.forEach((pathway: string) => {
                    const associated_cancers = pathway_cancer_mapping[pathway] || [];
                    associated_cancers.forEach((cancer: string) => {
                      if (!cancer_pathway_associations[cancer]) {
                        cancer_pathway_associations[cancer] = [];
                      }
                      cancer_pathway_associations[cancer].push(pathway);
                    });
                  });
                }
                
                // Create pathway analysis structure
                const constructedPathwayData: PathwayAnalysisData = {
                  disrupted_pathways: disrupted_pathways,
                  cancer_pathway_associations: cancer_pathway_associations,
                  summary: {
                    total_pathways_disrupted: disrupted_pathways.length,
                    highly_disrupted_pathways: disrupted_pathways.filter((p: DisruptedPathway) => p.significance > 50).length,
                    total_genes_affected: [...new Set(disrupted_pathways.flatMap((p: DisruptedPathway) => p.affected_genes))].length,
                    pathway_interaction_count: 0,
                    overall_burden_score: (pathwayBurdenSummary?.overall_burden_score ?? 0) as number,
                    high_burden_pathways: high_burden_pathways as string[]
                  }
                };
                pathwayData = constructedPathwayData;
               
                console.log('✅ Constructed pathway data:', pathwayData);
              }
              
              // Debug logging
              console.log('🔍 Final pathway data check:', {
                hasStructuredJson: !!pipelineResults?.structured_json?.pathway_analysis,
                hasPathwayBurdenResults: !!pipelineResults?.pathway_burden_results,
                finalPathwayData: !!pathwayData,
                disruptedPathways: pathwayData?.disrupted_pathways?.length || 0
              });
              
              if (!pathwayData || !pathwayData.disrupted_pathways || pathwayData.disrupted_pathways.length === 0) {
                // No pathway data available - show informative message
                return (
                  <div style={{
                    background: '#FFFFFF',
                    padding: '3rem',
                    borderRadius: '0.75rem',
                    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
                    border: '1px solid #E5E7EB',
                    textAlign: 'center'
                  }}>
                    <div style={{ 
                      color: '#6B7280',
                      fontSize: '1.125rem',
                      marginBottom: '1rem'
                    }}>
                      <svg 
                        style={{ width: '3rem', height: '3rem', marginBottom: '1rem' }}
                        fill="none" 
                        stroke="currentColor" 
                        viewBox="0 0 24 24"
                      >
                        <path 
                          strokeLinecap="round" 
                          strokeLinejoin="round" 
                          strokeWidth={1.5} 
                          d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" 
                        />
                      </svg>
                    </div>
                    <h3 style={{ 
                      color: '#111827',
                      fontSize: '1.25rem',
                      fontWeight: '600',
                      marginBottom: '0.5rem'
                    }}>
                      No Pathway Disruptions Found
                    </h3>
                    <p style={{ 
                      color: '#6B7280',
                      fontSize: '0.875rem',
                      maxWidth: '600px',
                      margin: '0 auto'
                    }}>
                      The analysis did not identify any significant pathway disruptions in your genetic data.
                      This suggests that your variants are not significantly affecting major cancer-related pathways.
                    </p>
                  </div>
                );
              }
              
              // If pathway data exists, use it
              return (
                <>
                  {/* Pathway Enrichment */}
                  <div id="pathway-disruption-analysis" style={{
                    background: '#FFFFFF',
                    padding: '2rem',
                    borderRadius: '0.75rem',
                    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
                    border: '1px solid #E5E7EB',
                    marginBottom: '2rem'
                  }}>
                    <div style={{ 
                      display: 'flex',
                      alignItems: 'center',
                      gap: '0.5rem',
                      marginBottom: '1rem'
                    }}>
                      <h3 style={{ 
                        color: '#111827',
                        fontSize: '1.125rem',
                        fontWeight: '600',
                        margin: 0
                      }}>
                        Affected Biological Pathways
                      </h3>
                      <div 
                        style={{ position: 'relative', display: 'inline-flex' }}
                        onMouseEnter={() => setHoveredTooltip('affected-pathways')}
                        onMouseLeave={() => setHoveredTooltip(null)}
                      >
                        <InformationCircleIcon 
                          style={{ width: '16px', height: '16px' }}
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
                        <SmartTooltip content="Biological pathways are cellular processes that work together to maintain normal cell function. When disrupted by genetic variants, these pathways can lead to cancer development. Each pathway shows the percentage of disruption based on the number and severity of variants affecting key genes in that pathway." isVisible={hoveredTooltip === 'affected-pathways'} triggerRef={null} />
                      </div>
                    </div>
                    
                    {/* Pathway Disruption Analysis */}
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem', marginBottom: '2rem' }}>
                      {pathwayData.disrupted_pathways.map((pathway: DisruptedPathway, index: number) => (
                        <div key={index} style={{
                          background: '#FFFFFF',
                          padding: '1.5rem',
                          border: `2px solid ${pathway.significance > 80 ? '#EF4444' : pathway.significance > 60 ? '#F59E0B' : '#3B82F6'}`,
                          borderRadius: '0.75rem',
                          boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
                          position: 'relative'
                        }}>
                          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                              <h4 style={{ 
                                color: '#111827',
                                fontSize: '1.125rem',
                                fontWeight: '600',
                                margin: 0
                              }}>
                                {pathway.name} Pathway
                              </h4>
                              <div 
                                style={{ position: 'relative', display: 'inline-flex' }}
                                onMouseEnter={() => setHoveredTooltip(`pathway-${index}`)}
                                onMouseLeave={() => setHoveredTooltip(null)}
                              >
                                <InformationCircleIcon 
                                  style={{ width: '16px', height: '16px' }}
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
                                <SmartTooltip 
                                  content={`${pathway.description}. This pathway contains ${pathway.genes_affected_ratio} genes with damaging variants. The ${pathway.significance}% disruption score indicates ${pathway.significance > 80 ? 'severe dysfunction requiring immediate medical attention' : pathway.significance > 60 ? 'moderate disruption with potential therapeutic implications' : 'mild disruption that should be monitored'}.`} 
                                  isVisible={hoveredTooltip === `pathway-${index}`} 
                                  triggerRef={null} 
                                />
                              </div>
                            </div>
                            <div style={{
                              background: pathway.significance > 80 ? '#EF4444' : pathway.significance > 60 ? '#F59E0B' : '#3B82F6',
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
                              background: pathway.significance > 80 ? '#EF4444' : pathway.significance > 60 ? '#F59E0B' : '#3B82F6',
                              height: '100%',
                              width: `${pathway.significance}%`,
                              transition: 'width 0.3s ease'
                            }}/>
                          </div>
                          
                          <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
                            {pathway.mutations.map((mutation: PathwayMutation, mutIndex: number) => (
                              <div key={mutIndex} style={{
                                display: 'flex',
                                alignItems: 'center',
                                padding: '0.75rem',
                                background: '#F9FAFB',
                                borderRadius: '0.5rem',
                                border: '1px solid #E5E7EB',
                                position: 'relative'
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
                                <div 
                                  style={{ position: 'relative', display: 'inline-flex' }}
                                  onMouseEnter={() => setHoveredTooltip(`mutation-${index}-${mutIndex}`)}
                                  onMouseLeave={() => setHoveredTooltip(null)}
                                >
                                  <InformationCircleIcon 
                                    style={{ width: '16px', height: '16px' }}
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
                                  <SmartTooltip 
                                    content={`${mutation.gene} is a key gene in the ${pathway.name} pathway. ${mutation.type === 'missense' ? 'Missense mutations change one amino acid and may alter protein function.' : mutation.type === 'frameshift' ? 'Frameshift mutations shift the reading frame and typically cause severe protein dysfunction.' : 'This mutation type affects protein structure and function.'} The variant contributes to the overall pathway disruption score.`} 
                                    isVisible={hoveredTooltip === `mutation-${index}-${mutIndex}`} 
                                    triggerRef={null} 
                                  />
                                </div>
                              </div>
                            ))}
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                  
                  {/* Cancer Risk by Pathway - NOW SEPARATE FROM pathway-disruption-analysis */}
                  <div id="cancer-risk-by-pathways" style={{
                    background: '#FFFFFF',
                    padding: '2rem',
                    borderRadius: '0.75rem',
                    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
                    border: '1px solid #E5E7EB',
                    marginBottom: '2rem'
                  }}>
                    <div style={{ 
                      display: 'flex',
                      alignItems: 'center',
                      gap: '0.5rem',
                      marginBottom: '1rem'
                    }}>
                      <h4 style={{ 
                        color: '#111827',
                        fontSize: '1.125rem',
                        fontWeight: '600',
                        margin: 0
                      }}>
                        Cancer Risk by Pathway Disruption
                      </h4>
                      <div 
                        style={{ position: 'relative', display: 'inline-flex' }}
                        onMouseEnter={() => setHoveredTooltip('cancer-risk-pathways')}
                        onMouseLeave={() => setHoveredTooltip(null)}
                      >
                        <InformationCircleIcon 
                          style={{ width: '16px', height: '16px' }}
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
                        <SmartTooltip content="Shows how pathway disruptions translate to cancer risk for different cancer types. Each cancer type is associated with specific pathways - when those pathways are disrupted, the risk for that cancer increases. The percentage represents the calculated risk based on the severity of pathway disruption and clinical evidence." isVisible={hoveredTooltip === 'cancer-risk-pathways'} triggerRef={null} />
                      </div>
                    </div>
                    
                    <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '1rem' }}>
                      {Object.entries(pathwayData.cancer_pathway_associations).map(([cancer, pathways]) => {
                        const riskFinding = genomicData.risk_findings.find(r => r.cancer_type === cancer);
                        if (!riskFinding) return null;
                        
                        return (
                          <div key={cancer} style={{
                            padding: '1rem',
                            border: `2px solid ${getRiskColor(riskFinding.risk_level)}`,
                            borderRadius: '0.5rem',
                            background: `${getRiskColor(riskFinding.risk_level)}10`,
                            position: 'relative'
                          }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginBottom: '0.5rem' }}>
                              <h5 style={{ 
                                color: '#111827',
                                fontSize: '1rem',
                                fontWeight: '600',
                                margin: 0,
                                textTransform: 'capitalize'
                              }}>
                                {cancer} Cancer
                              </h5>
                              <div 
                                style={{ position: 'relative', display: 'inline-flex' }}
                                onMouseEnter={() => setHoveredTooltip(`cancer-${cancer}`)}
                                onMouseLeave={() => setHoveredTooltip(null)}
                              >
                                <InformationCircleIcon 
                                  style={{ width: '14px', height: '14px' }}
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
                                <SmartTooltip 
                                  content={`${cancer.charAt(0).toUpperCase() + cancer.slice(1)} cancer risk is ${riskFinding.risk_percentage}% based on disruption in ${(pathways as string[]).length} pathway${(pathways as string[]).length > 1 ? 's' : ''}: ${(pathways as string[]).join(', ')}. This represents a ${riskFinding.risk_level} risk level requiring ${riskFinding.risk_level === 'high' ? 'immediate medical attention and enhanced screening' : riskFinding.risk_level === 'medium' ? 'regular monitoring and preventive measures' : 'standard screening protocols'}.`} 
                                  isVisible={hoveredTooltip === `cancer-${cancer}`} 
                                  triggerRef={null} 
                                />
                              </div>
                            </div>
                            <div style={{ 
                              fontSize: '1.5rem',
                              fontWeight: 'bold',
                              color: getRiskColor(riskFinding.risk_level),
                              marginBottom: '0.5rem'
                            }}>
                              {riskFinding.risk_percentage}%
                            </div>
                            <div style={{ display: 'flex', flexDirection: 'column', gap: '0.25rem' }}>
                              {(pathways as string[]).map((pathway: string, pathIndex: number) => (
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
                        );
                      })}
                    </div>
                  </div>
                </>
              );
            })()}
          </div>
        );
      
      case 'clinical':
        return (
          <div style={{ padding: '2rem' }}>
            <div style={{ 
              display: 'flex',
              justifyContent: 'space-between',
              alignItems: 'center',
              marginBottom: '1.5rem'
            }}>
              <div style={{ 
                display: 'flex',
                alignItems: 'center',
                gap: '0.5rem'
            }}>
              <h2 style={{ 
                color: '#111827',
                fontSize: '1.5rem',
                fontWeight: '600',
                margin: 0
              }}>
                In-Depth Report & Recommendations
              </h2>
              <div 
                style={{ position: 'relative', display: 'inline-flex' }}
                onMouseEnter={() => setHoveredTooltip('clinical-report')}
                onMouseLeave={() => setHoveredTooltip(null)}
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
                <SmartTooltip content="Comprehensive in-depth analysis combining survival predictions with personalized screening recommendations. Integrates genetic risk assessment with evidence-based research guidelines to provide actionable healthcare recommendations tailored to your specific genetic profile." isVisible={hoveredTooltip === 'clinical-report'} triggerRef={null} />
              </div>
              </div>
              <SubtabDownloadButton 
                subtabContent={subtabContents.clinical}
                isPDFGenerating={isPDFGenerating}
                setIsPDFGenerating={setIsPDFGenerating}
              />
            </div>
            
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
                <div style={{ 
                  display: 'flex',
                  alignItems: 'center',
                  gap: '0.5rem'
                }}>
                  <h3 style={{ 
                    color: '#111827',
                    fontSize: '1.125rem',
                    fontWeight: '600',
                    margin: 0
                  }}>
                    Survival Analysis vs Population Average
                  </h3>
                  <div 
                    style={{ position: 'relative', display: 'inline-flex' }}
                    onMouseEnter={() => setHoveredTooltip('survival-analysis')}
                    onMouseLeave={() => setHoveredTooltip(null)}
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
                    <SmartTooltip content="Comparative survival analysis showing how your genetic risk profile affects life expectancy compared to the general population. The red solid line represents your personalized risk profile, while the blue dashed line shows population averages. Early intervention and enhanced screening can significantly improve outcomes." isVisible={hoveredTooltip === 'survival-analysis'} triggerRef={null} />
                  </div>
                </div>
              </div>
              
              {/* Check if survival analysis data is available */}
              {(() => {
                const survivalData = pipelineResults?.structured_json?.survival_analysis;
                
                console.log('🔍 SURVIVAL DEBUG: survivalData exists:', !!survivalData);
                console.log('🔍 SURVIVAL DEBUG: patient_profile exists:', !!survivalData?.patient_profile);
                console.log('🔍 SURVIVAL DEBUG: estimated_survival exists:', !!survivalData?.patient_profile?.estimated_survival);
                console.log('🔍 SURVIVAL DEBUG: population_average exists:', !!survivalData?.population_average);
                
                if (!survivalData || 
                    !survivalData.patient_profile || 
                    !survivalData.patient_profile.estimated_survival || 
                    !Array.isArray(survivalData.patient_profile.estimated_survival) ||
                    !survivalData.population_average ||
                    !Array.isArray(survivalData.population_average)) {
                  // No proper survival data - show informative message
                  return (
                    <div style={{ 
                      padding: '3rem',
                      textAlign: 'center',
                      color: '#6B7280',
                      fontSize: '0.875rem'
                    }}>
                      <p style={{ marginBottom: '1rem' }}>
                        Survival analysis visualization is not yet available from the backend.
                      </p>
                      <p style={{ fontSize: '0.75rem', color: '#9CA3AF' }}>
                        This feature will provide personalized survival curves based on your genetic risk profile.
                      </p>
                      <div style={{ 
                        marginTop: '2rem',
                        padding: '1rem',
                        background: '#F3F4F6',
                        borderRadius: '0.5rem',
                        fontSize: '0.7rem',
                        color: '#6B7280'
                      }}>
                        <p style={{ marginBottom: '0.5rem' }}>
                          <strong>Expected Data Structure:</strong>
                        </p>
                        <ul style={{ 
                          textAlign: 'left',
                          margin: 0,
                          paddingLeft: '1.5rem',
                          fontSize: '0.65rem'
                        }}>
                          <li>survival_analysis.patient_profile.estimated_survival (array)</li>
                          <li>survival_analysis.population_average (array)</li>
                          <li>survival_analysis.patient_profile.risk_category (string)</li>
                        </ul>
                        <p style={{ marginTop: '0.5rem', fontStyle: 'italic' }}>
                          Backend implementation pending
                        </p>
                      </div>
                    </div>
                  );
                }
                
                // Use real survival data - with additional safety checks
                const patientSurvival = survivalData.patient_profile.estimated_survival || [];
                const populationAverage = survivalData.population_average || [];
                const riskCategory = survivalData.patient_profile.risk_category || 'Unknown';
                
                return (
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
                    {patientSurvival.map((point, i) => (
                      <g key={i}>
                        <line x1={80 + i * 80} y1="50" x2={80 + i * 80} y2="250" stroke="#E5E7EB" strokeWidth="1"/>
                        <text x={80 + i * 80} y="270" fill="#6B7280" fontSize="12" textAnchor="middle">{point.age || `Age ${i}`}</text>
                      </g>
                    ))}
                    
                    {/* Population average survival curve */}
                    {populationAverage.length > 0 && (
                      <path 
                        d={`M ${populationAverage.map((point, i) => 
                          `${80 + i * 80} ${250 - ((point.probability || 0) * 200)}`
                        ).join(' L ')}`}
                        stroke="#2563EB" 
                        strokeWidth="3" 
                        fill="none"
                        strokeDasharray="5,5"
                      />
                    )}
                    
                    {/* Patient risk profile survival curve */}
                    {patientSurvival.length > 0 && (
                      <path 
                        d={`M ${patientSurvival.map((point, i) => 
                          `${80 + i * 80} ${250 - ((point.probability || 0) * 200)}`
                        ).join(' L ')}`}
                        stroke="#EF4444" 
                        strokeWidth="3" 
                        fill="none"
                      />
                    )}
                    
                    {/* Legend */}
                    <g transform="translate(500, 80)">
                      <line x1="0" y1="0" x2="20" y2="0" stroke="#2563EB" strokeWidth="3" strokeDasharray="5,5"/>
                      <text x="25" y="5" fill="#2563EB" fontSize="12" fontWeight="600">Population Average</text>
                      <line x1="0" y1="20" x2="20" y2="20" stroke="#EF4444" strokeWidth="3"/>
                      <text x="25" y="25" fill="#EF4444" fontSize="12" fontWeight="600">{riskCategory} Risk Profile</text>
                    </g>
                    
                    {/* Axis labels */}
                    <text x="400" y="290" fill="#111827" fontSize="14" textAnchor="middle" fontWeight="600">Age</text>
                    <text x="30" y="150" fill="#111827" fontSize="14" textAnchor="middle" fontWeight="600" transform="rotate(-90 30 150)">Survival Probability</text>
                  </svg>
                );
              })()}
            </div>
            
            {/* Clinical Recommendations */}
            <div id="clinical-recommendations" style={{
              background: '#FFFFFF',
              padding: '2rem',
              borderRadius: '0.75rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)',
              border: '1px solid #E5E7EB'
            }}>
              <div style={{ 
                display: 'flex',
                alignItems: 'center',
                gap: '0.5rem',
                marginBottom: '1rem'
              }}>
                <h3 style={{ 
                  color: '#111827',
                  fontSize: '1.125rem',
                  fontWeight: '600',
                  margin: 0
                }}>
                  In-Depth Recommendations
                </h3>
                <div 
                  style={{ position: 'relative', display: 'inline-flex' }}
                  onMouseEnter={() => setHoveredTooltip('clinical-recommendations')}
                  onMouseLeave={() => setHoveredTooltip(null)}
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
                  <SmartTooltip content="Personalized screening and prevention recommendations based on your genetic risk profile. Each recommendation includes specific cancer types, risk percentages, screening protocols, and prevention strategies. High-risk findings require enhanced surveillance while medium-risk findings may benefit from earlier or more frequent screening." isVisible={hoveredTooltip === 'clinical-recommendations'} triggerRef={null} />
                </div>
              </div>
              
              {/* Check if clinical recommendations are available from backend */}
              {(() => {
                const clinicalRecs = pipelineResults?.structured_json?.clinical_recommendations;
                
                if (clinicalRecs && Array.isArray(clinicalRecs) && clinicalRecs.length > 0) {
                  // Use real clinical recommendations
                  return (
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                      {clinicalRecs.map((rec: ClinicalRecommendation, index: number) => (
                        <div key={index} style={{
                          padding: '1.5rem',
                          border: `2px solid ${getRiskColor(rec.risk_level)}`,
                          borderRadius: '0.75rem',
                          background: rec.risk_level === 'high' ? '#FEF2F2' : 
                                      rec.risk_level === 'medium' ? '#FFFBEB' : '#F0FDF4'
                        }}>
                          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '1rem' }}>
                            <h4 style={{ 
                              color: '#111827',
                              fontSize: '1rem',
                              fontWeight: '600',
                              textTransform: 'capitalize'
                            }}>
                              {rec.cancer_type} Cancer Screening
                            </h4>
                            <div style={{
                              background: getRiskColor(rec.risk_level),
                              color: '#FFFFFF',
                              padding: '0.25rem 0.5rem',
                              borderRadius: '0.375rem',
                              fontSize: '0.75rem',
                              fontWeight: '600'
                            }}>
                              {rec.risk_percentage}% Risk
                            </div>
                          </div>
                          <p style={{ 
                            color: '#4B5563',
                            fontSize: '0.875rem',
                            marginBottom: '1rem'
                          }}>
                            {rec.recommendation}
                          </p>
                          {rec.screening_protocol && (
                            <div style={{
                              background: '#FFFFFF',
                              padding: '1rem',
                              borderRadius: '0.5rem',
                              border: '1px solid #E5E7EB',
                              marginBottom: '1rem'
                            }}>
                              <h5 style={{ 
                                color: '#111827',
                                fontSize: '0.875rem',
                                fontWeight: '600',
                                marginBottom: '0.5rem'
                              }}>
                                Screening Protocol:
                              </h5>
                              <ul style={{ 
                                margin: 0,
                                paddingLeft: '1.5rem',
                                fontSize: '0.75rem',
                                color: '#4B5563'
                              }}>
                                <li>Test: {rec.screening_protocol.test}</li>
                                <li>Frequency: {rec.screening_protocol.frequency}</li>
                                <li>Start: {rec.screening_protocol.start_age}</li>
                              </ul>
                            </div>
                          )}
                          {rec.prevention_options && rec.prevention_options.length > 0 && (
                            <div style={{
                              background: '#F9FAFB',
                              padding: '1rem',
                              borderRadius: '0.5rem'
                            }}>
                              <h5 style={{ 
                                color: '#111827',
                                fontSize: '0.875rem',
                                fontWeight: '600',
                                marginBottom: '0.5rem'
                              }}>
                                Prevention Options:
                              </h5>
                              <ul style={{ 
                                margin: 0,
                                paddingLeft: '1.5rem',
                                fontSize: '0.75rem',
                                color: '#4B5563'
                              }}>
                                {rec.prevention_options.map((option: string, i: number) => (
                                  <li key={i}>{option}</li>
                                ))}
                              </ul>
                            </div>
                          )}
                        </div>
                      ))}
                    </div>
                  );
                } else {
                  // Use basic recommendations based on risk scores
                  return (
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                      {genomicData.risk_findings.filter(r => r.risk_level === 'high' || r.risk_level === 'medium').map((risk, index) => (
                        <div key={index} style={{
                          padding: '1.5rem',
                          border: `2px solid ${getRiskColor(risk.risk_level)}`,
                          borderRadius: '0.75rem',
                          background: risk.risk_level === 'high' ? '#FEF2F2' : '#FFFBEB'
                        }}>
                          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: '1rem' }}>
                            <div style={{ 
                              display: 'flex',
                              alignItems: 'center',
                              gap: '0.5rem'
                            }}>
                              <h4 style={{ 
                                color: '#111827',
                                fontSize: '1rem',
                                fontWeight: '600',
                                textTransform: 'capitalize',
                                margin: 0
                              }}>
                                {risk.cancer_type} Cancer Screening
                              </h4>
                              <div 
                                style={{ position: 'relative', display: 'inline-flex' }}
                                onMouseEnter={() => setHoveredTooltip(`screening-${risk.cancer_type}`)}
                                onMouseLeave={() => setHoveredTooltip(null)}
                              >
                                <InformationCircleIcon 
                                  style={{ width: '16px', height: '16px' }}
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
                                <SmartTooltip content={`Personalized screening recommendations for ${risk.cancer_type} cancer based on your genetic risk profile. Risk level: ${risk.risk_level}. Enhanced screening protocols may include earlier start ages, more frequent testing, or additional screening modalities compared to standard population guidelines.`} isVisible={hoveredTooltip === `screening-${risk.cancer_type}`} triggerRef={null} />
                              </div>
                            </div>
                            <div 
                              style={{ position: 'relative', display: 'inline-flex' }}
                              onMouseEnter={() => setHoveredTooltip(`risk-percentage-${risk.cancer_type}`)}
                              onMouseLeave={() => setHoveredTooltip(null)}
                            >
                              <div style={{
                                background: getRiskColor(risk.risk_level),
                                color: '#FFFFFF',
                                padding: '0.25rem 0.5rem',
                                borderRadius: '0.375rem',
                                fontSize: '0.75rem',
                                fontWeight: '600'
                              }}>
                                {risk.risk_percentage}% Risk
                              </div>
                              <SmartTooltip content={`${risk.risk_percentage}% lifetime risk of developing ${risk.cancer_type} cancer based on your genetic variants. This is ${risk.risk_level === 'high' ? 'significantly higher' : risk.risk_level === 'medium' ? 'moderately higher' : 'similar to or lower'} than the general population average. Risk is calculated using machine learning models trained on clinical outcomes and population data.`} isVisible={hoveredTooltip === `risk-percentage-${risk.cancer_type}`} triggerRef={null} />
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
                            border: '1px solid #E5E7EB'
                          }}>
                            <div style={{ 
                              display: 'flex',
                              alignItems: 'center',
                              gap: '0.5rem',
                              marginBottom: '0.5rem'
                            }}>
                              <h5 style={{ 
                                color: '#111827',
                                fontSize: '0.875rem',
                                fontWeight: '600',
                                margin: 0
                              }}>
                                Affected Genes:
                              </h5>
                              <div 
                                style={{ position: 'relative', display: 'inline-flex' }}
                                onMouseEnter={() => setHoveredTooltip(`affected-genes-${risk.cancer_type}`)}
                                onMouseLeave={() => setHoveredTooltip(null)}
                              >
                                <InformationCircleIcon 
                                  style={{ width: '14px', height: '14px' }}
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
                                <SmartTooltip content={`Genes containing variants that contribute to ${risk.cancer_type} cancer risk. These genes are involved in critical cellular pathways like DNA repair, cell cycle control, and tumor suppression. Variants in these genes can increase cancer susceptibility and influence screening recommendations.`} isVisible={hoveredTooltip === `affected-genes-${risk.cancer_type}`} triggerRef={null} />
                              </div>
                            </div>
                            <div style={{ display: 'flex', gap: '0.5rem', flexWrap: 'wrap' }}>
                              {risk.affected_genes.length > 0 ? risk.affected_genes.map((gene, geneIndex) => (
                                <span key={geneIndex} style={{
                                  background: getRiskColor(risk.risk_level),
                                  color: '#FFFFFF',
                                  padding: '0.25rem 0.5rem',
                                  borderRadius: '0.375rem',
                                  fontSize: '0.75rem',
                                  fontWeight: '600'
                                }}>
                                  {gene}
                                </span>
                              )) : (
                                <span style={{ color: '#6B7280', fontSize: '0.75rem' }}>
                                  No specific genes identified
                                </span>
                              )}
                            </div>
                          </div>
                        </div>
                      ))}
                      {genomicData.risk_findings.filter(r => r.risk_level === 'high' || r.risk_level === 'medium').length === 0 && (
                        <div style={{ 
                          padding: '2rem',
                          textAlign: 'center',
                          color: '#6B7280',
                          fontSize: '0.875rem'
                        }}>
                          No high or medium risk findings requiring special screening recommendations.
                          Standard population screening guidelines apply.
                        </div>
                      )}
                    </div>
                  );
                }
              })()}
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
        minHeight: '100vh',
        padding: '1rem 0'
      }}>
        <div style={{ 
          maxWidth: '1400px',
          margin: '0 auto',
          padding: '0 1rem'
        }}>
          {/* Header */}
          <div style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            marginBottom: '1rem'
          }}>
            <div>
              <h1 style={{
                fontSize: '1.5rem',
                fontWeight: 'bold',
                color: '#111827',
                marginBottom: '0.25rem'
              }}>
                In-Depth Analysis
              </h1>
              <p style={{
                color: '#4B5563',
                fontSize: '0.875rem'
              }}>
                Comprehensive genomic analysis • {sidebarData.riskLevel}
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
                  padding: '0.5rem 1rem',
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
            gridTemplateColumns: 'minmax(260px, 280px) 1fr', 
            gap: '1rem',
            alignItems: 'start'
          }}>
            {/* Sidebar */}
            <div style={{
              background: '#FFFFFF',
              borderRadius: '0.5rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB',
              position: 'sticky',
              top: '1rem',
              alignSelf: 'start'
            }}>
              {/* Risk Info */}
              <div style={{
                background: '#F9FAFB',
                padding: '1rem',
                borderBottom: '1px solid #E5E7EB'
              }}>
                <div style={{ 
                  fontSize: '1rem',
                  fontWeight: '600',
                  marginBottom: '0.5rem',
                  color: '#111827'
                }}>
                  {sidebarData.riskLevel}
                </div>
                <div style={{ 
                  fontSize: '0.8rem',
                  lineHeight: '1.3',
                  color: '#4B5563'
                }}>
                  <span dangerouslySetInnerHTML={{ __html: sidebarData.details }} />
                </div>
                <div style={{
                  background: '#DBEAFE',
                  color: '#2563EB',
                  padding: '0.5rem',
                  borderRadius: '0.375rem',
                  marginTop: '0.75rem',
                  fontSize: '0.8rem',
                  fontWeight: '500',
                  textAlign: 'center'
                }}>
                  Risk Score: {sidebarData.riskScore}
                </div>
              </div>

              {/* Navigation */}
              <div style={{ padding: '0.75rem' }}>
                <div style={{
                  fontSize: '0.75rem',
                  fontWeight: '600',
                  textTransform: 'uppercase',
                  letterSpacing: '0.05em',
                  color: '#4B5563',
                  marginBottom: '0.75rem'
                }}>
                  Analysis Workflow
                </div>
                
                {[
                  { id: 'analysis', label: 'Genomic Analysis' },
                  { id: 'variants', label: 'Variant Heatmap' },
                  { id: 'pathways', label: 'Pathway Analysis' },
                  { id: 'clinical', label: 'In-Depth Report' }
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
                  In-Depth Alerts
                </div>
                
                {sidebarData.alerts.map((alert: Alert, index: number) => (
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
                  validation={confidenceCheckValidation} 
                  isDetailed={true} 
                />
              </div>
            </div>

            {/* Main Content */}
            <div style={{
              background: '#FFFFFF',
              borderRadius: '0.5rem',
              boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1), 0 1px 2px 0 rgba(0, 0, 0, 0.06)',
              border: '1px solid #E5E7EB'
            }}>
              {renderTabContent()}
            </div>
          </div>
          
          {/* Analysis Disclaimer */}
          <DetailedDisclaimer />
        </div>
      </section>
    </Layout>
  );
};

export default ClinicalViewPage; 