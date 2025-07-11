"""
Report formatter for generating Markdown output with glossary support.
Includes template-based fallback for non-LLM mode.
"""

import logging
import os
import re
from datetime import datetime
from typing import Dict, Any, List, Optional, Set
from .config import ReportConfig

logger = logging.getLogger(__name__)


class MedicalGlossary:
    """Medical and genomic terminology glossary."""
    
    def __init__(self):
        self.terms = {
            "allele frequency": "The proportion of chromosomes in a population that carry a specific variant of a gene.",
            "CADD score": "Combined Annotation Dependent Depletion score - predicts the deleteriousness of genetic variants.",
            "clinical significance": "The medical importance of a genetic variant, ranging from benign to pathogenic.",
            "frameshift": "A genetic mutation that causes a shift in the reading frame of DNA, often resulting in a nonfunctional protein.",
            "heterozygous": "Having two different versions (alleles) of a gene.",
            "homozygous": "Having two identical versions (alleles) of a gene.",
            "missense variant": "A genetic change that results in a different amino acid being incorporated into the protein.",
            "nonsense variant": "A genetic change that creates a premature stop codon, resulting in a truncated protein.",
            "pathogenic": "A genetic variant that is known to cause disease.",
            "polygenic risk score": "A score that estimates disease risk based on the combined effect of many genetic variants.",
            "TCGA": "The Cancer Genome Atlas - a comprehensive database of cancer genomics data.",
            "variant of uncertain significance": "A genetic variant whose impact on health is not yet known.",
            "HGVS nomenclature": "Human Genome Variation Society standard notation for describing genetic variants.",
            "penetrance": "The likelihood that a person with a genetic variant will develop the associated condition.",
            "loss of function": "A genetic variant that reduces or eliminates the normal function of a gene.",
            "gain of function": "A genetic variant that increases or creates new gene function.",
            "splice site": "DNA sequences that mark the boundaries between exons and introns in genes.",
            "copy number variation": "Differences in the number of copies of a particular gene or chromosomal region.",
            "somatic mutation": "A genetic change that occurs in body cells (not inherited).",
            "germline mutation": "A genetic change present in reproductive cells that can be passed to offspring.",
            "tumor suppressor gene": "A gene that normally prevents cancer by controlling cell growth and division.",
            "oncogene": "A gene that, when mutated or overexpressed, can contribute to cancer development."
        }
    
    def get_definition(self, term: str) -> Optional[str]:
        """Get definition for a medical term."""
        return self.terms.get(term.lower())
    
    def find_terms_in_text(self, text: str) -> Set[str]:
        """Find medical terms that appear in the given text."""
        found_terms = set()
        text_lower = text.lower()
        
        for term in self.terms.keys():
            if term in text_lower:
                found_terms.add(term)
        
        return found_terms
    
    def add_term(self, term: str, definition: str):
        """Add a new term to the glossary."""
        self.terms[term.lower()] = definition


class ReportFormatter:
    """Formats genomic reports with LLM content or template fallback."""
    
    def __init__(self, config: ReportConfig, output_dir: str = "./reports"):
        self.config = config
        self.output_dir = output_dir
        self.glossary = MedicalGlossary()
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
    
    def format_report(self, 
                     data: Dict[str, Any], 
                     llm_content: Optional[str] = None,
                     report_id: Optional[str] = None) -> Dict[str, str]:
        """
        Format a complete report with LLM content or fallback templates.
        
        Args:
            data: Structured JSON data from pipeline
            llm_content: Generated content from LLM (if available)
            report_id: Unique identifier for the report
            
        Returns:
            Dictionary with paths to generated files
        """
        if report_id is None:
            report_id = f"report_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        # Generate markdown content
        if llm_content:
            markdown_content = self._format_llm_report(data, llm_content, report_id)
        else:
            markdown_content = self._format_template_report(data, report_id)
        
        # Save markdown file
        markdown_path = os.path.join(self.output_dir, f"{report_id}.md")
        with open(markdown_path, 'w', encoding='utf-8') as f:
            f.write(markdown_content)
        
        logger.info(f"Generated markdown report: {markdown_path}")
        
        result_paths = {"markdown": markdown_path}
        
        # Generate additional formats if requested
        if "html" in self.config.output_formats:
            html_path = self._convert_to_html(markdown_content, report_id)
            if html_path:
                result_paths["html"] = html_path
        
        if "pdf" in self.config.output_formats:
            pdf_path = self._convert_to_pdf(markdown_content, report_id)
            if pdf_path:
                result_paths["pdf"] = pdf_path
        
        if "txt" in self.config.output_formats:
            txt_path = self._convert_to_txt(markdown_content, report_id)
            if txt_path:
                result_paths["txt"] = txt_path
        
        return result_paths
    
    def _format_llm_report(self, data: Dict[str, Any], llm_content: str, report_id: str) -> str:
        """Format report using LLM-generated content."""
        
        # Extract metadata
        report_metadata = data.get("report_metadata", {})
        
        # Build header
        markdown = self._build_header(data, report_id, llm_enhanced=True)
        
        # Add LLM-generated content
        markdown += "\n## Clinical Assessment\n\n"
        markdown += llm_content
        
        # Add technical appendix if requested
        if self.config.include_technical_appendix:
            markdown += "\n\n"
            markdown += self._build_technical_appendix(data)
        
        # Add glossary if requested
        if self.config.include_glossary:
            markdown += "\n\n"
            markdown += self._build_glossary(llm_content)
        
        # Add footer
        markdown += "\n\n"
        markdown += self._build_footer(data, report_id)
        
        return markdown
    
    def _format_template_report(self, data: Dict[str, Any], report_id: str) -> str:
        """Format report using template-based fallback."""
        
        # Build header with fallback indicator
        markdown = self._build_header(data, report_id, llm_enhanced=False)
        
        # Add fallback content warning
        markdown += f"\n> **{self.config.dev_mode_indicator}**: {self.config.fallback_mode_indicator}\n\n"
        
        # Generate template-based sections
        markdown += self._build_template_overview(data)
        markdown += self._build_template_risk_summary(data)
        markdown += self._build_template_variants(data)
        markdown += self._build_template_recommendations(data)
        
        # Add technical appendix
        if self.config.include_technical_appendix:
            markdown += self._build_technical_appendix(data)
        
        # Add basic glossary
        if self.config.include_glossary:
            markdown += self._build_basic_glossary()
        
        # Add footer
        markdown += self._build_footer(data, report_id)
        
        return markdown
    
    def _build_header(self, data: Dict[str, Any], report_id: str, llm_enhanced: bool = False) -> str:
        """Build report header section."""
        
        header = f"""# Genomic Risk Assessment Report

**Generated:** {datetime.now().strftime('%B %d, %Y at %H:%M')}  
"""
        
        if llm_enhanced:
            header += "**Analysis Type:** AI-Enhanced Clinical Report  \n"
        else:
            header += "**Analysis Type:** Template-Based Report  \n"
        
        header += "\n---\n"
        
        return header
    
    def _build_template_overview(self, data: Dict[str, Any]) -> str:
        """Build template-based overview section."""
        
        summary = data.get("summary", {})
        risk_assessment = data.get("risk_assessment", {})
        
        total_variants = summary.get("total_variants_found", 0)
        qc_variants = summary.get("variants_passed_qc", 0)
        
        # Count high-risk findings
        high_risk_count = 0
        if risk_assessment and "scores" in risk_assessment:
            high_risk_count = sum(1 for score in risk_assessment["scores"].values() 
                                 if score > self.config.risk_threshold)
        
        overview = f"""## Case Overview

This genomic risk assessment analyzed **{total_variants}** genetic variants, of which **{qc_variants}** passed quality control filters. The analysis identified **{high_risk_count}** cancer types with elevated risk levels above the {self.config.risk_threshold}% baseline threshold.

The assessment utilized multiple computational approaches including population frequency analysis, tumor database comparisons, and machine learning risk models to evaluate cancer predisposition based on the genetic profile.

"""
        
        return overview
    
    def _build_template_risk_summary(self, data: Dict[str, Any]) -> str:
        """Build template-based risk summary section."""
        
        risk_assessment = data.get("risk_assessment", {})
        if not risk_assessment or "scores" not in risk_assessment:
            return "## Risk Summary\n\nNo risk assessment data available.\n\n"
        
        scores = risk_assessment["scores"]
        risk_genes = risk_assessment.get("risk_genes", {})
        
        # Filter for high-risk findings
        high_risk_findings = []
        baseline_findings = []
        
        for cancer_type, score in scores.items():
            if score > self.config.risk_threshold:
                high_risk_findings.append((cancer_type, score, risk_genes.get(cancer_type, [])))
            else:
                baseline_findings.append((cancer_type, score))
        
        # Sort by risk score
        high_risk_findings.sort(key=lambda x: x[1], reverse=True)
        
        risk_summary = "## Cancer Risk Assessment\n\n"
        
        if high_risk_findings:
            risk_summary += f"### Elevated Risk Findings (>{self.config.risk_threshold}% baseline)\n\n"
            
            for cancer_type, score, genes in high_risk_findings:
                risk_level = "High" if score >= 30 else "Moderate" if score >= 15 else "Elevated"
                risk_summary += f"**{cancer_type.title()} Cancer:** {score:.1f}% ({risk_level} Risk)\n"
                if genes:
                    risk_summary += f"- Associated genes: {', '.join(genes)}\n"
                risk_summary += "\n"
        else:
            risk_summary += "### No Elevated Risk Findings\n\n"
            risk_summary += "All analyzed cancer types show risk levels within normal baseline ranges. "
            risk_summary += "This suggests an average population risk profile.\n\n"
        
        if baseline_findings:
            risk_summary += f"### Baseline Risk Cancers (≤{self.config.risk_threshold}%)\n\n"
            baseline_cancers = [cancer.title() for cancer, _ in baseline_findings]
            risk_summary += f"{', '.join(baseline_cancers)}\n\n"
        
        return risk_summary
    
    def _build_template_variants(self, data: Dict[str, Any]) -> str:
        """Build template-based variants section."""
        
        variant_details = data.get("variant_details", [])
        if not variant_details:
            return ""
        
        variants_section = "## Key Genetic Variants\n\n"
        
        # Show top 5 variants
        key_variants = variant_details[:5]
        
        for i, variant in enumerate(key_variants, 1):
            quality = variant.get("quality_metrics", {})
            
            variants_section += f"### Variant {i}: {variant.get('gene', 'Unknown')}\n\n"
            variants_section += f"- **Location:** {variant.get('variant', 'Unknown')}\n"
            variants_section += f"- **Consequence:** {variant.get('consequence', 'Unknown')}\n"
            
            if variant.get('hgvs_c'):
                variants_section += f"- **HGVS (coding):** {variant['hgvs_c']}\n"
            if variant.get('hgvs_p'):
                variants_section += f"- **HGVS (protein):** {variant['hgvs_p']}\n"
            
            variants_section += f"- **Quality Score:** {quality.get('quality', 0)}\n"
            variants_section += f"- **Read Depth:** {quality.get('depth', 0)}\n"
            variants_section += f"- **Allele Frequency:** {quality.get('allele_freq', 0):.3f}\n"
            
            # Add CADD score if available
            cadd_scores = variant.get("cadd_scores", {})
            if cadd_scores and cadd_scores.get("phred"):
                score = cadd_scores["phred"]
                interpretation = "High impact" if score > 20 else "Moderate impact" if score > 10 else "Low impact"
                variants_section += f"- **CADD Score:** {score:.1f} ({interpretation})\n"
            
            # Add TCGA matches if available
            tcga_matches = variant.get("tcga_matches", {})
            if tcga_matches:
                variants_section += "- **TCGA Database Matches:**\n"
                for cancer_type, match_data in tcga_matches.items():
                    if match_data and match_data.get("frequency_percent", 0) > 0:
                        freq = match_data["frequency_percent"]
                        variants_section += f"  - {cancer_type.title()}: {freq:.1f}% frequency\n"
            
            variants_section += "\n"
        
        return variants_section
    
    def _build_template_recommendations(self, data: Dict[str, Any]) -> str:
        """Build template-based recommendations section."""
        
        risk_assessment = data.get("risk_assessment", {})
        
        recommendations = "## Clinical Recommendations\n\n"
        
        # Check for high-risk findings
        has_high_risk = False
        if risk_assessment and "scores" in risk_assessment:
            has_high_risk = any(score > self.config.risk_threshold 
                              for score in risk_assessment["scores"].values())
        
        if has_high_risk:
            recommendations += """### Immediate Actions

1. **Genetic Counseling:** Consultation with a certified genetic counselor is recommended to discuss these findings and their implications.

2. **Clinical Correlation:** These results should be interpreted in the context of personal and family medical history.

3. **Enhanced Screening:** Consider more frequent or earlier cancer screening based on identified risks.

### Follow-up Considerations

- Discuss findings with your healthcare provider
- Consider family testing if appropriate
- Review screening recommendations with specialists
- Stay informed about advances in genetic medicine

"""
        else:
            recommendations += """### Standard Care Recommendations

1. **Routine Screening:** Continue with standard cancer screening protocols appropriate for your age and risk factors.

2. **Healthy Lifestyle:** Maintain healthy lifestyle choices that can reduce cancer risk regardless of genetic profile.

3. **Stay Informed:** Keep up with evolving genetic research and screening guidelines.

### Important Notes

- These results do not eliminate all cancer risk
- Environmental and lifestyle factors remain important
- Genetic testing technology continues to evolve

"""
        
        return recommendations
    
    def _build_technical_appendix(self, data: Dict[str, Any]) -> str:
        """Build technical appendix section."""
        
        appendix = "## Technical Appendix\n\n"
        
        # Analysis methods
        appendix += "### Analysis Methods\n\n"
        appendix += "This report was generated using the GeneKnow genomic analysis pipeline, which includes:\n\n"
        appendix += "- **Quality Control:** Variant filtering based on quality scores, read depth, and allele frequency\n"
        appendix += "- **Population Analysis:** Comparison with reference populations and allele frequencies\n"
        appendix += "- **TCGA Integration:** Matching variants against The Cancer Genome Atlas database\n"
        appendix += "- **CADD Scoring:** Combined Annotation Dependent Depletion pathogenicity prediction\n"
        appendix += "- **Polygenic Risk Scoring:** Integration of multiple genetic variants for risk assessment\n"
        appendix += "- **Machine Learning:** Advanced risk models trained on cancer genomics data\n\n"
        
        # Quality metrics
        summary = data.get("summary", {})
        appendix += "### Quality Metrics\n\n"
        appendix += f"- **Total Variants Identified:** {summary.get('total_variants_found', 0)}\n"
        appendix += f"- **Variants Passing QC:** {summary.get('variants_passed_qc', 0)}\n"
        
        if summary.get('total_variants_found', 0) > 0:
            pass_rate = (summary.get('variants_passed_qc', 0) / summary['total_variants_found']) * 100
            appendix += f"- **QC Pass Rate:** {pass_rate:.1f}%\n"
        
        appendix += "\n"
        
        # Database information
        tcga_summary = data.get("tcga_summary", {})
        if tcga_summary:
            appendix += "### Database Coverage\n\n"
            cancer_types = tcga_summary.get("cancer_types_analyzed", [])
            if cancer_types:
                appendix += f"- **TCGA Cancer Types:** {len(cancer_types)} ({', '.join(cancer_types)})\n"
            
            variants_with_data = tcga_summary.get("variants_with_tcga_data", 0)
            appendix += f"- **Variants with TCGA Data:** {variants_with_data}\n\n"
        
        # Limitations
        appendix += "### Limitations\n\n"
        appendix += "- This analysis is based on current scientific knowledge and databases\n"
        appendix += "- Risk estimates are population-based and may not apply to all individuals\n"
        appendix += "- Not all genetic variants associated with cancer risk are included\n"
        appendix += "- Results should be interpreted by qualified healthcare professionals\n"
        appendix += "- This analysis is for research purposes and not for clinical diagnosis\n\n"
        
        return appendix
    
    def _build_glossary(self, content: str) -> str:
        """Build glossary based on terms found in content."""
        
        found_terms = self.glossary.find_terms_in_text(content)
        
        if not found_terms:
            return ""
        
        glossary_section = "## Glossary\n\n"
        
        # Sort terms alphabetically
        sorted_terms = sorted(found_terms)
        
        for term in sorted_terms:
            definition = self.glossary.get_definition(term)
            if definition:
                glossary_section += f"**{term.title()}:** {definition}\n\n"
        
        return glossary_section
    
    def _build_basic_glossary(self) -> str:
        """Build basic glossary with common terms."""
        
        common_terms = [
            "allele frequency",
            "CADD score", 
            "clinical significance",
            "pathogenic",
            "polygenic risk score",
            "TCGA",
            "variant of uncertain significance"
        ]
        
        glossary_section = "## Glossary\n\n"
        
        for term in common_terms:
            definition = self.glossary.get_definition(term)
            if definition:
                glossary_section += f"**{term.title()}:** {definition}\n\n"
        
        return glossary_section
    
    def _build_footer(self, data: Dict[str, Any], report_id: str) -> str:
        """Build report footer."""
        
        footer = """---

## Important Disclaimers

- This report is for research and educational purposes only
- Results should not be used for clinical decision-making without consultation with qualified healthcare providers
- Genetic risk assessment is based on current scientific understanding and may change as research evolves
- Not all genetic variants associated with disease risk are included in this analysis
- Environmental and lifestyle factors also contribute significantly to cancer risk

## Contact Information

For questions about this report or genetic counseling resources, please consult with your healthcare provider or a certified genetic counselor.

"""
        
        footer += f"*Report generated by GeneKnow Pipeline - Report ID: {report_id}*\n"
        
        return footer
    
    def _convert_to_html(self, markdown_content: str, report_id: str) -> Optional[str]:
        """Convert markdown to HTML (placeholder implementation)."""
        try:
            # Basic markdown to HTML conversion
            # In a real implementation, you might use a library like markdown or mistune
            html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Genomic Risk Assessment Report - {report_id}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
        h1, h2, h3 {{ color: #333; }}
        .highlight {{ background-color: #fff3cd; padding: 10px; border-left: 4px solid #ffc107; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
<pre>{markdown_content}</pre>
</body>
</html>"""
            
            html_path = os.path.join(self.output_dir, f"{report_id}.html")
            with open(html_path, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            logger.info(f"Generated HTML report: {html_path}")
            return html_path
            
        except Exception as e:
            logger.error(f"Failed to convert to HTML: {e}")
            return None
    
    def _convert_to_pdf(self, markdown_content: str, report_id: str) -> Optional[str]:
        """Convert markdown to PDF (placeholder implementation)."""
        # This would require a library like weasyprint, reportlab, or calling pandoc
        logger.info("PDF conversion not implemented in this version")
        return None
    
    def _convert_to_txt(self, markdown_content: str, report_id: str) -> Optional[str]:
        """Convert markdown to plain text."""
        try:
            # Simple markdown stripping
            txt_content = re.sub(r'[#*`_\[\]()]', '', markdown_content)
            txt_content = re.sub(r'\n\n+', '\n\n', txt_content)
            
            txt_path = os.path.join(self.output_dir, f"{report_id}.txt")
            with open(txt_path, 'w', encoding='utf-8') as f:
                f.write(txt_content)
            
            logger.info(f"Generated text report: {txt_path}")
            return txt_path
            
        except Exception as e:
            logger.error(f"Failed to convert to text: {e}")
            return None 