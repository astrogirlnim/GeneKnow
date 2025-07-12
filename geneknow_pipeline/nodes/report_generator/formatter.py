"""
Report formatter for generating Markdown output with glossary support.
Includes template-based fallback for non-LLM mode.
"""

import logging
import os
import re
from datetime import datetime
from typing import Dict, Any, Optional, Set
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
            "oncogene": "A gene that, when mutated or overexpressed, can contribute to cancer development.",
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

    def format_report(
        self, data: Dict[str, Any], llm_content: Optional[str] = None, report_id: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Format a complete report with LLM content or fallback templates.

        Args:
            data: Structured JSON data from pipeline
            llm_content: Generated content from LLM (if available)
            report_id: Unique identifier for the report

        Returns:
            Dictionary with report content (not file paths for HIPAA compliance)
        """
        if report_id is None:
            report_id = f"report_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

        # Generate markdown content in-memory
        if llm_content:
            markdown_content = self._format_llm_report(data, llm_content, report_id)
        else:
            markdown_content = self._format_template_report(data, report_id)

        logger.info(f"Generated markdown report in-memory (ID: {report_id})")

        # Return content directly instead of file paths for HIPAA compliance
        result_content = {"markdown": markdown_content}

        # Generate additional formats if requested (also in-memory)
        if "html" in self.config.output_formats:
            html_content = self._convert_to_html_content(markdown_content, report_id)
            if html_content:
                result_content["html"] = html_content

        if "pdf" in self.config.output_formats:
            # Note: PDF generation would need to be implemented in-memory
            # For now, we'll skip PDF to maintain HIPAA compliance
            logger.info("PDF generation skipped for HIPAA compliance")

        if "txt" in self.config.output_formats:
            txt_content = self._convert_to_txt_content(markdown_content, report_id)
            if txt_content:
                result_content["txt"] = txt_content

        return result_content

    def _format_llm_report(self, data: Dict[str, Any], llm_content: str, report_id: str) -> str:
        """Format report using LLM-generated content."""

        # Extract metadata
        data.get("report_metadata", {})

        # Build header
        markdown = self._build_header(data, report_id, llm_enhanced=True)

        # Add LLM-generated content under Clinical Assessment
        markdown += "\n## Clinical Assessment\n\n"
        markdown += llm_content

        # Always add standardized technical appendix
        markdown += "\n\n"
        markdown += self._build_standardized_technical_appendix(data)

        # Always add standardized glossary
        markdown += "\n\n"
        markdown += self._build_standardized_glossary(llm_content)

        # Always add standardized footer
        markdown += "\n\n"
        markdown += self._build_standardized_footer(data, report_id)

        return markdown

    def _format_template_report(self, data: Dict[str, Any], report_id: str) -> str:
        """Format report using template-based fallback."""

        # Build header with fallback indicator
        markdown = self._build_header(data, report_id, llm_enhanced=False)

        # Add fallback content warning
        markdown += f"\n> **{self.config.dev_mode_indicator}**: {self.config.fallback_mode_indicator}\n\n"

        # Add Clinical Assessment section with template content
        markdown += "\n## Clinical Assessment\n\n"

        # Generate template-based sections using the same structure as LLM reports
        markdown += self._build_template_summary(data)
        markdown += "\n\n"
        markdown += self._build_template_key_variants(data)
        markdown += "\n\n"
        markdown += self._build_template_risk_summary(data)
        markdown += "\n\n"
        markdown += self._build_template_clinical_interpretation(data)
        markdown += "\n\n"
        markdown += self._build_template_recommendations(data)

        # Always add standardized technical appendix
        markdown += "\n\n"
        markdown += self._build_standardized_technical_appendix(data)

        # Always add standardized glossary
        markdown += "\n\n"
        markdown += self._build_standardized_glossary("")

        # Always add standardized footer
        markdown += "\n\n"
        markdown += self._build_standardized_footer(data, report_id)

        return markdown

    def _build_header(self, data: Dict[str, Any], report_id: str, llm_enhanced: bool = False) -> str:
        """Build report header section."""

        header = """# Genomic Risk Assessment Report

**Generated:** {datetime.now().strftime('%B %d, %Y at %H:%M')}
"""

        if llm_enhanced:
            header += "**Analysis Type:** AI-Enhanced Clinical Report  \n"
        else:
            header += "**Analysis Type:** Template-Based Report  \n"

        header += "\n---\n"

        return header

    def _build_template_summary(self, data: Dict[str, Any]) -> str:
        """Build template-based summary section."""
        summary = data.get("summary", {})
        risk_assessment = data.get("risk_assessment", {})

        total_variants = summary.get("total_variants_found", 0)
        qc_variants = summary.get("variants_passed_qc", 0)
        qc_rate = (qc_variants / max(total_variants, 1)) * 100

        high_risk_count = 0
        if risk_assessment and "scores" in risk_assessment:
            high_risk_count = sum(1 for score in risk_assessment["scores"].values() if score > 5.0)

        content = "**Summary**\n\n"
        content += f"This genomic risk assessment analyzed {total_variants} genetic variants, with {qc_variants} variants passing quality control filters ({qc_rate:.1f}% pass rate). "

        if high_risk_count > 0:
            content += f"The analysis identified {high_risk_count} high-risk findings requiring clinical attention.\n\n"
            content += "The individual's genetic profile indicates elevated risk for certain cancer types that warrant further clinical evaluation and genetic counseling."
        else:
            content += "No high-risk findings were identified.\n\n"
            content += "The individual's genetic profile shows cancer risk levels within normal baseline ranges for all analyzed cancer types (<5%)."

        return content

    def _build_template_key_variants(self, data: Dict[str, Any]) -> str:
        """Build template-based key variants section."""
        variant_details = data.get("variant_details", [])

        content = "**Key Variants**\n\n"

        if not variant_details:
            content += (
                "No key pathogenic variants associated with elevated cancer risk were identified in this analysis."
            )
        else:
            # Show top 3 variants in simple format
            for i, variant in enumerate(variant_details[:3], 1):
                gene = variant.get("gene", "Unknown")
                consequence = variant.get("consequence", "Unknown")
                quality = variant.get("quality_metrics", {}).get("quality", 0)

                content += f"{i}. **{gene}**: {consequence} variant"
                if quality > 0:
                    content += f" (quality score: {quality})"
                content += "\n\n"

        return content

    def _build_template_risk_summary(self, data: Dict[str, Any]) -> str:
        """Build template-based risk summary section."""
        risk_assessment = data.get("risk_assessment", {})
        scores = risk_assessment.get("scores", {}) if risk_assessment else {}

        content = "**Risk Summary**\n\n"

        if not scores:
            content += "| Cancer Type | Risk (%) |\n|-------------|----------|\n| No data | N/A |\n\n"
            content += "Risk assessment data not available for this analysis."
        else:
            # Create table
            content += "| Cancer Type | Risk (%) |\n|-------------|----------|\n"
            for cancer_type, score in sorted(scores.items(), key=lambda x: x[1], reverse=True):
                content += f"| {cancer_type.title()} | {score:.1f} |\n"

            high_risk_count = sum(1 for score in scores.values() if score > 5.0)
            content += "\n"
            if high_risk_count > 0:
                content += "Risks above 5% are considered elevated and warrant further clinical attention."
            else:
                content += "All risks are within baseline population levels (<5%), indicating no elevated genetic predisposition was detected."

        return content

    def _build_template_clinical_interpretation(self, data: Dict[str, Any]) -> str:
        """Build template-based clinical interpretation section."""
        tcga_summary = data.get("tcga_summary", {})
        cadd_summary = data.get("cadd_summary", {})

        content = "**Clinical Interpretation**\n\n"
        content += "This genomic risk assessment utilized multiple computational approaches including population frequency analysis, "
        content += (
            f"TCGA tumor database comparison ({tcga_summary.get('variants_with_tcga_data', 0)} variants matched), "
        )
        content += "pathogenicity prediction algorithms, polygenic risk scoring, and machine learning risk models. "

        if cadd_summary and cadd_summary.get("enabled", False):
            content += f"CADD scoring was performed on {cadd_summary.get('variants_scored', 0)} variants. "

        content += "\n\nThe analysis provides research-grade insights that should be interpreted within the context of current scientific understanding and clinical guidelines. "
        content += "These findings should be correlated with family history, lifestyle factors, and clinical presentation for comprehensive risk assessment."

        return content

    def _build_template_recommendations(self, data: Dict[str, Any]) -> str:
        """Build template-based recommendations section."""
        risk_assessment = data.get("risk_assessment", {})
        scores = risk_assessment.get("scores", {}) if risk_assessment else {}
        high_risk_count = sum(1 for score in scores.values() if score > 5.0) if scores else 0

        content = "**Recommendations**\n\n"

        if high_risk_count > 0:
            content += "- Genetic counseling is recommended to discuss these findings and develop a personalized risk management plan\n"
            content += "- Enhanced screening protocols may be appropriate for elevated cancer risks\n"
            content += "- Consider consultation with oncology specialists for high-risk findings\n"
            content += "- Maintain regular follow-up care and health monitoring\n"
        else:
            content += "- Adhere to standard age-appropriate cancer screening protocols\n"
            content += "- Maintain healthy lifestyle including regular exercise and balanced diet\n"
            content += "- Continue routine preventive care and health maintenance\n"
            content += "- Consider genetic counseling if strong family history of cancer emerges\n"

        content += "- Consult a qualified healthcare provider for personalized guidance."

        return content

    def _build_standardized_technical_appendix(self, data: Dict[str, Any]) -> str:
        """Build standardized technical appendix that's the same for every report."""
        summary = data.get("summary", {})
        tcga_summary = data.get("tcga_summary", {})

        total_variants = summary.get("total_variants_found", 0)
        qc_variants = summary.get("variants_passed_qc", 0)
        qc_rate = (qc_variants / max(total_variants, 1)) * 100 if total_variants > 0 else 0

        appendix = """## Technical Appendix

### Analysis Methods

This report was generated using the GeneKnow genomic analysis pipeline, which includes:

- **Quality Control:** Variant filtering based on quality scores, read depth, and allele frequency
- **Population Analysis:** Comparison with reference populations and allele frequencies
- **TCGA Integration:** Matching variants against The Cancer Genome Atlas database
- **CADD Scoring:** Combined Annotation Dependent Depletion pathogenicity prediction
- **Polygenic Risk Scoring:** Integration of multiple genetic variants for risk assessment
- **Machine Learning:** Advanced risk models trained on cancer genomics data

### Quality Metrics

- **Total Variants Identified:** {total_variants}
- **Variants Passing QC:** {qc_variants}
- **QC Pass Rate:** {qc_rate:.1f}%

### Database Coverage

- **TCGA Cancer Types:** 5 (breast, colon, lung, prostate, blood)
- **Variants with TCGA Data:** {tcga_summary.get('variants_with_tcga_data', 0)}

### Limitations

- This analysis is based on current scientific knowledge and databases
- Risk estimates are population-based and may not apply to all individuals
- Not all genetic variants associated with cancer risk are included
- Results should be interpreted by qualified healthcare professionals
- This analysis is for research purposes and not for clinical diagnosis"""

        return appendix

    def _build_standardized_glossary(self, content: str) -> str:
        """Build standardized glossary with terms relevant to the report content."""
        glossary = "## Glossary\n\n"

        # Always include these basic terms
        glossary += "**Allele Frequency:** The proportion of chromosomes in a population that carry a specific variant of a gene.\n\n"
        glossary += "**Pathogenic:** A genetic variant that is known to cause disease.\n\n"

        # Add additional terms based on content
        if "frameshift" in content.lower():
            glossary += "**Frameshift:** A genetic mutation that causes a shift in the reading frame of DNA, often resulting in a nonfunctional protein.\n\n"

        if "missense" in content.lower():
            glossary += "**Missense Variant:** A genetic change that results in a different amino acid being incorporated into the protein.\n\n"

        # Remove trailing newlines
        return glossary.rstrip()

    def _build_standardized_footer(self, data: Dict[str, Any], report_id: str) -> str:
        """Build standardized footer that's the same for every report."""
        footer = """---

## Important Disclaimers

- This report is for research and educational purposes only
- Results should not be used for clinical decision-making without consultation with qualified healthcare providers
- Genetic risk assessment is based on current scientific understanding and may change as research evolves
- Not all genetic variants associated with disease risk are included in this analysis
- Environmental and lifestyle factors also contribute significantly to cancer risk

## Contact Information

For questions about this report or genetic counseling resources, please consult with your healthcare provider or a certified genetic counselor.

*Report generated by GeneKnow Pipeline - Report ID: {report_id}*"""

        return footer

    def _convert_to_html_content(self, markdown_content: str, report_id: str) -> Optional[str]:
        """Convert markdown to HTML (placeholder implementation)."""
        try:
            # Basic markdown to HTML conversion
            # In a real implementation, you might use a library like markdown or mistune
            html_content = """<!DOCTYPE html>
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

            return html_content

        except Exception as e:
            logger.error(f"Failed to convert to HTML: {e}")
            return None

    def _convert_to_pdf(self, markdown_content: str, report_id: str) -> Optional[str]:
        """Convert markdown to PDF (placeholder implementation)."""
        # This would require a library like weasyprint, reportlab, or calling pandoc
        logger.info("PDF conversion not implemented in this version")
        return None

    def _convert_to_txt_content(self, markdown_content: str, report_id: str) -> Optional[str]:
        """Convert markdown to plain text."""
        try:
            # Simple markdown stripping
            txt_content = re.sub(r"[#*`_\[\]()]", "", markdown_content)
            txt_content = re.sub(r"\n\n+", "\n\n", txt_content)

            return txt_content

        except Exception as e:
            logger.error(f"Failed to convert to text: {e}")
            return None
