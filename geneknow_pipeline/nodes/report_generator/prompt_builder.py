"""
Prompt builder for assembling LLM prompts from structured JSON data.
Focuses on high-risk findings (>5%) and clinician-appropriate language.
Each section is generated separately for maximum consistency.
"""

import logging
from typing import Dict, Any
from .config import ReportStyle

logger = logging.getLogger(__name__)


class PromptBuilder:
    """Builds structured prompts for LLM generation from genomic data."""

    def __init__(self, style: ReportStyle = ReportStyle.CLINICIAN):
        self.style = style
        self.risk_threshold = 5.0  # Only include >5% risk findings

    def build_summary_section_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating the Summary section only."""

        # Extract key information
        report_metadata = data.get("report_metadata", {})
        patient_data = data.get("patient_data", {})
        summary = data.get("summary", {})
        risk_assessment = data.get("risk_assessment", {})

        # Get high-risk findings
        high_risk_count = 0
        top_risks = []
        if risk_assessment and "scores" in risk_assessment:
            for cancer_type, score in risk_assessment["scores"].items():
                if score > self.risk_threshold:
                    high_risk_count += 1
                    top_risks.append(f"{cancer_type.title()}: {score:.1f}%")

        # Sort and limit to top 3
        top_risks = sorted(
            top_risks,
            key=lambda x: float(x.split(": ")[1].replace("%", "")),
            reverse=True,
        )[:3]

        prompt = self._get_style_prefix()
        prompt += """
Generate ONLY the Summary section content for a genomic risk assessment report.

INPUT DATA:
- File type: {patient_data.get('file_type', 'Unknown')}
- Total variants found: {summary.get('total_variants_found', 0)}
- Variants passed QC: {summary.get('variants_passed_qc', 0)}
- QC pass rate: {(summary.get('variants_passed_qc', 0) / max(summary.get('total_variants_found', 1), 1) * 100):.1f}%
- Processing time: {report_metadata.get('processing_time_seconds', 0):.1f} seconds
- High-risk findings: {high_risk_count}
- Top cancer risks: {', '.join(top_risks) if top_risks else 'All within baseline levels (<5%)'}

TASK:
Write exactly 1-2 paragraphs for the Summary section that:
1. First paragraph: Briefly describes the analysis performed (variants analyzed, QC results)
2. Second paragraph: States the overall risk profile and most significant findings
3. Uses professional medical language appropriate for clinicians
4. Is concise but informative
5. For high-risk cases: Mentions the elevated risks found
6. For low-risk cases: Emphasizes reassuring baseline results

OUTPUT FORMAT:
Return ONLY the paragraph text. Do not include section headers, bullets, or formatting.
Start directly with the content.

EXAMPLE OUTPUT (high-risk case):
"This genomic risk assessment analyzed [X] genetic variants from [file type] data, with [Y] variants passing quality control filters (Z% pass rate). The analysis was completed in [time] seconds and identified [count] high-risk findings requiring clinical attention.

The individual's genetic profile indicates significantly elevated risk for [top cancers], with the highest risk being [specific cancer] at [percentage]. These findings warrant further clinical evaluation and genetic counseling to discuss preventive measures and screening protocols."

EXAMPLE OUTPUT (low-risk case):
"This genomic risk assessment analyzed [X] genetic variants from [file type] data, with [Y] variants passing quality control filters (Z% pass rate). The analysis was completed in [time] seconds with no high-risk findings identified.

The individual's genetic profile shows cancer risk levels within normal baseline ranges for all analyzed cancer types (<5%). This reassuring result indicates that no known high-risk genetic variants were detected in major cancer susceptibility genes."
"""

        return prompt

    def build_key_variants_section_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating the Key Variants section only."""

        variant_details = data.get("variant_details", [])
        risk_assessment = data.get("risk_assessment", {})

        # Filter for high-impact variants
        key_variants = []
        for variant in variant_details[:10]:  # Limit to top 10
            # Check if variant is associated with high-risk cancer
            is_high_risk = False
            gene = variant.get("gene", "")

            if risk_assessment and "risk_genes" in risk_assessment:
                for cancer_type, score in risk_assessment.get("scores", {}).items():
                    if score > self.risk_threshold:
                        cancer_genes = risk_assessment["risk_genes"].get(
                            cancer_type, []
                        )
                        if gene in cancer_genes:
                            is_high_risk = True
                            break

            if (
                is_high_risk
                or variant.get("quality_metrics", {}).get("quality", 0) > 30
            ):
                key_variants.append(variant)

        prompt = self._get_style_prefix()

        if not key_variants:
            prompt += """
Generate ONLY the Key Variants section content for a case with no significant variants.

TASK:
Write exactly one sentence stating that no key pathogenic variants were identified.

OUTPUT FORMAT:
Return ONLY the sentence. Do not include section headers or formatting.

REQUIRED OUTPUT:
"No key pathogenic variants associated with elevated cancer risk were identified in this analysis."
"""
        else:
            prompt += """
Generate ONLY the Key Variants section content for a genomic risk assessment report.

KEY VARIANTS TO ANALYZE (top {min(len(key_variants), 5)}):
"""

            for i, variant in enumerate(key_variants[:5], 1):  # Top 5 variants
                quality = variant.get("quality_metrics", {})
                tcga_matches = variant.get("tcga_matches", {})
                cadd_scores = variant.get("cadd_scores", {})

                prompt += """
Variant {i}:
- Gene: {variant.get('gene', 'Unknown')}
- Variant ID: {variant.get('variant', 'Unknown')}
- Consequence: {variant.get('consequence', 'Unknown')}
- HGVS (coding): {variant.get('hgvs_c', 'N/A')}
- HGVS (protein): {variant.get('hgvs_p', 'N/A')}
- Quality score: {quality.get('quality', 0)}
- Read depth: {quality.get('depth', 0)}
- Allele frequency: {quality.get('allele_freq', 0):.3f}
"""

                if cadd_scores:
                    prompt += f"- CADD PHRED score: {cadd_scores.get('phred', 'N/A')}\n"

                if tcga_matches:
                    prompt += "- TCGA matches:\n"
                    for cancer_type, match_data in tcga_matches.items():
                        if match_data:
                            prompt += f"  * {cancer_type}: {match_data.get('frequency_percent', 0):.1f}% frequency\n"

                prompt += "\n"

            prompt += """
TASK:
Write a numbered list (1., 2., etc.) of the key variants with clinical interpretation.
For each variant, write 1-2 sentences that:
1. Identify the gene and mutation type (e.g., "BRCA1 frameshift mutation")
2. Explain the clinical significance and cancer risk association
3. Mention quality metrics if confidence is high (quality >90) or low (<50)
4. Use proper genetic nomenclature
5. Avoid overstating conclusions

OUTPUT FORMAT:
Return ONLY the numbered list. Do not include section headers or additional formatting.
Start directly with "1. **[GENE]**: [description]"

EXAMPLE OUTPUT:
"1. **BRCA1**: A frameshift mutation (c.5266dupC) resulting in a premature stop codon (p.Gln1756Profs*74). This pathogenic variant significantly increases the risk for breast and ovarian cancer. The quality score of 99 and read depth of 45 indicate high confidence in this finding.

2. **TP53**: A missense mutation (c.743G>A) causing an arginine to glutamine substitution at position 248 (p.Arg248Gln). This variant may contribute to increased cancer risk across multiple cancer types. The quality score of 85 indicates moderate confidence in this finding."
"""

        return prompt

    def build_risk_summary_section_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating the Risk Summary section only."""

        risk_assessment = data.get("risk_assessment", {})
        if not risk_assessment or "scores" not in risk_assessment:
            scores = {}
        else:
            scores = risk_assessment.get("scores", {})

        # Sort scores by descending risk
        sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)

        prompt = self._get_style_prefix()
        prompt += """
Generate ONLY the Risk Summary section content for a genomic risk assessment report.

CANCER RISK SCORES:
"""

        if sorted_scores:
            for cancer_type, score in sorted_scores:
                prompt += f"- {cancer_type.title()}: {score:.1f}%\n"
        else:
            prompt += "- No risk scores available\n"

        # Count high-risk findings
        high_risk_count = sum(
            1 for _, score in sorted_scores if score > self.risk_threshold
        )

        prompt += """
High-risk findings (>5%): {high_risk_count}

TASK:
Generate the Risk Summary section with:
1. A markdown table with columns "Cancer Type" and "Risk (%)"
2. All cancer types sorted by descending risk percentage
3. One sentence after the table explaining the significance

OUTPUT FORMAT:
Return the markdown table followed by exactly one explanatory sentence.
Do not include section headers or additional formatting.

REQUIRED TABLE FORMAT:
| Cancer Type | Risk (%) |
|-------------|----------|
| [Type]      | [X.X]    |
| [Type]      | [X.X]    |

REQUIRED SENTENCE (choose based on findings):
- If high-risk findings: "Risks above 5% are considered elevated and warrant further clinical attention."
- If no high-risk findings: "All risks are within baseline population levels (<5%), indicating no elevated genetic predisposition was detected."
"""

        return prompt

    def build_clinical_interpretation_section_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating the Clinical Interpretation section only."""

        summary = data.get("summary", {})
        risk_assessment = data.get("risk_assessment", {})
        tcga_summary = data.get("tcga_summary", {})
        cadd_summary = data.get("cadd_summary", {})

        # Count high-risk findings
        high_risk_count = 0
        if risk_assessment and "scores" in risk_assessment:
            high_risk_count = sum(
                1
                for score in risk_assessment["scores"].values()
                if score > self.risk_threshold
            )

        prompt = self._get_style_prefix()
        prompt += """
Generate ONLY the Clinical Interpretation section content for a genomic risk assessment report.

ANALYSIS SUMMARY:
- Total variants analyzed: {summary.get('total_variants_found', 0)}
- Quality-filtered variants: {summary.get('variants_passed_qc', 0)}
- High-risk cancer findings: {high_risk_count}
- TCGA database coverage: {len(tcga_summary.get('cancer_types_analyzed', []))} cancer types
- TCGA variants matched: {tcga_summary.get('variants_with_tcga_data', 0)}
"""

        if cadd_summary and cadd_summary.get("enabled", False):
            prompt += f"- CADD scoring: {cadd_summary.get('variants_scored', 0)} variants scored\n"

        prompt += """
TASK:
Write exactly 1-2 paragraphs for the Clinical Interpretation section that:
1. First paragraph: Explains the methodology used (TCGA matching, CADD scoring, ML models, etc.)
2. Second paragraph: Discusses reliability, limitations, and clinical context
3. Always includes the standard caveat about correlating with family history and clinical presentation
4. Uses professional medical language
5. Maintains balanced, evidence-based tone
6. Emphasizes research nature of analysis

OUTPUT FORMAT:
Return ONLY the paragraph text. Do not include section headers or formatting.
Start directly with the content.

REQUIRED ELEMENTS TO INCLUDE:
- Mention multiple computational approaches used
- Reference TCGA tumor database comparison
- Discuss pathogenicity prediction algorithms
- Note polygenic risk scoring and machine learning
- Include reliability/limitations discussion
- End with caveat about clinical correlation

REQUIRED ENDING SENTENCE:
"These findings should be correlated with family history, lifestyle factors, and clinical presentation for comprehensive risk assessment."
"""

        return prompt

    def build_recommendations_section_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating the Recommendations section only."""

        risk_assessment = data.get("risk_assessment", {})
        scores = risk_assessment.get("scores", {}) if risk_assessment else {}

        # Get high-risk cancers
        high_risk_cancers = []
        if scores:
            for cancer_type, score in scores.items():
                if score > self.risk_threshold:
                    high_risk_cancers.append((cancer_type, score))

        # Sort by risk score
        high_risk_cancers.sort(key=lambda x: x[1], reverse=True)

        prompt = self._get_style_prefix()
        prompt += """
Generate ONLY the Recommendations section content for a genomic risk assessment report.

RISK PROFILE:
- High-risk cancers found: {len(high_risk_cancers)}
"""

        if high_risk_cancers:
            prompt += "- Elevated risks:\n"
            for cancer_type, score in high_risk_cancers[:3]:  # Top 3
                prompt += f"  * {cancer_type.title()}: {score:.1f}%\n"
        else:
            prompt += "- All cancer risks within baseline levels\n"

        prompt += """
TASK:
Generate a bulleted list of 3-5 actionable clinical recommendations.
Tailor recommendations to the risk profile:

FOR HIGH-RISK CASES:
- Include genetic counseling recommendation
- Suggest specific screening protocols for elevated cancer types
- Mention consideration of preventive measures if appropriate
- Include standard follow-up advice

FOR LOW-RISK CASES:
- Recommend standard age-appropriate screening
- Mention continued adherence to preventive care
- Include lifestyle recommendations
- Still suggest genetic counseling if family history warrants

OUTPUT FORMAT:
Return ONLY the bulleted list using "- " format.
Do not include section headers or additional formatting.
Always end with: "- Consult a qualified healthcare provider for personalized guidance."

EXAMPLE OUTPUT (high-risk):
"- Genetic counseling is strongly recommended to discuss these findings and develop a personalized risk management plan
- Enhanced breast cancer screening including annual mammography and consideration of MRI screening
- Ovarian cancer surveillance with transvaginal ultrasound and CA-125 testing
- Consider prophylactic measures in consultation with oncology specialists
- Consult a qualified healthcare provider for personalized guidance."

EXAMPLE OUTPUT (low-risk):
"- Adhere to standard age-appropriate cancer screening protocols
- Maintain healthy lifestyle including regular exercise and balanced diet
- Continue routine preventive care and health maintenance
- Consider genetic counseling if strong family history of cancer emerges
- Consult a qualified healthcare provider for personalized guidance."
"""

        return prompt

    def _get_style_prefix(self) -> str:
        """Get style-specific prefix for prompts."""

        if self.style == ReportStyle.CLINICIAN:
            return """You are a clinical geneticist writing a professional genomic risk assessment report for healthcare providers. Use appropriate medical terminology, maintain scientific accuracy, and provide clinically relevant insights. """

        elif self.style == ReportStyle.TECHNICAL:
            return """You are a bioinformatician writing a detailed technical report for researchers and clinical scientists. Use precise scientific terminology, include methodological details, and provide comprehensive technical analysis. """

        elif self.style == ReportStyle.PATIENT:
            return """You are a genetic counselor writing a report for patients and their families. Use clear, accessible language while maintaining accuracy. Explain technical terms and provide context for understanding genetic risk. """

        return ""

    def build_full_report_prompt(self, data: Dict[str, Any]) -> str:
        """Build a comprehensive prompt for generating the entire report (legacy method)."""

        # This method is kept for backward compatibility but should not be used
        # The new approach generates each section separately
        logger.warning(
            "build_full_report_prompt is deprecated. Use individual section methods instead."
        )

        return self.build_summary_section_prompt(data)
