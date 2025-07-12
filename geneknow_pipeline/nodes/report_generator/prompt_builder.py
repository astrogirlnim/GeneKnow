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
        # Format the actual values
        file_type = patient_data.get('file_type', 'Unknown')
        
        # Handle both new and old summary format
        if 'total_variants_found' in summary:
            # New format with variant data
            total_variants = summary.get('total_variants_found', 0)
            qc_variants = summary.get('variants_passed_qc', 0)
        else:
            # Old format or missing data - try to get from file metadata
            file_metadata = patient_data.get('file_metadata', {})
            qc_stats = file_metadata.get('qc_stats', {})
            if qc_stats:
                total_variants = qc_stats.get('total_variants', 0)
                qc_variants = qc_stats.get('passed_qc', 0)
            else:
                total_variants = 0
                qc_variants = 0
        
        qc_pass_rate = (qc_variants / max(total_variants, 1)) * 100
        processing_time = report_metadata.get('processing_time_seconds', 0)
        top_risks_text = ', '.join(top_risks) if top_risks else 'All within baseline levels (<5%)'
        
        # Add debug logging
        print(f"\n=== PROMPT BUILDER DEBUG ===")
        print(f"summary data: {summary}")
        print(f"total_variants: {total_variants}")
        print(f"qc_variants: {qc_variants}")
        print(f"qc_pass_rate: {qc_pass_rate}")
        print(f"=== END PROMPT BUILDER DEBUG ===\n")

        prompt += f"""
Generate ONLY the Summary section content for a genomic risk assessment report.

INPUT DATA:
- File type: {file_type}
- Total variants found: {total_variants}
- Variants passed QC: {qc_variants}
- QC pass rate: {qc_pass_rate:.1f}%
- Processing time: {processing_time:.1f} seconds
- High-risk findings: {high_risk_count}
- Top cancer risks: {top_risks_text}

TASK:
Write exactly 2-3 detailed paragraphs for the Summary section that:
1. First paragraph: Describe the analysis performed including methodology overview (e.g., variant calling, QC filters, risk modeling), key statistics like total variants and QC pass rate.
2. Second paragraph: Detail the overall risk profile, highlighting top risks with brief explanations of contributing factors (e.g., specific genes or scores).
3. Third paragraph (optional): Discuss general implications and next steps.
4. Use professional medical language appropriate for healthcare providers.
5. Be informative and evidence-based, including brief references to analysis methods where relevant.
6. For high-risk cases: Elaborate on elevated risks and potential significance.
7. For low-risk cases: Emphasize reassuring aspects and standard recommendations.

OUTPUT FORMAT:
Return ONLY the paragraph text. ABSOLUTELY NO section headers, bullets, markdown, or extra formatting. Start directly with the content text.

GENERATE YOUR OWN TEXT based on the actual input data provided above. Do not use template language or placeholders.
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
            # Format the actual values
            total_variants = len(variant_details)
            high_risk_variants = len(key_variants)
            top_genes = list(set(v.get("gene", "") for v in key_variants if v.get("gene")))[:5]
            top_genes_text = ', '.join(top_genes) if top_genes else 'None identified'

            prompt += f"""
Generate ONLY the Key Variants section content for a genomic risk assessment report.

INPUT DATA:
- Total variants analyzed: {total_variants}
- High-risk variants: {high_risk_variants}
- Genes with elevated risk: {top_genes_text}

TASK:
Write a numbered list (1., 2., etc.) of the key variants with detailed interpretation.
1. Focus on variants that contribute most to cancer risk (pathogenic/likely pathogenic).
2. Explain the detailed significance, associated cancers, and evidence from databases like ClinVar/TCGA.
3. Include specific risk contributions and functional impact where known.
4. Use professional medical language.
5. Be specific about evidence and database sources.
6. Limit to the top 3-5 most significant variants if many exist.

OUTPUT FORMAT:
Return ONLY the numbered list. Do not include section headers or extra formatting.
Start directly with "1." for the first variant.

GENERATE YOUR OWN NUMBERED LIST based on the actual variant data provided. Each entry should describe a real variant from the analysis, not example variants.
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

        prompt += f"""
High-risk findings (>5%): {high_risk_count}

TASK:
Generate the Risk Summary section with:
1. A markdown table with columns "Cancer Type" and "Risk (%)"
2. All cancer types sorted by descending risk percentage
3. After the table, 1-2 detailed paragraphs explaining the scores, comparisons to population averages (assume 1-2% baseline), health implications, and factors influencing the risks.

OUTPUT FORMAT:
Return the markdown table followed by the explanatory paragraphs. ABSOLUTELY NO section headers. Start with the table.

REQUIRED TABLE FORMAT:
| Cancer Type | Risk (%) |
|-------------|----------|
| [Type]      | [X.X]    |

EXPLANATORY TEXT:
Explain scores relative to baselines, highlight elevations, discuss potential genetic factors.

REQUIRED SENTENCE (choose based on findings):
- If high-risk findings: "Risks above 5% are considered elevated and warrant further medical attention."
- If no high-risk findings: "All risks are within baseline population levels (<5%), indicating no elevated genetic predisposition was detected."
"""

        return prompt

    def build_clinical_interpretation_section_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating the Interpretation section only."""

        # Extract relevant data
        summary = data.get("summary", {})
        tcga_summary = data.get("tcga_summary", {})
        cadd_summary = data.get("cadd_summary", {})
        risk_assessment = data.get("risk_assessment", {})
        
        # Check if high-risk findings exist
        high_risk_count = self._count_high_risk_findings(data)
        
        # Build section prompt
        prompt = self._get_style_prefix()
        
        prompt += f"""
Generate ONLY the Interpretation section content for a genomic risk assessment report.

INPUT DATA:
- Total variants analyzed: {summary.get('total_variants_found', 0)}
- High-risk findings: {high_risk_count}
- TCGA variants matched: {tcga_summary.get('variants_with_tcga_data', 0)}
- CADD scoring enabled: {cadd_summary.get('enabled', False)}
- Risk assessment available: {bool(risk_assessment)}

TASK:
Write exactly 2-3 paragraphs for the Interpretation section that:
1. First paragraph: Explains the methodology used (TCGA matching, CADD scoring, ML models, etc.)
2. Second paragraph: Discusses reliability, limitations, and context
3. Always includes the standard caveat about correlating with family history and presentation
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
- End with caveat about correlation

REQUIRED ENDING SENTENCE:
"These findings should be correlated with family history, lifestyle factors, and presentation for comprehensive risk assessment."
"""

        return prompt

    def build_recommendations_section_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating the Recommendations section only."""

        risk_assessment = data.get("risk_assessment", {})
        scores = risk_assessment.get("scores", {}) if risk_assessment else {}
        summary = data.get("summary", {})

        # Get high-risk cancers
        high_risk_cancers = []
        if scores:
            for cancer_type, score in scores.items():
                if score > self.risk_threshold:
                    high_risk_cancers.append((cancer_type, score))

        # Sort by risk score
        high_risk_cancers.sort(key=lambda x: x[1], reverse=True)
        
        # Format the actual values
        high_risk_count = len(high_risk_cancers)
        top_risks = [f"{cancer}: {score:.1f}%" for cancer, score in high_risk_cancers[:3]]
        top_risks_text = ', '.join(top_risks) if top_risks else 'All within baseline levels'
        total_variants = summary.get('total_variants_found', 0)

        prompt = self._get_style_prefix()
        prompt += f"""
Generate ONLY the Recommendations section content for a genomic risk assessment report.

INPUT DATA:
- High-risk findings: {high_risk_count}
- Top cancer risks: {top_risks_text}
- Total variants analyzed: {total_variants}

TASK:
Generate a bulleted list of 4-6 actionable recommendations, each with 1-2 explanatory sentences.
Focus on:
1. Immediate actions based on risk level (genetic counseling, enhanced screening)
2. Preventive measures and lifestyle modifications
3. Family considerations and cascade screening when appropriate
4. Follow-up care and monitoring
5. Risk management plan based on current guidelines."

OUTPUT FORMAT:
Return ONLY the bulleted list. Use "- " for bullets. Do not include section headers or extra formatting.

GENERATE YOUR OWN RECOMMENDATIONS based on the actual risk levels and findings in the data. Do not use example recommendations.
"""

        return prompt

    def _get_style_prefix(self) -> str:
        """Get style-specific prefix for prompts."""

        if self.style == ReportStyle.CLINICIAN:
            return """You are a geneticist writing a professional genomic risk assessment report for healthcare providers. Use appropriate medical terminology, maintain scientific accuracy, and provide relevant insights. """

        elif self.style == ReportStyle.TECHNICAL:
            return """You are a bioinformatician writing a detailed technical report for researchers and scientists. Use precise scientific terminology, include methodological details, and provide comprehensive technical analysis. """

        elif self.style == ReportStyle.PATIENT:
            return """You are a genetic counselor writing a report for patients and their families. Use clear, accessible language while maintaining accuracy. Explain technical terms and provide context for understanding genetic risk. """

        return ""

    def _count_high_risk_findings(self, data: Dict[str, Any]) -> int:
        """Count the number of high-risk findings in the data."""
        risk_assessment = data.get("risk_assessment", {})
        if not risk_assessment or "scores" not in risk_assessment:
            return 0
        
        scores = risk_assessment.get("scores", {})
        return sum(1 for score in scores.values() if score > self.risk_threshold)

    def build_full_report_prompt(self, data: Dict[str, Any]) -> str:
        """Build a comprehensive prompt for generating the entire report (legacy method)."""

        # This method is kept for backward compatibility but should not be used
        # The new approach generates each section separately
        logger.warning(
            "build_full_report_prompt is deprecated. Use individual section methods instead."
        )

        return self.build_summary_section_prompt(data)
