"""
Prompt builder for assembling LLM prompts from structured JSON data.
Focuses on high-risk findings (>5%) and clinician-appropriate language.
"""

import logging
from typing import Dict, Any, List, Optional
from .config import ReportStyle

logger = logging.getLogger(__name__)


class PromptBuilder:
    """Builds structured prompts for LLM generation from genomic data."""
    
    def __init__(self, style: ReportStyle = ReportStyle.CLINICIAN):
        self.style = style
        self.risk_threshold = 5.0  # Only include >5% risk findings
    
    def build_case_overview_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating case overview section."""
        
        # Extract key information
        report_metadata = data.get("report_metadata", {})
        patient_data = data.get("patient_data", {})
        summary = data.get("summary", {})
        risk_assessment = data.get("risk_assessment", {})
        
        # Get high-risk findings
        high_risk_findings = []
        if risk_assessment and "scores" in risk_assessment:
            for cancer_type, score in risk_assessment["scores"].items():
                if score > self.risk_threshold:
                    high_risk_findings.append({
                        "cancer_type": cancer_type,
                        "score": score
                    })
        
        # Sort by risk score
        high_risk_findings.sort(key=lambda x: x["score"], reverse=True)
        
        prompt = self._get_style_prefix()
        prompt += f"""
Generate a professional case overview for a genomic risk assessment report.

INPUT DATA:
- File type: {patient_data.get('file_type', 'Unknown')}
- Total variants found: {summary.get('total_variants_found', 0)}
- Variants passed QC: {summary.get('variants_passed_qc', 0)}
- Processing time: {report_metadata.get('processing_time_seconds', 0):.1f} seconds
- High-risk findings: {len(high_risk_findings)}

HIGH-RISK CANCER TYPES (>5% risk):
"""
        
        if high_risk_findings:
            for finding in high_risk_findings:
                prompt += f"- {finding['cancer_type'].title()}: {finding['score']:.1f}%\n"
        else:
            prompt += "- None detected (all cancer risks within baseline levels)\n"
        
        prompt += f"""
TASK:
Write a concise case overview (2-3 paragraphs) that:
1. Summarizes the genomic analysis performed
2. Highlights the most significant findings
3. Uses appropriate medical terminology
4. Maintains a professional, clinical tone
5. Focuses only on findings above 5% risk threshold

Do not include specific variant details or technical parameters in this overview.
"""
        
        return prompt
    
    def build_key_variants_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating key variants section."""
        
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
                        cancer_genes = risk_assessment["risk_genes"].get(cancer_type, [])
                        if gene in cancer_genes:
                            is_high_risk = True
                            break
            
            if is_high_risk or variant.get("quality_metrics", {}).get("quality", 0) > 30:
                key_variants.append(variant)
        
        if not key_variants:
            return ""
        
        prompt = self._get_style_prefix()
        prompt += f"""
Generate a professional analysis of key genetic variants identified in this genomic assessment.

KEY VARIANTS TO ANALYZE:
"""
        
        for i, variant in enumerate(key_variants[:5], 1):  # Top 5 variants
            quality = variant.get("quality_metrics", {})
            tcga_matches = variant.get("tcga_matches", {})
            cadd_scores = variant.get("cadd_scores", {})
            
            prompt += f"""
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
        
        prompt += f"""
TASK:
Write a clinical interpretation of these key variants that:
1. Explains the significance of each variant in medical terms
2. Discusses the quality and reliability of the findings
3. Relates variants to cancer risk when applicable
4. Uses proper genetic nomenclature
5. Maintains scientific accuracy and clinical relevance
6. Avoids overstating conclusions

Focus on variants that contribute to elevated cancer risk (>5%).
"""
        
        return prompt
    
    def build_risk_summary_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating risk summary section."""
        
        risk_assessment = data.get("risk_assessment", {})
        if not risk_assessment or "scores" in risk_assessment:
            return ""
        
        scores = risk_assessment.get("scores", {})
        risk_genes = risk_assessment.get("risk_genes", {})
        
        # Filter for high-risk findings
        high_risk_findings = []
        for cancer_type, score in scores.items():
            if score > self.risk_threshold:
                genes = risk_genes.get(cancer_type, [])
                high_risk_findings.append({
                    "cancer_type": cancer_type,
                    "score": score,
                    "genes": genes
                })
        
        if not high_risk_findings:
            return self._build_low_risk_prompt()
        
        # Sort by risk score
        high_risk_findings.sort(key=lambda x: x["score"], reverse=True)
        
        prompt = self._get_style_prefix()
        prompt += f"""
Generate a comprehensive cancer risk assessment summary based on genomic analysis.

HIGH-RISK FINDINGS (>5% baseline risk):
"""
        
        for finding in high_risk_findings:
            prompt += f"""
{finding['cancer_type'].title()} Cancer:
- Risk score: {finding['score']:.1f}%
- Associated genes: {', '.join(finding['genes']) if finding['genes'] else 'None specified'}
- Risk level: {'High' if finding['score'] >= 30 else 'Moderate' if finding['score'] >= 15 else 'Elevated'}
"""
        
        prompt += f"""
BASELINE RISK CANCERS (<5% risk):
"""
        baseline_cancers = [cancer for cancer, score in scores.items() if score <= self.risk_threshold]
        if baseline_cancers:
            prompt += f"- {', '.join([cancer.title() for cancer in baseline_cancers])}\n"
        else:
            prompt += "- None (all analyzed cancers show elevated risk)\n"
        
        prompt += f"""
TASK:
Write a clinical risk assessment summary that:
1. Clearly communicates the cancer risk levels identified
2. Explains what these risk percentages mean in clinical context
3. Discusses the genetic basis for elevated risks
4. Provides appropriate caveats about risk interpretation
5. Uses language appropriate for healthcare providers
6. Emphasizes the need for clinical correlation
7. Mentions that genetic counseling may be beneficial

Focus on actionable findings and avoid speculation.
"""
        
        return prompt
    
    def build_interpretation_prompt(self, data: Dict[str, Any]) -> str:
        """Build prompt for generating overall interpretation section."""
        
        summary = data.get("summary", {})
        risk_assessment = data.get("risk_assessment", {})
        tcga_summary = data.get("tcga_summary", {})
        cadd_summary = data.get("cadd_summary", {})
        
        # Count high-risk findings
        high_risk_count = 0
        if risk_assessment and "scores" in risk_assessment:
            high_risk_count = sum(1 for score in risk_assessment["scores"].values() if score > self.risk_threshold)
        
        prompt = self._get_style_prefix()
        prompt += f"""
Generate a comprehensive clinical interpretation of this genomic risk assessment.

ANALYSIS SUMMARY:
- Total variants analyzed: {summary.get('total_variants_found', 0)}
- Quality-filtered variants: {summary.get('variants_passed_qc', 0)}
- High-risk cancer findings: {high_risk_count}
- TCGA database coverage: {len(tcga_summary.get('cancer_types_analyzed', []))} cancer types
"""
        
        if cadd_summary and cadd_summary.get("enabled", False):
            prompt += f"- CADD scoring: {cadd_summary.get('variants_scored', 0)} variants scored\n"
        
        prompt += f"""
CLINICAL CONTEXT:
This analysis used multiple computational approaches including:
- Population frequency analysis
- TCGA tumor database comparison ({tcga_summary.get('variants_with_tcga_data', 0)} variants matched)
- Pathogenicity prediction algorithms
- Polygenic risk scoring
- Machine learning risk models

TASK:
Write a comprehensive clinical interpretation that:
1. Synthesizes all findings into a coherent assessment
2. Discusses the reliability and limitations of the analysis
3. Provides clinical recommendations based on findings
4. Explains the methodology in accessible terms
5. Addresses appropriate follow-up actions
6. Emphasizes the research nature of this analysis
7. Recommends genetic counseling for significant findings

Maintain a balanced, evidence-based tone throughout.
"""
        
        return prompt
    
    def _build_low_risk_prompt(self) -> str:
        """Build prompt for cases with no high-risk findings."""
        
        prompt = self._get_style_prefix()
        prompt += f"""
Generate a reassuring but informative risk summary for a genomic analysis with no elevated cancer risks.

FINDINGS:
- All analyzed cancer types show risk levels within normal baseline ranges (<5%)
- No pathogenic variants identified in major cancer susceptibility genes
- Genetic profile suggests average population risk

TASK:
Write a clinical summary that:
1. Clearly communicates the reassuring nature of the findings
2. Explains what "baseline risk" means
3. Emphasizes that this doesn't eliminate all cancer risk
4. Mentions limitations of current genetic testing
5. Recommends standard screening protocols
6. Uses reassuring but medically accurate language

Avoid overconfidence while providing appropriate reassurance.
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
        """Build a comprehensive prompt for generating the entire report."""
        
        # Combine all sections into one comprehensive prompt
        sections = [
            "Generate a comprehensive genomic risk assessment report with the following sections:",
            "",
            "1. CASE OVERVIEW",
            self.build_case_overview_prompt(data),
            "",
            "2. KEY VARIANTS", 
            self.build_key_variants_prompt(data),
            "",
            "3. RISK SUMMARY",
            self.build_risk_summary_prompt(data),
            "",
            "4. CLINICAL INTERPRETATION",
            self.build_interpretation_prompt(data),
            "",
            "Format the output as a structured report with clear section headers."
        ]
        
        return "\n".join(sections) 