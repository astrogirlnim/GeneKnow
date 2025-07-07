"""
GenePredict ML Engine for Genomic Risk Assessment

This package provides AI-powered genomic risk assessment capabilities
with privacy-first local processing.
"""

__version__ = "0.1.0"
__author__ = "GenePredict Team"
__email__ = "dev@genepredict.ai"

from .models import RiskAssessmentModel, BreastCancerRiskModel
from .processors import VCFProcessor, BAMProcessor, FASTQProcessor
from .utils import SecurityManager, ConfigManager

__all__ = [
    "RiskAssessmentModel",
    "BreastCancerRiskModel", 
    "VCFProcessor",
    "BAMProcessor",
    "FASTQProcessor",
    "SecurityManager",
    "ConfigManager",
] 