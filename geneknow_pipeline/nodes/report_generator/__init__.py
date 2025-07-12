"""
GeneKnow Report Generator Module

Transforms structured variant-level JSON output from the local genomics pipeline
into professional, readable clinical reports using local LLMs or fallback templates.
"""

from .report_generator import process
from .model_interface import ModelInterface
from .prompt_builder import PromptBuilder
from .formatter import ReportFormatter
from .config import ReportConfig

__all__ = ["process", "ModelInterface", "PromptBuilder", "ReportFormatter", "ReportConfig"]
