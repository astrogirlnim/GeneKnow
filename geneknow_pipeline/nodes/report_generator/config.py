"""
Configuration management for the Report Generator module.
"""

import os
import yaml
from dataclasses import dataclass
from typing import List, Optional, Dict, Any
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class LLMBackend(Enum):
    """Supported LLM backends."""
    OLLAMA = "ollama"
    HUGGINGFACE = "huggingface"
    NONE = "none"  # Fallback mode


class ReportStyle(Enum):
    """Report writing styles."""
    CLINICIAN = "clinician"
    TECHNICAL = "technical"
    PATIENT = "patient"


@dataclass
class ReportConfig:
    """Configuration for report generation."""
    
    # LLM Configuration
    backend: LLMBackend = LLMBackend.NONE
    model_name: Optional[str] = None
    temperature: float = 0.3
    max_tokens: int = 2000
    enable_streaming: bool = True
    
    # Report Configuration
    style: ReportStyle = ReportStyle.CLINICIAN
    output_formats: List[str] = None
    include_glossary: bool = True
    include_technical_appendix: bool = True
    risk_threshold: float = 5.0  # Only include >5% risk findings
    
    # Fallback Configuration
    fallback_mode_indicator: str = "Generated without LLM assistance"
    dev_mode_indicator: str = "NON-LLM MODE"
    
    def __post_init__(self):
        """Set default values after initialization."""
        if self.output_formats is None:
            self.output_formats = ["markdown"]


def load_config(config_path: Optional[str] = None) -> ReportConfig:
    """
    Load configuration from YAML file or use defaults.
    
    Args:
        config_path: Path to config file. If None, looks for config.yaml in project root.
        
    Returns:
        ReportConfig instance
    """
    if config_path is None:
        # Look for config.yaml in project root
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        config_path = os.path.join(project_root, "config.yaml")
    
    config_data = {}
    
    if os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                yaml_data = yaml.safe_load(f)
                config_data = yaml_data.get("report_generator", {})
                logger.info(f"Loaded report generator config from {config_path}")
        except Exception as e:
            logger.warning(f"Failed to load config from {config_path}: {e}")
    else:
        logger.info(f"No config file found at {config_path}, using defaults")
    
    # Convert string enums to enum instances
    if "backend" in config_data:
        try:
            config_data["backend"] = LLMBackend(config_data["backend"])
        except ValueError:
            logger.warning(f"Invalid backend '{config_data['backend']}', using default")
            config_data["backend"] = LLMBackend.NONE
    
    if "style" in config_data:
        try:
            config_data["style"] = ReportStyle(config_data["style"])
        except ValueError:
            logger.warning(f"Invalid style '{config_data['style']}', using default")
            config_data["style"] = ReportStyle.CLINICIAN
    
    return ReportConfig(**config_data)


def save_config(config: ReportConfig, config_path: Optional[str] = None) -> None:
    """
    Save configuration to YAML file.
    
    Args:
        config: ReportConfig instance to save
        config_path: Path to save config file. If None, saves to project root.
    """
    if config_path is None:
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        config_path = os.path.join(project_root, "config.yaml")
    
    # Load existing config or create new
    existing_config = {}
    if os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                existing_config = yaml.safe_load(f) or {}
        except Exception as e:
            logger.warning(f"Failed to load existing config: {e}")
    
    # Convert config to dict
    config_dict = {
        "backend": config.backend.value,
        "model_name": config.model_name,
        "temperature": config.temperature,
        "max_tokens": config.max_tokens,
        "enable_streaming": config.enable_streaming,
        "style": config.style.value,
        "output_formats": config.output_formats,
        "include_glossary": config.include_glossary,
        "include_technical_appendix": config.include_technical_appendix,
        "risk_threshold": config.risk_threshold,
        "fallback_mode_indicator": config.fallback_mode_indicator,
        "dev_mode_indicator": config.dev_mode_indicator
    }
    
    # Update existing config
    existing_config["report_generator"] = config_dict
    
    try:
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, 'w') as f:
            yaml.dump(existing_config, f, default_flow_style=False, indent=2)
        logger.info(f"Saved report generator config to {config_path}")
    except Exception as e:
        logger.error(f"Failed to save config to {config_path}: {e}")


# Example config.yaml structure for reference
EXAMPLE_CONFIG = """
report_generator:
  backend: "ollama"  # ollama, huggingface, or none
  model_name: "llama3"  # null for auto-detect
  temperature: 0.3
  max_tokens: 2000
  enable_streaming: true
  style: "clinician"  # clinician, technical, or patient
  output_formats: ["markdown", "pdf"]
  include_glossary: true
  include_technical_appendix: true
  risk_threshold: 5.0
  fallback_mode_indicator: "Generated without LLM assistance"
  dev_mode_indicator: "NON-LLM MODE"
""" 