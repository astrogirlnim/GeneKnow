"""
Model interface for LLM backends (Ollama, HuggingFace) with auto-detection and fallback.
"""

import logging
import requests
import subprocess
import json
from typing import Optional, Dict, Any, List, Generator
from .config import ReportConfig, LLMBackend

logger = logging.getLogger(__name__)


class OllamaBackend:
    """Ollama LLM backend interface."""
    
    def __init__(self, base_url: str = "http://localhost:11434"):
        self.base_url = base_url
        self.available_models = []
        self._check_availability()
    
    def _check_availability(self) -> bool:
        """Check if Ollama is running and get available models."""
        try:
            response = requests.get(f"{self.base_url}/api/tags", timeout=5)
            if response.status_code == 200:
                data = response.json()
                self.available_models = [model["name"] for model in data.get("models", [])]
                logger.info(f"Ollama available with models: {self.available_models}")
                return True
        except Exception as e:
            logger.debug(f"Ollama not available: {e}")
        return False
    
    def is_available(self) -> bool:
        """Check if Ollama backend is available."""
        return len(self.available_models) > 0
    
    def get_best_model(self, preferred: Optional[str] = None) -> Optional[str]:
        """Get the best available model for report generation."""
        if not self.available_models:
            return None
        
        # If preferred model is specified and available, use it
        if preferred and preferred in self.available_models:
            return preferred
        
        # Priority order for medical/scientific writing
        priority_models = [
            "llama3",
            "llama3:8b", 
            "llama3:70b",
            "mistral",
            "mistral:7b",
            "codellama",
            "phi3",
            "gemma"
        ]
        
        for model in priority_models:
            if model in self.available_models:
                return model
        
        # Fall back to first available model
        return self.available_models[0]
    
    def generate(self, prompt: str, model: str, temperature: float = 0.3, max_tokens: int = 2000) -> str:
        """Generate text using Ollama."""
        try:
            payload = {
                "model": model,
                "prompt": prompt,
                "stream": False,
                "options": {
                    "temperature": temperature,
                    "num_predict": max_tokens
                }
            }
            
            response = requests.post(
                f"{self.base_url}/api/generate",
                json=payload,
                timeout=120
            )
            
            if response.status_code == 200:
                data = response.json()
                return data.get("response", "")
            else:
                logger.error(f"Ollama generation failed: {response.status_code}")
                return ""
                
        except Exception as e:
            logger.error(f"Ollama generation error: {e}")
            return ""
    
    def stream_generate(self, prompt: str, model: str, temperature: float = 0.3, max_tokens: int = 2000) -> Generator[str, None, None]:
        """Generate text using Ollama with streaming."""
        try:
            payload = {
                "model": model,
                "prompt": prompt,
                "stream": True,
                "options": {
                    "temperature": temperature,
                    "num_predict": max_tokens
                }
            }
            
            response = requests.post(
                f"{self.base_url}/api/generate",
                json=payload,
                stream=True,
                timeout=120
            )
            
            if response.status_code == 200:
                for line in response.iter_lines():
                    if line:
                        try:
                            data = json.loads(line.decode('utf-8'))
                            if 'response' in data:
                                yield data['response']
                        except json.JSONDecodeError:
                            continue
            else:
                logger.error(f"Ollama streaming failed: {response.status_code}")
                
        except Exception as e:
            logger.error(f"Ollama streaming error: {e}")


class HuggingFaceBackend:
    """HuggingFace transformers backend interface."""
    
    def __init__(self):
        self.pipeline = None
        self.available_models = []
        self._check_availability()
    
    def _check_availability(self) -> bool:
        """Check if HuggingFace transformers is available."""
        try:
            import transformers
            import torch
            
            # Check for common medical/scientific models
            self.available_models = [
                "microsoft/DialoGPT-medium",
                "microsoft/phi-2",
                "google/flan-t5-base",
                "allenai/scibert_scivocab_uncased"
            ]
            
            logger.info(f"HuggingFace transformers available")
            return True
            
        except ImportError:
            logger.debug("HuggingFace transformers not available")
            return False
    
    def is_available(self) -> bool:
        """Check if HuggingFace backend is available."""
        return len(self.available_models) > 0
    
    def get_best_model(self, preferred: Optional[str] = None) -> Optional[str]:
        """Get the best available model for report generation."""
        if not self.available_models:
            return None
        
        if preferred and preferred in self.available_models:
            return preferred
        
        # Priority for medical/scientific writing
        priority_models = [
            "microsoft/phi-2",
            "google/flan-t5-base", 
            "microsoft/DialoGPT-medium"
        ]
        
        for model in priority_models:
            if model in self.available_models:
                return model
        
        return self.available_models[0]
    
    def _load_pipeline(self, model_name: str):
        """Load the HuggingFace pipeline for a model."""
        if self.pipeline is None:
            try:
                from transformers import pipeline
                self.pipeline = pipeline("text-generation", model=model_name, max_length=2000)
                logger.info(f"Loaded HuggingFace model: {model_name}")
            except Exception as e:
                logger.error(f"Failed to load HuggingFace model {model_name}: {e}")
                return False
        return True
    
    def generate(self, prompt: str, model: str, temperature: float = 0.3, max_tokens: int = 2000) -> str:
        """Generate text using HuggingFace."""
        try:
            if not self._load_pipeline(model):
                return ""
            
            results = self.pipeline(
                prompt,
                max_length=len(prompt.split()) + max_tokens,
                temperature=temperature,
                do_sample=True,
                pad_token_id=self.pipeline.tokenizer.eos_token_id
            )
            
            if results and len(results) > 0:
                generated_text = results[0]['generated_text']
                # Remove the original prompt from the response
                if generated_text.startswith(prompt):
                    generated_text = generated_text[len(prompt):].strip()
                return generated_text
            
            return ""
            
        except Exception as e:
            logger.error(f"HuggingFace generation error: {e}")
            return ""
    
    def stream_generate(self, prompt: str, model: str, temperature: float = 0.3, max_tokens: int = 2000) -> Generator[str, None, None]:
        """Generate text using HuggingFace with pseudo-streaming."""
        # HuggingFace doesn't have native streaming, so we simulate it
        result = self.generate(prompt, model, temperature, max_tokens)
        if result:
            # Split into chunks for streaming effect
            words = result.split()
            chunk_size = max(1, len(words) // 20)  # ~20 chunks
            for i in range(0, len(words), chunk_size):
                chunk = " ".join(words[i:i + chunk_size])
                if i + chunk_size < len(words):
                    chunk += " "
                yield chunk


class ModelInterface:
    """Main interface for managing LLM backends with auto-detection."""
    
    def __init__(self, config: ReportConfig):
        self.config = config
        self.backend = None
        self.current_model = None
        self._initialize_backend()
    
    def _initialize_backend(self):
        """Initialize the appropriate LLM backend based on config and availability."""
        if self.config.backend == LLMBackend.NONE:
            logger.info("LLM backend disabled, using fallback mode")
            return
        
        # Try to initialize the configured backend
        if self.config.backend == LLMBackend.OLLAMA:
            ollama = OllamaBackend()
            if ollama.is_available():
                self.backend = ollama
                self.current_model = ollama.get_best_model(self.config.model_name)
                logger.info(f"Using Ollama backend with model: {self.current_model}")
                return
            else:
                logger.warning("Ollama backend requested but not available")
        
        elif self.config.backend == LLMBackend.HUGGINGFACE:
            hf = HuggingFaceBackend()
            if hf.is_available():
                self.backend = hf
                self.current_model = hf.get_best_model(self.config.model_name)
                logger.info(f"Using HuggingFace backend with model: {self.current_model}")
                return
            else:
                logger.warning("HuggingFace backend requested but not available")
        
        # Auto-detect available backends if configured backend failed
        logger.info("Auto-detecting available LLM backends...")
        
        # Try Ollama first
        ollama = OllamaBackend()
        if ollama.is_available():
            self.backend = ollama
            self.current_model = ollama.get_best_model()
            logger.info(f"Auto-detected Ollama backend with model: {self.current_model}")
            return
        
        # Try HuggingFace
        hf = HuggingFaceBackend()
        if hf.is_available():
            self.backend = hf
            self.current_model = hf.get_best_model()
            logger.info(f"Auto-detected HuggingFace backend with model: {self.current_model}")
            return
        
        # No LLM available, use fallback mode
        logger.info("No LLM backends available, using fallback mode")
    
    def is_available(self) -> bool:
        """Check if any LLM backend is available."""
        return self.backend is not None and self.current_model is not None
    
    def get_backend_info(self) -> Dict[str, Any]:
        """Get information about the current backend."""
        if not self.is_available():
            return {
                "backend": "none",
                "model": None,
                "available": False,
                "fallback_mode": True
            }
        
        backend_type = "ollama" if isinstance(self.backend, OllamaBackend) else "huggingface"
        return {
            "backend": backend_type,
            "model": self.current_model,
            "available": True,
            "fallback_mode": False
        }
    
    def generate(self, prompt: str, stream_callback=None) -> str:
        """Generate text, optionally with streaming."""
        if not self.is_available():
            logger.error("No LLM backend available")
            return ""
        
        if self.config.enable_streaming and stream_callback:
            # Stream generation
            full_response = ""
            for chunk in self.backend.stream_generate(
                prompt, 
                self.current_model,
                temperature=self.config.temperature,
                max_tokens=self.config.max_tokens
            ):
                full_response += chunk
                # Wrap chunk in structured format for callbacks
                stream_callback({
                    "section": "generation",
                    "content": chunk,
                    "total_length": len(full_response)
                })
            return full_response
        else:
            # Regular generation
            return self.backend.generate(
                prompt,
                self.current_model,
                temperature=self.config.temperature,
                max_tokens=self.config.max_tokens
            ) 