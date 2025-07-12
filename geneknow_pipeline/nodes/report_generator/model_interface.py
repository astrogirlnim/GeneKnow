"""
Model interface for LLM backends (Ollama) with auto-detection and fallback.
"""

import logging
import requests
import json
import os
from typing import Optional, Dict, Any, List, Generator
from .config import ReportConfig, LLMBackend

logger = logging.getLogger(__name__)


class OllamaBackend:
    """Ollama LLM backend interface."""

    def __init__(self, base_url: str = "http://localhost:11434"):
        self.base_url = base_url
        self.available_models = []
        # Use a session for connection pooling and better performance
        self.session = requests.Session()
        self.session.headers.update({'Content-Type': 'application/json'})
        self._check_availability()

    def _check_availability(self) -> bool:
        """Check if Ollama is running and get available models."""
        try:
            response = self.session.get(f"{self.base_url}/api/tags", timeout=5)
            if response.status_code == 200:
                data = response.json()
                self.available_models = [
                    model["name"] for model in data.get("models", [])
                ]
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
            "llama3.1:8b",  # Updated to prefer 3.1
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
        """Generate text using Ollama with optimized settings."""
        try:
            payload = {
                "model": model,
                "prompt": prompt,
                "stream": False,
                "options": {
                    "temperature": temperature,
                    "num_predict": max_tokens,
                    "num_ctx": 4096,  # Increased context window
                    "num_thread": 4,  # Optimize thread usage
                }
            }
            
            response = self.session.post(
                f"{self.base_url}/api/generate",
                json=payload,
                timeout=180  # Increased timeout for parallel requests
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

    def stream_generate(
        self, prompt: str, model: str, temperature: float = 0.3, max_tokens: int = 2000
    ) -> Generator[str, None, None]:
        """Generate text using Ollama with streaming."""
        try:
            payload = {
                "model": model,
                "prompt": prompt,
                "stream": True,
                "options": {
                    "temperature": temperature,
                    "num_predict": max_tokens,
                    "num_ctx": 4096,
                    "num_thread": 4,
                }
            }
            
            response = self.session.post(
                f"{self.base_url}/api/generate",
                json=payload,
                stream=True,
                timeout=180
            )
            if response.status_code == 200:
                for line in response.iter_lines():
                    if line:
                        try:
                            data = json.loads(line)
                            if "response" in data:
                                yield data["response"]
                        except json.JSONDecodeError:
                            continue
            else:
                logger.error(f"Ollama streaming failed: {response.status_code}")

        except Exception as e:
            logger.error(f"Ollama streaming error: {e}")



class ModelInterface:
    """Unified interface for LLM backends with auto-detection."""

    def __init__(self, config: ReportConfig):
        self.config = config
        self.backend = None
        self.current_model = None
        self._initialize_backend()

    def _initialize_backend(self):
        """Initialize the appropriate backend based on config or auto-detection."""
        if self.config.backend == LLMBackend.OLLAMA:
            ollama = OllamaBackend()
            if ollama.is_available():
                self.backend = ollama
                self.current_model = self._select_model(ollama)
                logger.info(f"Using Ollama backend with model: {self.current_model}")
            else:
                logger.warning("Ollama backend requested but not available")
        elif self.config.backend == LLMBackend.NONE:
            logger.info("Using fallback mode (no LLM)")
        else:
            # Auto-detect
            # Try Ollama first
            ollama = OllamaBackend()
            if ollama.is_available():
                self.backend = ollama
                self.current_model = self._select_model(ollama)
                logger.info(f"Auto-detected Ollama backend with model: {self.current_model}")
                return

            # No backend available
            logger.info("No LLM backend available, using fallback mode")

    def _select_model(self, backend):
        """Select the best available model for a backend."""
        if self.config.model_name and self.config.model_name != "auto":
            return self.config.model_name

        # Auto-select based on backend type
        if isinstance(backend, OllamaBackend):
            return backend.get_best_model()
        
        return None

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
                "fallback_mode": True,
            }

        backend_type = "ollama" if isinstance(self.backend, OllamaBackend) else "none"
        return {"backend": backend_type, "model": self.current_model, "available": True, "fallback_mode": False}

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
                max_tokens=self.config.max_tokens,
            ):
                full_response += chunk
                # Wrap chunk in structured format for callbacks
                stream_callback(
                    {
                        "section": "generation",
                        "content": chunk,
                        "total_length": len(full_response),
                    }
                )
            return full_response
        else:
            # Regular generation
            return self.backend.generate(
                prompt,
                self.current_model,
                temperature=self.config.temperature,
                max_tokens=self.config.max_tokens,
            )
