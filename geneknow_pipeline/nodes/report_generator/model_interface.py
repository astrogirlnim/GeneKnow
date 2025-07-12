"""
Model interface for LLM backends (Ollama, HuggingFace) with auto-detection and fallback.
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


class HuggingFaceBackend:
    """HuggingFace transformers backend interface with optimized loading."""
    
    def __init__(self):
        self.pipeline = None
        self.current_model = None
        self.available_models = []
        self._model_cache = {}  # Cache loaded models
        self._check_availability()

    def _check_availability(self) -> bool:
        """Check if HuggingFace transformers is available."""
        try:
            pass

            # Check for common medical/scientific models
            # Prioritize smaller, faster models for better user experience
            self.available_models = [
                "microsoft/DialoGPT-medium",
                "microsoft/phi-2", 
                "google/flan-t5-base",
                "distilgpt2",  # Smaller, faster alternative
                "gpt2",  # Fallback option
            ]

            logger.info("HuggingFace transformers available")
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
        
        # Priority for medical/scientific writing (optimized for speed)
        priority_models = [
            "microsoft/phi-2",  # Good balance of capability and speed
            "distilgpt2",  # Fast and efficient
            "google/flan-t5-base",  # Good for instruction following
            "microsoft/DialoGPT-medium",  # Larger but more capable
            "gpt2"  # Fallback
        ]
        for model in priority_models:
            if model in self.available_models:
                return model

        return self.available_models[0]

    def _load_pipeline(self, model_name: str):
        """Load the HuggingFace pipeline for a model with caching and optimization."""
        # Check if model is already cached
        if model_name in self._model_cache:
            self.pipeline = self._model_cache[model_name]
            self.current_model = model_name
            logger.info(f"Using cached HuggingFace model: {model_name}")
            return True
        
        # Load new model
        if self.pipeline is None or self.current_model != model_name:
            try:
                logger.info(f"Loading HuggingFace model: {model_name} (this may take a moment...)")
                
                from transformers import pipeline
                import torch
                
                # Optimize model loading parameters
                device = 0 if torch.cuda.is_available() else -1  # Use GPU if available
                
                # Create pipeline with optimized settings
                self.pipeline = pipeline(
                    "text-generation", 
                    model=model_name, 
                    device=device,
                    torch_dtype=torch.float16 if torch.cuda.is_available() else torch.float32,  # Use half precision on GPU
                    trust_remote_code=True,  # For some models
                    pad_token_id=50256,  # Common pad token
                    do_sample=True,
                    temperature=0.7,
                    top_p=0.9,
                    repetition_penalty=1.1,
                    return_full_text=False,  # Don't return the input prompt
                    clean_up_tokenization_spaces=True,  # Clean up tokenization
                )
                
                # Cache the loaded model
                self._model_cache[model_name] = self.pipeline
                self.current_model = model_name
                
                logger.info(f"Successfully loaded HuggingFace model: {model_name}")
                return True
            except Exception as e:
                logger.error(f"Failed to load HuggingFace model {model_name}: {e}")
                # Try to fall back to a simpler model
                if model_name != "distilgpt2":
                    logger.info("Attempting fallback to distilgpt2...")
                    return self._load_pipeline("distilgpt2")
                return False
        return True

    def generate(self, prompt: str, model: str, temperature: float = 0.3, max_tokens: int = 2000) -> str:
        """Generate text using HuggingFace with optimized parameters."""
        try:
            if not self._load_pipeline(model):
                return ""
            
            # Optimize generation parameters for longer, more coherent output
            # Remove conflicting parameters that cause warnings/hangs
            generation_params = {
                "max_new_tokens": max_tokens,  # Use max_new_tokens instead of max_length
                "temperature": max(temperature, 0.1),  # Ensure minimum temperature
                "do_sample": True,
                "top_p": 0.9,
                "top_k": 50,
                "repetition_penalty": 1.1,
                "pad_token_id": self.pipeline.tokenizer.eos_token_id,
                "eos_token_id": self.pipeline.tokenizer.eos_token_id,
                # Remove early_stopping and num_beams conflict
                "num_beams": 1,  # Use greedy search with sampling
                "return_full_text": False,  # Only return generated text, not prompt
            }
            
            # Add specific prompt formatting for better results
            formatted_prompt = self._format_prompt_for_generation(prompt, model)
            
            logger.info(f"Starting HuggingFace generation with {model}")
            
            # Use threading-based timeout for cross-platform compatibility
            import threading
            import queue
            
            result_queue = queue.Queue()
            exception_queue = queue.Queue()
            
            def run_generation():
                try:
                    results = self.pipeline(
                        formatted_prompt,
                        **generation_params
                    )
                    result_queue.put(results)
                except Exception as e:
                    exception_queue.put(e)
            
            # Start generation in a separate thread
            generation_thread = threading.Thread(target=run_generation)
            generation_thread.daemon = True
            generation_thread.start()
            
            # Wait for completion with timeout
            generation_thread.join(timeout=180)  # 3 minutes
            
            if generation_thread.is_alive():
                logger.error("HuggingFace generation timed out after 3 minutes")
                return ""
            
            # Check for exceptions
            if not exception_queue.empty():
                raise exception_queue.get()
            
            # Get results
            if not result_queue.empty():
                results = result_queue.get()
                
                if results and len(results) > 0:
                    generated_text = results[0]['generated_text']
                    
                    # Clean up the response - handle both cases where prompt might be included
                    if generated_text.startswith(formatted_prompt):
                        generated_text = generated_text[len(formatted_prompt):].strip()
                    
                    # Post-process to improve quality
                    generated_text = self._post_process_text(generated_text)
                    
                    logger.info(f"HuggingFace generation completed successfully")
                    return generated_text
                
                logger.warning("HuggingFace generation returned empty results")
                return ""
            else:
                logger.warning("HuggingFace generation completed but no results found")
                return ""
            
        except Exception as e:
            logger.error(f"HuggingFace generation error: {e}")
            return ""
    
    def _format_prompt_for_generation(self, prompt: str, model_name: str) -> str:
        """Format prompt for better generation based on model type."""
        if "phi" in model_name.lower():
            # Phi models work better with instruction format
            return f"Instruct: {prompt}\nOutput:"
        elif "flan" in model_name.lower():
            # Flan models are instruction-tuned
            return f"Please complete the following task:\n{prompt}\n\nResponse:"
        elif "gpt" in model_name.lower():
            # GPT models work well with direct prompts
            return f"{prompt}\n\n"
        else:
            # Default formatting
            return f"{prompt}\n\n"
    
    def _post_process_text(self, text: str) -> str:
        """Post-process generated text for better quality."""
        # Remove common artifacts
        text = text.strip()
        
        # Remove incomplete sentences at the end
        sentences = text.split('.')
        if len(sentences) > 1 and len(sentences[-1].strip()) < 10:
            text = '.'.join(sentences[:-1]) + '.'
        
        # Remove excessive whitespace
        import re
        text = re.sub(r'\n\s*\n', '\n\n', text)
        text = re.sub(r' +', ' ', text)
        
        return text
    
    def stream_generate(self, prompt: str, model: str, temperature: float = 0.3, max_tokens: int = 2000) -> Generator[str, None, None]:
        """Generate text using HuggingFace with pseudo-streaming (optimized)."""
        # HuggingFace doesn't have native streaming, so we simulate it
        result = self.generate(prompt, model, temperature, max_tokens)
        if result:
            # Split into chunks for streaming effect (more granular)
            words = result.split()
            chunk_size = max(1, len(words) // 30)  # ~30 chunks for smoother streaming
            for i in range(0, len(words), chunk_size):
                chunk = " ".join(words[i : i + chunk_size])
                if i + chunk_size < len(words):
                    chunk += " "
                yield chunk
                
                # Add small delay for realistic streaming effect
                import time
                time.sleep(0.1)


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
                logger.info(
                    f"Using HuggingFace backend with model: {self.current_model}"
                )
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
            logger.info(
                f"Auto-detected Ollama backend with model: {self.current_model}"
            )
            return

        # Try HuggingFace
        hf = HuggingFaceBackend()
        if hf.is_available():
            self.backend = hf
            self.current_model = hf.get_best_model()
            logger.info(
                f"Auto-detected HuggingFace backend with model: {self.current_model}"
            )
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
                "fallback_mode": True,
            }

        backend_type = (
            "ollama" if isinstance(self.backend, OllamaBackend) else "huggingface"
        )
        return {
            "backend": backend_type,
            "model": self.current_model,
            "available": True,
            "fallback_mode": False,
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
