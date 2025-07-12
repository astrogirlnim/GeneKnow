#!/usr/bin/env python3
"""
Simple test for LLM detection only.
"""

import sys
from pathlib import Path

# Add the parent directory to the path so we can import the module
sys.path.insert(0, str(Path(__file__).parent))

from nodes.report_generator.config import load_config
from nodes.report_generator.model_interface import ModelInterface


def test_detection():
    print("=== Testing LLM Detection ===")

    # Load current config without modifying it
    config = load_config()
    print(f"Config backend: {config.backend}")

    # Test detection
    interface = ModelInterface(config)
    print(f"Detected backend: {interface.backend}")
    print(f"Available: {interface.is_available()}")

    if interface.is_available():
        print(f"Current model: {interface.current_model}")
        backend_info = interface.get_backend_info()
        print(f"Backend info: {backend_info}")

        # Test generation
        print("\n=== Testing Generation ===")
        result = interface.generate("Write a brief medical summary about BRCA1 mutations.")
        print(f"Generated {len(result)} characters")
        if result:
            print(f"First 200 chars: {result[:200]}...")
    else:
        print("No LLM backend available")


if __name__ == "__main__":
    test_detection()
