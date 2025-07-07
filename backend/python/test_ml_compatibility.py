#!/usr/bin/env python3
"""
GenePredict ML Compatibility Test Script
Tests the availability and compatibility of ML libraries for genomic analysis.
"""

import sys
import importlib
import platform
from typing import Dict, List, Tuple

def test_package_import(package_name: str, display_name: str = None) -> Tuple[bool, str]:
    """Test if a package can be imported and return version info."""
    if display_name is None:
        display_name = package_name
    
    try:
        module = importlib.import_module(package_name)
        version = getattr(module, '__version__', 'unknown')
        return True, f"{display_name} {version}"
    except ImportError as e:
        return False, f"{display_name} - NOT AVAILABLE ({str(e)})"

def check_python_version() -> Dict[str, str]:
    """Check Python version and compatibility."""
    python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    platform_info = platform.platform()
    
    return {
        "python_version": python_version,
        "platform": platform_info,
        "supports_tensorflow": sys.version_info < (3, 13),
        "supports_jax": True,  # JAX supports Python 3.9+
    }

def main():
    """Main compatibility test function."""
    print("🧬 GenePredict ML Compatibility Test")
    print("=" * 50)
    
    # Check Python environment
    env_info = check_python_version()
    print(f"\n📋 Environment Information:")
    print(f"  Python Version: {env_info['python_version']}")
    print(f"  Platform: {env_info['platform']}")
    print(f"  TensorFlow Compatible: {'✅ Yes' if env_info['supports_tensorflow'] else '❌ No (Python 3.13+)'}")
    print(f"  JAX Compatible: {'✅ Yes' if env_info['supports_jax'] else '❌ No'}")
    
    # Test core ML libraries
    print(f"\n🔬 Core ML Libraries:")
    core_ml_packages = [
        ("numpy", "NumPy"),
        ("pandas", "Pandas"),
        ("sklearn", "Scikit-learn"),
        ("torch", "PyTorch"),
    ]
    
    for package, display_name in core_ml_packages:
        success, info = test_package_import(package, display_name)
        status = "✅" if success else "❌"
        print(f"  {status} {info}")
    
    # Test TensorFlow vs JAX
    print(f"\n🤖 Deep Learning Frameworks:")
    
    # TensorFlow
    tf_success, tf_info = test_package_import("tensorflow", "TensorFlow")
    tf_status = "✅" if tf_success else "❌"
    print(f"  {tf_status} {tf_info}")
    
    # JAX ecosystem
    jax_packages = [
        ("jax", "JAX"),
        ("jaxlib", "JAXlib"),
        ("flax", "Flax"),
        ("optax", "Optax"),
    ]
    
    jax_available = 0
    for package, display_name in jax_packages:
        success, info = test_package_import(package, display_name)
        status = "✅" if success else "❌"
        print(f"  {status} {info}")
        if success:
            jax_available += 1
    
    # Test genomic analysis libraries
    print(f"\n🧬 Genomic Analysis Libraries:")
    genomic_packages = [
        ("pysam", "PySAM"),
        ("cyvcf2", "cyvcf2"),
        ("Bio", "BioPython"),
    ]
    
    for package, display_name in genomic_packages:
        success, info = test_package_import(package, display_name)
        status = "✅" if success else "❌"
        print(f"  {status} {info}")
    
    # Test web framework
    print(f"\n🌐 Web Framework:")
    web_packages = [
        ("fastapi", "FastAPI"),
        ("uvicorn", "Uvicorn"),
        ("pydantic", "Pydantic"),
    ]
    
    for package, display_name in web_packages:
        success, info = test_package_import(package, display_name)
        status = "✅" if success else "❌"
        print(f"  {status} {info}")
    
    # Provide recommendations
    print(f"\n💡 Recommendations:")
    
    if not env_info['supports_tensorflow']:
        print("  📝 TensorFlow is not compatible with Python 3.13+")
        print("     ➡️  Using JAX as the primary ML framework")
        if jax_available >= 2:
            print("     ✅ JAX ecosystem is properly installed")
        else:
            print("     ❌ JAX ecosystem needs installation")
            print("        Run: pip install jax jaxlib flax optax")
    else:
        if not tf_success:
            print("  📝 TensorFlow is not installed but is compatible with your Python version")
            print("     ➡️  Run: pip install tensorflow>=2.13.0")
        
        if jax_available >= 2:
            print("  📝 Both TensorFlow and JAX are available - excellent for flexibility!")
    
    # Test simple ML operation
    print(f"\n🧪 Testing Basic ML Operations:")
    
    try:
        import numpy as np
        test_array = np.random.random((10, 5))
        mean_val = np.mean(test_array)
        print(f"  ✅ NumPy operations working (test mean: {mean_val:.4f})")
    except Exception as e:
        print(f"  ❌ NumPy operations failed: {e}")
    
    try:
        import pandas as pd
        test_df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        print(f"  ✅ Pandas operations working (test shape: {test_df.shape})")
    except Exception as e:
        print(f"  ❌ Pandas operations failed: {e}")
    
    # Test JAX if available
    try:
        import jax.numpy as jnp
        test_jax = jnp.array([1, 2, 3, 4, 5])
        result = jnp.mean(test_jax)
        print(f"  ✅ JAX operations working (test mean: {result:.4f})")
    except Exception as e:
        print(f"  ❌ JAX operations not available: {e}")
    
    print(f"\n🎯 Summary:")
    print("  Your GenePredict environment is being configured for optimal")
    print("  compatibility with your Python version and available libraries.")
    print("  The application will use the best available ML framework.")

if __name__ == "__main__":
    main() 