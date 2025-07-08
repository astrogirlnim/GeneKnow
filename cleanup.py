#!/usr/bin/env python3
"""
Cleanup script to remove temporary files while keeping test infrastructure.
"""
import os
import shutil
import glob

def cleanup():
    """Remove temporary files and directories."""
    
    print("🧹 Cleaning up temporary files...")
    
    # Files to remove in root directory
    root_files_to_remove = [
        "test_R1_aligned_sorted.bam",
        "test_R1_aligned_sorted.bam.bai",
        "test_R1_aligned.sam",
        "*.pyc",
        ".DS_Store"
    ]
    
    # Files to remove in langgraph directory
    langgraph_files_to_remove = [
        "langgraph/test_parallelization_output.txt",
        "langgraph/test_pipeline_results.json", 
        "langgraph/simple_server.log",
        "langgraph/server.log",
        "langgraph/test_variants.vcf",  # Temporary test file
        "langgraph/*.pyc"
    ]
    
    # Directories to remove
    dirs_to_remove = [
        "__pycache__",
        "langgraph/__pycache__",
        "langgraph/nodes/__pycache__",
        ".pytest_cache",
        "langgraph/.pytest_cache"
    ]
    
    # Clean up root files
    for pattern in root_files_to_remove:
        for file in glob.glob(pattern):
            try:
                os.remove(file)
                print(f"  ✓ Removed: {file}")
            except FileNotFoundError:
                pass
            except Exception as e:
                print(f"  ✗ Error removing {file}: {e}")
    
    # Clean up langgraph files
    for pattern in langgraph_files_to_remove:
        for file in glob.glob(pattern):
            try:
                os.remove(file)
                print(f"  ✓ Removed: {file}")
            except FileNotFoundError:
                pass
            except Exception as e:
                print(f"  ✗ Error removing {file}: {e}")
    
    # Clean up directories
    for dir_path in dirs_to_remove:
        if os.path.exists(dir_path):
            try:
                shutil.rmtree(dir_path)
                print(f"  ✓ Removed directory: {dir_path}")
            except Exception as e:
                print(f"  ✗ Error removing directory {dir_path}: {e}")
    
    print("\n✅ Cleanup complete!")
    print("\n📁 Kept:")
    print("  - Test scripts (test_pipeline.py, test_parallelization.py, etc.)")
    print("  - Test data files (test_R1.fastq.gz, test_R2.fastq.gz)")
    print("  - Core pipeline code")
    print("  - Configuration and setup files")
    print("  - Database files (SQLite, ML models)")
    

if __name__ == "__main__":
    cleanup() 