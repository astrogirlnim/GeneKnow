#!/usr/bin/env python3
"""
Clean up .tbi index files that bcftools downloads to the root directory
"""

import os
import glob
import shutil

def clean_index_files():
    """Move .tbi files from root to cache directory"""
    # Create cache directory if it doesn't exist
    cache_dir = "remote_index_cache"
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    
    # Find all .tbi files in root
    tbi_files = glob.glob("*.tbi")
    
    if not tbi_files:
        print("✅ No .tbi files in root directory")
        return
    
    print(f"🧹 Found {len(tbi_files)} .tbi files in root directory")
    
    # Move each file
    for tbi_file in tbi_files:
        dest = os.path.join(cache_dir, tbi_file)
        shutil.move(tbi_file, dest)
        print(f"  → Moved {tbi_file} to {cache_dir}/")
    
    print(f"✅ Cleaned up {len(tbi_files)} index files")

if __name__ == "__main__":
    clean_index_files() 