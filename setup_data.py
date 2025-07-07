#!/usr/bin/env python3
"""
Smart genomic data setup - downloads files once, shares across team
"""
import os
import subprocess
import hashlib
from pathlib import Path

# Shared data locations (in priority order)
DATA_LOCATIONS = [
    "/shared/genomics-data",  # Lab shared storage
    "~/genomics-data",        # User home directory  
    "./data-cache",           # Project local cache
]

# Required files with checksums for verification
REQUIRED_FILES = {
    "chr1_phase1": {
        "url": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz",
        "filename": "ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz",
        "size_gb": 11,
        "sha256": None  # Add if available
    },
    "chr22_phase3": {
        "url": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", 
        "filename": "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        "size_gb": 0.2,
        "sha256": None
    }
}

def find_or_create_data_dir():
    """Find the best location for shared data"""
    for location in DATA_LOCATIONS:
        path = Path(location).expanduser()
        
        # Check if location exists and is writable
        if path.exists() and os.access(path, os.W_OK):
            print(f"✅ Using existing data directory: {path}")
            return path
            
        # Try to create location
        try:
            path.mkdir(parents=True, exist_ok=True)
            print(f"✅ Created new data directory: {path}")
            return path
        except (PermissionError, OSError):
            print(f"❌ Cannot create directory: {path}")
            continue
    
    raise RuntimeError("No suitable data directory found!")

def check_file_exists(data_dir, filename):
    """Check if file exists and is complete"""
    file_path = data_dir / filename
    if not file_path.exists():
        return False
        
    # Quick size check (more sophisticated checks possible)
    size_mb = file_path.stat().st_size / (1024 * 1024)
    if size_mb < 100:  # Probably incomplete
        print(f"⚠️  File {filename} exists but seems incomplete ({size_mb:.1f}MB)")
        return False
        
    print(f"✅ Found existing file: {filename} ({size_mb:.1f}MB)")
    return True

def download_file(url, destination):
    """Download with resume capability"""
    print(f"📥 Downloading {destination.name} ({REQUIRED_FILES[list(REQUIRED_FILES.keys())[0]]['size_gb']}GB)")
    print("☕ This will take a while... grab coffee!")
    
    # Use wget with resume capability
    cmd = ["wget", "-c", "-O", str(destination), url]
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ Download completed: {destination}")
        return True
    except subprocess.CalledProcessError:
        print(f"❌ Download failed: {destination}")
        return False

def create_symlinks(data_dir):
    """Create symlinks in project directory"""
    test_data_dir = Path("test-data")
    test_data_dir.mkdir(exist_ok=True)
    
    for file_info in REQUIRED_FILES.values():
        source = data_dir / file_info["filename"]
        target = test_data_dir / file_info["filename"]
        
        # Remove existing symlink/file
        if target.exists() or target.is_symlink():
            target.unlink()
            
        # Create symlink
        try:
            target.symlink_to(source.absolute())
            print(f"🔗 Created symlink: {target} -> {source}")
        except OSError:
            # Fallback: copy file (Windows compatibility)
            import shutil
            shutil.copy2(source, target)
            print(f"📋 Copied file: {target}")

def main():
    """Main setup process"""
    print("🧬 Setting up genomic data for the team...")
    
    # Find best data location
    data_dir = find_or_create_data_dir()
    
    # Check/download each required file
    for key, file_info in REQUIRED_FILES.items():
        filename = file_info["filename"]
        file_path = data_dir / filename
        
        if not check_file_exists(data_dir, filename):
            print(f"📥 Need to download: {filename}")
            if not download_file(file_info["url"], file_path):
                continue
                
        # Create index if missing
        index_path = file_path.with_suffix(file_path.suffix + ".tbi")
        if not index_path.exists():
            print(f"🔍 Creating index for {filename}...")
            subprocess.run(["bcftools", "index", str(file_path)], check=True)
    
    # Create symlinks in project
    create_symlinks(data_dir)
    
    print("\n🎉 Data setup complete!")
    print("💡 Team members just need to run: python3 setup_data.py")
    print("💡 Large files are shared, only downloaded once per machine")

if __name__ == "__main__":
    main() 