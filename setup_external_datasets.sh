#!/bin/bash
# Script to download and organize genomic datasets on external SSD

# Configuration - update this path to your external SSD
EXTERNAL_SSD="/Volumes/YourSSDName"
DATA_DIR="$EXTERNAL_SSD/genomic-datasets"

# Create directory structure
echo "🗂️  Setting up directory structure..."
mkdir -p "$DATA_DIR/1000genomes/phase3"
mkdir -p "$DATA_DIR/1000genomes/phase1"
mkdir -p "$DATA_DIR/gnomad"
mkdir -p "$DATA_DIR/clinvar"

# Function to download with progress
download_with_progress() {
    local url=$1
    local output=$2
    echo "📥 Downloading $(basename $output)..."
    curl -L --progress-bar "$url" -o "$output"
}

# Download 1000 Genomes Phase 3 (latest, most complete)
echo "🧬 Downloading 1000 Genomes Phase 3 data..."
echo "This will take a while - each file is 200MB-1.8GB"

# Priority cancer-related chromosomes
PRIORITY_CHROMOSOMES=(1 2 3 7 9 10 11 12 13 17 19 20 22 X)

for chr in "${PRIORITY_CHROMOSOMES[@]}"; do
    FILE="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${FILE}"
    
    if [ ! -f "$DATA_DIR/1000genomes/phase3/${FILE}" ]; then
        download_with_progress "$URL" "$DATA_DIR/1000genomes/phase3/${FILE}"
        
        # Download index file too
        download_with_progress "${URL}.tbi" "$DATA_DIR/1000genomes/phase3/${FILE}.tbi"
    else
        echo "✅ Already have chr${chr}"
    fi
done

# Create symlinks in your project
echo "🔗 Creating symlinks in your project..."
ln -sf "$DATA_DIR" "$HOME/Gauntlet-Projects/LiteratureGapper/external-data"

echo "✅ Setup complete! Your external datasets are at: $DATA_DIR"
echo "📊 Space used: $(du -sh "$DATA_DIR" | cut -f1)" 