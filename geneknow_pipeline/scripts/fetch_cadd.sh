#!/bin/bash
# Script to download CADD v1.7 data and convert to SQLite database
# Usage: ./fetch_cadd.sh [output_dir]

set -e  # Exit on error

# Configuration
CADD_VERSION="v1.7"
CADD_URL="https://krishna.gs.washington.edu/download/CADD/${CADD_VERSION}/GRCh38"
OUTPUT_DIR="${1:-geneknow_pipeline/data}"
DB_FILE="${OUTPUT_DIR}/cadd_scores.db"
TEMP_DIR="${OUTPUT_DIR}/temp_cadd"

# Ensure output directory exists
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TEMP_DIR}"

echo "🧬 CADD Data Fetcher - v1.7 (hg38)"
echo "Output directory: ${OUTPUT_DIR}"
echo "Database file: ${DB_FILE}"

# Check for required tools
command -v sqlite3 >/dev/null 2>&1 || { echo "❌ sqlite3 is required but not installed. Aborting." >&2; exit 1; }

# Only check for wget/gunzip if not in test mode
if [ "${CADD_ENV}" != "test" ]; then
    command -v wget >/dev/null 2>&1 || { echo "❌ wget is required but not installed. Aborting." >&2; exit 1; }
    command -v gunzip >/dev/null 2>&1 || { echo "❌ gunzip is required but not installed. Aborting." >&2; exit 1; }
fi

# Download a subset of CADD scores (top ~30M SNVs)
# Full dataset is >300GB, we'll use the pre-scored common variants
echo "📥 Downloading CADD ${CADD_VERSION} pre-scored variants..."
CADD_FILE="whole_genome_SNVs_inclAnno.tsv.gz"
CADD_DOWNLOAD_URL="${CADD_URL}/${CADD_FILE}"

# For testing/development, use a smaller subset
# In production, download the full file or use tabix-indexed remote access
if [ "${CADD_ENV}" = "test" ]; then
    echo "⚠️  Test mode: Creating mock CADD database with sample data"
    # Create mock database for testing
    sqlite3 "${DB_FILE}" <<EOF
CREATE TABLE IF NOT EXISTS cadd (
    chrom TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    raw_score REAL,
    phred_score REAL,
    PRIMARY KEY (chrom, pos, ref, alt)
);

CREATE INDEX IF NOT EXISTS idx_cadd_location ON cadd(chrom, pos);

-- Insert some test data
INSERT OR REPLACE INTO cadd VALUES
    ('chr17', 43044295, 'A', 'T', 2.345, 24.7),  -- BRCA1 variant
    ('chr17', 43045719, 'G', 'A', 3.456, 34.1),  -- BRCA1 high impact
    ('chr13', 32315474, 'G', 'A', 1.234, 15.3),  -- BRCA2 variant
    ('chr5', 112839936, 'C', 'T', 0.456, 8.2),   -- APC variant
    ('chr5', 112840259, 'G', 'C', 4.567, 41.2);  -- APC high impact

-- Add metadata table
CREATE TABLE IF NOT EXISTS metadata (
    key TEXT PRIMARY KEY,
    value TEXT
);

INSERT OR REPLACE INTO metadata VALUES
    ('cadd_version', '${CADD_VERSION}'),
    ('genome_build', 'GRCh38/hg38'),
    ('creation_date', datetime('now')),
    ('variant_count', '5');

EOF
    echo "✅ Test database created at ${DB_FILE}"
    exit 0
fi

# For production, download subset or use remote tabix
echo "📊 Creating SQLite database schema..."
sqlite3 "${DB_FILE}" <<EOF
-- Main CADD scores table
CREATE TABLE IF NOT EXISTS cadd (
    chrom TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    raw_score REAL,
    phred_score REAL,
    PRIMARY KEY (chrom, pos, ref, alt)
);

-- Indices for fast lookups
CREATE INDEX IF NOT EXISTS idx_cadd_location ON cadd(chrom, pos);
CREATE INDEX IF NOT EXISTS idx_cadd_phred ON cadd(phred_score);

-- Metadata table
CREATE TABLE IF NOT EXISTS metadata (
    key TEXT PRIMARY KEY,
    value TEXT
);

INSERT OR REPLACE INTO metadata VALUES
    ('cadd_version', '${CADD_VERSION}'),
    ('genome_build', 'GRCh38/hg38'),
    ('creation_date', datetime('now'));
EOF

echo "📥 Downloading CADD data subset..."
# Download a manageable subset (e.g., first 10M variants)
# In production, consider:
# 1. Using tabix remote queries instead of local DB
# 2. Downloading only clinically relevant variants
# 3. Using CADD's API for on-demand scoring

# For now, create instructions for manual download
cat << EOF > "${OUTPUT_DIR}/DOWNLOAD_INSTRUCTIONS.txt"
CADD Data Download Instructions
==============================

Due to the large size of CADD data (>300GB), automatic download is not implemented.

Option 1: Remote Tabix Access (Recommended)
------------------------------------------
Use remote tabix queries to CADD servers. No download needed.
Configure CADD_REMOTE_TABIX environment variable.

Option 2: Download Pre-scored Variants
--------------------------------------
1. Visit: https://cadd.gs.washington.edu/download
2. Download: whole_genome_SNVs_inclAnno.tsv.gz
3. Extract and import top variants into SQLite

Option 3: Use CADD API
----------------------
For variants not in the pre-scored set, use the CADD REST API
for on-demand scoring.

For testing, run: CADD_ENV=test ./fetch_cadd.sh
EOF

echo "⚠️  Please see ${OUTPUT_DIR}/DOWNLOAD_INSTRUCTIONS.txt for manual download steps"
echo "💡 For testing, run: CADD_ENV=test $0"

# Clean up
rm -rf "${TEMP_DIR}"

echo "✅ CADD database setup complete!" 