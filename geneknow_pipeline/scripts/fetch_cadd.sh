#!/bin/bash
# Script to set up CADD scores table in existing population_variants database
# Usage: ./fetch_cadd.sh [--test]

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/.."
DB_FILE="${DATA_DIR}/population_variants.db"
CADD_VERSION="v1.7"

echo "🧬 CADD Scores Table Setup for GeneKnow Pipeline"
echo "================================================"
echo "Database: ${DB_FILE}"
echo ""

# Check if database exists
if [ ! -f "${DB_FILE}" ]; then
    echo "❌ Error: population_variants.db not found!"
    echo "   Please run create_population_database.py first"
    exit 1
fi

# Test mode - add sample data
if [ "$1" == "--test" ]; then
    echo "🧪 Test mode - Adding CADD scores table with sample data..."
    
    sqlite3 "${DB_FILE}" <<EOF
-- Create CADD scores table
CREATE TABLE IF NOT EXISTS cadd_scores (
    -- Primary key matching population_variants
    chrom TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    
    -- CADD scores
    raw_score REAL,
    phred_score REAL,
    
    -- Job tracking
    job_id TEXT,
    source TEXT DEFAULT 'local',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    PRIMARY KEY (chrom, pos, ref, alt),
    FOREIGN KEY (chrom, pos, ref, alt) 
        REFERENCES population_variants(chrom, pos, ref, alt)
);

-- Create indices for performance
CREATE INDEX IF NOT EXISTS idx_cadd_location ON cadd_scores(chrom, pos);
CREATE INDEX IF NOT EXISTS idx_cadd_phred ON cadd_scores(phred_score);
CREATE INDEX IF NOT EXISTS idx_cadd_job ON cadd_scores(job_id);

-- Create job tracking table
CREATE TABLE IF NOT EXISTS cadd_jobs (
    job_id TEXT PRIMARY KEY,
    job_type TEXT NOT NULL, -- 'batch_import', 'api_lookup', 'test'
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP,
    variant_count INTEGER DEFAULT 0,
    status TEXT DEFAULT 'running', -- 'running', 'completed', 'failed'
    metadata TEXT -- JSON string with additional info
);

-- Insert test job
INSERT INTO cadd_jobs (job_id, job_type, status, variant_count, completed_at, metadata) 
VALUES ('test_job_001', 'test', 'completed', 5, datetime('now'), 
    '{"cadd_version": "${CADD_VERSION}", "description": "Test data for CADD integration"}');

-- Insert test CADD scores (matching some variants that might be in population_variants)
INSERT OR REPLACE INTO cadd_scores (chrom, pos, ref, alt, raw_score, phred_score, job_id) VALUES
    ('17', 43044295, 'A', 'T', 2.345, 24.7, 'test_job_001'),  -- BRCA1 variant
    ('17', 43045719, 'G', 'A', 3.456, 34.1, 'test_job_001'),  -- BRCA1 high impact
    ('13', 32315474, 'G', 'A', 1.234, 15.3, 'test_job_001'),  -- BRCA2 variant
    ('5', 112839936, 'C', 'T', 0.456, 8.2, 'test_job_001'),   -- APC variant
    ('17', 42201391, 'C', 'T', 1.789, 18.9, 'test_job_001');  -- Test variant

-- Create view for easy variant lookup with both population and CADD data
CREATE VIEW IF NOT EXISTS variant_annotations AS
SELECT 
    pv.chrom, pv.pos, pv.ref, pv.alt, pv.gene,
    pv.gnomad_af, pv.clinical_significance, pv.is_pathogenic,
    pv.consequence, pv.review_status,
    cs.raw_score as cadd_raw, cs.phred_score as cadd_phred,
    cs.job_id as cadd_job_id, cs.created_at as cadd_scored_at
FROM population_variants pv
LEFT JOIN cadd_scores cs ON 
    pv.chrom = cs.chrom AND 
    pv.pos = cs.pos AND 
    pv.ref = cs.ref AND 
    pv.alt = cs.alt;

-- Show statistics
SELECT 'Population variants:' as metric, COUNT(*) as count FROM population_variants
UNION ALL
SELECT 'CADD scores added:', COUNT(*) FROM cadd_scores
UNION ALL
SELECT 'Variants with CADD:', COUNT(*) FROM variant_annotations WHERE cadd_phred IS NOT NULL;

EOF
    echo "✅ Test CADD scores table created and populated"
    exit 0
fi

# Production setup
echo "📊 Creating CADD scores table in population database..."
sqlite3 "${DB_FILE}" <<EOF
-- Create CADD scores table if not exists
CREATE TABLE IF NOT EXISTS cadd_scores (
    chrom TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    raw_score REAL,
    phred_score REAL,
    job_id TEXT,
    source TEXT DEFAULT 'local',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (chrom, pos, ref, alt),
    FOREIGN KEY (chrom, pos, ref, alt) 
        REFERENCES population_variants(chrom, pos, ref, alt)
);

-- Create indices
CREATE INDEX IF NOT EXISTS idx_cadd_location ON cadd_scores(chrom, pos);
CREATE INDEX IF NOT EXISTS idx_cadd_phred ON cadd_scores(phred_score);
CREATE INDEX IF NOT EXISTS idx_cadd_job ON cadd_scores(job_id);

-- Create job tracking table
CREATE TABLE IF NOT EXISTS cadd_jobs (
    job_id TEXT PRIMARY KEY,
    job_type TEXT NOT NULL,
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP,
    variant_count INTEGER DEFAULT 0,
    status TEXT DEFAULT 'running',
    metadata TEXT
);

-- Create variant annotations view
CREATE VIEW IF NOT EXISTS variant_annotations AS
SELECT 
    pv.chrom, pv.pos, pv.ref, pv.alt, pv.gene,
    pv.gnomad_af, pv.clinical_significance, pv.is_pathogenic,
    pv.consequence, pv.review_status,
    cs.raw_score as cadd_raw, cs.phred_score as cadd_phred,
    cs.job_id as cadd_job_id, cs.created_at as cadd_scored_at
FROM population_variants pv
LEFT JOIN cadd_scores cs ON 
    pv.chrom = cs.chrom AND 
    pv.pos = cs.pos AND 
    pv.ref = cs.ref AND 
    pv.alt = cs.alt;

EOF

echo ""
echo "📥 CADD Data Import Instructions"
echo "================================"
echo ""
echo "The CADD scores table has been created in population_variants.db"
echo ""
echo "To import CADD data:"
echo ""
echo "1. Download CADD v1.7 whole genome scores (warning: >100GB):"
echo "   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
echo ""
echo "2. Use the provided import script (create_cadd_import.py) to load data:"
echo "   python create_cadd_import.py --input whole_genome_SNVs.tsv.gz --job-id cadd_v1.7_import"
echo ""
echo "3. For development, use remote Tabix queries instead of local storage"
echo ""
echo "✅ CADD table structure created successfully!" 