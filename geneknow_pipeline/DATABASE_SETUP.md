# Population Database Setup Guide

This guide explains how to create and maintain the `population_variants.db` database that powers the GeneKnow pipeline's genetic variant validation system.

## Overview

The GeneKnow pipeline relies on a population database containing clinical variant information from ClinVar. This database is essential for:
- Validating genetic variants against known population data
- Assessing clinical significance (pathogenic vs benign)
- Calculating accurate cancer risk scores
- Preventing false positives in risk assessment

## Database Structure

The database contains multiple tables that work together to provide comprehensive variant annotations:

### 1. Main Population Variants Table

```sql
CREATE TABLE population_variants (
    chrom TEXT,                    -- Chromosome (1-22, X, Y)
    pos INTEGER,                   -- Position on chromosome
    ref TEXT,                      -- Reference allele
    alt TEXT,                      -- Alternative allele
    gene TEXT,                     -- Gene symbol
    gnomad_af REAL DEFAULT 0.0,    -- Population frequency
    clinical_significance TEXT,     -- ClinVar clinical significance
    is_pathogenic INTEGER DEFAULT 0, -- 1 if pathogenic, 0 if benign
    consequence TEXT,              -- Variant consequence type
    review_status TEXT,            -- ClinVar review status
    PRIMARY KEY (chrom, pos, ref, alt)
);
```

### 2. CADD Scores Table

```sql
CREATE TABLE cadd_scores (
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
```

### 3. Job Tracking Table

```sql
CREATE TABLE cadd_jobs (
    job_id TEXT PRIMARY KEY,
    job_type TEXT NOT NULL,      -- 'batch_import', 'api_lookup', 'test'
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP,
    variant_count INTEGER DEFAULT 0,
    status TEXT DEFAULT 'running', -- 'running', 'completed', 'failed'
    metadata TEXT                 -- JSON string with additional info
);
```

### 4. Variant Annotations View

A convenient view that joins population and CADD data:

```sql
CREATE VIEW variant_annotations AS
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
```

## Quick Setup (Recommended)

For team members setting up the project:

### 1. Download Pre-built Database
If available, download the pre-built database from your team's shared storage:
```bash
# Replace with your team's storage location
curl -o population_variants.db "https://your-team-storage/population_variants.db"
```

### 2. Create Database from Source
If you need to create the database from scratch:

```bash
# Navigate to the pipeline directory
cd geneknow_pipeline

# Install required dependencies
pip install requests

# Create database (cancer genes only - faster)
python create_population_database.py --cancer-genes-only

# OR create full database (slower but more comprehensive)
python create_population_database.py
```

### 3. Add CADD Scores Table
After creating the population database, add the CADD scores table:

```bash
# Add test CADD data for development
./scripts/fetch_cadd.sh --test

# OR for production setup (creates table structure only)
./scripts/fetch_cadd.sh
```

### 4. Verify Setup
```bash
# Check database exists and is valid
python -c "
import sqlite3
conn = sqlite3.connect('population_variants.db')
cursor = conn.cursor()
cursor.execute('SELECT COUNT(*) FROM population_variants')
print(f'Population variants: {cursor.fetchone()[0]:,}')
cursor.execute('SELECT COUNT(*) FROM cadd_scores')
print(f'CADD scores: {cursor.fetchone()[0]:,}')
conn.close()
"
```

## Detailed Setup Instructions

### Prerequisites

1. **Python 3.8+** with pip
2. **Internet connection** for downloading ClinVar data
3. **~2GB disk space** for temporary files and database

### Step-by-Step Database Creation

#### 1. Install Dependencies
```bash
pip install requests sqlite3
```

#### 2. Choose Your Database Size

**Option A: Cancer Genes Only (Recommended for Development)**
- Faster to create (~10-30 minutes)
- Smaller database size (~10-50 MB)
- Contains variants in 70+ cancer-associated genes
- Perfect for testing and development

```bash
python create_population_database.py --cancer-genes-only
```

**Option B: Full ClinVar Database (Production)**
- Longer to create (~1-3 hours)
- Larger database size (~100-500 MB)
- Contains all ClinVar variants
- Recommended for production use

```bash
python create_population_database.py
```

#### 3. Monitor Progress
The script will show progress updates:
```
2024-01-15 10:00:00 - INFO - Downloading https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz...
Progress: 45.2% (12,345,678/27,345,678 bytes)
2024-01-15 10:05:00 - INFO - Parsing ClinVar VCF file...
2024-01-15 10:05:00 - INFO - Processed 50,000 variants, kept 12,345
2024-01-15 10:10:00 - INFO - Creating database: population_variants.db
2024-01-15 10:12:00 - INFO - Database created successfully!
```

#### 4. Validation
The script automatically validates the database:
```
2024-01-15 10:12:00 - INFO - Validating database...
2024-01-15 10:12:00 - INFO - ✓ Found 1,234 variants in BRCA1
2024-01-15 10:12:00 - INFO - ✓ Found 987 variants in BRCA2
2024-01-15 10:12:00 - INFO - ✓ Found 2,345 variants in TP53
2024-01-15 10:12:00 - INFO - ✅ Database creation completed successfully!
```

## CADD Scores Integration

### Adding CADD Scores

The CADD (Combined Annotation Dependent Depletion) scores provide deleteriousness predictions for variants:

1. **Test Data Setup** (for development):
   ```bash
   ./scripts/fetch_cadd.sh --test
   ```
   This adds 5 sample CADD scores for testing.

2. **Production Import** (requires CADD data files):
   ```bash
   # Download CADD data (warning: >100GB)
   wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
   
   # Import using custom script (to be implemented)
   python create_cadd_import.py --input whole_genome_SNVs.tsv.gz --job-id cadd_v1.7_import
   ```

3. **Remote Lookup Alternative**:
   For development, the pipeline can use remote CADD lookups instead of local storage.
   Set `USE_REMOTE_CADD=true` environment variable.

### Job Tracking

Each CADD scoring run creates a job record for tracking:

```sql
-- View recent CADD jobs
SELECT job_id, job_type, status, variant_count, completed_at
FROM cadd_jobs
ORDER BY started_at DESC
LIMIT 10;

-- View variants from a specific job
SELECT COUNT(*) FROM cadd_scores WHERE job_id = 'cadd_20240115_123456_abc123';
```

## Usage Options

### Command Line Arguments

```bash
python create_population_database.py [OPTIONS]

Options:
  --cancer-genes-only    Only process variants in known cancer genes (faster)
  --force               Overwrite existing database without prompting
  -h, --help            Show help message
```

### Examples

```bash
# Create cancer-focused database (recommended for most use cases)
python create_population_database.py --cancer-genes-only

# Force overwrite existing database
python create_population_database.py --force

# Create full database for production
python create_population_database.py

# Quick test run (modify script to set limit=1000 for testing)
python create_population_database.py --cancer-genes-only --force
```

## Troubleshooting

### Common Issues

#### 1. Download Failures
```
ERROR - Failed to download https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz: HTTPError
```
**Solution:** Check internet connection and try again. ClinVar servers may be temporarily unavailable.

#### 2. Disk Space Issues
```
ERROR - No space left on device
```
**Solution:** Ensure you have at least 2GB free space. Clean up temporary files:
```bash
rm -rf temp_clinvar_download/
```

#### 3. Permission Errors
```
ERROR - Permission denied: population_variants.db
```
**Solution:** Ensure write permissions in the directory:
```bash
chmod 755 .
rm -f population_variants.db  # Remove if exists
```

#### 4. Memory Issues
```
ERROR - MemoryError during parsing
```
**Solution:** Use the `--cancer-genes-only` flag to reduce memory usage.

### Database Issues

#### Verify Database Integrity
```bash
# Check table structure
sqlite3 population_variants.db ".schema"

# Check record count
sqlite3 population_variants.db "SELECT COUNT(*) FROM population_variants;"

# Check CADD scores
sqlite3 population_variants.db "SELECT COUNT(*) FROM cadd_scores;"

# Check cancer genes
sqlite3 population_variants.db "SELECT gene, COUNT(*) FROM population_variants WHERE gene IN ('BRCA1', 'BRCA2', 'TP53') GROUP BY gene;"
```

#### Rebuild Database
```bash
# Remove corrupted database
rm -f population_variants.db

# Recreate
python create_population_database.py --cancer-genes-only --force

# Re-add CADD scores
./scripts/fetch_cadd.sh --test
```

## Database Maintenance

### Regular Updates

ClinVar data is updated monthly. To keep your database current:

```bash
# Check current database date
sqlite3 population_variants.db "SELECT COUNT(*) FROM population_variants;"

# Update database (every 3-6 months)
python create_population_database.py --force
```

### Backup Strategy

```bash
# Create backup
cp population_variants.db population_variants_backup_$(date +%Y%m%d).db

# Compress for storage
gzip population_variants_backup_$(date +%Y%m%d).db
```

## Team Setup Instructions

### For Project Lead

1. Create the database once:
   ```bash
   python create_population_database.py --cancer-genes-only
   ./scripts/fetch_cadd.sh --test
   ```

2. Share the database file with team members via:
   - Cloud storage (Google Drive, Dropbox)
   - Internal file server
   - Git LFS (if configured)

### For Team Members

1. Download the shared database file
2. Place it in the `geneknow_pipeline/` directory
3. Verify it works:
   ```bash
   python test_enhanced_api.py
   ```

## Performance Notes

### Database Size vs Speed Trade-offs

| Option | Database Size | Creation Time | Variants | Genes | Use Case |
|--------|---------------|---------------|----------|--------|-----------|
| Cancer genes only | ~10-50 MB | 10-30 min | ~50K-200K | 70+ | Development, Testing |
| Full ClinVar | ~100-500 MB | 1-3 hours | ~1M+ | All | Production |

### Query Performance

The database includes optimized indexes:
- `idx_gene`: Fast gene-based lookups
- `idx_position`: Fast position-based lookups
- `idx_pathogenic`: Fast pathogenicity filtering
- `idx_cadd_location`: Fast CADD score lookups
- `idx_cadd_phred`: Fast CADD threshold queries
- `idx_cadd_job`: Fast job-based filtering

Typical query times:
- Single variant lookup: <1ms
- Gene-based queries: <10ms
- Population frequency queries: <5ms
- CADD score lookup: <1ms
- Combined annotation lookup: <5ms

## Security Considerations

### Data Privacy
- ClinVar data is public and does not contain personal information
- The database contains only variant annotations, not patient data
- Safe to share within research teams

### Access Control
- No authentication required for database access
- Consider file-level permissions for sensitive environments
- Database is read-only during pipeline execution

## Support

### Getting Help

1. **Check the logs** for detailed error messages
2. **Verify prerequisites** (Python version, disk space, internet)
3. **Try the --cancer-genes-only flag** for faster, smaller database
4. **Check ClinVar status** at https://www.ncbi.nlm.nih.gov/clinvar/

### Reporting Issues

When reporting issues, please include:
- Error message (full traceback)
- Python version (`python --version`)
- Operating system
- Available disk space
- Command used

## Technical Details

### Data Source
- **Primary**: ClinVar VCF files (GRCh38)
- **URL**: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
- **Update frequency**: Monthly
- **Format**: VCF 4.2 with ClinVar annotations

### Processing Pipeline
1. Download compressed VCF file
2. Parse variant records
3. Extract clinical annotations
4. Normalize chromosome formats
5. Filter by cancer genes (if requested)
6. Create SQLite database with indexes
7. Validate database integrity

### Quality Assurance
- Automated validation of database structure
- Verification of cancer gene content
- Test queries against known variants
- Performance benchmarking

This database setup ensures that the GeneKnow pipeline can accurately validate genetic variants and provide reliable risk assessments. The reproducible creation process enables team collaboration and consistent results across different environments. 