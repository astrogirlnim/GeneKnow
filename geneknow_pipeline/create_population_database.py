#!/usr/bin/env python3
"""
Create Population Database Script
Downloads ClinVar data and creates population_variants.db for the GeneKnow pipeline.
This script ensures reproducibility by recreating the database from source data.

Usage:
    python create_population_database.py [--cancer-genes-only] [--force]

Options:
    --cancer-genes-only: Only process variants in known cancer genes (faster)
    --force: Overwrite existing database without prompting
"""

import os
import sys
import sqlite3
import requests
import gzip
import argparse
import logging
from typing import Dict, List

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Database path
DB_PATH = "population_variants.db"

# ClinVar FTP URLs
CLINVAR_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
CLINVAR_VCF_URL = f"{CLINVAR_BASE_URL}/clinvar.vcf.gz"
CLINVAR_VCF_INDEX_URL = f"{CLINVAR_BASE_URL}/clinvar.vcf.gz.tbi"

# Cancer-associated genes (subset for faster processing)
CANCER_GENES = {
    "BRCA1",
    "BRCA2",
    "TP53",
    "KRAS",
    "APC",
    "EGFR",
    "PIK3CA",
    "PTEN",
    "ATM",
    "CHEK2",
    "PALB2",
    "NBN",
    "NF1",
    "RB1",
    "VHL",
    "MLH1",
    "MSH2",
    "MSH6",
    "PMS2",
    "EPCAM",
    "CDKN2A",
    "STK11",
    "SMAD4",
    "BMPR1A",
    "MUTYH",
    "CDH1",
    "PTCH1",
    "SUFU",
    "TSC1",
    "TSC2",
    "WT1",
    "MEN1",
    "RET",
    "SDHB",
    "SDHC",
    "SDHD",
    "MAX",
    "TMEM127",
    "FLCN",
    "FH",
    "BRAF",
    "NRAS",
    "KIT",
    "PDGFRA",
    "JAK2",
    "FLT3",
    "NPM1",
    "CEBPA",
    "RUNX1",
    "ASXL1",
    "DNMT3A",
    "TET2",
    "IDH1",
    "IDH2",
    "SF3B1",
    "SRSF2",
    "U2AF1",
    "ZRSR2",
    "BCOR",
    "STAG2",
    "RAD21",
    "SMC1A",
    "SMC3",
    "CTCF",
    "ARID1A",
    "ARID1B",
    "ARID2",
    "SMARCA4",
    "SMARCB1",
    "PBRM1",
    "SETD2",
    "BAP1",
    "CREBBP",
    "EP300",
    "KDM6A",
    "KMT2D",
    "EZH2",
    "SUZ12",
    "EED",
    "FBXW7",
    "NOTCH1",
    "NOTCH2",
    "NOTCH3",
    "NOTCH4",
}


def download_file(url: str, filename: str, chunk_size: int = 8192) -> str:
    """Download a file with progress indication."""
    logger.info(f"Downloading {url}...")

    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()

        total_size = int(response.headers.get("content-length", 0))
        downloaded = 0

        with open(filename, "wb") as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)

                    if total_size > 0:
                        progress = (downloaded / total_size) * 100
                        sys.stdout.write(
                            f"\rProgress: {progress:.1f}% ({downloaded:,}/{total_size:,} bytes)"
                        )
                        sys.stdout.flush()

        print()  # New line after progress
        logger.info(f"Downloaded {filename} successfully")
        return filename

    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to download {url}: {e}")
        raise


def parse_clinvar_vcf(vcf_file: str, cancer_genes_only: bool = False) -> List[Dict]:
    """Parse ClinVar VCF file and extract relevant variant information."""
    logger.info(f"Parsing ClinVar VCF file: {vcf_file}")

    variants = []
    processed_count = 0
    kept_count = 0

    try:
        with gzip.open(vcf_file, "rt") as f:
            for line in f:
                line = line.strip()

                # Skip header lines
                if line.startswith("#"):
                    continue

                processed_count += 1

                # Progress indication
                if processed_count % 10000 == 0:
                    logger.info(
                        f"Processed {processed_count:,} variants, kept {kept_count:,}"
                    )

                # Parse VCF line
                fields = line.split("\t")
                if len(fields) < 8:
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                info = fields[7]

                # Skip multi-allelic sites for simplicity
                if "," in alt:
                    continue

                # Parse INFO field
                info_dict = {}
                for item in info.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info_dict[key] = value

                # Extract relevant fields
                gene_symbol = (
                    info_dict.get("GENEINFO", "").split(":")[0]
                    if "GENEINFO" in info_dict
                    else None
                )
                clinical_significance = info_dict.get("CLNSIG", "")
                review_status = info_dict.get("CLNREVSTAT", "")
                consequence = (
                    info_dict.get("MC", "").split("|")[1]
                    if "MC" in info_dict and "|" in info_dict["MC"]
                    else ""
                )

                # Skip if no gene symbol
                if not gene_symbol:
                    continue

                # Filter by cancer genes if requested
                if cancer_genes_only and gene_symbol not in CANCER_GENES:
                    continue

                # Determine pathogenicity
                is_pathogenic = 0
                if any(
                    term in clinical_significance.lower()
                    for term in ["pathogenic", "likely_pathogenic"]
                ):
                    is_pathogenic = 1
                elif any(
                    term in clinical_significance.lower()
                    for term in ["benign", "likely_benign"]
                ):
                    is_pathogenic = 0
                else:
                    is_pathogenic = (
                        0  # Default to non-pathogenic for uncertain variants
                    )

                # Extract population frequency from gnomAD if available
                gnomad_af = 0.0
                if "AF_ESP" in info_dict:
                    try:
                        gnomad_af = float(info_dict["AF_ESP"])
                    except ValueError:
                        pass
                elif "AF_EXAC" in info_dict:
                    try:
                        gnomad_af = float(info_dict["AF_EXAC"])
                    except ValueError:
                        pass
                elif "AF_TGP" in info_dict:
                    try:
                        gnomad_af = float(info_dict["AF_TGP"])
                    except ValueError:
                        pass

                # Create variant record
                variant = {
                    "chrom": chrom.replace("chr", ""),  # Normalize chromosome format
                    "pos": pos,
                    "ref": ref if ref != "." else "na",
                    "alt": alt if alt != "." else "na",
                    "gene": gene_symbol,
                    "gnomad_af": gnomad_af,
                    "clinical_significance": clinical_significance,
                    "is_pathogenic": is_pathogenic,
                    "consequence": consequence,
                    "review_status": review_status,
                }

                variants.append(variant)
                kept_count += 1

                # Limit for testing (remove in production)
                if kept_count >= 200000:  # Process up to 200k variants
                    logger.info("Reached processing limit, stopping...")
                    break

    except Exception as e:
        logger.error(f"Error parsing VCF file: {e}")
        raise

    logger.info(
        f"Parsing complete. Processed {processed_count:,} variants, kept {kept_count:,}"
    )
    return variants


def create_database(variants: List[Dict], db_path: str = DB_PATH) -> None:
    """Create the population_variants.db database from parsed variants."""
    logger.info(f"Creating database: {db_path}")

    # Remove existing database if it exists
    if os.path.exists(db_path):
        os.remove(db_path)

    # Create database connection
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Create table
        cursor.execute(
            """
        CREATE TABLE population_variants (
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            gene TEXT,
            gnomad_af REAL DEFAULT 0.0,
            clinical_significance TEXT,
            is_pathogenic INTEGER DEFAULT 0,
            consequence TEXT,
            review_status TEXT,
            PRIMARY KEY (chrom, pos, ref, alt)
        )
        """
        )

        # Create indexes
        cursor.execute("CREATE INDEX idx_gene ON population_variants(gene)")
        cursor.execute("CREATE INDEX idx_position ON population_variants(chrom, pos)")
        cursor.execute(
            "CREATE INDEX idx_pathogenic ON population_variants(is_pathogenic)"
        )

        # Insert variants
        logger.info("Inserting variants into database...")
        insert_query = """
        INSERT OR REPLACE INTO population_variants
        (chrom, pos, ref, alt, gene, gnomad_af, clinical_significance, is_pathogenic, consequence, review_status)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        batch_size = 1000
        for i in range(0, len(variants), batch_size):
            batch = variants[i : i + batch_size]
            batch_data = [
                (
                    v["chrom"],
                    v["pos"],
                    v["ref"],
                    v["alt"],
                    v["gene"],
                    v["gnomad_af"],
                    v["clinical_significance"],
                    v["is_pathogenic"],
                    v["consequence"],
                    v["review_status"],
                )
                for v in batch
            ]
            cursor.executemany(insert_query, batch_data)

            if (i + batch_size) % 10000 == 0:
                logger.info(f"Inserted {i + batch_size:,} variants...")

        conn.commit()

        # Get final statistics
        cursor.execute("SELECT COUNT(*) FROM population_variants")
        total_variants = cursor.fetchone()[0]

        cursor.execute(
            "SELECT COUNT(*) FROM population_variants WHERE is_pathogenic = 1"
        )
        pathogenic_variants = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(DISTINCT gene) FROM population_variants")
        unique_genes = cursor.fetchone()[0]

        logger.info("Database created successfully!")
        logger.info(f"  Total variants: {total_variants:,}")
        logger.info(f"  Pathogenic variants: {pathogenic_variants:,}")
        logger.info(f"  Unique genes: {unique_genes:,}")

    except Exception as e:
        logger.error(f"Error creating database: {e}")
        raise
    finally:
        conn.close()


def validate_database(db_path: str = DB_PATH) -> bool:
    """Validate the created database."""
    logger.info("Validating database...")

    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Check table exists
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='population_variants'"
        )
        if not cursor.fetchone():
            logger.error("Table 'population_variants' not found")
            return False

        # Check for cancer genes
        test_genes = ["BRCA1", "BRCA2", "TP53", "KRAS", "APC"]
        for gene in test_genes:
            cursor.execute(
                "SELECT COUNT(*) FROM population_variants WHERE gene = ?", (gene,)
            )
            count = cursor.fetchone()[0]
            if count > 0:
                logger.info(f"✓ Found {count} variants in {gene}")
            else:
                logger.warning(f"⚠ No variants found in {gene}")

        # Test the population_mapper query
        cursor.execute(
            """
        SELECT gene, gnomad_af, clinical_significance, is_pathogenic, consequence
        FROM population_variants
        WHERE chrom = ? AND pos = ?
        LIMIT 1
        """,
            ("17", 41234451),
        )

        row = cursor.fetchone()
        if row:
            logger.info(f"✓ Database query test passed: {row}")
        else:
            logger.info("ℹ No specific test variant found (this is normal)")

        conn.close()
        return True

    except Exception as e:
        logger.error(f"Database validation failed: {e}")
        return False


def main():
    """Main function to create the population database."""
    parser = argparse.ArgumentParser(
        description="Create population variants database from ClinVar data"
    )
    parser.add_argument(
        "--cancer-genes-only",
        action="store_true",
        help="Only process variants in known cancer genes (faster)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing database without prompting",
    )

    args = parser.parse_args()

    # Check if database already exists
    if os.path.exists(DB_PATH) and not args.force:
        response = input(f"Database {DB_PATH} already exists. Overwrite? (y/N): ")
        if response.lower() != "y":
            logger.info("Operation cancelled")
            return

    # Create temporary directory for downloads
    temp_dir = "temp_clinvar_download"
    os.makedirs(temp_dir, exist_ok=True)

    try:
        # Download ClinVar VCF file
        vcf_file = os.path.join(temp_dir, "clinvar.vcf.gz")
        download_file(CLINVAR_VCF_URL, vcf_file)

        # Parse VCF file
        variants = parse_clinvar_vcf(vcf_file, args.cancer_genes_only)

        if not variants:
            logger.error("No variants found in ClinVar data")
            return

        # Create database
        create_database(variants)

        # Validate database
        if validate_database():
            logger.info("✅ Database creation completed successfully!")

            # Show usage instructions
            print("\n" + "=" * 60)
            print("Database Creation Complete!")
            print("=" * 60)
            print(f"Created: {DB_PATH}")
            print(f"Size: {os.path.getsize(DB_PATH) / (1024*1024):.1f} MB")
            print("\nThe database is ready for use with the GeneKnow pipeline.")
            print("You can now run the pipeline with confidence that it will")
            print("properly validate genetic variants against population data.")
            print("=" * 60)
        else:
            logger.error("❌ Database validation failed")
            return

    except Exception as e:
        logger.error(f"Failed to create database: {e}")
        return
    finally:
        # Clean up temporary files
        if os.path.exists(temp_dir):
            import shutil

            shutil.rmtree(temp_dir)
            logger.info("Cleaned up temporary files")


if __name__ == "__main__":
    main()
