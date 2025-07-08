# LiteratureGapper

A comprehensive genomic data processing project with desktop UI and data processing pipelines for TCGA cancer data downloads and genomic variant analysis.

## Overview

LiteratureGapper combines a modern desktop application built with Tauri and React/TypeScript with powerful genomic data processing capabilities. The project includes:

- **Desktop Application**: Modern UI for genomic data visualization and analysis
- **Data Processing Pipelines**: FASTQ to VCF conversion and genomic variant extraction
- **TCGA Integration**: Blood cancer data download and processing
- **Genomic Analysis**: Precise region-based variant extraction from multiple data sources

## Prerequisites

Before getting started, ensure you have the following installed:

### Required Software
- **Node.js** (v16 or later) and **pnpm**
- **Rust** and **Cargo** (latest stable)
- **Python 3** (3.8 or later)
- **System Dependencies** (Linux):
  ```bash
  sudo apt update
  sudo apt install -y libgtk-3-dev libwebkit2gtk-4.1-dev libappindicator3-dev librsvg2-dev patchelf
  ```

### Genomic Tools (for data processing)
```bash
sudo apt install -y bcftools samtools tabix
```

### Verify Installation
```bash
node --version && pnpm --version
rustc --version && cargo --version
python3 --version
bcftools --version
```

## Quick Start

### 1. Install Dependencies
```bash
cd desktop/ui
pnpm install
```

### 2. Development Server
To run the desktop application in development mode:

```bash
# From the desktop directory
cd /path/to/LiteratureGapper/desktop
npx @tauri-apps/cli dev
```

**Note**: The development server will:
- Start the React frontend on http://localhost:5173/
- Launch the Tauri desktop application
- Enable hot-reload for both frontend and backend changes

### 3. Build for Production
```bash
# From the desktop directory
cd /path/to/LiteratureGapper/desktop
npx @tauri-apps/cli build
```

The built application will be available in `desktop/src-tauri/target/release/bundle/`.

## Genomic Data Processing

### Extract Variants by Region

The project includes a powerful genomic variant extraction tool that can work with multiple data sources:

```bash
# From the project root
python3 extract_by_region.py
```

#### Features:
- **Multiple Data Sources**: 1000 Genomes, gnomAD, local files
- **Parallel Processing**: Uses all CPU cores for fast extraction
- **Precise Boundaries**: Exact coordinate filtering (not just overlapping)
- **Performance Monitoring**: Detailed timing and memory usage statistics
- **Error Resilience**: Automatic retries and continue-on-error
- **Flexible Configuration**: Easy switching between remote and local data

#### Configuration

The extraction tool uses `config_data_source.py` for data source configuration:

```bash
# Use remote 1000 Genomes data (default - no download required)
export GENOMIC_DATA_SOURCE=remote

# Use local test data (limited chromosomes)
export GENOMIC_DATA_SOURCE=local

# Use gnomAD v3.1.2 (large files - 4-30GB each)
export GENOMIC_DATA_SOURCE=gnomad
```

#### Input Format

The tool expects a `regions.bed` file with tab-separated values:
```
chr1    35691274    35801992    TP73
chr17   7565097     7590868     TP53
chr13   32889611    32973805    BRCA2
```

Format: `chromosome    start    end    gene_name`

#### Output

- Creates `output_chunks/` directory with compressed VCF files
- Each gene gets its own `.vcf.gz` file
- Detailed performance metrics and variant counts
- Automatic cleanup of temporary files

## Project Structure

```
LiteratureGapper/
├── desktop/                    # Tauri desktop application
│   ├── ui/                    # React/TypeScript frontend
│   │   ├── src/               # Frontend source code
│   │   ├── package.json       # Frontend dependencies
│   │   └── ...
│   └── src-tauri/             # Rust backend
│       ├── src/               # Rust source code
│       ├── Cargo.toml         # Rust dependencies
│       └── tauri.conf.json    # Tauri configuration
├── extract_by_region.py       # Genomic variant extraction tool
├── config_data_source.py      # Data source configuration
├── regions.bed               # Example genomic regions
└── README.md                 # This file
```

## Development

### Frontend Development
```bash
cd desktop/ui
pnpm dev  # React development server only
```

### Backend Development
```bash
cd desktop/src-tauri
cargo run  # Rust backend only
```

### Full Stack Development
```bash
cd desktop
npx @tauri-apps/cli dev  # Both frontend and backend with hot-reload
```

## Data Sources

### Default: Remote 1000 Genomes
- **Pros**: No local storage required, comprehensive coverage
- **Cons**: Requires internet connection, slower than local data
- **Coverage**: All chromosomes 1-22 + X
- **Format**: Standard VCF with GRCh37 coordinates

### gnomAD v3.1.2
- **Pros**: Larger variant database, population frequency data
- **Cons**: Very large files (4-30GB each), requires significant bandwidth
- **Coverage**: Selected chromosomes
- **Format**: VCF with GRCh38 coordinates

### Local Data
- **Pros**: Fastest processing, no internet required
- **Cons**: Requires manual download and setup
- **Coverage**: Limited test chromosomes (1, 22)

## Performance

The genomic extraction tool is optimized for performance:

- **Parallel Processing**: Utilizes all CPU cores
- **Memory Efficient**: Processes regions independently
- **Compression**: Configurable compression levels (1-9)
- **Caching**: Remote index file caching
- **Monitoring**: Real-time performance metrics

Example output:
```
✅ TP53 (chr17): 245 variants, 23.4KB | Time: 2.34s | Memory: +12.5MB
📊 Total variants extracted: 8,432
⚡ Variants per second: 3,604
```

## Troubleshooting

### Common Issues

1. **Tauri CLI not found**
   ```bash
   # Use npx instead of global installation
   npx @tauri-apps/cli dev
   ```

2. **System dependencies missing (Linux)**
   ```bash
   sudo apt install -y libgtk-3-dev libwebkit2gtk-4.1-dev libappindicator3-dev librsvg2-dev patchelf
   ```

3. **bcftools not found**
   ```bash
   sudo apt install -y bcftools samtools tabix
   ```

4. **Network issues with remote data**
   ```bash
   # Switch to local test data
   export GENOMIC_DATA_SOURCE=local
   ```

### Getting Help

- Check the console output for detailed error messages
- Verify all prerequisites are installed
- Ensure internet connectivity for remote data sources
- Check file permissions for output directories

## Contributing

This project follows comprehensive development practices:

- **Thorough Planning**: Extensive analysis before implementation
- **Clean Architecture**: Well-organized folder structure
- **Comprehensive Logging**: Detailed logging throughout
- **Production Ready**: No mock data, real functionality from start
- **Memory Bank System**: Project documentation and knowledge management

## License

[License information to be added]

---

**Note**: This project is designed for genomic research and analysis. Ensure you have appropriate permissions and follow ethical guidelines when working with genomic data.
