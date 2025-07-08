# 📁 GenePredict File Structure Documentation

## Project Overview
This document outlines the complete file structure for the GenePredict project, including all directories and their purposes.

## Root Directory Structure
```
LiteratureGapper/
├── .cursor/                    # Cursor IDE configuration and rules
├── .git/                       # Git repository data
├── .github/                    # GitHub Actions workflows and templates
├── desktop/                    # Desktop application (Tauri + React)
│   ├── src-tauri/             # Rust backend application
│   │   ├── src/               # Rust source files
│   │   ├── Cargo.toml         # Rust dependencies
│   │   └── tauri.conf.json    # Tauri configuration
│   ├── python_ml/             # Python ML scripts for genomic processing
│   │   ├── fastq_to_vcf_pipeline.py     # FASTQ to VCF conversion pipeline
│   │   ├── extract_by_region.py         # Genomic region extraction tool
│   │   ├── config_data_source.py        # Data source configuration
│   │   └── generate_test_fastq.py       # Test FASTQ file generator
│   └── ui/                    # React frontend application
│       ├── src/               # React source files
│       ├── public/            # Static assets
│       ├── package.json       # Node.js dependencies
│       └── vite.config.ts     # Vite build configuration
├── docs/                      # Project documentation
│   ├── CI_CD_Pipeline_Guide.md
│   ├── CURRENT_STATUS.md
│   ├── Phase1_Foundation_Tauri_Setup_Plan.md
│   └── file_structure.md      # This file
├── documentation/             # Legacy documentation and concepts
│   ├── GenePredict_BrainLift/ # Project planning documents
│   ├── p4-BrainLift-altoghether/ # Problem statements and goals
│   ├── GenePredict-Concept.markdown
│   ├── Literature_Gap_Mapper_Idea.md
│   └── Parallel-Genomic-Plan.md
├── memory_bank/               # AI memory bank for project context
│   ├── memory_bank_projectbrief.md
│   ├── memory_bank_progress.md
│   └── memory_bank_activeContext.md
├── data_processing/           # Data processing and download scripts
│   ├── tcga_download/         # TCGA data download scripts
│   │   ├── src/               # Source code
│   │   ├── config/            # Configuration files
│   │   ├── logs/              # Log files
│   │   └── downloads/         # Downloaded data storage
│   ├── genomic_processing/    # Genomic data processing scripts
│   └── utils/                 # Shared utilities
├── .gitignore                 # Git ignore patterns
├── README.md                  # Main project documentation
└── PRD_V2.md                  # Product Requirements Document
```

## Directory Purposes

### `/desktop/` - Desktop Application
- **Purpose**: Cross-platform desktop application built with Tauri
- **Contents**: React frontend (UI), Rust backend (src-tauri), and Python ML scripts (python_ml)
- **Status**: Phase 1 Complete (Foundation), Phase 2 In Progress (Python Integration)

### `/desktop/python_ml/` - Python ML Scripts
- **Purpose**: Python scripts for genomic data processing and machine learning
- **Contents**: FASTQ/VCF processing, region extraction, test data generation
- **Status**: Phase 2 Development (Rust ⇄ Python integration complete)

### `/docs/` - Active Documentation
- **Purpose**: Current project documentation and guides
- **Contents**: Status reports, setup guides, and structure documentation
- **Status**: Actively maintained

### `/documentation/` - Legacy Documentation
- **Purpose**: Historical documentation and project concepts
- **Contents**: Original project concepts, planning documents, and brainstorming
- **Status**: Reference material

### `/memory_bank/` - AI Memory Bank
- **Purpose**: AI assistant memory and project context
- **Contents**: Project brief, progress tracking, and active context
- **Status**: Actively maintained by AI assistant

### `/data_processing/` - Data Processing Scripts
- **Purpose**: Python scripts for genomic data processing and download
- **Contents**: TCGA download scripts, genomic processing utilities
- **Status**: Phase 2 Development (New)

## File Naming Conventions

### Scripts
- Use snake_case for Python files: `tcga_download_manager.py`
- Use descriptive names: `blood_cancer_filter.py`
- Include version in major scripts: `gdc_api_client_v1.py`

### Configuration
- Use lowercase with underscores: `tcga_config.yaml`
- Include environment suffix: `production.env`, `development.env`

### Documentation
- Use PascalCase for main docs: `FileStructure.md`
- Use descriptive titles: `TCGA_Download_Guide.md`

## Access Patterns

### Data Flow
1. **Input**: Raw genomic data files (FASTQ, BAM, VCF)
2. **Processing**: Python scripts in `/data_processing/`
3. **Analysis**: Tauri application in `/desktop/`
4. **Output**: Risk reports and visualizations

### Development Workflow
1. **Planning**: Update memory bank and documentation
2. **Implementation**: Add scripts to appropriate `/data_processing/` subdirectories
3. **Testing**: Use local test data, log extensively
4. **Integration**: Connect to desktop application via Tauri plugins

## Security Considerations

### Data Privacy
- All processing remains local (no external API calls for analysis)
- Downloaded data stored in `/data_processing/*/downloads/` with appropriate .gitignore
- Authentication tokens stored in environment variables, never committed

### File Permissions
- Scripts require execute permissions: `chmod +x script.py`
- Configuration files should be read-only: `chmod 444 config.yaml`
- Log files should be append-only: `chmod 644 logs/*.log`

## Version Control

### Tracked Files
- All source code and documentation
- Configuration templates (without sensitive data)
- Project structure and setup scripts

### Ignored Files
- Downloaded genomic data files (large binary files)
- Log files (generated at runtime)
- Environment files with secrets
- Cache and temporary files

This structure supports the project's privacy-first, local-processing architecture while maintaining clear organization and scalability. 