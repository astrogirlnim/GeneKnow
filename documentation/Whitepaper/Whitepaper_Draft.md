# Geneknow: A Privacy-First Local Genomic Risk Assessment Platform

> **Disclaimer:** Geneknow is a statistical tool for investigative and research purposes only. It does not provide medical or clinical advice, diagnosis, or treatment recommendations. All results are intended for informational use by qualified professionals. Users must not rely on Geneknow for clinical decision-making. 

## 1. Executive Summary

**Narrative Overview:**
Geneknow is a completely free and open-source desktop application designed for privacy-preserving genomic risk assessment, enabling clinicians and researchers to analyze genetic data entirely on local hardware without external data transmission. Built on a Tauri-based architecture with React frontend, Rust integration, and Python ML backend, it processes formats such as FASTQ, BAM, VCF, and MAF to generate interpretable risk reports for specific cancers including blood, breast, colon, lung, and prostate, with potential for future expansion.

Key differentiators include its HIPAA-friendly design (not formally HIPAA-compliant), cross-platform support (Windows, macOS, Ubuntu, macOS Intel), and use of LangGraph for orchestrating accurate, reproducible genomic predictions. By leveraging aggregate statistics databases without storing actual patient data, FASTQs, or any Protected Health Information (PHI), Geneknow ensures complete data privacy while providing clinically relevant insights. The application is entirely free to use with no licensing fees, subscription costs, or commercial restrictions.

This whitepaper outlines Geneknow's architecture, scientific foundation, and value in clinical and research settings.

**Technical/Implementation Details:**
- **Statistical/Investigative Use Only:** Not a diagnostic or clinical tool; see disclaimer above.
- **Privacy-First:** No PHI or raw patient data is stored; all processing is local and offline.
- **Cross-Platform:** Runs on Windows, macOS (Intel/Apple Silicon), and Ubuntu Linux.
- **LangGraph Orchestration:** Ensures reproducible, modular, and traceable workflows.
- **No Licensing Fees:** 100% free and open-source, with no commercial restrictions.
- **Extensible Plugin System:** Supports custom ML nodes and analysis modules (infrastructure implemented).

## 2. Introduction & Motivation

**Narrative Overview:**
Genomic data analysis has revolutionized personalized medicine, but traditional cloud-based tools raise significant privacy concerns, with data breaches affecting millions. Regulatory frameworks like HIPAA and GDPR demand stringent data protection, yet many solutions require internet connectivity and external servers, limiting accessibility in low-resource or secure environments.

Geneknow addresses these challenges with a fully local, offline-capable platform that performs end-to-end genomic risk assessment on the user's device. Motivated by the need for accessible, compliant tools, Geneknow empowers users to analyze genetic variants for cancer risk without compromising privacy.

**Technical/Implementation Details:**
- All analysis is performed on the user's device, with no data ever leaving the machine.
- Designed for researchers and clinicians who need rapid, reproducible, and private genomic risk assessment for investigative purposes.

## 3. System Overview

**Narrative Overview:**
Geneknow employs a hybrid architecture combining Tauri's lightweight Rust core for system integration, a React-based frontend for intuitive user interaction, and a Python backend for ML-driven analysis. The system uses LangGraph to manage complex workflows, ensuring accurate predictions through modular, traceable node-based processing.

- **Input Handling:** Supports FASTQ, BAM, VCF, and MAF via secure, local file uploads with compressed format support (.gz extensions).
- **Processing Pipeline:** LangGraph-orchestrated nodes for variant calling, annotation, and risk scoring through 15 specialized processing nodes.
- **Output:** Comprehensive reports with visualizations, exportable in PDF and JSON formats.
- **Privacy Architecture:** All data remains on-device; no network calls post-installation.

Cross-platform compatibility is achieved through Tauri's bundling, with Rust handling low-level operations for efficiency.

**Technical/Implementation Details:**
- **Tauri (Rust) Core:** Handles system integration, plugin management, and secure file operations.
- **React Frontend:** Modern, tabbed UI for intuitive navigation, visualization, and report export.
- **Python Backend:** ML-driven analysis pipeline orchestrated by LangGraph, with 15+ specialized nodes.
- **Workflow:**
  1. **Input Handling:** Secure, local file uploads (FASTQ, BAM, VCF, MAF, .gz supported).
  2. **Processing Pipeline:** LangGraph DAG with parallel and sequential nodes for variant calling, annotation, risk scoring, and report generation.
  3. **Output:** Structured JSON and PDF reports, with visualizations and export options.
  4. **Privacy:** No network calls post-installation; all data remains on-device.
- **State Variables:** See `geneknow_pipeline/state.py` for all pipeline state fields (e.g., `file_path`, `risk_scores`, `structured_json`).

## 4. Core Technologies & Architecture

**Narrative Overview:**
- LangGraph workflow management for deterministic, reproducible predictions.
- Ensemble ML pipeline (Random Forest, Gradient Boosting, Linear Regression) trained on aggregate data from TCGA and COSMIC, with SHAP interpretability.
- Five static models: PRS, ClinVar, CADD, TCGA Mapper, Pathway Burden.
- Plugin system (infrastructure implemented) for future extensibility.
- SQLite databases with only aggregate statistics, no PHI.

**Technical/Implementation Details:**
### 4.1 LangGraph Workflow Orchestration
- **File:** `geneknow_pipeline/graph.py`
- **Nodes:** `file_input`, `preprocess`, `variant_calling`, `qc_filter`, `population_mapper`, `tcga_mapper`, `cadd_scoring`, `clinvar_annotator`, `prs_calculator`, `pathway_burden`, `ml_fusion_node`, `risk_model`, `shap_validator`, `metrics_calculator`, `formatter`, `report_writer`, and more.
- **Parallelization:** Variant calling/QC and static model scoring run in parallel for efficiency.

### 4.2 Machine Learning Pipeline
- **Files:** `ml_fusion_integration.py`, `ml_feature_extractor.py`, `ml_feature_extractor_no_leakage.py`, `ml_models/`
- **Models:** Random Forest, Gradient Boosting, Linear Regression (ensemble, trained on TCGA/COSMIC, no PHI)
- **Interpretability:** SHAP values for model transparency (`shap_validator` node)
- **No Data Leakage:** Training excludes clinical significance features (see `ml_feature_extractor_no_leakage.py`)

### 4.3 Static Models and Metrics
- **PRS:** `nodes/prs_calculator.py` (GWAS-based, population-specific, confidence-adjusted)
- **ClinVar:** `nodes/clinvar_annotator.py`, `population_variants.db` (local, curated, no PHI)
- **CADD:** `nodes/cadd_scoring.py` (offline, PHRED-scaled, cancer gene multipliers)
- **TCGA Mapper:** `nodes/population_mapper.py`, `gdc_api_client.py` (2,828 patient cohort, frequency-based)
- **Pathway Burden:** `nodes/pathway_burden.py` (gene/pathway disruption, burden aggregation)
- **ML Fusion:** `ml_fusion_node.py`, `ml_fusion_integration.py` (combines all static model outputs)

### 4.4 Plugin System
- **Files:** `desktop/python_ml/plugins/`, `desktop/src-tauri/src/plugin_registry.rs`, `plugin.rs`, `python_script_plugin.rs`
- **Features:** Manifest-based configuration, trait-based plugin interface, Python/Rust plugin support (infrastructure implemented, extensibility planned)

### 4.5 Data Management
- **Databases:** `population_variants.db`, `prs_snps.db` (aggregate stats only, no raw patient data)
- **Scripts:** `create_population_database.py` (reproducible, public data only)
- **No PHI:** Confirmed in all code and documentation

## 5. Privacy & Security Design

**Narrative Overview:**
Geneknow's HIPAA-friendly architecture processes all data locally, eliminating transmission risks and ensuring zero PHI storage. Features include encrypted temporary storage and automatic data cleanup post-analysis. While not formally HIPAA-certified, the design aligns with HIPAA technical safeguards including access controls and audit logs. The application bundles a complete Python runtime and dependencies to ensure no external dependencies are required during analysis. As a completely free and open-source solution, Geneknow removes cost barriers to genomic analysis while maintaining the highest privacy standards.

**Technical/Implementation Details:**
- **Local-Only Processing:** All analysis is performed on-device; no data leaves the user's machine (`enhanced_api_server.py`, `README.md`)
- **No PHI Storage:** Only aggregate statistics and annotations are stored; no raw sequences or patient identifiers
- **Encrypted Temp Storage & Cleanup:** Temporary files are encrypted and deleted post-analysis (`enhanced_api_server.py`, `UploadPage.tsx`)
- **No External Dependencies Post-Install:** All required resources are bundled; no internet required after setup
- **HIPAA-Friendly:** Designed to align with HIPAA technical safeguards, but not formally certified
- **Open Source:** Transparent, auditable codebase

## 6. User Interface and Experience

**Narrative Overview:**
Geneknow provides an intuitive, modern UI built with React and Tailwind CSS, emphasizing usability for clinical workflows.

### 6.1 Dashboard
The dashboard serves as the central hub for analysis results and quick insights, featuring:
- **Analysis Overview:** Displays probability scores and hazard scores with confidence indicators and SHAP validation when available.
- **Cancer Risk Assessment:** Shows cancer types with elevated risk above baseline thresholds, with risk percentages and affected genes.
- **Headline Metrics:** Interactive cards showing total variants found, processing time, and key variant details with tooltips.
- **Report Generation:** Tab-based interface for viewing AI-enhanced reports with markdown rendering and PDF export capabilities.
- **Visualization Widgets:** Risk distribution charts and variant type breakdowns with real pipeline data.

### 6.2 Clinical View (In-Depth Analysis Tab)
The Clinical View provides detailed, tabbed analysis for comprehensive genomic assessment:
- **Genomic Analysis Tab:** High-level summary with cancer risk cards, gene significance Manhattan plots, mutational signature analysis, and quality metrics display.
- **Variant Heatmap Tab:** Interactive heatmap showing gene-cancer associations based on pathway burden analysis, with summary statistics and pathway burden data visualization.
- **Pathway Analysis Tab:** Comprehensive pathway disruption analysis with cancer-pathway associations, disrupted pathway listings, and burden score visualizations.
- **Clinical Report Tab:** Survival analysis curves, clinical recommendations based on detected variants, and targeted therapy suggestions with prevention strategies.

Each tab includes dedicated export functionality with PDF generation capabilities that capture visualizations, summaries, and technical details for clinical documentation.

### 6.3 Export and Visualization Features
- **PDF Export:** High-resolution PDF generation with embedded visualizations using html2canvas and jsPDF libraries, including analysis summaries, technical details, and clinical insights.
- **Visualization Capture:** Automatic capture of interactive charts, heatmaps, and graphs for inclusion in exported reports.
- **Data Export:** JSON format exports for integration with external analysis tools and electronic health records.
- **Report Customization:** Selectable sections and configurable detail levels for different clinical use cases.

**Technical/Implementation Details:**
- **Frontend:** `desktop/ui/` (React, Tailwind CSS, Tauri integration)
- **Dashboard:** `DashboardPage.tsx` (headline metrics, risk assessment, variant details, report generation, visualization widgets)
- **Clinical View:** `ClinicalViewPage.tsx` (tabbed analysis: Genomic Analysis, Variant Heatmap, Pathway Analysis, Clinical Report)
- **Export & Visualization:** PDF export (jsPDF, html2canvas), JSON export, interactive charts, and heatmaps
- **Tab Navigation:** Dedicated tabs for each analysis type, with export options per tab
- **No Mock Data:** All visualizations and reports use real pipeline output
- **Accessibility:** Modern, responsive design with tooltips and clear navigation

## 7. Scientific Validation & Performance

**Narrative Overview:**
Models were trained on TCGA data (n=10,000+ samples), achieving AUC >0.85 for cancer risk prediction. The pipeline processes FASTQ files in approximately 0.7-1.0 seconds for 500 reads, VCF files in ~0.02 seconds for direct variant loading, with total pipeline execution under 1 second for most inputs. Local execution ensures consistency across hardware, with optimizations for multi-core processing.

**Technical/Implementation Details:**
- **Training Data:** TCGA (10,000+ samples), COSMIC, gnomAD (aggregate only)
- **Performance:**
  - FASTQ: ~0.7-1.0s for 500 reads
  - VCF: ~0.02s for direct variant loading
  - Total pipeline: <1s for most inputs
- **Validation:** AUC >0.85 for cancer risk prediction (see `README.md`, `test_pipeline_comprehensive.py`)
- **Testing:** Comprehensive test scripts, no mock data, real-world validation

## 8. Applications & Case Studies

**Narrative Overview:**
In clinical settings, Geneknow aids risk stratification for patients with suspected hereditary cancer syndromes. Research applications include batch processing for cohort studies, all performed offline. The platform supports point-of-care analysis in rural clinics and resource-limited settings where internet connectivity may be unreliable.

**Technical/Implementation Details:**
- **Clinical Research:** Risk stratification for hereditary cancer syndromes (investigative use only)
- **Batch Processing:** Cohort studies, offline analysis in research settings
- **Point-of-Care:** Usable in rural clinics and low-resource environments (no internet required)
- **User Workflows:** See `documentation/GenePredict_BrainLift/user-workflows.md`

## 9. Future Directions

**Narrative Overview:**
Planned expansions include additional cancer types beyond the current five supported types, completion of the plugin ecosystem for extensible analysis modules, enhanced survival analysis features, and potential mobile support for field-based genetic counseling.

**Technical/Implementation Details:**
- **Additional Cancer Types:** Planned expansion beyond current five
- **Plugin Ecosystem:** Full support for custom analysis modules (infrastructure in place)
- **Enhanced Survival Analysis:** More advanced prognostic modeling
- **Mobile Support:** Planned for future releases
- **LLM-Driven Reporting:** Ongoing improvements to narrative report generation

## 10. Conclusion

**Narrative Overview:**
Geneknow democratizes genomic analysis with privacy at its core, bridging accessibility and scientific rigor through local-first processing, comprehensive risk assessment, and clinician-friendly reporting. As a completely free and open-source platform that stores no PHI, Geneknow removes both cost and privacy barriers to genomic analysis, enabling equitable access to advanced genomic risk assessment tools for clinicians and researchers worldwide.

**Technical/Implementation Details:**
Geneknow democratizes genomic analysis with privacy at its core. By combining local-first processing, comprehensive risk assessment, and clinician-friendly reporting, it bridges accessibility and scientific rigor. As a free, open-source platform that stores no PHI, Geneknow removes both cost and privacy barriers, enabling equitable access to advanced genomic risk assessment for research and investigative purposes worldwide.

> **Disclaimer (repeated):** Geneknow is a statistical tool for investigative and research use only. It does not provide medical or clinical advice, diagnosis, or treatment. All results must be interpreted by qualified professionals and are not intended for clinical decision-making.

## References

### 10.1 Glossary
- **CADD**: Combined Annotation Dependent Depletion - variant deleteriousness score (PHRED-scaled, 0-40+ range).
- **PRS**: Polygenic Risk Score - cumulative genetic risk metric from GWAS studies.
- **SHAP**: SHapley Additive exPlanations - model interpretability method for ML predictions.
- **LangGraph**: Workflow orchestration framework for building complex, stateful processing pipelines.

### 10.2 Technical Specifications
- **Supported Formats:** FASTQ (.fastq, .fq), BAM (.bam), VCF (.vcf), MAF (.maf), with gzip compression support
- **Cancer Types:** Blood, breast, colon, lung, prostate (with future expansion capability)
- **Platforms:** Windows, macOS (Intel/Apple Silicon), Ubuntu Linux
- **Architecture:** Tauri 2.x + React 19 + Rust 1.88 + Python 3.11

### 10.3 API Documentation
- **POST /api/process:** Processes genomic file, returns job ID for status tracking
- **GET /api/results/{job_id}:** Returns comprehensive analysis results in structured JSON format
- **WebSocket /socket.io:** Real-time progress updates during analysis execution

---

This whitepaper is based on the current, implemented codebase and documentation. All claims are verifiable and reflect the actual state of the project as of this writing. For the latest updates, see the project repository and documentation. 