# GenePredict - Project Brief Memory Bank

## Project Identity
**Name:** GenePredict  
**Type:** AI-Powered Genomic Risk Assessment Platform  
**Architecture:** Local-first desktop application (Tauri + React + Rust)  
**Privacy Model:** Complete local processing, zero external data transmission  

## Core Mission
Build a privacy-first genomic risk assessment platform that processes genetic data entirely on the user's local machine. No data ever leaves the device, ensuring complete privacy and HIPAA compliance for sensitive genetic information.

## Primary Goals
1. **Privacy-First Processing:** All genomic analysis runs locally with zero external dependencies
2. **AI-Powered Risk Assessment:** TensorFlow-based models for disease risk prediction (initial focus: breast cancer)
3. **Interpretable Reports:** Llama 3.1-generated summaries in multiple languages (English, Hindi, Spanish)
4. **Cross-Platform Accessibility:** Secure desktop app for macOS, Windows, and Linux
5. **Developer-Friendly:** Open-source, modular codebase for community contributions

## Success Criteria
- Process FASTQ/BAM/VCF files with 85%+ accuracy in risk prediction
- Generate human-readable reports in under 1 minute
- Maintain 100% local processing (no external API calls)
- Support cross-platform deployment with consistent UX
- Provide extensible architecture for additional disease models

## Target Users
- **Primary:** Clinicians in academic hospitals and rural/low-resource clinics
- **Secondary:** Bioinformatics researchers and ML scientists
- **Tertiary:** NGOs and global health organizations
- **Future:** Patients seeking personalized genomic insights

## Non-Goals
- Cloud-based processing or data storage
- Real-time collaboration features
- Mobile app versions
- Electronic Health Record (EHR) direct integration (Phase 1)
- Commercial licensing or closed-source components

## Project Scope Boundaries
**In Scope:**
- Desktop application development
- Local ML model inference
- File format support (FASTQ, BAM, VCF)
- Basic data visualization
- PDF report generation
- Multi-language support

**Out of Scope (Current Phase):**
- Web application deployment
- Cloud infrastructure
- Mobile applications
- Real-time data streaming
- Advanced clinical decision support
- Regulatory approval processes

## Key Constraints
- **Privacy:** Zero external data transmission
- **Performance:** Must run on standard consumer hardware (single GPU optional)
- **Compliance:** GDPR and HIPAA compatible by design
- **Technology:** Open-source tools and frameworks only
- **Accessibility:** Support for low-resource computing environments

## Definition of Done
A fully functional desktop application that can:
1. Accept genomic file uploads via drag-and-drop
2. Process variants using local ML models
3. Generate interpretable risk reports
4. Export results as PDF documents
5. Operate completely offline
6. Install and run on major desktop platforms

This project brief serves as the immutable foundation for all development decisions and feature prioritization. 