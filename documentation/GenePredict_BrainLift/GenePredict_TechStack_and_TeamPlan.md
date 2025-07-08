
# 🧪 GenePredict – Tech Stack & Dev Assignment

## Recommended Tech Stack

### Platform
- **Tauri**: For building a secure, cross-platform desktop app (macOS, Windows, Ubuntu).
- **React + Tailwind CSS**: Modern UI with fast iteration and good dev velocity.
- **Rust (via Tauri)**: For secure, performant native functionality and file handling (FASTQ, BAM).
- **Python (via Tauri plugin)**: For invoking ML workflows, TensorFlow, and Llama 3.1.
- **No External Backend**: All inference runs locally for privacy and compliance.

### Genomic Processing & AI
- **TensorFlow**: Deep learning for genomic risk scoring.
- **OpenMined PySyft**: Differential privacy on ML models.
- **Llama 3.1 (via HuggingFace or local model runner)**: For generating readable reports.
- **Samtools / BioPython / pysam**: For FASTQ/BAM parsing and variant extraction.
- **DeepVariant**: For variant calling if needed.

---

## 👥 Dev Assignments by Feature

### Phase 1: Foundation ✅ **COMPLETED**
**Dev 2 (Full Stack)** – Setup and core architecture
- [x] Set up Tauri environment with React, Tailwind CSS.
- [x] Wire Rust backend to invoke Python ML functions.
- [x] Add plugin scaffolding for file parsing and ML execution.

**🎯 Phase 1 Achievements (January 2025):**
- ✅ **Cross-Platform Desktop App**: Tauri 2.6.2 + React 19.1.0 + Tailwind CSS 4.1
- ✅ **Rust-Python Integration**: Unified `execute_python()` helper with JSON data contracts
- ✅ **Development Workflow**: Hot reload, comprehensive testing, and build automation
- ✅ **Code Organization**: Consolidated Python ML scripts in `desktop/python_ml/`
- ✅ **Testing Infrastructure**: 6-level testing strategy with automated validation
- ✅ **Documentation**: Complete memory bank system and development guides
- ✅ **Quality Assurance**: 0 TypeScript errors, 0 ESLint issues, successful Rust compilation

### Phase 2: Data Layer 🚧 **IN PROGRESS**
**Dev 5 (Backend Expert)**
- [ ] Build Python/Rust interface for processing FASTQ/BAM files.
- [ ] Wrap DeepVariant / BioPython / pysam tools.
- [ ] Build TensorFlow pipeline for breast cancer prediction (initial model).
- [ ] Pipe results into JSON output usable by frontend.

**Dev 4 (Full Stack, Frontend Focus)**
- [ ] Design data schema and file manager interface in frontend.
- [ ] Create abstract components for processed sample display.
- [ ] Sync with Dev 5 on JSON data shape.

### Phase 3: Interface Layer
**Dev 3 (Frontend/UI/UX)**
- [ ] Build drag-and-drop uploader with file validation.
- [ ] Build variant/risk table view with filter/sort/search.
- [ ] Build loading states and interaction flows.

**Dev 1 (UX/Planning)**
- [ ] Define user flows and visual hierarchy.
- [ ] Design risk heatmap and report layout in Figma.
- [ ] Collaborate with Dev 3 on frontend animations & transitions.

### Phase 4: Reporting & Compliance
**Dev 2 (Full Stack)**
- [ ] Integrate Llama 3.1 to generate text summaries from output JSON.
- [ ] Support multilingual outputs (start with 3: English, Hindi, Spanish).
- [ ] Enable exporting reports to PDF.

**Dev 5 (Backend)**
- [ ] Integrate PySyft for privacy-preserving training (optional, flag for later).
- [ ] Ensure all data is processed locally (audit logs, file deletion routines).

### Phase 5: Explorer Mode
**Dev 4 (Full Stack, Frontend Focus)**
- [ ] Create interactive mode for hypothetical variant simulation.
- [ ] Visualize BRCA1/2 impact in chart form.
- [ ] Build frontend toggle to compare sample variants.

**Dev 1 (UX)**
- [ ] Design UX for "what-if" simulation and toggle interactions.

---

## 📈 Development Progress

**Overall Project Status**: 40% Complete
- **Phase 1 (Foundation)**: ✅ 100% Complete
- **Phase 2 (Data Layer)**: 🚧 Ready to Begin
- **Phase 3 (Interface Layer)**: 📋 Planned
- **Phase 4 (Reporting & Compliance)**: 📋 Planned
- **Phase 5 (Explorer Mode)**: 📋 Planned

**Key Technical Achievements:**
- **Robust Architecture**: Command, Strategy, Data Contract, and Observer patterns implemented
- **Cross-Platform Support**: macOS, Windows, and Linux compatibility verified
- **Performance Optimized**: 189KB bundle size, <500ms hot reload time
- **Security First**: All processing remains local, no external API calls
- **Developer Experience**: Comprehensive documentation and testing infrastructure

---

## Deliverable: Tech + Dev Planning Markdown

Each task is siloed to ensure frontend/backend pairing and minimize blockers.

**Next Priority**: Begin Phase 2 Data Layer development with TensorFlow integration and genomic file processing capabilities.
