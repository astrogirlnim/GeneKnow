
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

### Phase 1: Foundation
**Dev 2 (Full Stack)** – Setup and core architecture
- [x] Set up Tauri environment with React, Tailwind CSS.
- [ ] Wire Rust backend to invoke Python ML functions.
- [ ] Add plugin scaffolding for file parsing and ML execution.

### Phase 2: Data Layer
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
- [ ] Design UX for “what-if” simulation and toggle interactions.

---

## Deliverable: Tech + Dev Planning Markdown

Each task is siloed to ensure frontend/backend pairing and minimize blockers.
