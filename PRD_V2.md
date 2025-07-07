√### MVP Development Checklist — **GenePredict**

Template Overview
This checklist structures development into four dependency-ordered phases. Each sub-feature is a single, independently testable task for an AI developer.

---

#### Phases Overview

* [ ] **Phase 1: Foundation**
* [ ] **Phase 2: Data Layer**
* [ ] **Phase 3: Interface Layer**
* [ ] **Phase 4: Implementation Layer**

---

### Phase 1 — Foundation

*Criteria: essential scaffolding; all sub-features have **zero dependencies**.*

* [ ] **Project Scaffolding & Tooling**

  * [ ] Initialize Git repository with modular backend + frontend folders (independent — no dependencies)
  * [ ] Configure Python environment with Poetry/Pipenv (independent — no dependencies)
  * [ ] Bootstrap Node/React workspace with Vite + Tailwind (independent — no dependencies)
  * [ ] Add Docker Compose for backend, frontend, and model services (independent — no dependencies)

* [ ] **Core Libraries & Datasets**

  * [ ] Add TensorFlow 2.x and install PySyft for differential privacy (independent — no dependencies)
  * [ ] Pull Llama 3.1 405B weights via Hugging Face cache (independent — no dependencies)
  * [ ] Download 1000 Genomes & ClinVar reference data for local use (independent — no dependencies)
  * [ ] Create environment-variable config for local data paths (independent — no dependencies)

* [ ] **CI/CD & Testing Baseline**

  * [ ] Set up GitHub Actions pipeline for backend unit tests (independent — no dependencies)
  * [ ] Add frontend CI with ESLint, Jest, React Testing Library (independent — no dependencies)
  * [ ] Configure pre-commit hooks (Black, isort, Prettier) (independent — no dependencies)
  * [ ] Integrate Codecov for coverage reports (independent — no dependencies)

* [ ] **Security & Compliance Baseline**

  * [ ] Commit GDPR/HIPAA compliance checklist document to repo (independent — no dependencies)
  * [ ] Lock Docker network to local-only by default (independent — no dependencies)
  * [ ] Enable at-rest encryption for temporary genomic files (independent — no dependencies)
  * [ ] Add OSS license and CONTRIBUTING guidelines (independent — no dependencies)

---

### Phase 2 — Data Layer

*Criteria: data ingestion, storage, and model infrastructure; each sub-feature **depends on Phase 1**.*

* [ ] **VCF Upload & Parsing Engine**

  * [ ] Build client-side validator for VCF syntax (depends on Phase 1)
  * [ ] Implement backend VCF parser with cyvcf2 (depends on Phase 1)
  * [ ] Store parsed variants in an in-memory object model (depends on Phase 1)
  * [ ] Write unit tests for malformed VCF edge cases (depends on Phase 1)

* [ ] **Risk Prediction Pipeline Infrastructure**

  * [ ] Define TensorFlow model architecture for variant risk scoring (depends on Phase 1)
  * [ ] Provide training script using 1000 Genomes + ClinVar (depends on Phase 1)
  * [ ] Expose inference service endpoint (depends on Phase 1)
  * [ ] Implement batch processing for multiple samples (depends on Phase 1)

* [ ] **Variant Impact Analysis Module**

  * [ ] Map variants to genes via offline Ensembl/ANNOVAR DBs (depends on Phase 1)
  * [ ] Compute severity scores (CADD/PolyPhen) offline (depends on Phase 1)
  * [ ] Aggregate per-gene risk metrics into JSON (depends on Phase 1)
  * [ ] Serialize results for frontend consumption (depends on Phase 1)

* [ ] **Privacy Enforcement Layer**

  * [ ] Wrap model training with PySyft differential-privacy primitives (depends on Phase 1)
  * [ ] Encrypt all transient genomic data on disk (depends on Phase 1)
  * [ ] Add offline-mode unit tests (no outbound HTTP) (depends on Phase 1)
  * [ ] Implement log-scrubbing middleware (depends on Phase 1)

---

### Phase 3 — Interface Layer

*Criteria: user-facing components; each sub-feature **depends on Phase 1 & 2**.*

* [ ] **Global UI Framework**

  * [ ] Build React app shell with Tailwind & shadcn/ui (depends on Ph 1 & 2)
  * [ ] Implement dark/light mode toggle stored in localStorage (depends on Ph 1 & 2)
  * [ ] Apply clinician-friendly font & WCAG-AA contrast (depends on Ph 1 & 2)
  * [ ] Audit ARIA roles for full accessibility (depends on Ph 1 & 2)

* [ ] **File Upload Flow**

  * [ ] Add drag-and-drop VCF uploader with progress indicator (depends on Ph 1 & 2)
  * [ ] Display validation results to user (depends on Ph 1 & 2)
  * [ ] Trigger backend parse endpoint on success (depends on Ph 1 & 2)
  * [ ] Show animated processing status (depends on Ph 1 & 2)

* [ ] **Risk Visualization Components**

  * [ ] Render risk heatmap per gene/variant (depends on Ph 1 & 2)
  * [ ] Build searchable, filterable variant table (depends on Ph 1 & 2)
  * [ ] Add severity score bar chart component (depends on Ph 1 & 2)
  * [ ] Ensure responsive layout for tablet view (depends on Ph 1 & 2)

* [ ] **Risk Explorer Mode**

  * [ ] Create dashboard page with sample-profile selector (depends on Ph 1 & 2)
  * [ ] Enable zoom/pan interaction on variant plots (depends on Ph 1 & 2)
  * [ ] Persist explorer state via URL params (depends on Ph 1 & 2)
  * [ ] Add tooltips explaining metrics (depends on Ph 1 & 2)

---

### Phase 4 — Implementation Layer

*Criteria: core value delivery; each sub-feature **depends on Phases 1, 2 & 3**.*

* [ ] **Interpretable Report Generation Service**

  * [ ] Craft Llama 3.1 prompt to output clinician-level narrative (depends on Ph 1–3)
  * [ ] Generate layperson summary section (depends on Ph 1–3)
  * [ ] Implement multilingual output (≥20 languages) (depends on Ph 1–3)
  * [ ] Export PDF via WeasyPrint with embedded charts (depends on Ph 1–3)

* [ ] **Report Delivery & Export**

  * [ ] Provide backend endpoint bundling JSON + PDF (depends on Ph 1–3)
  * [ ] Add “Download report” button with loading state (depends on Ph 1–3)
  * [ ] Apply print-ready CSS for hard-copy fidelity (depends on Ph 1–3)
  * [ ] Cache last report in IndexedDB for offline access (depends on Ph 1–3)

* [ ] **Developer & Research Extensibility**

  * [ ] Publish REST/OpenAPI spec for all endpoints (depends on Ph 1–3)
  * [ ] Create plugin interface to swap models/datasets (depends on Ph 1–3)
  * [ ] Add example Jupyter notebooks for custom analyses (depends on Ph 1–3)
  * [ ] Release core packages to PyPI & NPM (depends on Ph 1–3)

* [ ] **Demo Fictional Patient Journey**

  * [ ] Generate synthetic VCF mimicking real-world variants (depends on Ph 1–3)
  * [ ] Pre-compute risk predictions for demo (depends on Ph 1–3)
  * [ ] Build guided walkthrough page highlighting results (depends on Ph 1–3)
  * [ ] Draft screencast script for future promo (depends on Ph 1–3)

---

**Feature Independence Rules Recap**

* Each sub-feature can be implemented, tested, and rolled back in isolation.
* Dependencies flow strictly from Phase 1 → 2 → 3 → 4.

Use the status marks as you progress:
`[ ]` Not Started `[~]` In Progress `[!]` Blocked `[x]` Completed
