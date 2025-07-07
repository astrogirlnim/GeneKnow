
# 🧬 GenePredict – User Flow

This document outlines the user journey through the GenePredict desktop application. It defines each key interaction point and the transitions between core features. It is intended to guide both UI design and application architecture.

---

## 1. Application Launch

**Trigger**: User opens the GenePredict desktop app (Tauri-based).

**Outcome**:
- Main dashboard is displayed.
- User is presented with a CTA to upload genomic data (FASTQ/BAM).

---

## 2. File Upload and Preprocessing

**Trigger**: User drags and drops or selects a FASTQ or BAM file.

**Outcome**:
- Client-side validation of file type and structure.
- Parsing and normalization via Rust/Python backend.
- Once validated, user is taken to the processing state.

**Transitions to**:
- Sample processing view (Phase 2 backend functions triggered).

---

## 3. Variant Processing and Risk Scoring

**Trigger**: Successfully parsed genomic data.

**Outcome**:
- DeepVariant or equivalent pipeline processes raw data.
- TensorFlow model evaluates risk factors (initial focus: breast cancer).
- Output includes: top variant risks, severity scores, gene associations.

**Transitions to**:
- Display of processed data and risk summary.

---

## 4. Results Dashboard

**Trigger**: Successful model output.

**Outcome**:
- UI displays risk heatmap, variant table, and risk categorization.
- User can search, filter, and inspect variant-level details.

**User Options**:
- Generate Report (next step)
- Enter Risk Explorer Mode (optional)
- Return to upload screen

---

## 5. Report Generation

**Trigger**: User selects "Generate Report".

**Outcome**:
- Llama 3.1 generates a readable summary of variant impacts and risk explanation.
- User selects language (default: English).
- Report is rendered in UI.

**User Options**:
- Export report as PDF
- Regenerate (if desired)
- Return to dashboard

---

## 6. Risk Explorer Mode (Optional)

**Trigger**: User clicks "Explore Risk" or toggles simulation mode.

**Outcome**:
- User can simulate adding/removing variants.
- Risk prediction updates dynamically.
- UI displays delta in overall risk and explanatory notes.

**Transitions to**:
- Back to Results Dashboard
- Optional report generation from simulated state

---

## 7. Application Exit

**Trigger**: User closes the application.

**Outcome**:
- Optionally clear local cache (configurable).
- All processing completed locally — no data retention unless user exports.

---

## Future Expansion Hooks

> These flows are not implemented in MVP but are noted for architecture planning.

- Multi-disease risk support
- EHR Integration
- Clinician account management
- Feedback loop (clinician ratings, risk interpretation scoring)
