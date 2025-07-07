
# 🧬 GenePredict – User Workflows by Persona

This document defines end-to-end user workflows aligned with key user types. Each flow is based on real-world tasks, limited to the scope of features in the current MVP, and designed to ensure intuitive, complete interactions through the app.

---

## 👩‍⚕️ User Flow 1: Clinician – Single Patient Risk Assessment

**Goal**: Quickly assess a patient’s risk for breast cancer based on a newly received BAM file.

1. Open GenePredict on local desktop machine.
2. Upload patient’s BAM file via drag-and-drop or file picker.
3. File is validated and normalized automatically.
4. Processing begins – DeepVariant generates variants.
5. TensorFlow model evaluates risk and returns variant-level scores.
6. Risk dashboard displays:
   - Risk heatmap by gene
   - Top variants table
   - Severity filters
7. Clinician clicks “Generate Report”.
8. Llama 3.1 generates a risk summary in plain English.
9. Clinician downloads PDF and attaches to EHR manually.

---

## 🧑‍🔬 User Flow 2: Researcher – Compare Sample Variants

**Goal**: Analyze and compare results across two sample datasets.

1. Open GenePredict and upload two BAM or FASTQ files.
2. Each file is processed independently and scored.
3. Results appear in separate tabs or split view.
4. Researcher toggles between heatmaps and variant tables.
5. Notes differences in BRCA1 mutation between samples.
6. Opens “Risk Explorer” to simulate variant swaps.
7. Generates comparative report showing relative risk delta.

---

## 🧑‍💻 User Flow 3: Data Scientist – Validate Pipeline Locally

**Goal**: Run the full risk assessment process locally and confirm output structure.

1. Opens the application in test mode (CLI or GUI).
2. Uploads test VCF, BAM, or FASTQ from local dataset.
3. Watches console or log overlay for processing state.
4. Validates JSON output matches expected schema.
5. Optionally skips report generation.
6. Uses outputs for integration into external notebook/pipeline.

---

## 🌍 User Flow 4: Global Health Worker – Multilingual Report for Patient

**Goal**: Use local-only app to screen a rural patient and explain risk in native language.

1. Launches GenePredict in Spanish mode.
2. Uploads FASTQ sample from on-site genomic testing.
3. Processes file and views simplified dashboard.
4. Clicks “Generate Report”.
5. Llama 3.1 creates plain-language report in Spanish.
6. PDF report printed and explained to patient during visit.

---

## 🧪 User Flow 5: Developer – Try “Risk Explorer” for Variant Simulation

**Goal**: Test how changes in variant presence affect predicted risk.

1. Loads the app and uploads a BAM file.
2. Views baseline risk score and variant table.
3. Enters “Risk Explorer” mode.
4. Adds/removes synthetic variants manually.
5. Watches heatmap and score update dynamically.
6. Takes screenshot or exports simulation report.

---

## 🧑‍⚕️ User Flow 6: Medical Resident – Education & Interpretation

**Goal**: Learn to interpret genomic reports using guided interface.

1. Launches app and loads sample file included with GenePredict.
2. Walks through dashboard explanations:
   - Hover tooltips on gene regions
   - Info popups for variant types
3. Generates report to understand AI’s summary logic.
4. Cross-references with literature or ClinVar via tooltips.
5. Exports report and shares feedback for future iterations.

---

## ✅ Notes

- All workflows respect local-only processing and do not rely on server roundtrips.
- All roles access the same interface but may use features differently.
- No authentication or cloud syncing needed in MVP.
