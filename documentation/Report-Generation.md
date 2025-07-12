GeneKnow Report Generator Module

Purpose

The GeneKnow Report Generator module transforms structured variant-level JSON output from the local genomics pipeline into a professional, readable clinical report tailored for clinicians, researchers, and medically trained users. It leverages locally available LLMs (via Ollama) to craft a well-structured scientific narrative that remains understandable.

Design Goals
	•	Fully offline: All operations must work without an internet connection.
	•	Modular LLM Backend: Use local models from Ollama as the LLM writing assistant.
	•	Scientific Language: Use appropriate medical/genomic terminology.
	•	Readable & Trustworthy: Prioritize clarity, accuracy, and structure.

Inputs

Structured JSON with fields including:
	•	case_id
	•	filtered_variants: Array of variant objects with annotations
	•	model_outputs: Results from PRS, CADD, TCGA, ClinVar, etc.
	•	cancer_risk_scores: Dictionary with per-cancer risk percentages

Outputs
	•	Markdown (.md) report (intermediate)
	•	Optionally exported as:
	•	PDF (via Pandoc or html/pdf generator)
	•	Plaintext (.txt)
	•	HTML (.html)

Major Components

1. Prompt Builder (prompt_builder.py)
	•	Parses JSON inputs and prepares natural language descriptions
	•	Assembles contextual prompt for the LLM, including model outputs
	•	May use structured prompt templates

2. Model Interface (model_interface.py)

Supports:
	•	ollama backend (e.g., llama3, mistral, codellama)
	•	Selects backend from config

3. Formatter (formatter.py)
	•	Takes raw JSON + LLM narrative and outputs:
	•	Markdown file (base template)
	•	Calls optional export functions (pdf, html, txt)

4. Report Templates
	•	Markdown template with sections:
	•	Case Overview
	•	Key Variants
	•	Per-Cancer Risk Summary
	•	Interpretation Summary
	•	Technical Appendix

User Configuration
	•	Configuration file (config.yaml) lets user set:
	•	Preferred LLM backend
	•	Report style (technical, clinician-focused)
	•	Output format(s)

Model Usage
	•	LLM is used to:
	•	Interpret structured variant data into coherent narratives
	•	Generate readable summaries of per-model outputs
	•	Translate risk scores and statistics into clinical language

Fallback Mode

If no LLM is available, the system can generate a non-AI templated report using pre-defined phrases and logic.

Example Output Sections

Patient Case: 2025A-13892
--------------------------------------
Top Variant: BRCA1 (c.5266dupC, p.Gln1756Profs*74)
Consequence: Frameshift, Pathogenic (ClinVar)
PRS Score (Colon): 0.84 (High Risk)
TCGA Match (BRCA1): 15.0% tumor freq, 0.1% normal (Enrichment: 150x)

Interpretation:
The presence of a known pathogenic BRCA1 variant with high TCGA enrichment and an elevated PRS score suggests a strong predisposition to breast and colorectal cancers. Clinical validation is recommended.

Future Additions
	•	Allow user notes and comments to be inserted in appendices
	•	Add graph/heatmap image export support for visual summaries
	•	Plugin-style support for downstream EHR export or FHIR bundles

## Implementation Phases

We'll break this into iterative phases to keep things manageable. Each phase includes deliverables, estimated steps, and any tool usage (e.g., code edits).

- **Phase 1: Module Setup and LLM Interface (Backend Foundation)**
  - Create the new `report_generator` directory under `geneknow_pipeline/nodes/` with core files: `model_interface.py` (for LLM detection/backends/fallback), `prompt_builder.py` (for assembling prompts from JSON), `formatter.py` (for Markdown output with glossary), and a main `report_generator.py` (entry point to replace `report_writer.process`).
  - Implement LLM auto-detection: Check for Ollama (e.g., via API ping). Default to fallback if none found, with a clear dev indicator in the output.
  - Add basic config handling (e.g., read from a `config.yaml` in the project root, with defaults for style="clinician", output_formats=["markdown"]).
  - Deliverable: Working module that can take sample JSON, detect LLM, and generate a basic templated report (non-LLM fallback) with a glossary stub.

- **Phase 2: Prompt Building and Narrative Generation**
  - Flesh out `prompt_builder.py` to parse the input JSON (focusing on >5% risk items) and create structured prompts (e.g., "Generate a scientific narrative for this high-risk variant: [data] using clinician terminology.").
  - Implement generation in `model_interface.py`: Use LLM if available to create narratives for sections (e.g., Case Overview, Key Variants, Risk Summary, Interpretation). Add glossary generation (e.g., scan for terms like "Frameshift" and append definitions).
  - Handle streaming: Generate in chunks and yield them (to be hooked into the API later).
  - Deliverable: Module can produce full LLM-enhanced Markdown reports from JSON, with fallback mode clearly marked.

- **Phase 3: Pipeline Integration and Replacement**
  - Replace the `report_writer` node in `geneknow_pipeline/graph.py`: Update the graph to call the new `report_generator.process` instead, passing the JSON (e.g., `state["report_sections"]` or `state["structured_json"]`).
  - Modify `enhanced_api_server.py` to send the original JSON to the frontend immediately, then trigger report generation in parallel (e.g., via a background task). Add WebSocket endpoint for streaming report chunks to the frontend.
  - Ensure the original JSON is unchanged and sent first; attach report paths/streams as extras (e.g., `{"original_json": ..., "report_stream": websocket_id}`).
  - Deliverable: Pipeline runs end-to-end, generating and streaming the new report without disrupting existing frontend data.

- **Phase 4: Frontend Updates**
  - In `desktop/ui/src/`, add a new tab/component (e.g., `ReportPage.tsx`) in the dashboard to display the streamed Markdown report (using a Markdown renderer like `react-markdown`).
  - Add a settings page (e.g., `SettingsPage.tsx`) with a button/form for LLM config (e.g., dropdowns for backend/model, saved to local storage or a config file).
  - Implement PDF download button (e.g., using `jsPDF` or a Tauri plugin to convert Markdown to PDF).
  - Deliverable: Users can view/stream the report in a new tab and configure LLMs via settings.

- **Phase 5: Testing and Polish**
  - Add unit tests in `geneknow_pipeline/test_report_generator.py` (e.g., test LLM detection, prompt building, fallback mode).
  - End-to-end tests: Run pipeline with your sample files, verify >5% risk filtering, glossary presence, and streaming.
  - Polish: Add dev indicators (e.g., "NON-LLM MODE" header), error handling, and documentation updates.
  - Deliverable: Fully tested module, ready for your manual verification with known files.