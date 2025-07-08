# 🛠️ Phase 1 – Foundation, Step 2: Rust ⇄ Python ML Wiring Implementation Plan

## 1. Verification Summary
- Inspected `desktop/src-tauri/src/lib.rs` and confirmed five `#[command]` functions already spawn **Python 3** subprocesses via `std::process::Command`.
- Verified main entry `desktop/src-tauri/src/main.rs` delegates to `app_lib::run()`; no Python logic there.
- Cargo manifest (`desktop/src-tauri/Cargo.toml`) currently has **no** `PyO3` / `tauri-plugin-python` dependency – integration is pure CLI-level.
- Located **duplicate** Python scripts:
  - Root-level: `fastq_to_vcf_pipeline.py`, `extract_by_region.py`, `config_data_source.py`
  - UI copy: `desktop/ui/src/api/…`
  - Referenced path in Rust: `src/api/<script>.py` (relative to runtime CWD) – **mismatch** may break execution.

## 2. Objective & Scope
Wire the Rust backend (Tauri) to reliably invoke Python-based ML workflows across macOS / Windows / Linux **without duplication** or brittle paths, and expose the results to the React UI.

## 3. Related Code & Key Symbols
```
- desktop/src-tauri/src/lib.rs
  • struct FastqToVcfOptions { reference, fastq1, fastq2, output_prefix, threads, aligner }
  • struct FastqToVcfResult  { success, vcf_file, bam_file, variant_count, error, execution_time }
  • struct ExtractRegionOptions { bed_file, output_dir, max_processes }
  • struct ExtractRegionResult  { success, total_regions, successful_extractions, total_variants, total_size, execution_time, results, error }
  • #[command] convert_fastq_to_vcf()
  • #[command] extract_genomic_regions()
  • #[command] check_genomic_dependencies()
  • #[command] get_available_vcf_files()
  • #[command] generate_test_fastq()

- Duplicate Python scripts
  • /fastq_to_vcf_pipeline.py
  • /extract_by_region.py
  • /config_data_source.py
  • /desktop/ui/src/api/(same names)
```

## 4. Proposed Architecture
```
+----------------------+            spawn             +---------------------------+
|  React Frontend      |  Tauri invoke("convert…") →  |  Rust backend (Tauri cmd) |
+----------------------+                             |  • execute_python() util   |
                                                     |  • JSON stdout / stderr    |
                                                     +-------------┬-------------+
                                                                   | subprocess
                                                                   V
                                                     +---------------------------+
                                                     |  Python ML Scripts        |
                                                     |  • fastq_to_vcf_pipeline  |
                                                     |  • extract_by_region      |
                                                     +---------------------------+
```

## 5. Comprehensive Checklist

### 🔹 Phase 1 – Foundation
- [ ] **Step 2**: Rust ⇄ Python wiring complete

#### 5.1 Python Script Consolidation
- [ ] **Identify** authoritative copies of each ML script
- [ ] **Create** single folder `desktop/python_ml/` (referenced in `docs/file_structure.md`)
- [ ] **Move** `fastq_to_vcf_pipeline.py`, `extract_by_region.py`, `config_data_source.py` into that folder
- [ ] **Delete** duplicate copies in `desktop/ui/src/api/` after verification

#### 5.2 Rust Command Abstraction
- [ ] Add `utils::execute_python(script: &str, args: &[&str]) -> Result<Output>` helper in Rust
- [ ] Update all `#[command]` functions to call the helper (DRY pattern)
- [ ] Resolve script path at runtime using `tauri::api::path::resource_dir()` for cross-platform correctness

#### 5.3 Data Contract Standardisation
- [ ] Ensure every Python script prints **single-line JSON** to stdout (`print(json.dumps(result))`)
- [ ] Update Rust to parse with `serde_json` instead of ad-hoc regex where applicable
- [ ] Document JSON schema in `docs/api_contracts.md`

#### 5.4 Environment & Packaging
- [ ] Ship **embedded** Python environment using [`python-launcher`](https://github.com/pyenv) or instruct users to install Python 3.10+
- [ ] Add pre-flight checker `check_python_available()` in Rust – returns version
- [ ] Prepare `Makefile` target `make venv` to set up local virtualenv with required wheels (`tensorflow`, `pysam`, etc.)

#### 5.5 Logging & Error Surfacing
- [ ] Pipe `stderr` of Python into Rust and forward via `tauri-plugin-log`
- [ ] Add contextual logs before & after each invocation (script, args, working dir)

#### 5.6 Cross-Platform Path Handling
- [ ] Convert incoming paths (`fastq1`, `bed_file`, etc.) to OS-specific separators with `std::path::PathBuf`
- [ ] Escape spaces on Windows when spawning subprocesses

#### 5.7 Tests & Validation
- [ ] **CLI smoke test**: `cargo test --package app_lib -- test_convert_fastq_to_vcf`
- [ ] **End-to-end**: UI button triggers pipeline, UI receives populated JSON
- [ ] **Performance** benchmark on 1000 read pairs – record execution time

#### 5.8 Documentation & Linting
- [ ] Update `docs/file_structure.md` with new `/desktop/python_ml/` directory
- [ ] Inline ///doc comments on every new Rust helper
- [ ] ESLint + `ruff` for Python scripts

### 🔸 Firebase Considerations (Forward-Looking)
- Although GenePredict is **local-first**, future optional modules may use:
  - Firebase **Remote Config** to toggle experimental ML pipelines (requires network opt-in)
  - Firebase **Crashlytics** for anonymous crash reports of the desktop app
- Action items for this step:
  - [ ] Keep config keys (`FIREBASE_API_KEY`, etc.) behind `TAURI_ALLOW_NETWORK=false` flag
  - [ ] Document potential integration points in `docs/firebase_strategy.md` (placeholder)

---

## 6. Time & Resource Estimate
| Task Group | Effort |
|-----------|--------|
| Consolidation & Paths | 0.5 day |
| Rust Abstractions & Refactor | 1 day |
| JSON Contract & Parsing | 0.5 day |
| Env / Packaging Scripts | 0.5 day |
| Logging & Cross-platform quirks | 0.5 day |
| Tests & Benchmarks | 1 day |
| Docs & Cleanup | 0.5 day |
| **Total** | **~4.5 developer-days** |

> Once every item above is checked, Step 2 of the Foundation phase is **DONE** ✅ 