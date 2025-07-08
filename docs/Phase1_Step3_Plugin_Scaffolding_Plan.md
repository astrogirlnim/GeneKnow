# Phase 1 – Step 3: Plugin Scaffolding Implementation Plan

> Scope: Add dynamic plugin discovery, loading, and execution for genomic file parsing & ML workflows inside the cross-platform Tauri desktop application.

---

## ✅ Current State Verification

| Area | Observation |
|------|-------------|
| Rust backend (`desktop/src-tauri/src`) | `lib.rs` directly invokes hard-coded Python scripts via `utils::execute_python(…)`. No trait‐based abstraction or registry exists. |
| Utility layer (`utils.rs`) | Contains generic `execute_python` helper but **no plugin concept**. |
| Python ML scripts (`desktop/python_ml/*.py`) | Individual standalone scripts, no manifest or metadata. |
| Front-end (React) | Components call specific Tauri commands; no dynamic plugin UI. |
| Docs | `GenePredict_TechStack_and_TeamPlan.md` lists Step 3 as **incomplete**. |

Conclusion → **Plugin scaffolding has not been implemented yet.**

---

## 🚀 High-Level Goals

- Abstract *any* genomic/ML workflow into a **Plugin** implementing a common contract.
- Discover plugins automatically at runtime across macOS, Windows, Linux.
- Provide a **registry** exposed via Tauri commands & TypeScript typings for the UI.
- Allow per-plugin configuration via JSON (local first, Firebase Remote Config override optional).
- Maintain full local execution & privacy.

---

## 🗂️ Planned File/Module Changes

- **Rust**
  - `desktop/src-tauri/src/plugin.rs` (new)
  - `desktop/src-tauri/src/plugin_registry.rs` (new)
  - `desktop/src-tauri/src/lib.rs` (edit – wire registry & commands)
  - `desktop/src-tauri/src/utils.rs` (edit – generic helpers)
- **Python**
  - `desktop/python_ml/plugins/` (new folder for plugin entry scripts)
  - Each plugin has:
    - `<plugin_name>.py` (script)
    - `manifest.json`
- **TypeScript (React UI)**
  - `desktop/ui/src/api/plugins.ts` (new)
  - `desktop/ui/src/components/PluginRunner.tsx` (new)
- **Config & Docs**
  - `docs/Phase1_Step3_Plugin_Scaffolding_Plan.md` ← *this file*
  - Update `docs/file_structure.md` after implementation.

---

## 🧩 Core Architecture

1. **Trait Contract** (`plugin.rs`)
   ```rust
   pub trait GenomicPlugin {
       /// Unique identifier, eg. "fastq_to_vcf"
       fn id(&self) -> &'static str;
       /// Human-readable name
       fn name(&self) -> &'static str;
       /// Short description
       fn description(&self) -> &'static str;
       /// Execute the plugin with arbitrary args encoded as JSON
       fn run(&self, args: serde_json::Value) -> anyhow::Result<serde_json::Value>;
   }
   ```
   Key structs/enums:
   - `PluginManifest` (parsed from `manifest.json`)
   - `PluginError`

2. **Registry** (`plugin_registry.rs`)
   - Scans `desktop/python_ml/plugins/**/manifest.json` on app start.
   - Validates manifest → builds `Vec<Box<dyn GenomicPlugin>>`.
   - Exposes:
     - `fn list() -> Vec<PluginSummary>`
     - `fn run(id, args_json) -> Result<serde_json::Value>`

3. **Python Adapter**
   - Concrete struct `PythonScriptPlugin` that implements `GenomicPlugin` and calls `utils::execute_python` using manifest info.

4. **Tauri Commands** (add to `lib.rs`)
   - `list_plugins() -> Vec<PluginSummary>`
   - `run_plugin(id: String, args_json: String) -> Result<String>`

5. **Front-end Flow**
   - Fetch plugin list on load (`plugins.ts` → `invoke("list_plugins")`).
   - Display selectable plugins in UI (`PluginRunner.tsx`).
   - Send JSON-encoded arguments → `run_plugin`.

6. **Configuration Layer**
   - Default `manifest.json` per plugin holds config & schema.
   - Optional override: fetch remote config key `plugins/<id>` from Firebase RC when internet is available; merge with local manifest (keeping local execution only).

---

## 📋 Implementation Checklist (nested)

- [ ] **Phase 1 – Foundation**
  - [ ] **Step 3: Plugin Scaffolding**
    - [ ] Design `GenomicPlugin` trait & error types
    - [ ] Implement `PythonScriptPlugin` adapter
    - [ ] Build `plugin_registry.rs` with directory scan & caching
    - [ ] Expose `list_plugins` & `run_plugin` Tauri commands
    - [ ] Update `lib.rs` to initialize registry and register commands
    - [ ] Refactor existing hard-coded calls in `lib.rs` to use registry (backwards-compat alias)
    - [ ] Move existing Python scripts into `python_ml/plugins/<id>/`
    - [ ] Create `manifest.json` for each moved plugin (include schema)
    - [ ] Add TypeScript `plugins.ts` API hook
    - [ ] Build `PluginRunner.tsx` component (minimal MVP)
    - [ ] Update Tailwind styles as needed
    - [ ] Add extensive logging across Rust + Python
    - [ ] Update `docs/file_structure.md`
    - [ ] **Firebase**
      - [ ] Add optional Remote Config fetch util in Rust (respect offline mode)
      - [ ] Document Crashlytics integration path (future)
    - [ ] Manual cross-platform tests:
      - [ ] macOS
      - [ ] Windows
      - [ ] Linux
    - [ ] Update Memory Bank `activeContext` after completion

---

## 🔑 Key Variables & Functions Reference

| File | Identifier | Purpose |
|------|------------|---------|
| `plugin.rs` | `GenomicPlugin` | Core trait every plugin must implement |
| `` | `PluginManifest` | Struct representing `manifest.json` |
| `` | `PluginError` | Unified error enum |
| `plugin_registry.rs` | `PLUGIN_REGISTRY` | Lazy-static global registry |
| `` | `scan_plugins()` | Discovers & registers plugins |
| `PythonScriptPlugin` | `run()` | Bridges to `utils::execute_python` |
| `utils.rs` | `execute_python()` | Already existing helper, reused |
| `lib.rs` | `list_plugins`, `run_plugin` | Tauri commands exposed to UI |
| `plugins.ts` | `listPlugins()` | Front-end fetch util |
| `` | `runPlugin()` | Invoke plugin with args |

---

## 🌐 Firebase Considerations

1. **Remote Config**
   - Key namespace: `plugins/<plugin_id>` holds latest manifest JSON.
   - Fetch on app start (non-blocking); fallback to local manifest.

2. **Crashlytics / Analytics**
   - *Optional*: Send anonymized plugin execution metrics.
   - Respect user privacy toggle.

3. **Hosting / Update Distribution**
   - Plugin updates shipped via app auto-update; Firebase only supplies metadata.

---

## ⏭️ Next Steps

Upon approval of this document:
1. Create stubs for new Rust modules & move one existing script (`fastq_to_vcf_pipeline`) into plugin format as POC.
2. Iterate through checklist to MVP completion, committing after each logical chunk (no pushes, no slashes in messages).
3. Perform manual QA on all three OS targets.
4. Update Memory Bank & documentation. 