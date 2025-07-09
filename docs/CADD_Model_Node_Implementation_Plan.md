# 🧬 CADD Risk Model Node – Implementation Plan

> **Context**  
> We are migrating GeneKnow’s single-node `risk_model` design to a five-model static-annotation architecture (PRS, ClinVar, CADD, TCGA Frequency Match, Gene/Pathway Burden) that feeds a TensorFlow **Risk-Fusion** model.  
> This document details the **CADD node** implementation (Phase 2) while outlining the broader refactor roadmap.  
> *Location:* `docs/CADD_Model_Node_Implementation_Plan.md`

---

## 0. Executive Summary
1. Introduce a dedicated LangGraph node `cadd_scoring` that enriches `filtered_variants` with pre-computed CADD PHRED scores.
2. Insert the node into the pipeline **after** `population_mapper` and **before** the new `feature_vector_builder` (placeholder) – temporarily bypassing the existing `risk_model` until the Risk-Fusion phase.
3. Retain the current `risk_model` file (deprecated flag) to preserve backward compatibility.
4. Ship a local SQLite cache of CADD scores for top ~30 M SNVs (hg38); implement on-demand Tabix fallback for uncached variants.
5. Provide exhaustive logging, configuration flags, and CLI hooks per project rules.

---

## 1. Current Relevant Files (No Duplicates)
| Purpose | File | Key Variables / Functions |
|---------|------|---------------------------|
| LangGraph orchestration | `geneknow_pipeline/graph.py` | `create_genomic_pipeline()`, `route_after_preprocess`, `merge_variants` |
| State definition | `geneknow_pipeline/state.py` | `GenomicState.filtered_variants`, `GenomicState.variant_count`, etc. |
| Population frequencies | `geneknow_pipeline/nodes/population_mapper.py` | `process()`, `query_population_database()` |
| Existing ML risk node | `geneknow_pipeline/nodes/risk_model.py` | `_simple_risk_calculation()`, `process()` |
| QC logic | `geneknow_pipeline/nodes/qc_filter.py` | `apply_qc_filters()` |
| Test data & utils | `geneknow_pipeline/test_*`, `geneknow_pipeline/test_data/` |  |

**Verified** via code search – no existing `cadd`-related files; safe to create new node `geneknow_pipeline/nodes/cadd_scoring.py`.

---

## 2. New Files / Edits
1. `geneknow_pipeline/nodes/cadd_scoring.py`  
   – Node implementation, PHRED thresholding, logging.
2. `geneknow_pipeline/data/cadd_scores.db` *(git-ignored large binary, downloaded script)*.
3. `geneknow_pipeline/scripts/fetch_cadd.sh` – One-shot downloader / converter to SQLite.
4. `geneknow_pipeline/nodes/feature_vector_builder.py` *(stub, returns passthrough until remaining models are ready).*  
5. **Edits**:
   - `geneknow_pipeline/graph.py` – insert conditional edge: `population_mapper → cadd_scoring → feature_vector_builder` (temporarily redirect existing edge to `risk_model`).
   - `geneknow_pipeline/state.py` – add fields: `cadd_enriched_variants`, `cadd_stats`.
   - `geneknow_pipeline/enhanced_api_server.py` – expose new pipeline node names in `/api/pipeline-info`.

---

## 3. Phase Breakdown & Checklists

### Phase 1 – Design & Data Prep
- [x] Finalise CADD version (v1.7  hg38) & licensing.
- [x] Write `scripts/fetch_cadd.sh` to download tsv.gz ➜ SQLite (`cadd_scores.db`).
- [x] Document schema: `CREATE TABLE cadd (chrom TEXT, pos INT, ref TEXT, alt TEXT, phred REAL, raw FLOAT, PRIMARY KEY(chrom,pos,ref,alt));`
- [x] Add GHA cache step for DB artefact (skip in repo).

### Phase 2 – Node Implementation (`cadd_scoring.py`)
- [x] Function `lookup_cadd_score(chrom,pos,ref,alt) -> float|None`  
      • First SQLite lookup  
      • Fallback: remote Tabix query (`https://api.caddscores...`)  
      • Cache miss logging.
- [x] `process(state)`:
  1. Iterate over `state['filtered_variants']`.
  2. Annotate each with `cadd_phred`, `cadd_raw`, `cadd_risk_weight` (`min-max` scaled 0-1).
  3. Build summary stats (`mean`, `max`, `variants_gt20`).
  4. Append to `state['completed_nodes']`.
  5. Write extensive console logs per variant.
- [x] Update `state` keys:
```python
state['cadd_stats'] = {
    'mean_phred': 12.6,
    'max_phred': 34.1,
    'variants_gt20': 3,
    'lookup_missing': 5
}
```

### Phase 3 – Pipeline Wiring
- [x] Modify `graph.py`:
```python
workflow.add_node('cadd_scoring', cadd_scoring.process)
workflow.add_edge('population_mapper','cadd_scoring')
workflow.add_edge('cadd_scoring','feature_vector_builder')
# Temporarily: feature_vector_builder -> formatter (bypass risk_model)
```
- [x] Feature-toggle: CLI flag `--legacy-risk` to re-enable old path.

### Phase 4 – Testing & Validation
- [x] Unit tests: `test_cadd_scoring.py` with mock variants.
- [x] Integration: run full pipeline on sample VCF, assert `cadd_phred` present.
- [x] Performance benchmark (<200 ms per 10 k lookups local DB).
- [x] Memory profile (<200 MB for SQLite pagecache).

### Phase 5 – Documentation & Logging
- [x] Update `README.md` architecture diagram.
- [x] Extend `docs/TESTING_GUIDE.md` with CADD scenarios.
- [x] Ensure **hundreds** of `logger.debug()` lines across node for traceability per user rules.

### Phase 6 – Risk-Fusion Integration (Future)
- Placeholder `feature_vector_builder.py` collects outputs from **all five** static nodes (only CADD initially) → saves `state['feature_vector']`.
- `risk_fusion_model.py` (TensorFlow) will consume this vector; existing `risk_model.py` remains but is flagged `@deprecated`.

---

## 4. Variable & Config Inventory
| Variable | Type | Source | Description |
|----------|------|--------|-------------|
| `cadd_phred` | float | CADD DB | PHRED-scaled deleteriousness score (higher = worse) |
| `cadd_raw`   | float | CADD DB | Raw CADD score |
| `cadd_risk_weight` | float 0-1 | Node calc | Normalised risk contribution |
| `CADD_DB_PATH` | str | `cadd_scoring` | Path to local SQLite file |
| `CADD_REMOTE_TABIX` | str | env | Remote Tabix endpoint for misses |

---

## 5. Code Architecture Considerations
1. **Placement**: New node lives in `geneknow_pipeline/nodes/` to match existing structure.
2. **Caching**: Bundle `.db` in user cache dir (`~/.geneknow/cadd`) to keep desktop bundle small.
3. **Concurrency**: Use connection pool or one connection per node invocation (LangGraph runs nodes sequentially by default).
4. **Error Handling**: Any lookup failure downgrades to `risk_weight = 0.2` + warning, does **not** fail pipeline.
5. **Config Flags** (`desktop/src-tauri/tauri.conf.json`): expose `CADD_CACHE_DIR`, `USE_REMOTE_CADD`.
6. **UI Surface**: Update React hooks to display **CADD max PHRED** in variant table tooltip.

---

## 6. Milestone Timeline *(t-shirt estimates)*
| Phase | Effort | Owner |
|-------|--------|-------|
| 1 – Data Prep | 1 d | Dev A |
| 2 – Node Impl | 2 d | Dev A |
| 3 – Wiring    | 0.5 d | Dev A |
| 4 – Testing   | 1 d | QA |
| 5 – Docs      | 0.5 d | Dev A |
| **Total**     | **5 d** |  |

---

## 7. Open Questions
1. Target genome build – confirm **hg38** across pipeline.
2. Do we need indel support immediately? (CADD v1.7 covers many indels but DB size ↑3×.)
3. Desired PHRED threshold for “deleterious” flag in UI (>20?).
4. Multi-threaded lookup vs. sequential (depends on Tabix fallback latency).

---

*Prepared by Cursor AI – v0.1* 