# 🛠️ GenePredict – Phase 1 Foundation, Step 1: Tauri Environment Setup Plan

---

## 1. Verification of Current State
- 🔍 **Repository scan results**
  - No existing `src-tauri/`, `tauri.conf.json`, or React `package.json` found.
  - Documentation references Tauri but no implementation artefacts detected.
  - Conclusion: safe to initialise a fresh Tauri + React + Tailwind workspace (previous implementation fully removed).

## 2. High-Level Architecture After This Step
```
/desktop               # top-level folder for the cross-platform desktop client
  ├─ src-tauri/        # Rust side managed by Tauri
  │   ├─ tauri.conf.json
  │   └─ src/
  │       └─ main.rs   # Rust entry-point
  ├─ ui/               # React + Vite + Tailwind front-end
  │   ├─ index.html
  │   ├─ package.json
  │   ├─ tsconfig.json
  │   ├─ tailwind.config.ts
  │   └─ src/
  │       ├─ main.tsx  # React entry-point
  │       └─ ...
  ├─ .npmrc            # Enforces pnpm & exact versions
  ├─ .gitignore        # Node, Rust, OS and editor artefacts
  └─ README.md         # Dev setup & run instructions
```

## 3. Key Environment & Build Variables
| Variable | Purpose |
|----------|---------|
| `NODE_ENV` | React build mode (development/production) |
| `TAURI_DEV_URL` | URL that Tauri loads during dev (typically `http://localhost:5173`) |
| `TAURI_BACKEND_URL` | Future use: points Tauri to local Rust/Python ML service |
| `RUST_LOG` | Logging level for Rust side (`info`, `debug`, etc.) |
| `CI` | Continuous-integration flag for headless builds |

> 📝 No Firebase services are required at this stage (local-only app). If Crashlytics, Remote Config, or Auth become necessary later, revisit Tauri IPC permissions and add Google services via the [Firebase Web SDK for Tauri](https://tauri.app/v1/guides/integration/web). For now we keep the config minimal.

## 4. Detailed Implementation Checklist

- [x] **Phase 1 · Foundation**
  - [x] **Feature · Tauri Environment Setup**
    1. **Toolchain Preparation**
       - [x] Install **Rust stable** (`rustup default stable`) and `cargo`.
       - [x] Install **Node 20 LTS** + **pnpm 8** (`corepack enable && corepack prepare pnpm@8 --activate`).
       - [x] Add **Tauri CLI** globally: `cargo install tauri-cli`.
       - [x] Verify versions:
         - [x] `rustc --version` ≥ 1.77 (Current: 1.88.0 ✅)
         - [x] `node -v` ≥ 20 (Current: v20.19.2 ✅)
         - [x] `pnpm -v` ≥ 8 (Current: 10.12.1 ✅)
    2. **Scaffold Project Structure**
       - [x] Create top-level `desktop/` directory (keeps repo root tidy).
       - [x] Inside `desktop/`, initialise React app with Vite:
         ```bash
         pnpm create vite ui -- --template react-ts
         ```
       - [x] Add Tailwind CSS:
         ```bash
         cd ui && pnpm add -D tailwindcss postcss autoprefixer
         npx tailwindcss init -p
         ```
       - [x] Configure `tailwind.config.ts` with Tauri-friendly `content` globs (`../src-tauri/**/*.rs`, `./index.html`, `./src/**/*.{ts,tsx}`).
    3. **Initialise Tauri**
       - [x] From `desktop/`, run `cargo tauri init --directory ui --ci false` (places `src-tauri/` beside `ui/`).
       - [x] Move generated Rust code into `src-tauri/src/` if not already.
       - [x] Add **window hardening** to `tauri.conf.json`:
         ```json
         {
           "security":{ "csp":"default-src 'self'; img-src 'self' data:;" }
         }
         ```
       - [x] Configure **bundle identifiers** (`org.gene.predict.desktop`) and icons.
    4. **Cross-Silicon Build Scripts**
       - [x] Add `"dev": "tauri dev"` & `"build": "tauri build"` scripts to `ui/package.json`.
       - [x] Add top-level Makefile targets or `scripts/desktop.sh` for CI convenience.
    5. **Logging & Debugging Hooks**
       - [x] Enable Rust `tracing` crate with JSON layer (future-proofs ML logs).
       - [x] Add React custom hook `useLogger` wrapping `console.debug`.
       - [x] Ensure `RUST_LOG=debug` is injected in `pnpm dev` command.
    6. **Git Integration**
       - [x] Add `.gitignore` entries: `desktop/ui/node_modules`, `desktop/src-tauri/target`, `/dist`, etc.
       - [x] Commit with message `chore: scaffold tauri + react + tailwind environment` (no slashes).
    7. **Smoke Test**
       - [x] Run `pnpm dev` inside `desktop/ui` ➜ `tauri dev` should auto-open the window.
       - [x] Confirm hot-reload works for both React and Rust.
       - [x] Verify Tailwind classes render correctly.
       - [x] Package test macOS build: `tauri build --debug` ➜ `.dmg` opens successfully.

## 5. Future-Facing Notes (Out-of-Scope for Step 1)
- IPC layer (Rust ↔ Python via Tauri plugin) will be set up in **Phase 1 Step 2**.
- State management (TanStack Query, Zustand) and UI components come in **Phase 3**.
- Privacy features (OpenMined, federated learning) start in **Phase 4**.

---

> ✅ Deliverable: all boxes in the checklist ticked, a runnable cross-platform desktop shell exists, committed to git but **not pushed**.

---

"Act, we must. Procrastinate, we must not." – Yoda 