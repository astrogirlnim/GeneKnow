# 🧪 GenePredict Testing Guide

This guide explains how to test the Rust ⇄ Python ML integration implemented in Phase 1, Step 2.

## 🚀 Quick Test (Recommended)

Run the comprehensive test script from the `desktop/` directory:

```bash
cd desktop
./test_implementation.sh
```

This script tests:
- ✅ Python scripts with JSON output
- ✅ Rust compilation
- ✅ Frontend development server
- ✅ Integration architecture

## 🔍 Individual Component Testing

### 1. Python Scripts Testing

Test each Python script individually from the `desktop/` directory:

```bash
# Generate test FASTQ files
python3 python_ml/generate_test_fastq.py --output-dir test_output --num-reads 5 --json

# List available VCF files  
python3 python_ml/config_data_source.py --list-vcf-files --json

# Get data source info
python3 python_ml/config_data_source.py --get-data-info --json

# Test extract by region (requires BED file)
python3 python_ml/extract_by_region.py --help
```

### 2. Rust Compilation Testing

From the `desktop/src-tauri/` directory:

```bash
cd src-tauri

# Check compilation
cargo check

# Run with debug logging
RUST_LOG=debug cargo check

# Build optimized version
cargo build --release
```

### 3. Frontend Testing

From the `desktop/ui/` directory:

```bash
cd ui

# Install dependencies
pnpm install

# Start development server
pnpm dev

# Check for linting issues
pnpm run lint

# Build for production
pnpm build
```

## 🦀 Rust-Python Integration Testing

To test the actual Rust-Python integration:

### Option A: Through Tauri Commands

Start the application and test through the UI:

```bash
# Terminal 1: Start frontend
cd desktop/ui && pnpm dev

# Terminal 2: Start Tauri app  
cd desktop/src-tauri && cargo tauri dev --no-dev-server
```

### Option B: Direct Integration Test

Create a simple test binary to test the `execute_python` function:

```rust
// In src-tauri/src/bin/integration_test.rs
use app_lib::utils::execute_python;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let result = execute_python("config_data_source", &["--list-vcf-files", "--json"])?;
    
    if result.status.success() {
        println!("✅ Integration test passed!");
        let stdout = String::from_utf8_lossy(&result.stdout);
        if let Ok(json) = serde_json::from_str::<serde_json::Value>(&stdout) {
            println!("📊 JSON parsed successfully");
        }
    }
    
    Ok(())
}
```

Build and run:
```bash
cd src-tauri
cargo build --bin integration_test
./target/debug/integration_test
```

## 🎯 Expected Results

### Python Scripts
Each script should output valid JSON when called with `--json` flag:

```json
{
  "success": true,
  "data": "...",
  "error": null
}
```

### Rust Compilation
Should compile with only minor warnings about unused functions:

```
warning: function `execute_python_json` is never used
warning: function `normalize_path` is never used  
warning: function `escape_path_spaces` is never used
```

### Integration Test
Should output:
```
✅ Python script executed successfully
✅ JSON parsing successful - 23 VCF files listed
🎉 Rust-Python integration test PASSED!
```

## 🐛 Troubleshooting

### "Python script not found"
- Ensure you're running from the correct directory (`desktop/`)
- Check that Python scripts exist in `python_ml/` directory
- Verify script names don't include `.py` extension when calling `execute_python`

### "Command not found: tauri"
```bash
cd desktop/ui
pnpm install  # Reinstalls @tauri-apps/cli
```

### Port conflicts
If port 5173 is in use:
```bash
# Kill existing processes
pkill -f "vite|cargo|tauri"

# Or edit vite.config.ts to use different port
```

### JSON parsing errors
- Ensure Python scripts are called with `--json` flag
- Check that scripts output single-line JSON to stdout
- Verify no extra logging is printed to stdout in JSON mode

## 📊 Performance Benchmarks

Expected execution times on modern hardware:

| Test | Expected Time |
|------|---------------|
| Python script execution | < 1 second |
| Rust compilation (check) | < 5 seconds |
| Frontend dev server start | < 10 seconds |
| Full integration test | < 30 seconds |

## 🔄 Continuous Testing

For ongoing development, use these commands:

```bash
# Watch for Rust changes
cd desktop/src-tauri && cargo watch -x check

# Watch for frontend changes  
cd desktop/ui && pnpm dev

# Auto-test Python scripts on change
find python_ml -name "*.py" | entr python3 python_ml/config_data_source.py --json
```

## ✅ Test Checklist

Before submitting changes, verify:

- [ ] `./test_implementation.sh` passes completely
- [ ] All Python scripts work with `--json` flag
- [ ] Rust code compiles without errors
- [ ] Frontend development server starts successfully
- [ ] Integration test shows successful Rust → Python communication
- [ ] No temporary test files left behind
- [ ] README instructions are accurate and tested

---

**Last Updated**: January 2025  
**Implementation**: Phase 1, Step 2 - Rust ⇄ Python Integration 