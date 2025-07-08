# GenePredict - System Patterns Memory Bank

## Core Architecture Pattern

### 🏗️ **Local-First Privacy Architecture**
```
┌─────────────────────────────────────────────────────────────────┐
│                    Desktop Application                           │
├─────────────────────────────────────────────────────────────────┤
│  React Frontend (TypeScript)                                   │
│  ├── Component Layer (UI/UX)                                   │
│  ├── State Management (React Hooks)                            │
│  ├── API Integration (Tauri Commands)                          │
│  └── Logging Infrastructure (useLogger)                        │
├─────────────────────────────────────────────────────────────────┤
│  Tauri Bridge Layer                                             │
│  ├── Command Interface (#[command] functions)                  │
│  ├── Security (CSP, Local-only)                                │
│  └── Cross-Platform Compatibility                              │
├─────────────────────────────────────────────────────────────────┤
│  Rust Backend (Native Performance)                             │
│  ├── Python Execution Engine (execute_python)                 │
│  ├── File System Operations                                    │
│  ├── JSON Data Contracts                                       │
│  └── Error Handling & Logging                                  │
├─────────────────────────────────────────────────────────────────┤
│  Python ML Layer                                               │
│  ├── Genomic Processing Scripts                                │
│  ├── Data Validation & Parsing                                 │
│  ├── ML Model Inference (Future)                               │
│  └── Standardized JSON Output                                  │
└─────────────────────────────────────────────────────────────────┘
```

## Key Design Patterns

### 1. **Command Pattern (Tauri Integration)**
```rust
#[command]
pub async fn convert_fastq_to_vcf(options: FastqToVcfOptions) -> Result<FastqToVcfResult, String> {
    // Unified execution pattern
    let result = execute_python("fastq_to_vcf_pipeline", &args)?;
    // Standardized JSON parsing
    parse_json_output(&result)
}
```

**Benefits:**
- Consistent error handling across all commands
- Centralized logging and debugging
- Type-safe data contracts
- Cross-platform compatibility

### 2. **Strategy Pattern (Python Script Execution)**
```rust
pub fn execute_python(script: &str, args: &[&str]) -> Result<Output, Box<dyn Error>> {
    // Cross-platform path discovery
    let project_root = find_project_root()?;
    let script_path = project_root.join("python_ml").join(format!("{}.py", script));
    
    // Unified execution strategy
    Command::new("python3")
        .arg(&script_path)
        .args(args)
        .current_dir(&project_root)
        .output()
}
```

**Benefits:**
- DRY principle: Single point of Python execution
- Consistent path handling across platforms
- Centralized error handling and logging
- Easy to extend with new scripts

### 3. **Data Contract Pattern (JSON Standardization)**
```python
# Python side
def main():
    result = process_genomic_data()
    if args.json:
        print(json.dumps(result))
    else:
        print_human_readable(result)
```

```rust
// Rust side
fn parse_json_output(output: &Output) -> Result<T, String> {
    let stdout = String::from_utf8_lossy(&output.stdout);
    serde_json::from_str(&stdout)
        .or_else(|_| legacy_parsing(&stdout))  // Fallback pattern
}
```

**Benefits:**
- Consistent data exchange format
- Backward compatibility through fallback
- Type safety with serde_json
- Easy debugging with human-readable fallback

### 4. **Observer Pattern (Logging Infrastructure)**
```typescript
// React side
const logger = useLogger();
logger.info("Processing genomic data...");

// Rust side
log::info!("Executing Python script: {} with args: {:?}", script, args);
```

**Benefits:**
- Comprehensive debugging capabilities
- Consistent log format across layers
- Real-time development feedback
- Production troubleshooting support

## Architecture Principles

### 🔒 **Security First**
- **Zero External Dependencies**: All processing happens locally
- **CSP Policies**: Strict content security policies in Tauri
- **Input Validation**: Comprehensive validation at all layers
- **Secure File Handling**: Rust-powered native file operations

### 🚀 **Performance Optimization**
- **Native Performance**: Rust backend for computational tasks
- **Efficient Bundling**: Vite with tree-shaking (189KB bundle)
- **Hot Reload**: Sub-500ms reload times in development
- **Memory Management**: Proper cleanup of temporary files

### 🌐 **Cross-Platform Compatibility**
- **Path Handling**: `std::path::PathBuf` for OS-specific paths
- **Process Execution**: Platform-agnostic Python subprocess management
- **Build System**: Consistent builds across macOS, Windows, Linux
- **UI Responsiveness**: Tailwind CSS responsive design

### 📊 **Data Integrity**
- **Type Safety**: TypeScript frontend + Rust backend
- **Schema Validation**: JSON schema validation at boundaries
- **Error Propagation**: Comprehensive error handling chain
- **Fallback Mechanisms**: Graceful degradation patterns

## Integration Patterns

### 🔄 **Frontend ↔ Backend Communication**
```
React Component → Tauri Command → Rust Function → Python Script → JSON Output → Rust Parsing → React State
```

### 🐍 **Python Script Integration**
```
execute_python("script_name", &["--arg1", "value1", "--json"])
↓
Project Root Discovery
↓
python3 /path/to/project/python_ml/script_name.py --arg1 value1 --json
↓
JSON Output Parsing
↓
Structured Data Return
```

### 📁 **File System Organization**
```
desktop/
├── src-tauri/          # Rust backend
│   ├── src/
│   │   ├── lib.rs      # Command implementations
│   │   ├── utils.rs    # Shared utilities
│   │   └── main.rs     # Entry point
├── ui/                 # React frontend
│   ├── src/
│   │   ├── components/ # UI components
│   │   ├── hooks/      # React hooks
│   │   └── api/        # API integrations
├── python_ml/          # Python ML scripts
│   ├── config_data_source.py
│   ├── fastq_to_vcf_pipeline.py
│   ├── extract_by_region.py
│   └── generate_test_fastq.py
└── testing/            # Test infrastructure
    ├── test_implementation.sh
    └── test_rust_integration.sh
```

## Quality Assurance Patterns

### 🧪 **Testing Strategy**
- **Unit Tests**: Rust functions with `cargo test`
- **Integration Tests**: Python-Rust communication validation
- **End-to-End Tests**: Full application workflow testing
- **Performance Tests**: Benchmarking of critical operations

### 📋 **Code Quality**
- **Linting**: ESLint + TypeScript for frontend
- **Formatting**: Rustfmt for backend code
- **Documentation**: Comprehensive inline documentation
- **Type Safety**: Strict TypeScript + Rust type checking

### 🔍 **Debugging Infrastructure**
- **Structured Logging**: Consistent log format across layers
- **Error Propagation**: Clear error messages at all levels
- **Development Tools**: Hot reload + DevTools integration
- **Testing Scripts**: Automated testing and validation

## Future Extension Patterns

### 🔌 **Plugin Architecture (Planned)**
- **ML Model Plugins**: Extensible ML model integration
- **File Format Plugins**: Support for additional genomic formats
- **Report Plugins**: Customizable report generation
- **UI Plugins**: Extensible visualization components

### 🌍 **Internationalization Pattern**
- **Multi-Language Support**: English, Hindi, Spanish
- **Locale-Aware Formatting**: Date, number, currency formatting
- **Cultural Adaptation**: Region-specific medical terminologies
- **Accessibility**: Screen reader and keyboard navigation support

This system architecture provides a robust foundation for privacy-first genomic analysis while maintaining extensibility for future enhancements. 