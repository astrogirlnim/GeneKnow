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

# 🏗️ System Patterns

## Architecture Overview
GeneKnow follows a **LangGraph-based pipeline architecture** for genomic analysis, with a desktop application frontend (Tauri + React) connected to a Python ML backend.

### Core Patterns

#### 1. **LangGraph Pipeline Pattern**
```python
# Each node follows this pattern:
def process(state: Dict[str, Any]) -> Dict[str, Any]:
    logger.info(f"Starting {node_name}")
    state["current_node"] = node_name
    
    try:
        # Process logic
        state["result_key"] = process_data(state["input_key"])
        state["completed_nodes"].append(node_name)
    except Exception as e:
        state["errors"].append({
            "node": node_name,
            "error": str(e),
            "timestamp": datetime.now()
        })
    
    return state
```

#### 2. **Unified Database Architecture**
The system uses a single SQLite database (`population_variants.db`) with multiple related tables:

```sql
-- Main population data
CREATE TABLE population_variants (
    chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    gene TEXT, gnomad_af REAL, clinical_significance TEXT,
    is_pathogenic INTEGER, consequence TEXT, review_status TEXT,
    PRIMARY KEY (chrom, pos, ref, alt)
);

-- CADD scores with foreign key relationship
CREATE TABLE cadd_scores (
    chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    raw_score REAL, phred_score REAL,
    job_id TEXT, source TEXT, created_at TIMESTAMP,
    PRIMARY KEY (chrom, pos, ref, alt),
    FOREIGN KEY (chrom, pos, ref, alt) REFERENCES population_variants
);

-- Job tracking for audit trail
CREATE TABLE cadd_jobs (
    job_id TEXT PRIMARY KEY, job_type TEXT,
    started_at TIMESTAMP, completed_at TIMESTAMP,
    variant_count INTEGER, status TEXT, metadata TEXT
);

-- Convenient view joining all data
CREATE VIEW variant_annotations AS
SELECT pv.*, cs.raw_score as cadd_raw, cs.phred_score as cadd_phred
FROM population_variants pv
LEFT JOIN cadd_scores cs USING (chrom, pos, ref, alt);
```

**Benefits:**
- Single database file for all variant annotations
- Foreign key constraints ensure data integrity
- Job tracking provides audit trail
- Views simplify complex queries
- Easy to backup and share

#### 3. **Tauri Desktop Integration**
```typescript
// Frontend calls Tauri commands
const result = await invoke<PipelineResult>('run_genomic_pipeline', {
    filePath: selectedFile,
    options: pipelineOptions
});

// Backend processes through Python
#[tauri::command]
async fn run_genomic_pipeline(file_path: String, options: Value) -> Result<Value> {
    let output = Command::new("python")
        .arg("geneknow_pipeline/enhanced_api_server.py")
        .arg("--file").arg(file_path)
        .output()?;
    // ...
}
```

#### 4. **Error Handling Pattern**
Every component follows fail-safe principles:
- Errors are logged but don't crash the pipeline
- Missing data gets default values
- Each node can recover from previous failures

#### 5. **Risk Assessment Architecture**
Moving from single model to five-model ensemble:
```
Current: filtered_variants → risk_model → risk_scores
Future:  filtered_variants → [CADD, PRS, ClinVar, TCGA, Pathway] → feature_vector → risk_fusion → risk_scores
```

## File Organization

### Pipeline Components
```
geneknow_pipeline/
├── nodes/              # LangGraph nodes (one per processing step)
├── models/             # ML models and configs
├── data/               # Static data files
├── scripts/            # Setup and utility scripts  
├── test_*.py          # Test files for each component
└── enhanced_api_server.py  # Main API server
```

### Desktop Application
```
desktop/
├── src-tauri/         # Rust backend
├── ui/                # React frontend
└── python_ml/         # Python ML integration
```

## Data Flow Patterns

### 1. **Variant Processing Flow**
```
FASTQ/BAM → Alignment → Variant Calling → VCF → MAF → 
Population Mapping → CADD Scoring → Risk Assessment → Report
```

### 2. **State Management**
- Single `GenomicState` object flows through pipeline
- Each node adds its results to state
- Immutable updates (copy state, modify, return)

### 3. **Database Query Pattern**
```python
# Normalize chromosome format
norm_chrom = chrom.replace("chr", "")

# Try exact match first
cursor.execute("""
    SELECT * FROM variant_annotations 
    WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
""", (norm_chrom, pos, ref, alt))

# Fall back to position-based lookup
if not result:
    cursor.execute("""
        SELECT * FROM variant_annotations
        WHERE chrom = ? AND pos = ?
    """, (norm_chrom, pos))
```

## Configuration Management

### Environment Variables
```bash
# Database paths
POPULATION_DB_PATH=/path/to/population_variants.db

# Feature flags
USE_LEGACY_RISK=false
USE_REMOTE_CADD=true

# API configuration
GENEKNOW_API_PORT=8000
```

### Pipeline Configuration
- Node-specific configs in each module
- Global configs in environment
- Override capability for testing

## Testing Patterns

### 1. **Unit Tests**
Each node has dedicated test file:
```python
class TestNodeName(unittest.TestCase):
    def setUp(self):
        # Create test data
        
    def test_process(self):
        # Test main functionality
        
    def test_error_handling(self):
        # Test failure cases
```

### 2. **Integration Tests**
Full pipeline tests with real data:
```python
# Run complete pipeline
state = run_pipeline(test_file)
assert state["pipeline_status"] == "completed"
assert len(state["errors"]) == 0
```

### 3. **Database Tests**
Test with temporary database:
```python
def setUp(self):
    self.temp_db = tempfile.NamedTemporaryFile()
    # Create test schema and data
```

## Performance Patterns

### 1. **Database Optimization**
- Indexes on all lookup columns
- Prepared statements for repeated queries
- Connection pooling for concurrent access

### 2. **Memory Management**
- Stream large files instead of loading
- Process variants in batches
- Clear intermediate results

### 3. **Parallel Processing**
- Independent pipeline paths run in parallel
- Batch database operations
- Async I/O where possible

## Security Patterns

### 1. **Data Isolation**
- Each pipeline run gets unique job ID
- No cross-contamination between runs
- Audit trail for all operations

### 2. **Input Validation**
- File format validation before processing
- Bounds checking on all numeric inputs
- Sanitization of user-provided strings

### 3. **Error Information**
- Log detailed errors internally
- Return sanitized errors to users
- No sensitive paths in error messages

## Deployment Patterns

### 1. **Desktop Distribution**
- Tauri bundles for each platform
- Python environment embedded
- Auto-update capability

### 2. **Server Deployment**
- Docker containers for consistency
- Environment-based configuration
- Health check endpoints

### 3. **Database Distribution**
- Compressed database for downloads
- Incremental update scripts
- Version tracking in metadata

## Evolution Patterns

### 1. **Adding New Nodes**
1. Create node file in `nodes/`
2. Implement `process()` function
3. Add to pipeline graph
4. Create tests
5. Update documentation

### 2. **Adding Database Tables**
1. Add table creation to scripts
2. Update foreign key relationships
3. Create/update views
4. Add migration logic
5. Update documentation

### 3. **Feature Toggles**
- Environment variables for new features
- Backward compatibility maintained
- Gradual rollout capability 