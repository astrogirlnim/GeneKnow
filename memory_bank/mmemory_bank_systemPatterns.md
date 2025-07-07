# GenePredict System Patterns

## Architecture Overview

### Three-Layer Desktop Architecture
1. **Frontend Layer**: React + Tailwind CSS for user interface
2. **Native Layer**: Rust (Tauri) for secure file handling, system integration, and plugin orchestration
3. **AI/ML Layer**: Python for machine learning, genomic processing, and LLM integration

### Key Design Patterns

#### 1. Plugin-Based Architecture
- **Rust Trait System**: Defined plugin interfaces for extensible genomic processors
- **Python Interop**: Rust invokes Python ML workflows through secure plugin system
- **Modular Design**: Each genomic file type (FASTQ, BAM, VCF) has dedicated processor plugin

#### 2. Privacy-First Data Flow
- **Local-Only Processing**: All genomic data remains on user's filesystem
- **Temporary File Management**: Encrypted temp files automatically cleaned up after processing
- **No Network Calls**: TCGA reference data downloaded once, cached locally for offline analysis
- **Audit Logging**: Complete processing trail logged for transparency

#### 3. Async Processing Pipeline
```
File Upload → Rust File Validation → Python ML Plugin → Risk Calculation → Frontend Display
     ↓                ↓                     ↓               ↓              ↓
  Progress UI    File Type Detection   TCGA Reference   Confidence     Interactive
    Updates       & Metadata         Data Matching     Intervals        Graphs
```

## Technical Patterns

### 1. Error Handling Strategy
- **Rust Error Types**: Custom error enums for different failure modes
- **Python Exception Mapping**: Clean error translation between Rust and Python layers
- **User-Friendly Messages**: Technical errors converted to actionable user guidance
- **Graceful Degradation**: Partial results displayed when possible

### 2. Configuration Management
- **Environment Variables**: Local paths, cache settings, ML model parameters
- **Layered Config**: Default → Environment → User Preferences → Runtime Overrides
- **Secure Secrets**: API keys and sensitive config encrypted at rest

### 3. Data Validation Patterns
- **Schema Validation**: Pydantic models ensure data integrity across language boundaries
- **File Format Detection**: Magic number validation before expensive processing
- **Genomic Data Validation**: Variant format checking, reference genome matching

### 4. Performance Optimization
- **Lazy Loading**: TCGA reference data loaded on-demand per cancer type
- **Streaming Processing**: Large genomic files processed in chunks to manage memory
- **Parallel Processing**: CPU-intensive ML operations distributed across available cores
- **Smart Caching**: Computed risk scores cached with genomic file hashes

## Integration Patterns

### TCGA Data Integration
- **API Client Pattern**: Async HTTP client for GDC Data Portal with rate limiting
- **Local Data Mirror**: Reference datasets cached locally for offline analysis
- **Incremental Updates**: New TCGA releases integrated without disrupting existing cache
- **Data Versioning**: Track TCGA data release versions for reproducible results

### LLM Integration Pattern
- **Local Llama 3.1**: HuggingFace models run locally for report generation
- **Prompt Engineering**: Structured prompts ensure consistent, clinical-appropriate language
- **Multi-Language Support**: Language-specific prompts for localized report generation
- **Fallback Strategy**: English reports if target language model unavailable

### Frontend-Backend Communication
- **Tauri Commands**: Type-safe IPC between React frontend and Rust backend
- **Event Streaming**: Real-time progress updates during long-running genomic analysis
- **JSON Schema**: Standardized data formats ensure consistent frontend rendering

## Security Patterns

### File System Security
- **Sandboxed Processing**: All genomic files processed in isolated temporary directories
- **Automatic Cleanup**: Sensitive temp files securely deleted after processing
- **Permission Validation**: File access permissions checked before processing
- **Path Sanitization**: User file paths validated to prevent directory traversal

### Data Privacy Patterns
- **Zero Network Transmission**: Genetic data never sent over network
- **Local Encryption**: Temporary files encrypted with user-specific keys
- **Memory Safety**: Rust memory management prevents genetic data leaks
- **Audit Trail**: Complete processing history logged for compliance verification

## Scalability Patterns

### Resource Management
- **Memory Bounds**: Maximum memory usage capped based on system resources
- **Processing Queues**: Multiple genomic files queued for sequential processing
- **Background Tasks**: Non-critical operations (cache updates) run in background
- **Resource Monitoring**: Real-time tracking of CPU, memory, and disk usage

### Extensibility Patterns
- **Plugin Discovery**: Dynamic loading of new genomic processor plugins
- **Version Compatibility**: Backward compatibility maintained for plugin interfaces
- **Configuration Hot-Reload**: Settings changes applied without application restart 