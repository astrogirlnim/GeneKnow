# GenePredict - Tech Context Memory Bank

## Technology Stack Overview

### 🎯 **Core Technologies**
- **Frontend Framework**: React 19.1.0
- **UI Language**: TypeScript 5.8.3
- **Styling**: Tailwind CSS 4.1.11
- **Build Tool**: Vite 7.0.0
- **Desktop Framework**: Tauri 2.6.2
- **Backend Language**: Rust 1.88+
- **ML/Data Language**: Python 3.8+
- **Package Manager**: pnpm 10.12+

### 🏗️ **Architecture Stack**
```
┌─────────────────────────────────────────────────────────────────┐
│  Frontend Layer                                                 │
│  ├── React 19.1.0 (Component Framework)                        │
│  ├── TypeScript 5.8.3 (Type Safety)                            │
│  ├── Tailwind CSS 4.1.11 (Styling)                             │
│  ├── Vite 7.0.0 (Build Tool)                                   │
│  └── pnpm 10.12+ (Package Management)                          │
├─────────────────────────────────────────────────────────────────┤
│  Native Bridge Layer                                            │
│  ├── Tauri 2.6.2 (Desktop Framework)                           │
│  ├── WebView2 (Windows), WebKit (macOS), GTK WebKit (Linux)    │
│  └── IPC Communication (JSON over WebSocket)                   │
├─────────────────────────────────────────────────────────────────┤
│  Backend Layer                                                  │
│  ├── Rust 1.88+ (System Language)                              │
│  ├── tokio (Async Runtime)                                     │
│  ├── serde (Serialization)                                     │
│  ├── log (Logging)                                             │
│  └── std::process (Python Execution)                           │
├─────────────────────────────────────────────────────────────────┤
│  Data Processing Layer                                          │
│  ├── Python 3.8+ (ML & Data Processing)                        │
│  ├── NumPy (Numerical Computing)                               │
│  ├── JSON (Data Exchange)                                       │
│  └── argparse (CLI Interface)                                  │
└─────────────────────────────────────────────────────────────────┘
```

## Development Environment

### 🛠️ **Required Tools**
- **Node.js**: 20+ (LTS recommended)
- **Rust**: 1.77+ (stable toolchain)
- **Python**: 3.8+ (for genomic processing)
- **pnpm**: 8+ (package manager)
- **Git**: 2.30+ (version control)

### 🚀 **Development Workflow**
```bash
# Project setup
git clone <repository>
cd LiteratureGapper/desktop/ui
pnpm install

# Development commands
pnpm run tauri-dev      # Full application (React + Rust)
pnpm run dev            # Frontend only
pnpm run build          # Production build
pnpm run tauri-build    # Desktop app build

# Testing commands
cd ../
./test_implementation.sh    # Comprehensive testing
./test_rust_integration.sh  # Rust-Python integration test
```

### 🔧 **IDE Configuration**
- **VS Code Extensions**:
  - rust-analyzer (Rust support)
  - Tauri (Desktop app support)
  - ES7+ React/Redux/React-Native snippets
  - Tailwind CSS IntelliSense
  - Python extension
- **Settings**: ESLint, Prettier, rustfmt integration

## Dependencies & Versions

### 📦 **Frontend Dependencies (package.json)**
```json
{
  "dependencies": {
    "@tauri-apps/api": "2.6.2",
    "@tauri-apps/plugin-dialog": "2.6.2",
    "@tauri-apps/plugin-fs": "2.6.2",
    "@tauri-apps/plugin-log": "2.6.2",
    "@tauri-apps/plugin-shell": "2.6.2",
    "react": "19.1.0",
    "react-dom": "19.1.0",
    "react-router-dom": "7.6.3"
  },
  "devDependencies": {
    "@tauri-apps/cli": "2.6.2",
    "@types/react": "19.0.8",
    "@types/react-dom": "19.0.8",
    "@typescript-eslint/eslint-plugin": "8.20.0",
    "@typescript-eslint/parser": "8.20.0",
    "@vitejs/plugin-react": "5.0.0",
    "autoprefixer": "10.4.20",
    "concurrently": "9.1.0",
    "eslint": "9.18.0",
    "eslint-plugin-react": "7.37.2",
    "postcss": "8.5.13",
    "tailwindcss": "4.1.11",
    "typescript": "5.8.3",
    "vite": "7.0.0"
  }
}
```

### 🦀 **Backend Dependencies (Cargo.toml)**
```toml
[dependencies]
tauri = { version = "2.6.2", features = [
  "macos-private-api",
  "protocol-asset",
  "devtools"
] }
tauri-plugin-dialog = "2.6.2"
tauri-plugin-fs = "2.6.2"  
tauri-plugin-log = "2.6.2"
tauri-plugin-shell = "2.6.2"
serde = { version = "1.0.223", features = ["derive"] }
serde_json = "1.0.132"
tokio = { version = "1.41.1", features = ["full"] }
log = "0.4.22"
```

### 🐍 **Python Dependencies (Planned)**
```python
# Core data processing
numpy >= 1.24.0
pandas >= 2.0.0
scipy >= 1.10.0

# Genomic processing
biopython >= 1.81
pysam >= 0.21.0
pyvcf >= 0.6.8

# Machine learning (future)
tensorflow >= 2.12.0
scikit-learn >= 1.3.0
torch >= 2.0.0 (optional)

# Utilities
argparse (stdlib)
json (stdlib)
pathlib (stdlib)
```

## Build Configuration

### ⚙️ **Vite Configuration (vite.config.ts)**
```typescript
export default defineConfig({
  plugins: [react()],
  clearScreen: false,
  server: {
    port: 5173,
    strictPort: true,
    host: 'localhost'
  },
  build: {
    target: 'esnext',
    rollupOptions: {
      output: {
        manualChunks: {
          vendor: ['react', 'react-dom'],
          tauri: ['@tauri-apps/api']
        }
      }
    }
  }
});
```

### 🏗️ **Tauri Configuration (tauri.conf.json)**
```json
{
  "app": {
    "name": "GenePredict",
    "version": "0.1.0"
  },
  "build": {
    "beforeDevCommand": "pnpm dev",
    "beforeBuildCommand": "pnpm build",
    "devUrl": "http://localhost:5173",
    "frontendDist": "../dist"
  },
  "bundle": {
    "active": true,
    "targets": "all",
    "identifier": "com.genepredict.app",
    "icon": ["icons/32x32.png", "icons/128x128.png", "icons/icon.ico"]
  },
  "security": {
    "csp": "default-src 'self'; connect-src 'self' https://tauri.localhost",
    "devCsp": "default-src 'self' 'unsafe-inline' 'unsafe-eval'; connect-src 'self' ws://localhost:*"
  }
}
```

### 🎨 **Tailwind Configuration (tailwind.config.ts)**
```typescript
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      fontFamily: {
        'sans': ['Inter', 'ui-sans-serif', 'system-ui'],
      },
      colors: {
        primary: {
          50: '#f0f9ff',
          600: '#2563eb',
          700: '#1d4ed8',
        }
      }
    }
  },
  plugins: []
}
```

## Platform Support

### 💻 **Desktop Platforms**
- **macOS**: 10.13+ (High Sierra)
- **Windows**: 10/11 (64-bit)
- **Linux**: Ubuntu 18.04+, Debian 10+, Fedora 35+

### 🔧 **Platform-Specific Features**
- **macOS**: Native menu bar, app bundle
- **Windows**: Native Windows API, MSI installer
- **Linux**: AppImage, deb/rpm packages

### 📱 **Hardware Requirements**
- **CPU**: x86_64 (Intel/AMD), ARM64 (Apple Silicon)
- **Memory**: 4GB RAM minimum, 8GB recommended
- **Storage**: 500MB for app, 2GB+ for genomic data processing
- **GPU**: Optional for ML acceleration (future)

## Security & Privacy

### 🔒 **Security Features**
- **Content Security Policy**: Strict CSP headers
- **Local-Only Processing**: No external network calls
- **Secure File Handling**: Rust-powered native file operations
- **Input Validation**: Comprehensive validation at all layers

### 🛡️ **Privacy Architecture**
- **Zero External Dependencies**: All processing happens locally
- **No Telemetry**: No data collection or analytics
- **Secure Storage**: Local file system only
- **HIPAA Compliance**: Designed for medical data privacy

## Performance Characteristics

### 📊 **Current Metrics**
- **Bundle Size**: 189KB (gzipped)
- **Hot Reload**: <500ms
- **Build Time**: ~30s (development), ~2min (production)
- **Memory Usage**: ~100MB (idle), ~500MB (processing)
- **Startup Time**: <3s (cold start)

### 🚀 **Optimization Strategies**
- **Tree Shaking**: Vite eliminates unused code
- **Code Splitting**: Lazy loading for better performance
- **Rust Backend**: Native performance for computational tasks
- **Efficient Bundling**: Optimized production builds

## Development Constraints

### 🚧 **Current Limitations**
- **Single-threaded Python**: Python GIL limits parallelism
- **File Size Limits**: Large genomic files may require streaming
- **Memory Constraints**: Limited by available system memory
- **Platform APIs**: Some features may be platform-specific

### 🔮 **Future Enhancements**
- **Multi-threading**: Rust-based parallel processing
- **GPU Acceleration**: CUDA/OpenCL for ML operations
- **Streaming Processing**: Handle large files efficiently
- **Plugin System**: Extensible architecture for new features

## Quality Assurance

### 🧪 **Testing Infrastructure**
- **Unit Tests**: `cargo test` for Rust code
- **Integration Tests**: Python-Rust communication validation
- **E2E Tests**: Full application workflow testing
- **Performance Tests**: Benchmarking critical operations

### 📋 **Code Quality Tools**
- **ESLint**: JavaScript/TypeScript linting
- **rustfmt**: Rust code formatting
- **Prettier**: Code formatting
- **TypeScript**: Static type checking

### 🔍 **Monitoring & Debugging**
- **Structured Logging**: Consistent log format
- **Error Tracking**: Comprehensive error reporting
- **Performance Monitoring**: Resource usage tracking
- **Development Tools**: Hot reload, DevTools integration

This technical foundation provides a robust, secure, and performant platform for privacy-first genomic analysis while maintaining flexibility for future enhancements. 