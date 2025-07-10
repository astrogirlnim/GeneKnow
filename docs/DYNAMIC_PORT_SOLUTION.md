# 🔌 Dynamic Port Allocation for GenePredict

## Overview

GenePredict now supports dynamic port allocation for its local API server, ensuring the application works even when the default port (5000) is occupied. The server automatically finds the next available port and communicates it to the frontend.

## 🚀 How It Works

### 1. **Server Side (Python)**
The `enhanced_api_server.py` now includes:
- `find_available_port()` function that tries ports starting from 5000
- Writes the actual port to a file specified by `--port-file` argument
- Outputs `API_SERVER_PORT:XXXX` to stdout for process monitoring

```python
def find_available_port(start_port=5000, max_attempts=100):
    """Find an available port starting from start_port"""
    for port in range(start_port, start_port + max_attempts):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('', port))
                s.close()
                return port
        except OSError:
            continue
```

### 2. **Tauri Backend (Rust)**
The Rust backend manages the server lifecycle:
- Starts the Python server with a port file argument
- Monitors the port file creation
- Stores the port in a global state
- Provides `get_api_port` command for frontend

```rust
// Global state for API port
static API_SERVER_PORT: Lazy<Mutex<Option<u16>>> = Lazy::new(|| Mutex::new(None));

#[command]
async fn get_api_port() -> Result<u16, String> {
    let port_guard = API_SERVER_PORT.lock().unwrap();
    match *port_guard {
        Some(port) => Ok(port),
        None => Err("API server not started".to_string()),
    }
}
```

### 3. **Frontend (TypeScript)**
The frontend uses a centralized API configuration:
- `apiConfig.ts` manages dynamic port discovery
- Queries Tauri backend for the current port
- Falls back to environment variables in development

```typescript
class ApiConfig {
    async initialize(): Promise<void> {
        if (window.__TAURI__) {
            this.port = await invoke<number>('get_api_port');
            this.baseUrl = `http://localhost:${this.port}`;
        } else {
            // Development mode fallback
            this.baseUrl = import.meta.env.VITE_API_URL || 'http://localhost:5001';
        }
    }
}
```

## 🔧 Implementation Details

### Bundled Startup Script
The `start_api_server.sh` script in bundled resources:
```bash
# Create port file path
PORT_FILE="$SCRIPT_DIR/.api_server_port"

# Start with port file
exec "$PYTHON_EXE" -m geneknow_pipeline.enhanced_api_server --port-file "$PORT_FILE" "$@"
```

### Port Discovery Flow
1. User launches GenePredict app
2. Tauri backend starts Python server
3. Python finds available port (5000, 5001, 5002...)
4. Python writes port to `.api_server_port` file
5. Rust reads port from file and stores in memory
6. Frontend queries port via `get_api_port` command
7. All API calls use the dynamic port

## 🧪 Testing

### Manual Test
```bash
# Test dynamic port allocation
cd desktop
./scripts/test-dynamic-port.sh
```

### Test Script Features
- Starts server with port capture
- Waits for port file creation
- Tests health endpoint on dynamic port
- Cleans up after testing

## 🎯 Benefits

1. **No Port Conflicts**: App works even if port 5000 is busy
2. **Automatic Discovery**: Frontend automatically finds the correct port
3. **Seamless Experience**: Users don't need to configure ports
4. **Development Friendly**: Works in both dev and production modes

## 🚨 Error Handling

- If no port available in range 5000-5100: Error message
- If server fails to start: Cleanup and retry
- If port file not created: Timeout after 10 seconds
- Frontend fallback: Uses default port if Tauri command fails

## 📝 Configuration

No user configuration needed! The system handles everything automatically:
- Tries ports starting from 5000
- Finds next available port up to 5100
- Communicates port to all components
- Works on all platforms (macOS, Windows, Linux)

## 🔄 Migration Notes

For existing installations:
- No migration needed
- Old hardcoded port 5001 references updated
- Backward compatible with dev environments
- Environment variables still work for overrides 