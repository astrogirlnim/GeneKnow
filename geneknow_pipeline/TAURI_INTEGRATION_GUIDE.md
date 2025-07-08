# Tauri Integration Guide for GeneKnow Pipeline API

This guide explains how to integrate the GeneKnow LangGraph Pipeline API with your Tauri desktop application.

## Architecture Overview

```
┌─────────────────────┐     ┌─────────────────────┐     ┌─────────────────────┐
│  Tauri Frontend     │────▶│    Rust Backend     │────▶│   Python API        │
│  (React/TypeScript) │     │   (Tauri Commands)  │     │  (Flask/SocketIO)   │
└─────────────────────┘     └─────────────────────┘     └─────────────────────┘
         ▲                                                         │
         │                    WebSocket Updates                    │
         └─────────────────────────────────────────────────────────┘
```

## Setup Steps

### 1. Start the Python API Server

The API server must be running locally for the desktop app to work:

```bash
cd geneknow_pipeline
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
python enhanced_api_server.py
```

The server will run on `http://localhost:5001`

### 2. Add Rust Dependencies

In `desktop/src-tauri/Cargo.toml`, ensure you have:

```toml
[dependencies]
tauri = { version = "2", features = ["protocol-asset", "shell-open"] }
serde = { version = "1", features = ["derive"] }
serde_json = "1"
tokio = { version = "1", features = ["full"] }
reqwest = { version = "0.11", features = ["json"] }
```

### 3. Create Tauri Commands

Add these commands to `desktop/src-tauri/src/lib.rs`:

```rust
use serde::{Deserialize, Serialize};
use serde_json::json;

#[derive(Debug, Serialize, Deserialize)]
struct ProcessFileRequest {
    file_path: String,
    preferences: Option<serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
struct JobResponse {
    job_id: String,
    status: String,
    message: String,
}

#[tauri::command]
async fn process_genomic_file(
    file_path: String,
    preferences: Option<serde_json::Value>
) -> Result<String, String> {
    let client = reqwest::Client::new();
    
    let request_body = json!({
        "file_path": file_path,
        "preferences": preferences.unwrap_or(json!({}))
    });
    
    let response = client
        .post("http://localhost:5001/api/process")
        .json(&request_body)
        .send()
        .await
        .map_err(|e| format!("Request failed: {}", e))?;
    
    if response.status().is_success() {
        let job_data: JobResponse = response
            .json()
            .await
            .map_err(|e| format!("Failed to parse response: {}", e))?;
        
        Ok(job_data.job_id)
    } else {
        let error_data: serde_json::Value = response
            .json()
            .await
            .unwrap_or(json!({"error": "Unknown error"}));
        
        Err(error_data["error"].as_str().unwrap_or("Unknown error").to_string())
    }
}

#[tauri::command]
async fn get_job_status(job_id: String) -> Result<serde_json::Value, String> {
    let client = reqwest::Client::new();
    
    let response = client
        .get(format!("http://localhost:5001/api/status/{}", job_id))
        .send()
        .await
        .map_err(|e| format!("Request failed: {}", e))?;
    
    if response.status().is_success() {
        response
            .json()
            .await
            .map_err(|e| format!("Failed to parse response: {}", e))
    } else {
        Err("Failed to get job status".to_string())
    }
}

#[tauri::command]
async fn get_job_results(job_id: String) -> Result<serde_json::Value, String> {
    let client = reqwest::Client::new();
    
    let response = client
        .get(format!("http://localhost:5001/api/results/{}", job_id))
        .send()
        .await
        .map_err(|e| format!("Request failed: {}", e))?;
    
    if response.status().is_success() {
        response
            .json()
            .await
            .map_err(|e| format!("Failed to parse response: {}", e))
    } else {
        Err("Failed to get job results".to_string())
    }
}

#[tauri::command]
async fn check_api_health() -> Result<bool, String> {
    let client = reqwest::Client::new();
    
    let response = client
        .get("http://localhost:5001/api/health")
        .timeout(std::time::Duration::from_secs(5))
        .send()
        .await;
    
    match response {
        Ok(resp) => Ok(resp.status().is_success()),
        Err(_) => Ok(false),
    }
}
```

### 4. Register Commands

In the `run()` function, add the new commands:

```rust
.invoke_handler(tauri::generate_handler![
    // ... existing commands ...
    process_genomic_file,
    get_job_status,
    get_job_results,
    check_api_health,
])
```

### 5. Frontend Integration

Use the TypeScript API client in your React components:

```typescript
// Install socket.io-client first:
// pnpm add socket.io-client

import { invoke } from '@tauri-apps/api/core';
import { open } from '@tauri-apps/plugin-dialog';
import { useGeneKnowPipeline } from '../api/geneknowPipeline';

function FileProcessor() {
  const { client } = useGeneKnowPipeline();
  const [processing, setProcessing] = useState(false);
  const [progress, setProgress] = useState(0);
  const [results, setResults] = useState(null);

  const selectAndProcessFile = async () => {
    try {
      // Check if API is running
      const apiHealthy = await invoke<boolean>('check_api_health');
      if (!apiHealthy) {
        alert('GeneKnow API is not running. Please start it first.');
        return;
      }

      // Select file using Tauri dialog
      const selected = await open({
        multiple: false,
        filters: [{
          name: 'Genomic Files',
          extensions: ['fastq', 'fq', 'bam', 'vcf', 'maf']
        }]
      });

      if (!selected) return;

      setProcessing(true);
      setProgress(0);

      // Process file through Tauri command
      const jobId = await invoke<string>('process_genomic_file', {
        filePath: selected,
        preferences: {
          language: 'en',
          include_technical: true,
          patient_data: {
            age: 45,
            sex: 'F'
          }
        }
      });

      // Subscribe to progress updates
      client.subscribeToJobProgress(jobId, (update) => {
        setProgress(update.progress);
        console.log(`Progress: ${update.progress}% - ${update.current_step}`);
      });

      // Poll for completion (or use WebSocket)
      let completed = false;
      while (!completed) {
        const status = await invoke('get_job_status', { jobId });
        
        if (status.status === 'completed') {
          completed = true;
          const results = await invoke('get_job_results', { jobId });
          setResults(results);
          setProcessing(false);
        } else if (status.status === 'failed') {
          throw new Error(status.error || 'Processing failed');
        }

        await new Promise(resolve => setTimeout(resolve, 1000));
      }

    } catch (error) {
      console.error('Processing error:', error);
      setProcessing(false);
      alert(`Error: ${error.message}`);
    }
  };

  return (
    <div>
      <button onClick={selectAndProcessFile} disabled={processing}>
        {processing ? 'Processing...' : 'Select Genomic File'}
      </button>
      
      {processing && (
        <div>
          <progress value={progress} max={100} />
          <span>{progress}%</span>
        </div>
      )}
      
      {results && (
        <div>
          <h3>Results</h3>
          <p>Variants: {results.variant_count}</p>
          <h4>Risk Scores:</h4>
          <ul>
            {Object.entries(results.risk_scores).map(([type, score]) => (
              <li key={type}>{type}: {(score * 100).toFixed(1)}%</li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
}
```

### 6. Environment Configuration

Create `.env` files for different environments:

**`.env.development`**
```
VITE_API_URL=http://localhost:5001
VITE_WS_URL=http://localhost:5001
```

**`.env.production`**
```
VITE_API_URL=http://localhost:5001
VITE_WS_URL=http://localhost:5001
```

### 7. Auto-start API Server (Optional)

You can auto-start the Python API when the Tauri app launches:

```rust
use std::process::Command;

#[tauri::command]
async fn start_api_server() -> Result<bool, String> {
    // Check if already running
    let health_check = check_api_health().await?;
    if health_check {
        return Ok(true);
    }

    // Start the API server
    #[cfg(target_os = "windows")]
    let mut cmd = Command::new("cmd");
    #[cfg(target_os = "windows")]
    cmd.args(&["/C", "cd", "geneknow_pipeline", "&&", "python", "enhanced_api_server.py"]);

    #[cfg(not(target_os = "windows"))]
    let mut cmd = Command::new("sh");
    #[cfg(not(target_os = "windows"))]
    cmd.args(&["-c", "cd geneknow_pipeline && python enhanced_api_server.py"]);

    match cmd.spawn() {
        Ok(_) => {
            // Wait for server to start
            tokio::time::sleep(tokio::time::Duration::from_secs(3)).await;
            Ok(true)
        }
        Err(e) => Err(format!("Failed to start API server: {}", e))
    }
}
```

## Security Considerations

1. **Local Only**: The API runs on localhost only, no external access
2. **File Access**: Use Tauri's file dialog to ensure proper permissions
3. **Path Validation**: Always validate file paths on the Rust side
4. **API Key**: Consider adding a simple API key for local authentication

## Deployment

For production deployment:

1. **Bundle Python**: Include Python runtime and dependencies
2. **Service Management**: Use system services to manage the API
3. **Logging**: Implement proper logging for debugging
4. **Error Recovery**: Add retry logic and graceful degradation

## Testing

Test the integration:

```bash
# 1. Start the API server
cd geneknow_pipeline
python enhanced_api_server.py

# 2. Test with the test script
python test_enhanced_api.py

# 3. Run the Tauri app in dev mode
cd desktop
pnpm tauri dev
```

## Troubleshooting

### API Server Not Starting
- Check Python installation and virtual environment
- Verify all dependencies are installed
- Check port 5001 is not in use

### Connection Refused
- Ensure API server is running
- Check firewall settings
- Verify localhost connectivity

### CORS Issues
- API server includes proper CORS headers for Tauri
- Check the origin in browser dev tools

### File Path Issues
- Use absolute paths from Tauri file dialog
- Ensure file exists and is readable
- Check file size limits (5GB default)

## Next Steps

1. Implement progress persistence across app restarts
2. Add batch processing support
3. Create background service installer
4. Add user preferences storage
5. Implement result caching 