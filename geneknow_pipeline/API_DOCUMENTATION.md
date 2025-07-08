# GeneKnow Pipeline API Documentation

## Overview

The GeneKnow Pipeline API provides a RESTful interface with WebSocket support for processing genomic files through the LangGraph pipeline. It's designed for seamless integration with the Tauri desktop application and supports offline operation.

## Features

- **RESTful API** for file processing and job management
- **WebSocket support** for real-time progress updates
- **Async job processing** with background workers
- **File upload** and local file processing support
- **Comprehensive error handling** and validation
- **CORS support** for Tauri integration

## Getting Started

### Prerequisites

- Python 3.8+
- All genomic processing tools (BWA, SAMtools, etc.)
- Required Python packages (see requirements.txt)

### Installation

1. Navigate to the pipeline directory:
```bash
cd geneknow_pipeline
```

2. Create and activate a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Start the API server:
```bash
python enhanced_api_server.py
```

The server will start on `http://localhost:5001` by default.

### Configuration

Environment variables:
- `PORT`: Server port (default: 5001)
- `DEBUG`: Enable debug mode (default: False)

## API Endpoints

### Health Check

Check if the API server is running and healthy.

```http
GET /api/health
```

**Response:**
```json
{
  "status": "healthy",
  "timestamp": "2024-01-15T10:30:00Z",
  "service": "GeneKnow Pipeline API",
  "version": "2.0.0",
  "jobs_active": 2
}
```

### Pipeline Information

Get detailed information about pipeline capabilities.

```http
GET /api/pipeline-info
```

**Response:**
```json
{
  "name": "GeneKnow Genomic Risk Assessment",
  "version": "2.0.0",
  "capabilities": {
    "supported_formats": [".fastq", ".bam", ".vcf", ".maf"],
    "max_file_size_gb": 5.0,
    "pipeline_nodes": [
      {
        "id": "file_input",
        "name": "File Validation",
        "description": "Validates input file format"
      }
      // ... more nodes
    ],
    "cancer_types": ["breast", "colon", "lung", "prostate", "blood"],
    "output_formats": ["json", "pdf", "html"],
    "languages": ["en", "es", "hi"]
  }
}
```

### Supported Formats

Get detailed information about supported file formats.

```http
GET /api/supported-formats
```

**Response:**
```json
{
  "formats": [
    {
      "extension": ".fastq",
      "description": "Raw sequencing reads",
      "compressed": [".fastq.gz"],
      "paired_end_support": true
    }
    // ... more formats
  ]
}
```

### Upload and Process File

Upload a genomic file for processing.

```http
POST /api/upload
Content-Type: multipart/form-data
```

**Request Body:**
- `file`: The genomic file to process
- `preferences`: JSON string with processing preferences (optional)

**Example using curl:**
```bash
curl -X POST http://localhost:5001/api/upload \
  -F "file=@sample.vcf" \
  -F 'preferences={"language":"en","include_technical":true}'
```

**Response:**
```json
{
  "job_id": "123e4567-e89b-12d3-a456-426614174000",
  "filename": "sample.vcf",
  "file_size": 1024000,
  "file_type": "vcf",
  "status": "uploaded",
  "message": "File uploaded successfully. Processing started."
}
```

### Process Local File

Process a file that's already on the local system (ideal for Tauri integration).

```http
POST /api/process
Content-Type: application/json
```

**Request Body:**
```json
{
  "file_path": "/path/to/local/file.vcf",
  "preferences": {
    "language": "en",
    "include_technical": true,
    "patient_data": {
      "age": 45,
      "sex": "F",
      "family_history": true
    }
  }
}
```

**Response:**
```json
{
  "job_id": "123e4567-e89b-12d3-a456-426614174000",
  "status": "processing",
  "message": "Processing started"
}
```

### Get Job Status

Check the status of a processing job.

```http
GET /api/status/{job_id}
```

**Response:**
```json
{
  "id": "123e4567-e89b-12d3-a456-426614174000",
  "status": "processing",
  "file_type": "vcf",
  "created_at": "2024-01-15T10:30:00Z",
  "started_at": "2024-01-15T10:30:05Z",
  "progress": 45,
  "current_step": "variant_calling"
}
```

### Get Job Results

Retrieve the full results of a completed job.

```http
GET /api/results/{job_id}
```

**Response:**
```json
{
  "pipeline_status": "completed",
  "variant_count": 1523,
  "risk_scores": {
    "breast": 0.23,
    "colon": 0.15,
    "lung": 0.08
  },
  "report_sections": {
    "summary": "...",
    "detailed_findings": "...",
    "recommendations": "..."
  },
  "processing_time_seconds": 45.3
}
```

### Download Results

Download results as a JSON file.

```http
GET /api/results/{job_id}/download
```

**Response:** Binary file download

### Cancel Job

Cancel a running job.

```http
POST /api/cancel/{job_id}
```

**Response:**
```json
{
  "message": "Job cancelled",
  "job_id": "123e4567-e89b-12d3-a456-426614174000"
}
```

### List Jobs

List all jobs with optional filtering.

```http
GET /api/jobs?status=completed&limit=10
```

**Query Parameters:**
- `status` (optional): Filter by status (pending, processing, completed, failed, cancelled)
- `limit` (optional): Maximum number of jobs to return (default: 50)

**Response:**
```json
{
  "total": 25,
  "jobs": [
    {
      "id": "123e4567-e89b-12d3-a456-426614174000",
      "status": "completed",
      "created_at": "2024-01-15T10:30:00Z",
      "progress": 100
    }
    // ... more jobs
  ]
}
```

## WebSocket Events

### Connection

Connect to the WebSocket server at the same URL as the HTTP server.

```javascript
const socket = io('http://localhost:5001');
```

### Events

#### `connect`
Fired when connected to the server.

#### `disconnect`
Fired when disconnected from the server.

#### `subscribe_job`
Subscribe to job progress updates.

**Emit:**
```json
{
  "job_id": "123e4567-e89b-12d3-a456-426614174000"
}
```

#### `job_progress`
Receive job progress updates.

**Receive:**
```json
{
  "job_id": "123e4567-e89b-12d3-a456-426614174000",
  "status": "processing",
  "progress": 65,
  "current_step": "risk_model"
}
```

## Integration Examples

### JavaScript/TypeScript (Frontend)

```typescript
import { GeneKnowPipelineClient } from './api/geneknowPipeline';

const client = new GeneKnowPipelineClient();

// Process a file with progress tracking
async function processFile(file: File) {
  try {
    const { job_id } = await client.uploadFile(file, {
      language: 'en',
      include_technical: true
    });
    
    // Subscribe to progress updates
    client.subscribeToJobProgress(job_id, (progress) => {
      console.log(`Progress: ${progress.progress}% - ${progress.current_step}`);
    });
    
    // Wait for completion
    const job = await client.waitForJobCompletion(job_id);
    
    // Get results
    const results = await client.getJobResults(job_id);
    console.log('Risk scores:', results.risk_scores);
    
  } catch (error) {
    console.error('Processing failed:', error);
  }
}
```

### Python (Testing/Scripts)

```python
import requests
import json

# Process a local file
def process_genomic_file(file_path):
    # Start processing
    response = requests.post(
        'http://localhost:5001/api/process',
        json={
            'file_path': file_path,
            'preferences': {
                'language': 'en',
                'include_technical': True
            }
        }
    )
    
    job_data = response.json()
    job_id = job_data['job_id']
    
    # Poll for status
    while True:
        status_response = requests.get(
            f'http://localhost:5001/api/status/{job_id}'
        )
        status = status_response.json()
        
        print(f"Progress: {status['progress']}% - {status.get('current_step', 'N/A')}")
        
        if status['status'] in ['completed', 'failed', 'cancelled']:
            break
        
        time.sleep(1)
    
    # Get results if completed
    if status['status'] == 'completed':
        results_response = requests.get(
            f'http://localhost:5001/api/results/{job_id}'
        )
        return results_response.json()
    else:
        raise Exception(f"Job failed: {status.get('error', 'Unknown error')}")
```

### Tauri Integration

```rust
use serde_json::json;

#[tauri::command]
async fn process_genomic_file(file_path: String) -> Result<String, String> {
    // Call the Python API
    let client = reqwest::Client::new();
    let response = client
        .post("http://localhost:5001/api/process")
        .json(&json!({
            "file_path": file_path,
            "preferences": {
                "language": "en",
                "include_technical": true
            }
        }))
        .send()
        .await
        .map_err(|e| e.to_string())?;
    
    let job_data: serde_json::Value = response
        .json()
        .await
        .map_err(|e| e.to_string())?;
    
    Ok(job_data["job_id"].as_str().unwrap().to_string())
}
```

## Error Handling

All endpoints return appropriate HTTP status codes:

- `200 OK`: Successful request
- `202 Accepted`: Job accepted for processing
- `400 Bad Request`: Invalid request parameters
- `404 Not Found`: Resource not found
- `413 Payload Too Large`: File exceeds size limit
- `500 Internal Server Error`: Server error

Error responses include a JSON body:
```json
{
  "error": "Detailed error message"
}
```

## Security Considerations

1. **File Validation**: All uploaded files are validated for type and size
2. **Path Traversal Protection**: File paths are sanitized
3. **CORS**: Configured for Tauri and localhost only
4. **Resource Limits**: Maximum file size and concurrent jobs limited
5. **Temporary File Cleanup**: Automatic cleanup on server shutdown

## Performance Tips

1. Use WebSocket subscriptions for real-time updates instead of polling
2. Process files locally (via file path) when possible to avoid upload overhead
3. Implement client-side file validation before upload
4. Use the job listing endpoint with filters to manage job history

## Troubleshooting

### Common Issues

1. **Connection Refused**: Ensure the API server is running on the correct port
2. **CORS Errors**: Check that the origin is allowed in CORS configuration
3. **File Not Found**: Verify file paths are absolute when using local processing
4. **WebSocket Connection Failed**: Check firewall settings and WebSocket support

### Debug Mode

Enable debug mode for detailed logging:
```bash
DEBUG=True python enhanced_api_server.py
```

## Future Enhancements

- [ ] Authentication and user sessions
- [ ] Batch processing support
- [ ] Result caching
- [ ] Export to additional formats (PDF, HTML)
- [ ] Webhook notifications for job completion
- [ ] Rate limiting and quotas
- [ ] S3/cloud storage integration 