"""
Enhanced API server for GeneKnow LangGraph pipeline.
Designed for seamless integration with Tauri desktop app.
Supports file processing, progress tracking, and real-time updates.
"""
from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from flask_socketio import SocketIO, emit
import os
import tempfile
import shutil
import uuid
import threading
import time
from datetime import datetime
from typing import Dict, Any, Optional
import logging
import json
from pathlib import Path
from werkzeug.utils import secure_filename
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

# Import the pipeline
try:
    # Try relative import first (when running as module)
    from .graph import run_pipeline
except ImportError:
    # Fall back to absolute import (when running directly)
    from graph import run_pipeline

# Custom JSON encoder to handle numpy types and datetime
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if HAS_NUMPY:
            if isinstance(obj, (np.integer, np.int64)):
                return int(obj)
            elif isinstance(obj, (np.floating, np.float64)):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
        # Handle datetime objects
        if isinstance(obj, datetime):
            return obj.isoformat()
        # Handle pandas NaN/None
        if obj is None or (hasattr(obj, '__str__') and str(obj) == 'nan'):
            return None
        return super().default(obj)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize Flask app with SocketIO for real-time updates
app = Flask(__name__)
# Note: app.json_encoder is deprecated in newer Flask versions, we'll handle this in the response
CORS(app, origins=["tauri://localhost", "http://localhost:*"])  # Allow Tauri and local dev
socketio = SocketIO(app, cors_allowed_origins="*")

# Configuration
class Config:
    MAX_FILE_SIZE = 5 * 1024 * 1024 * 1024  # 5GB
    ALLOWED_EXTENSIONS = {'.fastq', '.fq', '.fastq.gz', '.fq.gz', '.bam', '.vcf', '.vcf.gz', '.maf', '.maf.gz'}
    UPLOAD_FOLDER = tempfile.mkdtemp(prefix='geneknow_uploads_')
    RESULTS_FOLDER = tempfile.mkdtemp(prefix='geneknow_results_')
    SESSION_TIMEOUT = 3600  # 1 hour
    
# In-memory storage for job tracking
jobs: Dict[str, Any] = {}
job_lock = threading.Lock()

# Utility functions
def allowed_file(filename: str) -> bool:
    """Check if file extension is allowed."""
    ext = os.path.splitext(filename.lower())[1]
    if filename.endswith('.gz'):
        ext = os.path.splitext(os.path.splitext(filename.lower())[0])[1] + '.gz'
    return ext in Config.ALLOWED_EXTENSIONS

def get_file_type(filename: str) -> str:
    """Determine file type from extension."""
    filename_lower = filename.lower()
    if filename_lower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        return 'fastq'
    elif filename_lower.endswith('.bam'):
        return 'bam'
    elif filename_lower.endswith(('.vcf', '.vcf.gz')):
        return 'vcf'
    elif filename_lower.endswith(('.maf', '.maf.gz')):
        return 'maf'
    return 'unknown'

def convert_numpy_types(obj):
    """Convert numpy types and other non-JSON-serializable types to Python native types."""
    if HAS_NUMPY:
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, np.str_):
            return str(obj)
    
    # Handle datetime objects
    if isinstance(obj, datetime):
        return obj.isoformat()
    
    # Handle pandas NaN/None and other special values
    if obj is None:
        return None
    if hasattr(obj, '__str__') and str(obj) == 'nan':
        return None
    
    # Handle dictionaries
    if isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    
    # Handle lists and tuples
    elif isinstance(obj, (list, tuple)):
        return [convert_numpy_types(item) for item in obj]
    
    # Handle sets
    elif isinstance(obj, set):
        return [convert_numpy_types(item) for item in obj]
    
    # Return as-is for other types
    return obj

def create_job(file_path: str, preferences: dict) -> str:
    """Create a new job entry."""
    job_id = str(uuid.uuid4())
    with job_lock:
        jobs[job_id] = {
            'id': job_id,
            'status': 'pending',
            'file_path': file_path,
            'file_type': get_file_type(file_path),
            'preferences': preferences,
            'created_at': datetime.now().isoformat(),
            'started_at': None,
            'completed_at': None,
            'progress': 0,
            'current_step': None,
            'result': None,
            'error': None
        }
    return job_id

def update_job(job_id: str, updates: dict):
    """Update job status and emit progress via WebSocket."""
    with job_lock:
        if job_id in jobs:
            # Convert any numpy types before updating
            updates = convert_numpy_types(updates)
            jobs[job_id].update(updates)
            # Emit progress update via WebSocket
            socketio.emit('job_progress', {
                'job_id': job_id,
                'status': jobs[job_id]['status'],
                'progress': jobs[job_id]['progress'],
                'current_step': jobs[job_id]['current_step']
            }, room=job_id)

def process_file_async(job_id: str):
    """Process file in background thread."""
    try:
        job = jobs[job_id]
        update_job(job_id, {
            'status': 'processing',
            'started_at': datetime.now().isoformat(),
            'progress': 5
        })
        
        # Create custom preferences with progress callback
        preferences = job['preferences'].copy()
        
        def progress_callback(node: str, progress: int):
            """Callback to update progress during pipeline execution."""
            update_job(job_id, {
                'current_step': node,
                'progress': min(progress, 95)  # Reserve last 5% for finalization
            })
        
        preferences['progress_callback'] = progress_callback
        
        # Run the pipeline
        result = run_pipeline(job['file_path'], preferences)
        
        # Process results
        if result.get('pipeline_status') == 'completed':
            # Convert all numpy types in the result before saving
            converted_result = convert_numpy_types(result)
            
            # Save results to file
            result_file = os.path.join(Config.RESULTS_FOLDER, f"{job_id}_result.json")
            with open(result_file, 'w') as f:
                json.dump(converted_result, f, indent=2)
            
            update_job(job_id, {
                'status': 'completed',
                'completed_at': datetime.now().isoformat(),
                'progress': 100,
                'result': convert_numpy_types({
                    'variant_count': result.get('variant_count', 0),
                    'risk_scores': result.get('risk_scores', {}),
                    'report_sections': result.get('report_sections', {}),
                    'processing_time': result.get('processing_time_seconds', 0),
                    'cadd_stats': result.get('cadd_stats', {}),
                    'tcga_matches': result.get('tcga_matches', {}),
                    'tcga_cohort_sizes': result.get('tcga_cohort_sizes', {}),
                    'prs_results': result.get('prs_results', {}),
                    'prs_summary': result.get('prs_summary', {}),
                    'structured_json': result.get('structured_json', {}),
                    'result_file': result_file
                })
            })
        else:
            raise Exception(result.get('errors', ['Unknown error'])[0])
            
    except Exception as e:
        logger.error(f"Job {job_id} failed: {str(e)}")
        update_job(job_id, {
            'status': 'failed',
            'completed_at': datetime.now().isoformat(),
            'error': str(e)
        })

# API Routes

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.now().isoformat(),
        'service': 'GeneKnow Pipeline API',
        'version': '2.0.0',
        'jobs_active': len([j for j in jobs.values() if j['status'] == 'processing'])
    })

@app.route('/api/pipeline-info', methods=['GET'])
def pipeline_info():
    """Get detailed pipeline capabilities."""
    return jsonify({
        'name': 'GeneKnow Genomic Risk Assessment',
        'version': '2.0.0',
        'capabilities': {
            'supported_formats': list(Config.ALLOWED_EXTENSIONS),
            'max_file_size_gb': Config.MAX_FILE_SIZE / (1024**3),
            'pipeline_nodes': [
                {'id': 'file_input', 'name': 'File Validation', 'description': 'Validates input file format'},
                {'id': 'preprocess', 'name': 'Preprocessing', 'description': 'Aligns FASTQ or loads variants'},
                {'id': 'variant_calling', 'name': 'Variant Calling', 'description': 'Identifies genetic variants'},
                {'id': 'qc_filter', 'name': 'Quality Control', 'description': 'Filters low-quality variants'},
                {'id': 'population_mapper', 'name': 'Population Mapping', 'description': 'Maps variants to population frequencies and ClinVar'},
                {'id': 'tcga_mapper', 'name': 'TCGA Frequency Analysis', 'description': 'Calculates tumor enrichment scores from TCGA database'},
                {'id': 'cadd_scoring', 'name': 'CADD Scoring', 'description': 'Enriches variants with deleteriousness scores'},
                {'id': 'feature_vector_builder', 'name': 'Feature Vector Builder', 'description': 'Combines static model outputs'},
                {'id': 'risk_model', 'name': 'Risk Assessment', 'description': 'ML-based risk prediction with multi-signal fusion'},
                {'id': 'formatter', 'name': 'Result Formatting', 'description': 'Formats results for frontend'},
                {'id': 'report_writer', 'name': 'Report Generation', 'description': 'Creates final report'}
            ],
            'cancer_types': ['breast', 'colon', 'lung', 'prostate', 'blood'],
            'output_formats': ['json', 'pdf', 'html'],
            'languages': ['en', 'es', 'hi']
        }
    })

@app.route('/api/supported-formats', methods=['GET'])
def supported_formats():
    """Get supported file formats with descriptions."""
    return jsonify({
        'formats': [
            {
                'extension': '.fastq',
                'description': 'Raw sequencing reads',
                'compressed': ['.fastq.gz'],
                'paired_end_support': True
            },
            {
                'extension': '.bam',
                'description': 'Aligned sequencing reads',
                'compressed': [],
                'paired_end_support': False
            },
            {
                'extension': '.vcf',
                'description': 'Variant Call Format',
                'compressed': ['.vcf.gz'],
                'paired_end_support': False
            },
            {
                'extension': '.maf',
                'description': 'Mutation Annotation Format',
                'compressed': ['.maf.gz'],
                'paired_end_support': False
            }
        ]
    })

@app.route('/api/upload', methods=['POST'])
def upload_file():
    """Upload a file for processing."""
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        
        if not allowed_file(file.filename):
            return jsonify({
                'error': f'Invalid file type. Allowed: {", ".join(Config.ALLOWED_EXTENSIONS)}'
            }), 400
        
        # Save uploaded file
        filename = secure_filename(file.filename)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        unique_filename = f"{timestamp}_{filename}"
        file_path = os.path.join(Config.UPLOAD_FOLDER, unique_filename)
        file.save(file_path)
        
        # Check file size
        file_size = os.path.getsize(file_path)
        if file_size > Config.MAX_FILE_SIZE:
            os.remove(file_path)
            return jsonify({
                'error': f'File too large. Maximum size: {Config.MAX_FILE_SIZE / (1024**3):.1f}GB'
            }), 413
        
        # Get preferences from form data or JSON
        preferences = {}
        if request.form.get('preferences'):
            preferences = json.loads(request.form.get('preferences'))
        
        # Create job
        job_id = create_job(file_path, preferences)
        
        # Start processing in background
        thread = threading.Thread(target=process_file_async, args=(job_id,))
        thread.start()
        
        return jsonify({
            'job_id': job_id,
            'filename': filename,
            'file_size': file_size,
            'file_type': get_file_type(filename),
            'status': 'uploaded',
            'message': 'File uploaded successfully. Processing started.'
        }), 202
        
    except Exception as e:
        logger.error(f"Upload error: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/process', methods=['POST'])
def process_file():
    """Process a file that's already on the system (for Tauri integration)."""
    try:
        data = request.get_json()
        if not data or 'file_path' not in data:
            return jsonify({'error': 'file_path is required'}), 400
        
        file_path = data['file_path']
        
        # Validate file exists
        if not os.path.exists(file_path):
            return jsonify({'error': f'File not found: {file_path}'}), 404
        
        # Validate file type
        if not allowed_file(file_path):
            return jsonify({
                'error': f'Invalid file type. Allowed: {", ".join(Config.ALLOWED_EXTENSIONS)}'
            }), 400
        
        # Get preferences
        preferences = data.get('preferences', {})
        
        # Create job
        job_id = create_job(file_path, preferences)
        
        # Start processing
        thread = threading.Thread(target=process_file_async, args=(job_id,))
        thread.start()
        
        return jsonify({
            'job_id': job_id,
            'status': 'processing',
            'message': 'Processing started'
        }), 202
        
    except Exception as e:
        logger.error(f"Process error: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/status/<job_id>', methods=['GET'])
def job_status(job_id: str):
    """Get job status."""
    with job_lock:
        if job_id not in jobs:
            return jsonify({'error': 'Job not found'}), 404
        
        job = jobs[job_id].copy()
        # Remove internal data
        job.pop('file_path', None)
        
        # Ensure all numpy types are converted before returning
        converted_job = convert_numpy_types(job)
        return jsonify(converted_job)

@app.route('/api/results/<job_id>', methods=['GET'])
def job_results(job_id: str):
    """Get job results."""
    with job_lock:
        if job_id not in jobs:
            return jsonify({'error': 'Job not found'}), 404
        
        job = jobs[job_id]
        if job['status'] != 'completed':
            return jsonify({
                'error': f'Job not completed. Current status: {job["status"]}'
            }), 400
        
        # Read full results from file
        result_file = job['result'].get('result_file')
        if result_file and os.path.exists(result_file):
            with open(result_file, 'r') as f:
                full_results = json.load(f)
            # Ensure all numpy types are converted before returning
            converted_results = convert_numpy_types(full_results)
            return jsonify(converted_results)
        else:
            # Fallback to job result, ensure it's converted
            converted_result = convert_numpy_types(job['result'])
            return jsonify(converted_result)

@app.route('/api/results/<job_id>/download', methods=['GET'])
def download_results(job_id: str):
    """Download results as JSON file."""
    with job_lock:
        if job_id not in jobs:
            return jsonify({'error': 'Job not found'}), 404
        
        job = jobs[job_id]
        if job['status'] != 'completed':
            return jsonify({
                'error': f'Job not completed. Current status: {job["status"]}'
            }), 400
        
        result_file = job['result'].get('result_file')
        if result_file and os.path.exists(result_file):
            return send_file(
                result_file,
                as_attachment=True,
                download_name=f'geneknow_results_{job_id}.json',
                mimetype='application/json'
            )
        else:
            return jsonify({'error': 'Results file not found'}), 404

@app.route('/api/cancel/<job_id>', methods=['POST'])
def cancel_job(job_id: str):
    """Cancel a running job."""
    with job_lock:
        if job_id not in jobs:
            return jsonify({'error': 'Job not found'}), 404
        
        job = jobs[job_id]
        if job['status'] not in ['pending', 'processing']:
            return jsonify({
                'error': f'Cannot cancel job with status: {job["status"]}'
            }), 400
        
        # Update job status
        job['status'] = 'cancelled'
        job['completed_at'] = datetime.now().isoformat()
        
        return jsonify({
            'message': 'Job cancelled',
            'job_id': job_id
        })

@app.route('/api/jobs', methods=['GET'])
def list_jobs():
    """List all jobs with optional filtering."""
    status_filter = request.args.get('status')
    limit = int(request.args.get('limit', 50))
    
    with job_lock:
        job_list = list(jobs.values())
        
        # Filter by status if provided
        if status_filter:
            job_list = [j for j in job_list if j['status'] == status_filter]
        
        # Sort by created_at descending
        job_list.sort(key=lambda x: x['created_at'], reverse=True)
        
        # Apply limit
        job_list = job_list[:limit]
        
        # Remove sensitive data
        for job in job_list:
            job.pop('file_path', None)
            if job['result']:
                job['result'].pop('result_file', None)
        
        return jsonify({
            'total': len(job_list),
            'jobs': job_list
        })

# WebSocket events for real-time updates

@socketio.on('connect')
def handle_connect():
    """Handle client connection."""
    logger.info(f"Client connected: {request.sid}")
    emit('connected', {'message': 'Connected to GeneKnow Pipeline'})

@socketio.on('disconnect')
def handle_disconnect():
    """Handle client disconnection."""
    logger.info(f"Client disconnected: {request.sid}")

@socketio.on('subscribe_job')
def handle_job_subscription(data):
    """Subscribe to job updates."""
    job_id = data.get('job_id')
    if job_id and job_id in jobs:
        # Join room for this job
        socketio.server.enter_room(request.sid, job_id)
        emit('subscribed', {'job_id': job_id})
        
        # Send current status
        with job_lock:
            job = jobs[job_id]
            emit('job_progress', {
                'job_id': job_id,
                'status': job['status'],
                'progress': job['progress'],
                'current_step': job['current_step']
            })
    else:
        emit('error', {'message': 'Invalid job ID'})

@socketio.on('unsubscribe_job')
def handle_job_unsubscription(data):
    """Unsubscribe from job updates."""
    job_id = data.get('job_id')
    if job_id:
        socketio.server.leave_room(request.sid, job_id)
        emit('unsubscribed', {'job_id': job_id})

# Cleanup temporary files on shutdown
import atexit

def cleanup_temp_files():
    """Clean up temporary directories."""
    try:
        shutil.rmtree(Config.UPLOAD_FOLDER, ignore_errors=True)
        shutil.rmtree(Config.RESULTS_FOLDER, ignore_errors=True)
        logger.info("Cleaned up temporary files")
    except Exception as e:
        logger.error(f"Error cleaning up temp files: {e}")

atexit.register(cleanup_temp_files)

def find_available_port(start_port=5000, max_attempts=100):
    """Find an available port starting from start_port"""
    import socket
    
    for port in range(start_port, start_port + max_attempts):
        try:
            # Try to bind to the port
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('', port))
                s.close()
                return port
        except OSError:
            # Port is in use, try the next one
            continue
    
    raise RuntimeError(f"No available ports found in range {start_port}-{start_port + max_attempts}")

# Main entry point
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='GeneKnow Pipeline API Server')
    parser.add_argument('--port', type=int, default=5000, help='Starting port to try (will find next available)')
    parser.add_argument('--port-file', help='File to write the actual port to')
    args = parser.parse_args()
    
    # Find an available port
    actual_port = find_available_port(args.port)
    
    # Write the port to a file if requested (for Tauri to read)
    if args.port_file:
        with open(args.port_file, 'w') as f:
            f.write(str(actual_port))
        logger.info(f"Wrote port {actual_port} to {args.port_file}")
    
    # Also write to stdout for the parent process to capture
    print(f"API_SERVER_PORT:{actual_port}", flush=True)
    
    debug = os.environ.get('DEBUG', 'False').lower() == 'true'
    
    logger.info(f"Starting Enhanced GeneKnow API Server on port {actual_port}")
    logger.info(f"Upload folder: {Config.UPLOAD_FOLDER}")
    logger.info(f"Results folder: {Config.RESULTS_FOLDER}")
    logger.info(f"WebSocket support enabled")
    
    socketio.run(app, host='0.0.0.0', port=actual_port, debug=debug) 