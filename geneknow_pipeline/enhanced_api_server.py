"""
Enhanced API server for GeneKnow LangGraph pipeline.
Designed for seamless integration with Tauri desktop app.
Supports file processing, progress tracking, and real-time updates.
"""

# Import standard library modules first
import tempfile
import uuid
import shutil
import argparse
import json

# Import third-party modules
from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from flask_socketio import SocketIO, emit, join_room, leave_room
from werkzeug.utils import secure_filename
from datetime import datetime
import threading
import os
import logging
from typing import Dict, Any

# Optional numpy import
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

# Import local modules
try:
    # Try relative import first (when run as module)
    from .graph import run_pipeline
except ImportError:
    # Fall back to absolute import (when run directly)
    import sys

    # Add the parent directory to the path
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
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
        if obj is None or (hasattr(obj, "__str__") and str(obj) == "nan"):
            return None
        return super().default(obj)

# Configure logging with proper formatting and unbuffered output
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler()
    ],  # Explicitly use StreamHandler for console output
)
logger = logging.getLogger(__name__)

# Force unbuffered output for Python (important when running as subprocess)
import sys

sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(line_buffering=True)

# Print startup message immediately (before any imports that might affect output)
print("Starting GeneKnow Enhanced API Server...", flush=True)

# Initialize Flask app with SocketIO for real-time updates
app = Flask(__name__)
# Note: app.json_encoder is deprecated in newer Flask versions, we'll handle this in the response

# Configure Flask for large file uploads
app.config['MAX_CONTENT_LENGTH'] = 5 * 1024 * 1024 * 1024  # 5GB
app.config['UPLOAD_FOLDER'] = tempfile.mkdtemp(prefix="geneknow_uploads_")
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

CORS(
    app, origins=["tauri://localhost", "http://localhost:*"]
)  # Allow Tauri and local dev

# Configure eventlet properly to avoid blocking issues
try:
    import eventlet
    # Patch standard library for async operation, but preserve stdout/stderr
    eventlet.monkey_patch(socket=True, select=True, thread=False)
    print("Eventlet monkey patching completed", flush=True)
except ImportError:
    print("WARNING: eventlet not installed, falling back to threading mode", flush=True)
    eventlet = None

# Initialize SocketIO with appropriate async mode
async_mode = 'eventlet' if eventlet else 'threading'
socketio = SocketIO(app, cors_allowed_origins="*", async_mode=async_mode, logger=True, engineio_logger=True)
print(f"SocketIO initialized with async_mode: {async_mode}", flush=True)

# Configuration
class Config:
    MAX_FILE_SIZE = 5 * 1024 * 1024 * 1024  # 5GB
    ALLOWED_EXTENSIONS = {
        ".fastq",
        ".fq",
        ".fastq.gz",
        ".fq.gz",
        ".bam",
        ".vcf",
        ".vcf.gz",
        ".maf",
        ".maf.gz",
    }
    UPLOAD_FOLDER = tempfile.mkdtemp(prefix="geneknow_uploads_")
    RESULTS_FOLDER = tempfile.mkdtemp(prefix="geneknow_results_")
    SESSION_TIMEOUT = 3600  # 1 hour


# In-memory storage for job tracking
jobs: Dict[str, Any] = {}
job_lock = threading.Lock()


# Utility functions
def allowed_file(filename: str) -> bool:
    """Check if file extension is allowed."""
    ext = os.path.splitext(filename.lower())[1]
    if filename.endswith(".gz"):
        ext = os.path.splitext(os.path.splitext(filename.lower())[0])[1] + ".gz"
    return ext in Config.ALLOWED_EXTENSIONS


def get_file_type(filename: str) -> str:
    """Determine file type from extension."""
    filename_lower = filename.lower()
    if filename_lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "fastq"
    elif filename_lower.endswith(".bam"):
        return "bam"
    elif filename_lower.endswith((".vcf", ".vcf.gz")):
        return "vcf"
    elif filename_lower.endswith((".maf", ".maf.gz")):
        return "maf"
    return "unknown"


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
    if hasattr(obj, "__str__") and str(obj) == "nan":
        return None

    # Handle dictionaries
    if isinstance(obj, dict):
        converted = {}
        for k, v in obj.items():
            # Skip non-serializable keys like ML model instances
            if k in ["ml_fusion_model_instance", "ml_fusion_feature_matrix"]:
                continue  # Skip these ML objects that can't be serialized
            converted[k] = convert_numpy_types(v)
        return converted

    # Handle lists and tuples
    elif isinstance(obj, (list, tuple)):
        return [convert_numpy_types(item) for item in obj]

    # Handle sets
    elif isinstance(obj, set):
        return [convert_numpy_types(item) for item in obj]

    # Return as-is for other types
    return obj


def format_pipeline_results_for_frontend(pipeline_state: dict) -> dict:
    """
    Transform full pipeline state into frontend-compatible PipelineResult format.
    
    This function extracts the key fields that the frontend expects and formats
    them according to the TypeScript PipelineResult interface.
    """
    # Extract core fields that frontend expects
    formatted_result = {
        "pipeline_status": pipeline_state.get("pipeline_status", "completed"),
        "variant_count": pipeline_state.get("variant_count", 0),
        "risk_scores": pipeline_state.get("risk_scores", {}),
        "risk_genes": pipeline_state.get("risk_genes", {}),
        "processing_time_seconds": pipeline_state.get("processing_time_seconds", 0),
        
        # Report sections - transform if needed
        "report_sections": pipeline_state.get("report_sections", {}),
        
        # Enhanced report content
        "enhanced_report_content": pipeline_state.get("enhanced_report_content", {}),
        
        # Report generator info
        "report_generator_info": pipeline_state.get("report_generator_info", {}),
        
        # TCGA and analysis data
        "tcga_matches": pipeline_state.get("tcga_matches", {}),
        "cadd_stats": pipeline_state.get("cadd_stats", {}),
        
        # Pathway burden results
        "pathway_burden_results": pipeline_state.get("pathway_burden_results", {}),
        "pathway_burden_summary": pipeline_state.get("pathway_burden_summary", {}),
        
        # Structured JSON for detailed frontend components
        "structured_json": pipeline_state.get("structured_json", {}),
        
        # Variant details for tables
        "variants": format_variants_for_frontend(pipeline_state.get("variant_details", []))
    }
    
    return formatted_result


def format_variants_for_frontend(variant_details: list) -> list:
    """Format variant details for frontend consumption."""
    if not variant_details:
        return []
    
    formatted_variants = []
    for variant in variant_details[:10]:  # Limit to top 10 variants
        formatted_variant = {
            "gene": variant.get("gene", "Unknown"),
            "position": variant.get("position", 0),
            "type": variant.get("mutation_type", variant.get("consequence", "Unknown")),
            "impact": variant.get("functional_impact", "Unknown"),
            "quality_score": variant.get("quality_metrics", {}).get("quality", 0),
            "clinical_significance": variant.get("clinical_significance", "Unknown")
        }
        formatted_variants.append(formatted_variant)
    
    return formatted_variants


def create_job(file_path: str, preferences: dict) -> str:
    """Create a new job entry."""
    job_id = str(uuid.uuid4())
    with job_lock:
        jobs[job_id] = {
            "id": job_id,
            "status": "pending",
            "file_path": file_path,
            "file_type": get_file_type(file_path),
            "preferences": preferences,
            "created_at": datetime.now().isoformat(),
            "started_at": None,
            "completed_at": None,
            "progress": 0,
            "current_step": None,
            "result": None,
            "error": None,
        }
    return job_id


def update_job(job_id: str, updates: Dict[str, Any]):
    """Update job information."""
    with job_lock:
        if job_id in jobs:
            # Convert any numpy types before updating
            updates = convert_numpy_types(updates)
            jobs[job_id].update(updates)
            # Track last activity
            jobs[job_id]["last_activity"] = datetime.now().isoformat()
            # Emit update via WebSocket
            socketio.emit(
                "job_progress",
                {
                    "job_id": job_id,
                    "status": jobs[job_id]["status"],
                    "progress": jobs[job_id]["progress"],
                    "current_step": jobs[job_id]["current_step"],
                },
                room=job_id,
            )


def cleanup_job_files(job_id: str, cleanup_results: bool = False):
    """Clean up temporary files associated with a job."""
    with job_lock:
        job = jobs.get(job_id)
        if not job:
            return

        # Clean up the uploaded file if it exists in temp directory
        file_path = job.get("file_path")
        if file_path and os.path.exists(file_path):
            # Only delete if it's in a temp directory
            temp_markers = [
                "/tmp/",
                "/var/folders/",
                Config.UPLOAD_FOLDER,
                "geneknow_temp",
            ]
            if any(marker in file_path for marker in temp_markers):
                try:
                    os.remove(file_path)
                    logger.info(f"Cleaned up temporary file: {file_path}")
                except Exception as e:
                    logger.error(f"Failed to clean up file {file_path}: {e}")

        # Only clean up result files if explicitly requested (e.g., on server shutdown)
        if cleanup_results:
            result = job.get("result", {})
            if isinstance(result, dict):
                result_file = result.get("result_file")
                if result_file and os.path.exists(result_file):
                    try:
                        os.remove(result_file)
                        logger.info(f"Cleaned up result file: {result_file}")
                    except Exception as e:
                        logger.error(
                            f"Failed to clean up result file {result_file}: {e}"
                        )


def process_file_async(job_id: str):
    """Process file in background thread."""
    try:
        job = jobs[job_id]
        update_job(
            job_id,
            {
                "status": "processing",
                "started_at": datetime.now().isoformat(),
                "progress": 5,
            },
        )

        # Create custom preferences with progress callback
        preferences = job["preferences"].copy()

        def progress_callback(node: str, progress: int):
            """Callback to update progress during pipeline execution."""
            update_job(
                job_id,
                {
                    "current_step": node,
                    "progress": min(progress, 95),
                },  # Reserve last 5% for finalization
            )

        preferences["progress_callback"] = progress_callback

        # Run the pipeline
        result = run_pipeline(job["file_path"], preferences)

        # Process results
        if result.get("pipeline_status") == "completed":
            # Convert all numpy types in the result before saving
            converted_result = convert_numpy_types(result)

            # Save results to file
            result_file = os.path.join(Config.RESULTS_FOLDER, f"{job_id}_result.json")
            with open(result_file, "w") as f:
                json.dump(converted_result, f, indent=2)

            update_job(
                job_id,
                {
                    "status": "completed",
                    "completed_at": datetime.now().isoformat(),
                    "progress": 100,
                    "result": convert_numpy_types(
                        {
                            "variant_count": result.get("variant_count", 0),
                            "risk_scores": result.get("risk_scores", {}),
                            "report_sections": result.get("report_sections", {}),
                            "processing_time_seconds": result.get(
                                "processing_time_seconds", 0
                            ),
                            "cadd_stats": result.get("cadd_stats", {}),
                            "tcga_matches": result.get("tcga_matches", {}),
                            "tcga_cohort_sizes": result.get("tcga_cohort_sizes", {}),
                            "prs_results": result.get("prs_results", {}),
                            "prs_summary": result.get("prs_summary", {}),
                            "structured_json": result.get("structured_json", {}),
                            "ml_fusion_results": result.get("ml_fusion_results", {}),
                            "ml_risk_assessment": result.get("ml_risk_assessment", {}),
                            "metrics": result.get("metrics", {}),
                            "metrics_summary": result.get("metrics_summary", {}),
                            "completed_nodes": result.get("completed_nodes", []),
                            "warnings": result.get("warnings", []),
                            "enhanced_report_content": result.get(
                                "enhanced_report_content", {}
                            ),
                            "report_generator_info": result.get(
                                "report_generator_info", {}
                            ),
                            "result_file": result_file,
                        }
                    ),
                },
            )
            cleanup_job_files(
                job_id, cleanup_results=False
            )  # Clean up temp files but preserve results
        else:
            raise Exception(result.get("errors", ["Unknown error"])[0])

    except Exception as e:
        logger.error(f"Job {job_id} failed: {str(e)}")
        update_job(
            job_id,
            {
                "status": "failed",
                "completed_at": datetime.now().isoformat(),
                "error": str(e),
            },
        )
        cleanup_job_files(job_id, cleanup_results=True)  # Clean up all files on failure


# API Routes


@app.route("/api/health", methods=["GET"])
def health_check():
    """Health check endpoint."""
    return jsonify(
        {
            "status": "healthy",
            "timestamp": datetime.now().isoformat(),
            "service": "GeneKnow Pipeline API",
            "version": "2.0.0",
            "jobs_active": len(
                [j for j in jobs.values() if j["status"] == "processing"]
            ),
        }
    )


@app.route("/api/pipeline-info", methods=["GET"])
def pipeline_info():
    """Get detailed pipeline capabilities."""
    return jsonify(
        {
            "name": "GeneKnow Genomic Risk Assessment",
            "version": "2.0.0",
            "capabilities": {
                "supported_formats": list(Config.ALLOWED_EXTENSIONS),
                "max_file_size_gb": Config.MAX_FILE_SIZE / (1024**3),
                "pipeline_nodes": [
                    {
                        "id": "file_input",
                        "name": "File Validation",
                        "description": "Validates input file format",
                    },
                    {
                        "id": "preprocess",
                        "name": "Preprocessing",
                        "description": "Aligns FASTQ or loads variants",
                    },
                    {
                        "id": "variant_calling",
                        "name": "Variant Calling",
                        "description": "Identifies genetic variants",
                    },
                    {
                        "id": "qc_filter",
                        "name": "Quality Control",
                        "description": "Filters low-quality variants",
                    },
                    {
                        "id": "population_mapper",
                        "name": "Population Mapping",
                        "description": "Maps variants to population frequencies and ClinVar",
                    },
                    {
                        "id": "tcga_mapper",
                        "name": "TCGA Frequency Analysis",
                        "description": "Calculates tumor enrichment scores from TCGA database",
                    },
                    {
                        "id": "cadd_scoring",
                        "name": "CADD Scoring",
                        "description": "Enriches variants with deleteriousness scores",
                    },
                    {
                        "id": "feature_vector_builder",
                        "name": "Feature Vector Builder",
                        "description": "Combines static model outputs",
                    },
                    {
                        "id": "risk_model",
                        "name": "Risk Assessment",
                        "description": "ML-based risk prediction with multi-signal fusion",
                    },
                    {
                        "id": "formatter",
                        "name": "Result Formatting",
                        "description": "Formats results for frontend",
                    },
                    {
                        "id": "report_writer",
                        "name": "Report Generation",
                        "description": "Creates final report",
                    },
                ],
                "cancer_types": ["breast", "colon", "lung", "prostate", "blood"],
                "output_formats": ["json", "pd", "html"],
                "languages": ["en", "es", "hi"],
            },
        }
    )


@app.route("/api/supported-formats", methods=["GET"])
def supported_formats():
    """Get supported file formats with descriptions."""
    return jsonify(
        {
            "formats": [
                {
                    "extension": ".fastq",
                    "description": "Raw sequencing reads",
                    "compressed": [".fastq.gz"],
                    "paired_end_support": True,
                },
                {
                    "extension": ".bam",
                    "description": "Aligned sequencing reads",
                    "compressed": [],
                    "paired_end_support": False,
                },
                {
                    "extension": ".vcf",
                    "description": "Variant Call Format",
                    "compressed": [".vcf.gz"],
                    "paired_end_support": False,
                },
                {
                    "extension": ".maf",
                    "description": "Mutation Annotation Format",
                    "compressed": [".maf.gz"],
                    "paired_end_support": False,
                },
            ]
        }
    )


@app.route("/api/upload", methods=["POST"])
def upload_file():
    """Upload a file for processing with streaming support for large files."""
    try:
        if "file" not in request.files:
            return jsonify({"error": "No file provided"}), 400

        file = request.files["file"]
        if not file.filename or file.filename == "":
            return jsonify({"error": "No file selected"}), 400

        if not allowed_file(file.filename):
            return (
                jsonify(
                    {
                        "error": f'Invalid file type. Allowed: {", ".join(Config.ALLOWED_EXTENSIONS)}'
                    }
                ),
                400,
            )

        # Save uploaded file with streaming to handle large files
        filename = secure_filename(file.filename)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        unique_filename = f"{timestamp}_{filename}"
        file_path = os.path.join(Config.UPLOAD_FOLDER, unique_filename)
        
        # Stream the file to disk instead of loading into memory
        logger.info(f"Starting streaming upload for {filename}")
        bytes_written = 0
        chunk_size = 64 * 1024  # 64KB chunks for better performance
        
        try:
            with open(file_path, 'wb') as f:
                while True:
                    chunk = file.stream.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    bytes_written += len(chunk)
                    
                    # Log progress for very large files
                    if bytes_written % (100 * 1024 * 1024) == 0:  # Every 100MB
                        logger.info(f"Uploaded {bytes_written / (1024*1024):.1f}MB of {filename}")
            
            logger.info(f"Upload completed: {filename} ({bytes_written / (1024*1024):.1f}MB)")
            
        except Exception as e:
            # Clean up partial file on error
            if os.path.exists(file_path):
                os.remove(file_path)
            raise e

        # Check file size
        file_size = os.path.getsize(file_path)
        if file_size > Config.MAX_FILE_SIZE:
            os.remove(file_path)
            return (
                jsonify(
                    {
                        "error": f"File too large. Maximum size: {Config.MAX_FILE_SIZE / (1024**3):.1f}GB"
                    }
                ),
                413,
            )

        # Get preferences from form data or JSON
        preferences = {}
        preferences_str = request.form.get("preferences")
        if preferences_str:
            preferences = json.loads(preferences_str)

        # Create job
        job_id = create_job(file_path, preferences)

        # Start processing in background
        thread = threading.Thread(target=process_file_async, args=(job_id,))
        thread.start()

        return (
            jsonify(
                {
                    "job_id": job_id,
                    "filename": filename,
                    "file_size": file_size,
                    "file_type": get_file_type(filename),
                    "status": "uploaded",
                    "message": "File uploaded successfully. Processing started.",
                }
            ),
            202,
        )

    except Exception as e:
        logger.error(f"Upload error: {str(e)}")
        return jsonify({"error": str(e)}), 500


@app.route("/api/process", methods=["POST"])
def process_file():
    """Process a file that's already on the system (for Tauri integration)."""
    try:
        data = request.get_json()
        if not data or "file_path" not in data:
            return jsonify({"error": "file_path is required"}), 400

        file_path = data["file_path"]

        # Validate file exists
        if not os.path.exists(file_path):
            return jsonify({"error": f"File not found: {file_path}"}), 404

        # Validate file type
        if not allowed_file(file_path):
            return (
                jsonify(
                    {
                        "error": f'Invalid file type. Allowed: {", ".join(Config.ALLOWED_EXTENSIONS)}'
                    }
                ),
                400,
            )

        # Get preferences
        preferences = data.get("preferences", {})

        # Create job
        job_id = create_job(file_path, preferences)

        # Start processing
        thread = threading.Thread(target=process_file_async, args=(job_id,))
        thread.start()

        return (
            jsonify(
                {
                    "job_id": job_id,
                    "status": "processing",
                    "message": "Processing started",
                }
            ),
            202,
        )

    except Exception as e:
        logger.error(f"Process error: {str(e)}")
        return jsonify({"error": str(e)}), 500


@app.route("/api/status/<job_id>", methods=["GET"])
def job_status(job_id: str):
    """Get job status."""
    with job_lock:
        if job_id not in jobs:
            return jsonify({"error": "Job not found"}), 404

        job = jobs[job_id].copy()
        # Remove internal data
        job.pop("file_path", None)

        # Ensure all numpy types are converted before returning
        converted_job = convert_numpy_types(job)
        return jsonify(converted_job)


@app.route("/api/results/<job_id>", methods=["GET"])
def job_results(job_id: str):
    """Get job results."""
    with job_lock:
        if job_id not in jobs:
            return jsonify({"error": "Job not found"}), 404

        job = jobs[job_id]
        if job["status"] != "completed":
            return (
                jsonify(
                    {"error": f'Job not completed. Current status: {job["status"]}'}
                ),
                400,
            )

        # Read full results from file
        result_file = job["result"].get("result_file")
        logger.info(f"Looking for result file: {result_file}")
        logger.info(f"Result file exists: {result_file and os.path.exists(result_file)}")
        
        if result_file and os.path.exists(result_file):
            logger.info("Reading full results from file and formatting for frontend")
            with open(result_file, "r") as f:
                full_results = json.load(f)
            
            # Transform full pipeline state into frontend-compatible format
            formatted_results = format_pipeline_results_for_frontend(full_results)
            logger.info(f"Formatted results keys: {list(formatted_results.keys())}")
            
            # Ensure all numpy types are converted before returning
            converted_results = convert_numpy_types(formatted_results)
            return jsonify(converted_results)
        else:
            logger.info("Using fallback job result")
            # Fallback to job result, ensure it's converted
            converted_result = convert_numpy_types(job["result"])
            return jsonify(converted_result)


@app.route("/api/results/<job_id>/download", methods=["GET"])
def download_results(job_id: str):
    """Download results as JSON file."""
    with job_lock:
        if job_id not in jobs:
            return jsonify({"error": "Job not found"}), 404

        job = jobs[job_id]
        if job["status"] != "completed":
            return (
                jsonify(
                    {"error": f'Job not completed. Current status: {job["status"]}'}
                ),
                400,
            )

        result_file = job["result"].get("result_file")
        if result_file and os.path.exists(result_file):
            return send_file(
                result_file,
                as_attachment=True,
                download_name=f"geneknow_results_{job_id}.json",
                mimetype="application/json",
            )
        else:
            return jsonify({"error": "Results file not found"}), 404


@app.route("/api/cancel/<job_id>", methods=["POST"])
def cancel_job(job_id: str):
    """Cancel a running job."""
    with job_lock:
        if job_id not in jobs:
            return jsonify({"error": "Job not found"}), 404

        job = jobs[job_id]
        if job["status"] not in ["pending", "processing"]:
            return (
                jsonify({"error": f'Cannot cancel job with status: {job["status"]}'}),
                400,
            )

        # Update job status
        job["status"] = "cancelled"
        job["completed_at"] = datetime.now().isoformat()

        return jsonify({"message": "Job cancelled", "job_id": job_id})


@app.route("/api/jobs", methods=["GET"])
def list_jobs():
    """List all jobs with optional filtering."""
    status_filter = request.args.get("status")
    limit = int(request.args.get("limit", 50))

    with job_lock:
        job_list = list(jobs.values())

        # Filter by status if provided
        if status_filter:
            job_list = [j for j in job_list if j["status"] == status_filter]

        # Sort by created_at descending
        job_list.sort(key=lambda x: x["created_at"], reverse=True)

        # Apply limit
        job_list = job_list[:limit]

        # Remove sensitive data
        for job in job_list:
            job.pop("file_path", None)
            if job["result"]:
                job["result"].pop("result_file", None)

        return jsonify({"total": len(job_list), "jobs": job_list})


# Report Generator Configuration Endpoints

@app.route('/api/report-generator/config', methods=['GET'])
def get_report_config():
    """Get current report generator configuration."""
    try:
        from nodes.report_generator.config import load_config
        config = load_config()
        
        return jsonify({
            'backend': config.backend.value,
            'model_name': config.model_name,
            'temperature': config.temperature,
            'max_tokens': config.max_tokens,
            'style': config.style.value,
            'enable_streaming': config.enable_streaming,
            'enable_parallel_generation': config.enable_parallel_generation,
            'max_parallel_workers': config.max_parallel_workers,
            'risk_threshold': config.risk_threshold,
            'include_glossary': config.include_glossary,
            'include_technical_appendix': config.include_technical_appendix,
            'output_formats': config.output_formats
        })
    except Exception as e:
        logger.error(f"Error loading report config: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/report-generator/config', methods=['POST'])
def save_report_config():
    """Save report generator configuration."""
    try:
        from nodes.report_generator.config import save_config, ReportConfig, LLMBackend, ReportStyle
        
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No data provided'}), 400
        
        # Validate and convert backend
        backend_str = data.get('backend', 'none')
        try:
            backend = LLMBackend(backend_str)
        except ValueError:
            return jsonify({'error': f'Invalid backend: {backend_str}'}), 400
        
        # Validate and convert style
        style_str = data.get('style', 'clinician')
        try:
            style = ReportStyle(style_str)
        except ValueError:
            return jsonify({'error': f'Invalid style: {style_str}'}), 400
        
        # Create config object
        config = ReportConfig(
            backend=backend,
            model_name=data.get('model_name'),
            temperature=float(data.get('temperature', 0.3)),
            max_tokens=int(data.get('max_tokens', 2000)),
            style=style,
            enable_streaming=bool(data.get('enable_streaming', True)),
            enable_parallel_generation=bool(data.get('enable_parallel_generation', True)),
            max_parallel_workers=int(data.get('max_parallel_workers', 5)),
            risk_threshold=float(data.get('risk_threshold', 5.0)),
            include_glossary=bool(data.get('include_glossary', True)),
            include_technical_appendix=bool(data.get('include_technical_appendix', True)),
            output_formats=data.get('output_formats', ['markdown'])
        )
        
        # Save configuration
        save_config(config)
        
        return jsonify({'message': 'Configuration saved successfully'})
    
    except Exception as e:
        logger.error(f"Error saving report config: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/report-generator/available-models', methods=['GET'])
def get_available_models():
    """Get available LLM models for report generation."""
    try:
        from nodes.report_generator.model_interface import OllamaBackend
        
        # Check Ollama
        ollama = OllamaBackend()
        ollama_available = ollama.is_available()
        ollama_models = ollama.available_models if ollama_available else []
        
        return jsonify({
            'status': {
                'ollama': ollama_available
            },
            'models': {
                'ollama': ollama_models
            },
            'recommended': {
                'ollama': ['llama3', 'mistral', 'codellama']
            }
        })
    except Exception as e:
        logger.error(f"Error getting available models: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/report-generator/warm-model', methods=['POST'])
def warm_model():
    """Warm up a model by loading it into memory (not needed for Ollama)."""
    try:
        data = request.get_json()
        model_name = data.get('model_name', 'auto')
        backend = data.get('backend', 'ollama')
        
        if backend != 'ollama':
            return jsonify({'message': 'Model warming only supported for Ollama models'}), 200
        
        # Ollama models are loaded on-demand, no warming needed
        return jsonify({
            'success': True,
            'message': 'Ollama models are loaded on-demand',
            'model': model_name,
            'backend': backend
        })
        
    except Exception as e:
        logger.error(f"Error warming model: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/report-generator/model-status', methods=['GET'])
def get_model_status():
    """Get the current status of loaded models."""
    try:
        from nodes.report_generator.model_interface import OllamaBackend
        
        ollama = OllamaBackend()
        
        return jsonify({
            'ollama': {
                'available': ollama.is_available(),
                'models': ollama.available_models if ollama.is_available() else []
            }
        })
        
    except Exception as e:
        logger.error(f"Error getting model status: {e}")
        return jsonify({'error': str(e)}), 500

# WebSocket events for real-time updates


@socketio.on("connect")
def handle_connect():
    """Handle client connection."""
    logger.info("Client connected")
    emit("connected", {"message": "Connected to GeneKnow Pipeline"})


@socketio.on("disconnect")
def handle_disconnect():
    """Handle client disconnection."""
    logger.info("Client disconnected")


@socketio.on("subscribe_job")
def handle_job_subscription(data):
    """Subscribe to job updates."""
    job_id = data.get("job_id")
    if job_id and job_id in jobs:
        # Join room for this job
        join_room(job_id)
        emit("subscribed", {"job_id": job_id})

        # Send current status
        with job_lock:
            job = jobs[job_id]
            emit(
                "job_progress",
                {
                    "job_id": job_id,
                    "status": job["status"],
                    "progress": job["progress"],
                    "current_step": job["current_step"],
                },
            )
    else:
        emit("error", {"message": "Invalid job ID"})


@socketio.on("unsubscribe_job")
def handle_job_unsubscription(data):
    """Unsubscribe from job updates."""
    job_id = data.get("job_id")
    if job_id:
        leave_room(job_id)
        emit("unsubscribed", {"job_id": job_id})


# Cleanup temporary files on shutdown
import atexit


def cleanup_temp_files():
    """Clean up temporary directories and result files."""
    try:
        # Clean up result files for all jobs
        with job_lock:
            for job_id in list(jobs.keys()):
                cleanup_job_files(job_id, cleanup_results=True)

        # Clean up temporary directories
        shutil.rmtree(Config.UPLOAD_FOLDER, ignore_errors=True)
        shutil.rmtree(Config.RESULTS_FOLDER, ignore_errors=True)
        logger.info("Cleaned up temporary files and result files")
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
                s.bind(("", port))
                s.close()
                return port
        except OSError:
            # Port is in use, try the next one
            continue

    raise RuntimeError(
        f"No available ports found in range {start_port}-{start_port + max_attempts}"
    )


# Main entry point
if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="GeneKnow Pipeline API Server")
    parser.add_argument("--port", type=int, help="Port to run the server on")
    parser.add_argument("--debug", action="store_true", help="Run in debug mode")

    args = parser.parse_args()

    # Find available port
    port = args.port if args.port else find_available_port()
    actual_port = port

    # Always announce the port to stdout for Rust backend to capture
    print(f"API_SERVER_PORT:{actual_port}", flush=True)

    # Start the server (Flask development server - for Gunicorn, this block won't run)
    print(f"Starting GeneKnow Pipeline API Server on port {actual_port}", flush=True)
    logger.info(f"Starting GeneKnow Pipeline API Server on port {actual_port}")
    debug = args.debug or os.environ.get("DEBUG", "").lower() == "true"

    # Print configuration for debugging
    print(f"Debug mode: {debug}", flush=True)
    print(f"Async mode: {socketio.async_mode}", flush=True)

    # Note: allow_unsafe_werkzeug is only needed for production mode with Flask dev server
    # Gunicorn doesn't need this flag
    socketio.run(
        app, host="0.0.0.0", port=actual_port, debug=debug, allow_unsafe_werkzeug=True
    )

    print("Server stopped", flush=True)

# Expose app for Gunicorn
application = app  # Some WSGI servers look for 'application'
