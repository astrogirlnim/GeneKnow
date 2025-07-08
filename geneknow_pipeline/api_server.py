"""
API server for the genomic risk assessment pipeline.
Provides REST endpoints for processing genomic files.
"""
from flask import Flask, request, jsonify
import os
import tempfile
from datetime import datetime
from graph import run_pipeline
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)

# Configuration
MAX_FILE_SIZE = 5 * 1024 * 1024 * 1024  # 5GB max file size
ALLOWED_EXTENSIONS = {'.fastq', '.fq', '.fastq.gz', '.fq.gz', '.bam', '.vcf', '.vcf.gz', '.maf'}


@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "service": "genomic-risk-assessment"
    })


@app.route('/process', methods=['POST'])
def process_file():
    """
    Process a genomic file through the pipeline.
    
    Expected JSON payload:
    {
        "file_path": "/path/to/file.vcf",
        "file_type": "vcf",  # Optional, will be auto-detected
        "preferences": {
            "language": "en",
            "include_technical": true
        }
    }
    """
    try:
        # Get request data
        data = request.get_json()
        if not data:
            return jsonify({
                "status": "error",
                "message": "No JSON data provided"
            }), 400
        
        # Validate required fields
        file_path = data.get('file_path')
        if not file_path:
            return jsonify({
                "status": "error", 
                "message": "file_path is required"
            }), 400
        
        # Check file exists
        if not os.path.exists(file_path):
            return jsonify({
                "status": "error",
                "message": f"File not found: {file_path}"
            }), 404
        
        # Check file extension
        file_ext = os.path.splitext(file_path.lower())[1]
        if file_ext not in ALLOWED_EXTENSIONS:
            return jsonify({
                "status": "error",
                "message": f"Unsupported file type. Allowed: {', '.join(ALLOWED_EXTENSIONS)}"
            }), 400
        
        # Check file size
        file_size = os.path.getsize(file_path)
        if file_size > MAX_FILE_SIZE:
            return jsonify({
                "status": "error",
                "message": f"File too large. Maximum size: {MAX_FILE_SIZE / (1024**3):.1f}GB"
            }), 413
        
        # Get preferences
        preferences = data.get('preferences', {})
        
        # Log the request
        logger.info(f"Processing file: {file_path}")
        logger.info(f"File size: {file_size / (1024**2):.2f}MB")
        logger.info(f"Preferences: {preferences}")
        
        # Run the pipeline
        result = run_pipeline(file_path, preferences)
        
        # Check pipeline status
        if result.get('pipeline_status') == 'completed':
            # Successful completion
            response = {
                "status": "success",
                "processing_time": result.get('processing_time_seconds', 0),
                "results": {
                    "variant_count": result.get('variant_count', 0),
                    "risk_scores": result.get('risk_scores', {}),
                    "report_sections": result.get('report_sections', {}),
                    "file_metadata": result.get('file_metadata', {}),
                    "completed_nodes": result.get('completed_nodes', [])
                }
            }
            
            # Add warnings if any
            if result.get('warnings'):
                response['warnings'] = result['warnings']
            
            return jsonify(response), 200
            
        else:
            # Pipeline failed
            errors = result.get('errors', [])
            error_message = "Pipeline failed"
            if errors:
                error_message = errors[-1].get('error', error_message)
            
            return jsonify({
                "status": "error",
                "message": error_message,
                "errors": errors,
                "completed_nodes": result.get('completed_nodes', [])
            }), 500
            
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        return jsonify({
            "status": "error",
            "message": f"Internal server error: {str(e)}"
        }), 500


@app.route('/info', methods=['GET'])
def pipeline_info():
    """Get information about the pipeline capabilities."""
    return jsonify({
        "name": "Genomic Risk Assessment Pipeline",
        "version": "1.0.0",
        "supported_formats": list(ALLOWED_EXTENSIONS),
        "max_file_size_gb": MAX_FILE_SIZE / (1024**3),
        "pipeline_nodes": [
            "file_input",
            "preprocess", 
            "variant_calling",
            "qc_filter",
            "tcga_mapper",
            "risk_model",
            "formatter",
            "report_writer"
        ],
        "supported_languages": ["en", "es", "hi"],
        "cancer_types": ["breast", "colon", "lung", "prostate", "blood"]
    })


if __name__ == '__main__':
    # Run the server
    port = int(os.environ.get('PORT', 5001))
    debug = os.environ.get('DEBUG', 'False').lower() == 'true'
    
    logger.info(f"Starting API server on port {port}")
    logger.info(f"Debug mode: {debug}")
    logger.info(f"Supported file types: {', '.join(ALLOWED_EXTENSIONS)}")
    
    app.run(host='0.0.0.0', port=port, debug=debug) 