"""
Simple API server for GeneKnow LangGraph pipeline.
Provides REST endpoints for Tauri to call.
"""
from flask import Flask, request, jsonify
from flask_cors import CORS
import os
import tempfile
from werkzeug.utils import secure_filename
from graph import run_pipeline
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)  # Enable CORS for Tauri

# Configure upload settings
UPLOAD_FOLDER = tempfile.gettempdir()
ALLOWED_EXTENSIONS = {'fastq', 'fq', 'bam', 'sam', 'gz'}

def allowed_file(filename):
    return '.' in filename and \
           any(filename.lower().endswith(ext) for ext in ALLOWED_EXTENSIONS)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({"status": "healthy", "service": "GeneKnow Pipeline API"})

@app.route('/analyze', methods=['POST'])
def analyze_genomic_file():
    """
    Main endpoint for genomic analysis.
    Accepts a file upload and user preferences.
    """
    try:
        # Check if file was uploaded
        if 'file' not in request.files:
            return jsonify({"error": "No file provided"}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({"error": "No file selected"}), 400
        
        if not allowed_file(file.filename):
            return jsonify({"error": "Invalid file type. Allowed: FASTQ, BAM, SAM"}), 400
        
        # Save uploaded file temporarily
        filename = secure_filename(file.filename)
        filepath = os.path.join(UPLOAD_FOLDER, filename)
        file.save(filepath)
        
        # Get user preferences from request
        preferences = request.form.get('preferences', {})
        if isinstance(preferences, str):
            import json
            preferences = json.loads(preferences)
        
        logger.info(f"Processing file: {filename}")
        
        # Run the pipeline
        result = run_pipeline(filepath, preferences)
        
        # Clean up temporary file
        os.remove(filepath)
        
        # Return results
        if result['pipeline_status'] == 'completed':
            return jsonify({
                "status": "success",
                "data": result['structured_json'],
                "report": result.get('report_markdown', ''),
                "report_sections": result.get('report_sections', {}),
                "processing_time": result.get('processing_time_seconds', 0)
            })
        else:
            return jsonify({
                "status": "failed",
                "errors": result.get('errors', []),
                "warnings": result.get('warnings', []),
                "completed_nodes": result.get('completed_nodes', [])
            }), 500
            
    except Exception as e:
        logger.error(f"API error: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/analyze_path', methods=['POST'])
def analyze_genomic_path():
    """
    Alternative endpoint that accepts a file path instead of upload.
    Useful for testing or when file is already on the server.
    """
    try:
        data = request.get_json()
        if not data or 'file_path' not in data:
            return jsonify({"error": "No file_path provided"}), 400
        
        file_path = data['file_path']
        preferences = data.get('preferences', {})
        
        if not os.path.exists(file_path):
            return jsonify({"error": f"File not found: {file_path}"}), 404
        
        logger.info(f"Processing file path: {file_path}")
        
        # Run the pipeline
        result = run_pipeline(file_path, preferences)
        
        # Return results
        if result['pipeline_status'] == 'completed':
            return jsonify({
                "status": "success",
                "data": result['structured_json'],
                "report": result.get('report_markdown', ''),
                "report_sections": result.get('report_sections', {}),
                "processing_time": result.get('processing_time_seconds', 0)
            })
        else:
            return jsonify({
                "status": "failed", 
                "errors": result.get('errors', []),
                "warnings": result.get('warnings', []),
                "completed_nodes": result.get('completed_nodes', [])
            }), 500
            
    except Exception as e:
        logger.error(f"API error: {str(e)}")
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    print("🧬 GeneKnow Pipeline API Server")
    print("=" * 40)
    print("Endpoints:")
    print("  GET  /health         - Health check")
    print("  POST /analyze        - Upload file for analysis")  
    print("  POST /analyze_path   - Analyze file by path")
    print("=" * 40)
    print("Starting server on http://localhost:5001")
    print("Press Ctrl+C to stop")
    
    # Run without debug mode to avoid hanging issues
    app.run(host='0.0.0.0', port=5001, debug=False) 