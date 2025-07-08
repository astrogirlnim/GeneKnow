"""
Simplified API server for testing.
"""
from flask import Flask, jsonify, request
import json

app = Flask(__name__)

@app.route('/health', methods=['GET'])
def health():
    return jsonify({"status": "healthy", "service": "GeneKnow Pipeline API"})

@app.route('/test', methods=['POST'])
def test_endpoint():
    data = request.get_json()
    return jsonify({
        "message": "Test successful",
        "received": data
    })

if __name__ == '__main__':
    print("Starting simple server on http://localhost:5001")
    app.run(host='0.0.0.0', port=5001, debug=False) 