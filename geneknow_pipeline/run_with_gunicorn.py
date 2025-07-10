#!/usr/bin/env python3
"""
Run the GeneKnow API server with Gunicorn
This wrapper ensures proper initialization and WebSocket support
"""

import os
import sys
import subprocess

def main():
    # Get port file from command line if provided
    port_file = None
    if '--port-file' in sys.argv:
        idx = sys.argv.index('--port-file')
        if idx + 1 < len(sys.argv):
            port_file = sys.argv[idx + 1]
            os.environ['PORT_FILE'] = port_file
    
    # Set up Gunicorn command
    gunicorn_cmd = [
        sys.executable, '-m', 'gunicorn',
        '--config', 'gunicorn_config.py',
        'enhanced_api_server:app',  # Note: using app, not socketio
        '--worker-class', 'eventlet',
        '-w', '2',  # 2 workers for desktop app
        '--timeout', '120',
        '--log-level', 'info'
    ]
    
    # Run Gunicorn
    try:
        subprocess.run(gunicorn_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Gunicorn: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main() 