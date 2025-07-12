#!/usr/bin/env python3
"""
Run the GeneKnow API server with Gunicorn
This wrapper ensures proper initialization and WebSocket support
"""

import sys
import subprocess


def main():
    # Set up Gunicorn command
    gunicorn_cmd = [
        sys.executable,
        "-m",
        "gunicorn",
        "--config",
        "gunicorn_config.py",
        "enhanced_api_server:app",  # Note: using app, not socketio
        "--worker-class",
        "eventlet",
        "-w",
        "2",  # 2 workers for desktop app
        "--timeout",
        "120",
        "--log-level",
        "info",
    ]

    # Run Gunicorn
    try:
        subprocess.run(gunicorn_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Gunicorn: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
