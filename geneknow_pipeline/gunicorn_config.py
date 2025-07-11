"""
Gunicorn configuration for GeneKnow API server
Optimized for single-user desktop application
"""

import multiprocessing
import os
import socket

# Dynamic port allocation
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
    raise RuntimeError("No available ports found")

# Server socket
# Try to get port from environment or find available
if os.environ.get('API_PORT'):
    port = int(os.environ.get('API_PORT'))
else:
    port = find_available_port()
    
bind = f"127.0.0.1:{port}"  # Only localhost for security
backlog = 64

# Always announce the port to stdout for Rust backend to capture
print(f"API_SERVER_PORT:{port}", flush=True)

# Worker processes
# For desktop app, 2-4 workers is optimal
workers = min(4, multiprocessing.cpu_count())
worker_class = "eventlet"  # Required for Flask-SocketIO
worker_connections = 100  # Low for single user
timeout = 120  # Longer timeout for genomic processing
keepalive = 2

# Logging
accesslog = "-"  # stdout
errorlog = "-"   # stderr
loglevel = "info"
access_log_format = '%(h)s %(l)s %(u)s %(t)s "%(r)s" %(s)s %(b)s "%(f)s" "%(a)s" %(D)s'

# Process naming
proc_name = 'geneknow-api'

# Server mechanics
daemon = False
pidfile = None
umask = 0
user = None
group = None
tmp_upload_dir = None

# SSL (disabled for local desktop app)
keyfile = None
certfile = None

# Disable stats
statsd_host = None

# Development helpers
reload = False
reload_engine = "auto"
reload_extra_files = []
spew = False

# Server hooks
def when_ready(server):
    server.log.info("GeneKnow API Server is ready. Listening at: %s", server.address)

def worker_int(worker):
    worker.log.info("Worker interrupted!")

def pre_fork(server, worker):
    server.log.info("Worker spawned (pid: %s)", worker.pid) 