# 🦄 Gunicorn Deployment for GenePredict

## Overview

GenePredict now supports **Gunicorn** as a production-ready WSGI server alternative to Flask's development server. This provides better performance, stability, and proper worker management for the desktop application.

## 🚀 Benefits of Gunicorn

1. **Production-Ready**: Designed for production use, unlike Flask's development server
2. **Better Performance**: Multi-worker architecture for handling concurrent requests
3. **Stable Memory Usage**: Proper worker recycling prevents memory leaks
4. **WebSocket Support**: Works with eventlet for real-time progress updates
5. **No More Warnings**: Eliminates the "unsafe Werkzeug" warning

## 📋 Implementation Details

### Files Added/Modified

1. **`geneknow_pipeline/requirements-lite.txt`**
   - Added `gunicorn>=21.2.0` 
   - Added `eventlet>=0.33.0` for WebSocket support

2. **`geneknow_pipeline/gunicorn_config.py`**
   - Production-optimized configuration
   - Dynamic port allocation
   - 2-4 workers for desktop app
   - Eventlet worker class for WebSocket support

3. **`desktop/bundled_resources/start_api_server.sh`**
   - Updated to support both Flask and Gunicorn
   - Defaults to Flask for compatibility
   - Can enable Gunicorn with `USE_GUNICORN=true`

4. **`geneknow_pipeline/enhanced_api_server.py`**
   - Exposed `app` object for Gunicorn
   - Added `application = app` alias
   - Specified eventlet async mode

## 🔧 Configuration

### Gunicorn Settings (Optimized for Desktop App)

```python
workers = min(4, multiprocessing.cpu_count())  # 2-4 workers
worker_class = "eventlet"                       # WebSocket support
worker_connections = 100                        # Low for single user
timeout = 120                                   # Long for genomic processing
bind = "127.0.0.1:PORT"                        # Localhost only
```

### Dynamic Port Allocation

Gunicorn config automatically:
1. Finds an available port starting from 5000
2. Writes the port to `.api_server_port` file
3. Outputs `API_SERVER_PORT:XXXX` for Tauri to capture

## 🏃 Usage

### Enable Gunicorn in Production

Set the environment variable before building:
```bash
export USE_GUNICORN=true
./bundle-python-optimized.sh
```

### Manual Testing

```bash
# Test with Gunicorn
cd desktop/bundled_resources
USE_GUNICORN=true ./start_api_server.sh

# Test with Flask (default)
./start_api_server.sh
```

### Verify It's Working

When Gunicorn is running, you'll see:
```
🦄 Using Gunicorn production server...
[2025-01-10 16:30:00 +0000] [12345] [INFO] Starting gunicorn 21.2.0
[2025-01-10 16:30:00 +0000] [12345] [INFO] Listening at: http://127.0.0.1:5000 (12345)
[2025-01-10 16:30:00 +0000] [12345] [INFO] Using worker: eventlet
[2025-01-10 16:30:00 +0000] [12346] [INFO] Booting worker with pid: 12346
[2025-01-10 16:30:00 +0000] [12347] [INFO] Booting worker with pid: 12347
```

## 🔄 Migration Path

### Phase 1: Current State (Default Flask)
- Flask development server with `allow_unsafe_werkzeug=True`
- Works but shows warnings in production
- Single process, adequate for desktop app

### Phase 2: Optional Gunicorn (Now Available)
- Set `USE_GUNICORN=true` to enable
- Automatic fallback to Flask if Gunicorn not installed
- Fully backward compatible

### Phase 3: Default Gunicorn (Future)
- After thorough testing, make Gunicorn the default
- Keep Flask as fallback for development

## 🧪 Performance Comparison

| Metric | Flask Dev Server | Gunicorn |
|--------|-----------------|----------|
| Concurrent Requests | 1 (threaded) | 2-4 (workers) |
| Memory Stability | May grow over time | Stable with recycling |
| WebSocket Support | Basic | Optimized with eventlet |
| Production Ready | No (with warnings) | Yes |
| Startup Time | ~1s | ~2-3s |

## 🐛 Troubleshooting

### Gunicorn Not Starting
```bash
# Check if installed
python3 -c "import gunicorn; print(gunicorn.__version__)"

# If not, install it
pip install gunicorn eventlet
```

### Port Conflicts
```bash
# Gunicorn automatically finds next available port
# Check the log file for actual port:
cat bundled_resources/logs/api_server_*.log | grep "Listening at"
```

### WebSocket Issues
- Ensure eventlet is installed
- Check that `async_mode='eventlet'` is set in SocketIO initialization
- Verify `worker_class = "eventlet"` in gunicorn config

## 🔒 Security Notes

- Gunicorn binds to `127.0.0.1` only (localhost)
- No external access possible
- Workers run as the current user
- No SSL/TLS needed for local desktop app

## 📝 Summary

Gunicorn provides a production-ready alternative to Flask's development server, eliminating warnings and improving stability. The implementation is:

1. **Backward Compatible**: Defaults to Flask if not enabled
2. **Easy to Enable**: Just set `USE_GUNICORN=true`
3. **Desktop Optimized**: Configuration tuned for single-user app
4. **Future Proof**: Ready to become the default server 