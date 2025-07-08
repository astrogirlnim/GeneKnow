# 🚀 GeneKnow Server Startup Guide

## 🎯 Quick Start

Get GeneKnow running in 3 steps:

```bash
# 1. Start the API Server
cd geneknow_pipeline
python enhanced_api_server.py

# 2. Start the Desktop App (in new terminal)
cd desktop/ui
npm run dev

# 3. Open your browser to http://localhost:5173
```

---

## 📋 Prerequisites

### **System Requirements**
- **Python 3.8+** with pip
- **Node.js 16+** with npm
- **Git** (for cloning repository)
- **4GB+ RAM** for genomic processing

### **Check Your Environment**
```bash
# Verify Python version
python --version  # Should show 3.8+

# Verify Node.js version
node --version    # Should show 16+

# Verify npm version
npm --version     # Should show 8+
```

---

## 🛠️ Installation & Setup

### **1. Clone Repository**
```bash
git clone <repository-url>
cd LiteratureGapper
```

### **2. Setup Python Environment**
```bash
# Navigate to pipeline directory
cd geneknow_pipeline

# Create virtual environment
python -m venv venv

# Activate virtual environment
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install Python dependencies
pip install -r requirements.txt
```

### **3. Setup Node.js Environment**
```bash
# Navigate to UI directory
cd desktop/ui

# Install Node.js dependencies
npm install
```

### **4. Verify Installation**
```bash
# Test Python environment
cd geneknow_pipeline
python -c "import flask, flask_cors, flask_socketio; print('Python deps OK')"

# Test Node.js environment
cd desktop/ui
npm run build  # Should complete without errors
```

---

## 🚀 Starting the Servers

### **Method 1: Development Mode (Recommended)**

#### **Terminal 1: Start API Server**
```bash
cd geneknow_pipeline

# Activate Python virtual environment
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Start the enhanced API server
python enhanced_api_server.py

# You should see:
# INFO:__main__:Starting Enhanced GeneKnow API Server on port 5001
# INFO:__main__:Upload folder: /tmp/geneknow_uploads_xxx
# INFO:__main__:Results folder: /tmp/geneknow_results_xxx
# INFO:__main__:WebSocket support enabled
# * Running on all addresses (0.0.0.0)
# * Running on http://127.0.0.1:5001
# * Running on http://[::1]:5001
```

#### **Terminal 2: Start Desktop App**
```bash
cd desktop/ui

# Start development server
npm run dev

# You should see:
# ➜  Local:   http://localhost:5173/
# ➜  Network: use --host to expose
# ➜  press h + enter to show help
```

#### **Terminal 3: Start Tauri App (Optional)**
```bash
cd desktop/ui

# Start Tauri development app
npm run tauri:dev

# This will open the native desktop application
```

### **Method 2: Production Mode**

#### **API Server Production**
```bash
cd geneknow_pipeline
source venv/bin/activate

# Set production environment
export DEBUG=false
export PORT=5001

# Start with gunicorn (install if needed)
pip install gunicorn
gunicorn -w 4 -b 0.0.0.0:5001 enhanced_api_server:app
```

#### **Desktop App Production**
```bash
cd desktop/ui

# Build production version
npm run build

# Build Tauri app
npm run tauri:build

# The built app will be in src-tauri/target/release/
```

---

## 🔧 Configuration

### **Environment Variables**

#### **API Server Configuration**
```bash
# Port (default: 5001)
export PORT=5001

# Debug mode (default: false)
export DEBUG=true

# CORS origins (default: tauri://localhost)
export CORS_ORIGINS="http://localhost:5173,tauri://localhost"
```

#### **Desktop App Configuration**
```bash
# API URL (default: http://localhost:5001)
export VITE_API_URL=http://localhost:5001

# WebSocket URL (default: http://localhost:5001)
export VITE_WS_URL=http://localhost:5001
```

### **Configuration Files**

#### **API Server: `geneknow_pipeline/enhanced_api_server.py`**
```python
# Configuration class (lines 58-64)
class Config:
    MAX_FILE_SIZE = 5 * 1024 * 1024 * 1024  # 5GB
    ALLOWED_EXTENSIONS = {'.fastq', '.fq', '.fastq.gz', '.fq.gz', '.bam', '.vcf', '.vcf.gz', '.maf', '.maf.gz'}
    UPLOAD_FOLDER = tempfile.mkdtemp(prefix='geneknow_uploads_')
    RESULTS_FOLDER = tempfile.mkdtemp(prefix='geneknow_results_')
    SESSION_TIMEOUT = 3600  # 1 hour
```

#### **Desktop App: `desktop/ui/vite.config.ts`**
```typescript
// Development server configuration
export default defineConfig({
  server: {
    port: 5173,
    host: true
  },
  // ... other config
})
```

---

## 🧪 Testing the Setup

### **1. Test API Server**
```bash
# Health check
curl http://localhost:5001/api/health

# Expected response:
# {
#   "status": "healthy",
#   "timestamp": "2024-01-15T10:30:00Z",
#   "service": "GeneKnow Pipeline API",
#   "version": "2.0.0",
#   "jobs_active": 0
# }
```

### **2. Test Pipeline Info**
```bash
# Get pipeline capabilities
curl http://localhost:5001/api/pipeline-info

# Should return pipeline nodes and capabilities
```

### **3. Test Desktop App**
```bash
# Open browser to
http://localhost:5173

# You should see the GeneKnow interface
```

### **4. Test File Processing**
```bash
# Generate demo files first
python setup_demo_files.py

# Test with sample VCF file
curl -X POST http://localhost:5001/api/process \
  -H "Content-Type: application/json" \
  -d '{"file_path": "demo_files/patient_variants.vcf"}'

# Should return a job_id
```

---

## 🐛 Troubleshooting

### **Common Issues**

#### **🔴 "Port 5001 already in use"**
```bash
# Check what's using the port
lsof -i :5001

# Kill the process
kill -9 <PID>

# Or use a different port
PORT=5002 python enhanced_api_server.py
```

#### **🔴 "Module not found" errors**
```bash
# Make sure virtual environment is activated
source venv/bin/activate

# Reinstall dependencies
pip install -r requirements.txt

# For Node.js modules
cd desktop/ui
rm -rf node_modules
npm install
```

#### **🔴 "Permission denied" errors**
```bash
# Fix Python permissions
chmod +x enhanced_api_server.py

# Fix Node.js permissions
sudo chown -R $(whoami) node_modules
```

#### **🔴 "CORS errors" in browser**
```bash
# Check API server CORS configuration
# Make sure CORS_ORIGINS includes your frontend URL

# Or disable CORS temporarily for testing
CORS_ORIGINS="*" python enhanced_api_server.py
```

#### **🔴 "Connection refused" errors**
```bash
# Check if API server is running
curl http://localhost:5001/api/health

# Check firewall settings
sudo ufw status  # Linux
# Make sure port 5001 is open

# Check localhost resolution
ping localhost
```

#### **🔴 "File upload fails"**
```bash
# Check upload directory permissions
ls -la /tmp/geneknow_uploads_*

# Check file size limits
# Default is 5GB, see Config.MAX_FILE_SIZE

# Check file format
file demo_files/patient_variants.vcf
```

### **Debug Mode**

#### **Enable API Debug Mode**
```bash
cd geneknow_pipeline
DEBUG=true python enhanced_api_server.py

# More verbose logging will appear
```

#### **Enable Frontend Debug Mode**
```bash
cd desktop/ui
npm run dev

# Open browser dev tools (F12)
# Check console for error messages
```

### **Log Files**

#### **API Server Logs**
```bash
# View real-time logs
tail -f geneknow_pipeline/logs/api.log

# Check for errors
grep ERROR geneknow_pipeline/logs/api.log
```

#### **Desktop App Logs**
```bash
# Browser console (F12)
# Check for JavaScript errors

# Tauri app logs (if using)
# Check terminal output when running npm run tauri:dev
```

---

## 📊 Monitoring & Health Checks

### **API Server Monitoring**
```bash
# Health endpoint
curl http://localhost:5001/api/health

# Check active jobs
curl http://localhost:5001/api/jobs

# Server statistics
curl http://localhost:5001/api/pipeline-info
```

### **System Resources**
```bash
# Check memory usage
free -h

# Check CPU usage
top -p $(pgrep -f enhanced_api_server.py)

# Check disk space
df -h
```

### **Network Monitoring**
```bash
# Check listening ports
netstat -tuln | grep :5001

# Monitor network traffic
sudo tcpdump -i lo port 5001
```

---

## 🔄 Restarting Services

### **Graceful Restart**
```bash
# Stop API server (Ctrl+C in terminal)
# Then restart
python enhanced_api_server.py

# Stop desktop app (Ctrl+C in terminal)
# Then restart
npm run dev
```

### **Force Restart**
```bash
# Kill API server process
pkill -f enhanced_api_server.py

# Kill Node.js processes
pkill -f "node.*vite"

# Then restart normally
```

---

## 🚀 Production Deployment

### **API Server Production**
```bash
# Install production server
pip install gunicorn

# Create systemd service (Linux)
sudo tee /etc/systemd/system/geneknow-api.service > /dev/null <<EOF
[Unit]
Description=GeneKnow API Server
After=network.target

[Service]
Type=simple
User=geneknow
WorkingDirectory=/path/to/geneknow_pipeline
Environment=PATH=/path/to/venv/bin
ExecStart=/path/to/venv/bin/gunicorn -w 4 -b 0.0.0.0:5001 enhanced_api_server:app
Restart=always

[Install]
WantedBy=multi-user.target
EOF

# Start service
sudo systemctl enable geneknow-api
sudo systemctl start geneknow-api
```

### **Desktop App Production**
```bash
# Build production version
npm run build

# Build Tauri app
npm run tauri:build

# Distribute the built application
ls src-tauri/target/release/
```

---

## 📝 Quick Reference

### **Essential Commands**
```bash
# Start everything
cd geneknow_pipeline && python enhanced_api_server.py &
cd desktop/ui && npm run dev &

# Stop everything
pkill -f enhanced_api_server.py
pkill -f "node.*vite"

# Health check
curl http://localhost:5001/api/health

# View logs
tail -f geneknow_pipeline/logs/api.log
```

### **Default URLs**
- **API Server**: http://localhost:5001
- **Desktop App**: http://localhost:5173
- **API Documentation**: http://localhost:5001/api/pipeline-info
- **Health Check**: http://localhost:5001/api/health

### **Default Ports**
- **API Server**: 5001
- **Desktop App**: 5173
- **WebSocket**: 5001 (same as API)

---

## ✅ Success Indicators

When everything is working correctly, you should see:

1. **API Server**: 
   - ✅ "Starting Enhanced GeneKnow API Server on port 5001"
   - ✅ Health check returns `{"status": "healthy"}`

2. **Desktop App**: 
   - ✅ "Local: http://localhost:5173/"
   - ✅ Browser shows GeneKnow interface

3. **Integration**: 
   - ✅ File uploads work
   - ✅ Processing completes
   - ✅ Results display correctly

**🎉 Ready to analyze genomes locally and privately!**

---

## 🆘 Need Help?

If you're still having issues:

1. **Check Prerequisites**: Verify Python 3.8+ and Node.js 16+
2. **Follow Setup Steps**: Don't skip the virtual environment
3. **Check Logs**: Look for specific error messages
4. **Test Components**: Verify each piece works individually
5. **Ask for Help**: Share specific error messages and steps taken

**Remember**: GeneKnow processes genetic data locally for maximum privacy and security! 🧬🔒 