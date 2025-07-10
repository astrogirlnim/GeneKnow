#!/bin/bash
# test-dynamic-port.sh - Test the dynamic port allocation for GeneKnow API server

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLE_DIR="$SCRIPT_DIR/../bundled_resources"

echo "🧪 Testing Dynamic Port Allocation"
echo "================================="

# Check if bundle exists
if [ ! -d "$BUNDLE_DIR" ]; then
    echo "❌ Bundle directory not found at: $BUNDLE_DIR"
    echo "   Please run ./bundle-python-optimized.sh first"
    exit 1
fi

# Test 1: Start server and capture port
echo ""
echo "1️⃣ Testing server startup with port capture..."

cd "$BUNDLE_DIR"
PYTHON_EXE="$BUNDLE_DIR/python_runtime/bin/python3"

# Start the server and capture output
echo "Starting server..."
PORT_FILE="$BUNDLE_DIR/.test_api_port"
rm -f "$PORT_FILE"

# Run the server in background and capture output
"$PYTHON_EXE" -m geneknow_pipeline.enhanced_api_server --port-file "$PORT_FILE" > server.log 2>&1 &
SERVER_PID=$!

echo "Server PID: $SERVER_PID"
echo "Waiting for port file..."

# Wait for port file to be created
TIMEOUT=10
COUNTER=0
while [ ! -f "$PORT_FILE" ] && [ $COUNTER -lt $TIMEOUT ]; do
    sleep 1
    COUNTER=$((COUNTER + 1))
    echo "  Waiting... ($COUNTER/$TIMEOUT)"
done

if [ -f "$PORT_FILE" ]; then
    PORT=$(cat "$PORT_FILE")
    echo "✅ Server started on port: $PORT"
    
    # Test 2: Check health endpoint
    echo ""
    echo "2️⃣ Testing health endpoint..."
    sleep 2  # Give server time to fully start
    
    if curl -s "http://localhost:$PORT/api/health" | python3 -m json.tool; then
        echo "✅ Health check passed!"
    else
        echo "❌ Health check failed"
    fi
else
    echo "❌ Port file not created within timeout"
    cat server.log
fi

# Cleanup
echo ""
echo "🧹 Cleaning up..."
kill $SERVER_PID 2>/dev/null || true
rm -f "$PORT_FILE" server.log

echo ""
echo "✅ Test complete!" 