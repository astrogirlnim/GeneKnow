#!/bin/bash

# 🦀 Rust-Python Integration Test
# Tests the actual execute_python helper function

set -e

echo "🦀 Testing Rust execute_python Integration"
echo "=========================================="

cd src-tauri

# Create a simple test in Rust that calls our execute_python function
cat > src/test_integration.rs << 'EOF'
use crate::utils::execute_python;
use std::path::PathBuf;

pub fn test_python_integration() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔍 Testing execute_python with generate_test_fastq.py...");
    
    // Test calling generate_test_fastq.py
    let args = vec![
        "--output-dir", "test_rust_output",
        "--num-reads", "3",
        "--json"
    ];
    
    let result = execute_python("generate_test_fastq.py", &args)?;
    
    if result.status.success() {
        let stdout = String::from_utf8_lossy(&result.stdout);
        println!("✅ Python script executed successfully");
        println!("📄 Output: {}", stdout.trim());
        
        // Try to parse as JSON
        if let Ok(json) = serde_json::from_str::<serde_json::Value>(&stdout) {
            println!("✅ JSON parsing successful");
            if let Some(success) = json.get("success") {
                if success.as_bool() == Some(true) {
                    println!("✅ Python script reported success");
                    return Ok(());
                }
            }
        }
        
        Err("JSON parsing failed or success=false".into())
    } else {
        let stderr = String::from_utf8_lossy(&result.stderr);
        eprintln!("❌ Python script failed: {}", stderr);
        Err("Python execution failed".into())
    }
}
EOF

# Add the test module to lib.rs
if ! grep -q "mod test_integration;" src/lib.rs; then
    echo "" >> src/lib.rs
    echo "#[cfg(test)]" >> src/lib.rs
    echo "mod test_integration;" >> src/lib.rs
fi

# Add a simple test binary
cat > src/bin/test_integration.rs << 'EOF'
use app::utils::execute_python;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🦀 Rust-Python Integration Test");
    println!("===============================");
    
    println!("🔍 Testing execute_python with config_data_source.py...");
    
    let result = execute_python("config_data_source.py", &["--list-vcf-files", "--json"])?;
    
    if result.status.success() {
        let stdout = String::from_utf8_lossy(&result.stdout);
        println!("✅ Python script executed successfully");
        
        // Parse JSON and count entries
        if let Ok(json) = serde_json::from_str::<serde_json::Value>(&stdout) {
            if let Some(obj) = json.as_object() {
                println!("✅ JSON parsing successful - {} VCF files listed", obj.len());
            }
        }
    } else {
        let stderr = String::from_utf8_lossy(&result.stderr);
        eprintln!("❌ Python script failed: {}", stderr);
        return Err("Python execution failed".into());
    }
    
    println!("🎉 Rust-Python integration test PASSED!");
    Ok(())
}
EOF

echo "🔨 Building test binary..."
cargo build --bin test_integration

echo "🚀 Running integration test..."
./target/debug/test_integration

echo ""
echo "✅ Rust-Python integration verified!"

# Cleanup test files
rm -f src/test_integration.rs src/bin/test_integration.rs
echo "🧹 Cleaned up test files" 