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
