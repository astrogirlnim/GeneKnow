// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(
  all(not(debug_assertions), target_os = "windows"),
  windows_subsystem = "windows"
)]

// Learn more about Tauri commands at https://tauri.app/develop/calling-rust/
use tauri::{command, Emitter};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::process::Command;
use std::fs;
use std::thread;
use std::time::Duration;
use std::path::PathBuf;

pub mod utils;
pub mod plugin;
pub mod python_script_plugin;
pub mod plugin_registry;

use utils::execute_python;
use plugin::PluginSummary;
use plugin_registry::PluginRegistryApi;

// Add these imports for managing the API server
use std::process::{Child, Stdio};
use std::sync::Mutex;
use once_cell::sync::Lazy;

// Global API server process handle and port
static API_SERVER_PROCESS: Lazy<Mutex<Option<Child>>> = Lazy::new(|| Mutex::new(None));
static API_SERVER_PORT: Lazy<Mutex<Option<u16>>> = Lazy::new(|| Mutex::new(None));

#[derive(Debug, Serialize, Deserialize)]
pub struct FastqToVcfOptions {
    reference: String,
    fastq1: String,
    fastq2: Option<String>,
    output_prefix: String,
    threads: Option<u32>,
    aligner: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FastqToVcfResult {
    success: bool,
    vcf_file: Option<String>,
    bam_file: Option<String>,
    variant_count: Option<u32>,
    error: Option<String>,
    execution_time: Option<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ExtractRegionOptions {
    bed_file: String,
    output_dir: String,
    max_processes: Option<u32>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ExtractRegionResult {
    success: bool,
    total_regions: u32,
    successful_extractions: u32,
    total_variants: u32,
    total_size: u32,
    execution_time: f64,
    results: Vec<String>,
    error: Option<String>,
}

#[command]
async fn convert_fastq_to_vcf(options: FastqToVcfOptions) -> Result<FastqToVcfResult, String> {
    // Build arguments for the Python script
    let mut args = vec![
        "-r", &options.reference,
        "-1", &options.fastq1,
        "-o", &options.output_prefix,
        "--json", // Request JSON output
    ];
    
    // Add optional arguments
    let fastq2_str;
    if let Some(fastq2) = &options.fastq2 {
        fastq2_str = fastq2.clone();
        args.push("-2");
        args.push(&fastq2_str);
    }
    
    let threads_str;
    if let Some(threads) = options.threads {
        threads_str = threads.to_string();
        args.push("-t");
        args.push(&threads_str);
    }
    
    let aligner_str;
    if let Some(aligner) = &options.aligner {
        aligner_str = aligner.clone();
        args.push("-a");
        args.push(&aligner_str);
    }
    
    // Execute the Python script using our utility function
    let start_time = std::time::Instant::now();
    let output = execute_python("fastq_to_vcf_pipeline", &args)
        .map_err(|e| e.to_string())?;
    let execution_time = start_time.elapsed().as_secs_f64();
    
    if output.status.success() {
        // Parse JSON output
        let stdout = String::from_utf8_lossy(&output.stdout);
        
        // Try to parse as JSON, fall back to legacy parsing if needed
        if let Ok(json_result) = serde_json::from_str::<serde_json::Value>(&stdout) {
            // Parse JSON result
            let success = json_result.get("success").and_then(|v| v.as_bool()).unwrap_or(false);
            let vcf_file = json_result.get("vcf_file").and_then(|v| v.as_str()).map(|s| s.to_string());
            let bam_file = json_result.get("bam_file").and_then(|v| v.as_str()).map(|s| s.to_string());
            let variant_count = json_result.get("variant_count").and_then(|v| v.as_u64()).map(|v| v as u32);
            let error = json_result.get("error").and_then(|v| v.as_str()).map(|s| s.to_string());
            
            Ok(FastqToVcfResult {
                success,
                vcf_file,
                bam_file,
                variant_count,
                error,
                execution_time: Some(execution_time),
            })
        } else {
            // Legacy parsing fallback
            let variant_count = stdout
                .lines()
                .find(|line| line.contains("Variants found:"))
                .and_then(|line| line.split(':').last())
                .and_then(|num| num.trim().parse::<u32>().ok());
            
            Ok(FastqToVcfResult {
                success: true,
                vcf_file: Some(format!("{}.vcf", options.output_prefix)),
                bam_file: Some(format!("{}.bam", options.output_prefix)),
                variant_count,
                error: None,
                execution_time: Some(execution_time),
            })
        }
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        Err(format!("Pipeline failed: {}", stderr))
    }
}

#[command]
async fn extract_genomic_regions(options: ExtractRegionOptions) -> Result<ExtractRegionResult, String> {
    // Build arguments for the Python script
    let mut args = vec![
        "--bed-file", &options.bed_file,
        "--output-dir", &options.output_dir,
        "--json", // Request JSON output
    ];
    
    // Add optional max_processes argument
    let max_processes_str;
    if let Some(max_processes) = options.max_processes {
        max_processes_str = max_processes.to_string();
        args.push("--max-processes");
        args.push(&max_processes_str);
    }
    
    // Execute the Python script using our utility function
    let start_time = std::time::Instant::now();
    let output = execute_python("extract_by_region", &args)
        .map_err(|e| e.to_string())?;
    let execution_time = start_time.elapsed().as_secs_f64();
    
    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        
        // Try to parse as JSON, fall back to legacy parsing if needed
        if let Ok(json_result) = serde_json::from_str::<serde_json::Value>(&stdout) {
            // Parse JSON result
            let success = json_result.get("success").and_then(|v| v.as_bool()).unwrap_or(false);
            let total_regions = json_result.get("total_regions").and_then(|v| v.as_u64()).unwrap_or(0) as u32;
            let successful_extractions = json_result.get("successful_extractions").and_then(|v| v.as_u64()).unwrap_or(0) as u32;
            let total_variants = json_result.get("total_variants").and_then(|v| v.as_u64()).unwrap_or(0) as u32;
            let total_size = json_result.get("total_size").and_then(|v| v.as_u64()).unwrap_or(0) as u32;
            let results = json_result.get("results")
                .and_then(|v| v.as_array())
                .map(|arr| arr.iter().filter_map(|v| v.as_str()).map(|s| s.to_string()).collect())
                .unwrap_or_default();
            let error = json_result.get("error").and_then(|v| v.as_str()).map(|s| s.to_string());
            
            Ok(ExtractRegionResult {
                success,
                total_regions,
                successful_extractions,
                total_variants,
                total_size,
                execution_time,
                results,
                error,
            })
        } else {
            // Legacy parsing fallback
            let lines: Vec<String> = stdout.lines().map(|s| s.to_string()).collect();
            let results: Vec<String> = lines.iter()
                .filter(|line| line.contains("✅") || line.contains("❌"))
                .cloned()
                .collect();
            
            // Extract summary statistics
            let successful_extractions = results.iter().filter(|r| r.contains("✅")).count() as u32;
            let total_regions = results.len() as u32;
            
            // Extract total variants from summary
            let total_variants = stdout
                .lines()
                .find(|line| line.contains("Total variants extracted:"))
                .and_then(|line| line.split(':').last())
                .and_then(|num| num.trim().replace(",", "").parse::<u32>().ok())
                .unwrap_or(0);
            
            // Extract total size
            let total_size = stdout
                .lines()
                .find(|line| line.contains("Total output size:"))
                .and_then(|line| line.split(':').last())
                .and_then(|size_str| {
                    // Extract the KB value
                    size_str.split("KB").next()
                        .and_then(|num| num.trim().parse::<f32>().ok())
                        .map(|kb| (kb * 1024.0) as u32) // Convert to bytes
                })
                .unwrap_or(0);
            
            Ok(ExtractRegionResult {
                success: true,
                total_regions,
                successful_extractions,
                total_variants,
                total_size,
                execution_time,
                results,
                error: None,
            })
        }
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        Err(format!("Extraction failed: {}", stderr))
    }
}

#[command]
async fn check_genomic_dependencies() -> Result<HashMap<String, bool>, String> {
    let tools = vec!["samtools", "bcftools", "bwa", "minimap2", "bowtie2"];
    let mut results = HashMap::new();
    
    for tool in tools {
        let output = Command::new("which")
            .arg(tool)
            .output()
            .map_err(|e| e.to_string())?;
        
        results.insert(tool.to_string(), output.status.success());
    }
    
    Ok(results)
}

#[command]
async fn get_available_vcf_files() -> Result<HashMap<String, String>, String> {
    // Execute the Python config script to get available VCF files
    let args = vec!["--list-vcf-files", "--json"];
    let output = execute_python("config_data_source", &args)
        .map_err(|e| e.to_string())?;
    
    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        serde_json::from_str(&stdout).map_err(|e| e.to_string())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        Err(format!("Failed to get VCF files: {}", stderr))
    }
}

#[command]
async fn generate_test_fastq(output_dir: String, num_reads: u32) -> Result<HashMap<String, String>, String> {
    // Execute the Python script to generate test FASTQ files
    let num_reads_str = num_reads.to_string();
    let args = vec![
        "--output-dir", &output_dir,
        "--num-reads", &num_reads_str,
        "--json"
    ];
    
    let output = execute_python("generate_test_fastq", &args)
        .map_err(|e| e.to_string())?;
    
    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        serde_json::from_str(&stdout).map_err(|e| e.to_string())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        Err(format!("Failed to generate test FASTQ: {}", stderr))
    }
}

// ========================================
// Plugin Management Commands
// ========================================

/// List all available plugins
#[command]
async fn list_plugins() -> Result<Vec<PluginSummary>, String> {
    log::info!("Listing available plugins via Tauri command");
    
    PluginRegistryApi::list_plugins()
        .map_err(|e| {
            log::error!("Failed to list plugins: {}", e);
            e.to_string()
        })
}

/// Execute a plugin with the given arguments
#[command]
async fn run_plugin(plugin_id: String, args_json: String) -> Result<String, String> {
    log::info!("Running plugin '{}' with args: {}", plugin_id, args_json);
    
    // Parse JSON arguments
    let args: serde_json::Value = serde_json::from_str(&args_json)
        .map_err(|e| format!("Invalid JSON arguments: {}", e))?;
    
    // Execute plugin
    let result = PluginRegistryApi::run_plugin(&plugin_id, args)
        .map_err(|e| {
            log::error!("Plugin execution failed for '{}': {}", plugin_id, e);
            e.to_string()
        })?;
    
    // Return result as JSON string
    serde_json::to_string(&result)
        .map_err(|e| format!("Failed to serialize plugin result: {}", e))
}

/// Check if a plugin exists
#[command]
async fn has_plugin(plugin_id: String) -> Result<bool, String> {
    log::debug!("Checking if plugin '{}' exists", plugin_id);
    
    PluginRegistryApi::has_plugin(&plugin_id)
        .map_err(|e| e.to_string())
}

/// Get plugin metadata
#[command]
async fn get_plugin_metadata(plugin_id: String) -> Result<Option<String>, String> {
    log::debug!("Getting metadata for plugin '{}'", plugin_id);
    
    let metadata = PluginRegistryApi::get_plugin_metadata(&plugin_id)
        .map_err(|e| e.to_string())?;
    
    if let Some(metadata) = metadata {
        serde_json::to_string(&metadata)
            .map(Some)
            .map_err(|e| format!("Failed to serialize metadata: {}", e))
    } else {
        Ok(None)
    }
}

/// Reload all plugins
#[command]
async fn reload_plugins() -> Result<(), String> {
    log::info!("Reloading all plugins via Tauri command");
    
    PluginRegistryApi::reload_plugins()
        .map_err(|e| {
            log::error!("Failed to reload plugins: {}", e);
            e.to_string()
        })
}

/// Get plugin registry statistics
#[command]
async fn get_plugin_registry_stats() -> Result<String, String> {
    log::debug!("Getting plugin registry statistics");
    
    let stats = PluginRegistryApi::get_stats()
        .map_err(|e| e.to_string())?;
    
    serde_json::to_string(&stats)
        .map_err(|e| format!("Failed to serialize stats: {}", e))
}

// ========================================
// File Processing Helpers
// ========================================

#[command]
async fn save_temp_file(file_name: String, file_content: Vec<u8>) -> Result<String, String> {
    use std::fs;
    
    // Create a temporary directory for genomic files
    let temp_dir = std::env::temp_dir().join("geneknow_temp");
    fs::create_dir_all(&temp_dir).map_err(|e| e.to_string())?;
    
    // Save the file
    let file_path = temp_dir.join(&file_name);
    fs::write(&file_path, file_content).map_err(|e| e.to_string())?;
    
    log::info!("Saved temporary file: {:?}", file_path);
    Ok(file_path.to_string_lossy().to_string())
}

#[command]
async fn delete_temp_file(file_path: String) -> Result<(), String> {
    use std::fs;
    use std::path::Path;
    
    let path = Path::new(&file_path);
    
    // Security check: only delete files in temp directories
    let temp_markers = ["geneknow_temp", "tmp", "temp", "T/"];
    let path_str = file_path.to_lowercase();
    
    if !temp_markers.iter().any(|marker| path_str.contains(marker)) {
        return Err("Can only delete files in temporary directories".to_string());
    }
    
    if path.exists() {
        fs::remove_file(path).map_err(|e| e.to_string())?;
        log::info!("Deleted temporary file: {}", file_path);
    }
    
    Ok(())
}

// ========================================
// GeneKnow Pipeline API Server Management
// ========================================

#[command]
async fn get_api_port() -> Result<u16, String> {
    let port_guard = API_SERVER_PORT.lock().unwrap();
    match *port_guard {
        Some(port) => Ok(port),
        None => Err("API server not started".to_string()),
    }
}

#[command]
async fn check_api_health() -> Result<bool, String> {
    // Get the current port
    let port = {
        let port_guard = API_SERVER_PORT.lock().unwrap();
        match *port_guard {
            Some(port) => port,
            None => return Ok(false), // Server not started
        }
    };
    
    let client = reqwest::Client::new();
    
    let response = client
        .get(format!("http://localhost:{}/api/health", port))
        .timeout(std::time::Duration::from_secs(2))
        .send()
        .await;
    
    match response {
        Ok(resp) => Ok(resp.status().is_success()),
        Err(_) => Ok(false),
    }
}

#[command]
async fn start_api_server() -> Result<bool, String> {
    // Check if already running
    if check_api_health().await.unwrap_or(false) {
        return Ok(true);
    }

    // Determine if we're running in production (bundled) or development mode
    let (startup_script, working_dir) = if cfg!(debug_assertions) {
        // Development mode: use local Python and script
        log::info!("Running in development mode");
        
        // Get the path to the Python script relative to the workspace root
        let mut api_script = std::env::current_exe()
            .map_err(|e| format!("Failed to get executable path: {}", e))?;
        
        // Go up to find the workspace root
        for _ in 0..5 {
            api_script.pop();
            if api_script.join("geneknow_pipeline/enhanced_api_server.py").exists() {
                break;
            }
        }
        
        let api_script = api_script.join("geneknow_pipeline/enhanced_api_server.py");
        
        // Verify the script exists
        if !api_script.exists() {
            return Err(format!("API server script not found at: {:?}", api_script));
        }
        
        // In development, use system Python
        #[cfg(target_os = "windows")]
        let python_cmd = PathBuf::from("python");
        #[cfg(not(target_os = "windows"))]
        let python_cmd = PathBuf::from("python3");
        
        (python_cmd, api_script.parent().unwrap().to_path_buf())
    } else {
        // Production mode: use bundled Python and resources
        log::info!("Running in production mode");
        
        // Get the bundled resources directory from the app bundle
        let resource_dir = std::env::current_exe()
            .map_err(|e| format!("Failed to get executable path: {}", e))?
            .parent()
            .ok_or("Failed to get executable directory")?
            .parent()
            .ok_or("Failed to get app directory")?
            .join("Resources")
            .join("_up_")
            .join("bundled_resources");
        
        log::info!("Bundled resources directory: {:?}", resource_dir);
        
        // Use the appropriate startup script for the platform
        #[cfg(target_os = "windows")]
        let startup_script = resource_dir.join("start_api_server.bat");
        #[cfg(not(target_os = "windows"))]
        let startup_script = resource_dir.join("start_api_server.sh");
        
        if !startup_script.exists() {
            return Err(format!("Startup script not found at: {:?}", startup_script));
        }
        
        (startup_script.clone(), resource_dir)
    };

    let child = if cfg!(debug_assertions) {
        // Development mode: run Python directly
        log::info!("Starting Python API server in development mode");
        Command::new(&startup_script)
            .arg("enhanced_api_server.py")
            .current_dir(&working_dir)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Failed to start API server: {}", e))?
    } else {
        // Production mode: run the startup script
        log::info!("Starting server with script: {:?}", startup_script);
        log::info!("Working directory: {:?}", working_dir);
        
        // Make sure the script is executable
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mut perms = fs::metadata(&startup_script)
                .map_err(|e| format!("Failed to get script metadata: {}", e))?
                .permissions();
            perms.set_mode(0o755);
            fs::set_permissions(&startup_script, perms)
                .map_err(|e| format!("Failed to set script permissions: {}", e))?;
        }
        
        let mut cmd = Command::new(&startup_script);
        cmd.current_dir(&working_dir)
           .stdout(Stdio::piped())
           .stderr(Stdio::piped());
        
        // Set environment variables to help Python find its modules
        cmd.env("PYTHONPATH", format!("{}:{}/geneknow_pipeline", 
                                     working_dir.display(), 
                                     working_dir.display()));
        
        let child = cmd.spawn()
            .map_err(|e| format!("Failed to start API server: {} (script: {:?})", e, startup_script))?;
        
        child
    };
    
    // Capture stdout to get the port announcement
    let mut captured_port: Option<u16> = None;
    let mut last_error: Option<String> = None;
    
    // Store the process handle first
    {
        let mut process_guard = API_SERVER_PROCESS.lock().unwrap();
        *process_guard = Some(child);
    }
    
    // Handle stdout output in a separate thread
    let stdout = {
        let mut process_guard = API_SERVER_PROCESS.lock().unwrap();
        if let Some(ref mut child) = *process_guard {
            child.stdout.take()
        } else {
            None
        }
    };
    
    if let Some(stdout) = stdout {
        use std::sync::mpsc;
        let (port_tx, port_rx) = mpsc::channel();
        
        thread::spawn(move || {
            use std::io::{BufRead, BufReader};
            let reader = BufReader::new(stdout);
            for line in reader.lines() {
                if let Ok(line) = line {
                    println!("[Python API] {}", line);
                    
                    // Look for port announcement
                    if line.starts_with("API_SERVER_PORT:") {
                        if let Some(port_str) = line.strip_prefix("API_SERVER_PORT:") {
                            if let Ok(port) = port_str.trim().parse::<u16>() {
                                let _ = port_tx.send(port);
                            }
                        }
                    }
                }
            }
        });
        
        // Wait for port announcement with timeout
        let start_time = std::time::Instant::now();
        let timeout = std::time::Duration::from_secs(20);
        
        while start_time.elapsed() < timeout {
            // Check if process is still running
            {
                let mut process_guard = API_SERVER_PROCESS.lock().unwrap();
                if let Some(ref mut child) = *process_guard {
                    match child.try_wait() {
                        Ok(Some(status)) => {
                            last_error = Some(format!("API server process exited with status: {:?}", status));
                            log::error!("Server startup failed: {}", last_error.as_ref().unwrap());
                            break;
                        }
                        Ok(None) => {
                            // Still running, continue
                        }
                        Err(e) => {
                            last_error = Some(format!("Failed to check process status: {}", e));
                            break;
                        }
                    }
                }
            }
            
            // Check for port announcement
            if let Ok(port) = port_rx.try_recv() {
                captured_port = Some(port);
                log::info!("Captured API server port from stdout: {}", port);
                break;
            }
            
            thread::sleep(Duration::from_millis(100));
        }
    }
    
    // Handle stderr output in a separate thread
    let stderr = {
        let mut process_guard = API_SERVER_PROCESS.lock().unwrap();
        if let Some(ref mut child) = *process_guard {
            child.stderr.take()
        } else {
            None
        }
    };
    
    if let Some(stderr) = stderr {
        thread::spawn(move || {
            use std::io::{BufRead, BufReader};
            let reader = BufReader::new(stderr);
            for line in reader.lines() {
                if let Ok(line) = line {
                    eprintln!("[Python API ERROR] {}", line);
                }
            }
        });
    }
    
    // Store the port
    if let Some(port) = captured_port {
        let mut port_guard = API_SERVER_PORT.lock().unwrap();
        *port_guard = Some(port);
    } else {
        // Clean up the process if we failed to get the port
        stop_api_server().await?;
        
        let error_msg = if let Some(err) = last_error {
            format!("Failed to capture API server port: {}", err)
        } else {
            "Failed to capture API server port: timeout".to_string()
        };
        
        return Err(error_msg);
    }

    // Wait for server to start (with timeout)
    let start_time = std::time::Instant::now();
    let timeout = std::time::Duration::from_secs(10);
    
    while start_time.elapsed() < timeout {
        if check_api_health().await.unwrap_or(false) {
            log::info!("GeneKnow API server started successfully");
            return Ok(true);
        }
        tokio::time::sleep(tokio::time::Duration::from_millis(500)).await;
    }

    // If we get here, server didn't start in time
    stop_api_server().await?;
    Err("API server failed to start within timeout".to_string())
}

#[command]
async fn stop_api_server() -> Result<(), String> {
    let mut process_guard = API_SERVER_PROCESS.lock().unwrap();
    
    if let Some(mut child) = process_guard.take() {
        child.kill().map_err(|e| format!("Failed to stop API server: {}", e))?;
        log::info!("GeneKnow API server stopped");
    }
    
    // Clear the port
    let mut port_guard = API_SERVER_PORT.lock().unwrap();
    *port_guard = None;
    
    Ok(())
}

#[command]
async fn get_api_server_status() -> Result<serde_json::Value, String> {
    let is_healthy = check_api_health().await.unwrap_or(false);
    
    if is_healthy {
        // Get the current port
        let port = {
            let port_guard = API_SERVER_PORT.lock().unwrap();
            match *port_guard {
                Some(port) => port,
                None => return Ok(serde_json::json!({
                    "status": "stopped",
                    "service": "GeneKnow Pipeline API"
                })),
            }
        };
        
        // Get detailed health info
        let client = reqwest::Client::new();
        let response = client
            .get(format!("http://localhost:{}/api/health", port))
            .send()
            .await
            .map_err(|e| e.to_string())?;
            
        let health_data = response.json::<serde_json::Value>().await
            .unwrap_or_else(|_| serde_json::json!({
                "status": "running",
                "details": "Unable to get detailed health info"
            }));
            
        Ok(health_data)
    } else {
        Ok(serde_json::json!({
            "status": "stopped",
            "service": "GeneKnow Pipeline API"
        }))
    }
}

#[command]
async fn process_genomic_file(
    file_path: String,
    preferences: Option<serde_json::Value>
) -> Result<String, String> {
    // Ensure API server is running
    if !check_api_health().await.unwrap_or(false) {
        start_api_server().await?;
    }

    // Get the current port
    let port = get_api_port().await?;

    let client = reqwest::Client::new();
    
    let request_body = serde_json::json!({
        "file_path": file_path,
        "preferences": preferences.unwrap_or(serde_json::json!({}))
    });
    
    let response = client
        .post(format!("http://localhost:{}/api/process", port))
        .json(&request_body)
        .send()
        .await
        .map_err(|e| format!("Request failed: {}", e))?;
    
    if response.status().is_success() {
        let job_data: serde_json::Value = response
            .json()
            .await
            .map_err(|e| format!("Failed to parse response: {}", e))?;
        
        Ok(job_data["job_id"].as_str().unwrap_or("").to_string())
    } else {
        let error_data: serde_json::Value = response
            .json()
            .await
            .unwrap_or(serde_json::json!({"error": "Unknown error"}));
        
        Err(error_data["error"].as_str().unwrap_or("Unknown error").to_string())
    }
}

#[command]
async fn get_job_status(job_id: String) -> Result<serde_json::Value, String> {
    // Get the current port
    let port = get_api_port().await?;
    
    let client = reqwest::Client::new();
    
    let response = client
        .get(format!("http://localhost:{}/api/status/{}", port, job_id))
        .send()
        .await
        .map_err(|e| format!("Request failed: {}", e))?;
    
    if response.status().is_success() {
        response
            .json()
            .await
            .map_err(|e| format!("Failed to parse response: {}", e))
    } else {
        Err("Failed to get job status".to_string())
    }
}

#[command]
async fn check_environment() -> Result<bool, String> {
    // Check if bundled Python runtime exists
    if cfg!(debug_assertions) {
        // Development mode: check system Python
        let output = Command::new("python3")
            .arg("--version")
            .output();
        
        match output {
            Ok(output) => Ok(output.status.success()),
            Err(_) => Ok(false),
        }
    } else {
        // Production mode: check bundled Python
        let resource_dir = std::env::current_exe()
            .map_err(|e| format!("Failed to get executable path: {}", e))?
            .parent()
            .ok_or("Failed to get executable directory")?
            .parent()
            .ok_or("Failed to get app directory")?
            .join("Resources")
            .join("_up_")
            .join("bundled_resources");
        
        #[cfg(target_os = "windows")]
        let python_exe = resource_dir.join("python_runtime").join("python.exe");
        #[cfg(not(target_os = "windows"))]
        let python_exe = resource_dir.join("python_runtime").join("bin").join("python3");
        
        Ok(python_exe.exists())
    }
}

#[command]
async fn check_database_exists() -> Result<bool, String> {
    if cfg!(debug_assertions) {
        // Development mode: check in geneknow_pipeline directory
        let mut db_path = std::env::current_exe()
            .map_err(|e| format!("Failed to get executable path: {}", e))?;
        
        for _ in 0..5 {
            db_path.pop();
            if db_path.join("geneknow_pipeline").exists() {
                break;
            }
        }
        
        let db_path = db_path.join("geneknow_pipeline").join("population_variants.db");
        Ok(db_path.exists())
    } else {
        // Production mode: check in bundled resources
        let resource_dir = std::env::current_exe()
            .map_err(|e| format!("Failed to get executable path: {}", e))?
            .parent()
            .ok_or("Failed to get executable directory")?
            .parent()
            .ok_or("Failed to get app directory")?
            .join("Resources")
            .join("_up_")
            .join("bundled_resources");
        
        let db_path = resource_dir.join("geneknow_pipeline").join("population_variants.db");
        Ok(db_path.exists())
    }
}

#[command]
async fn initialize_database() -> Result<bool, String> {
    // In our case, the database should already be created during bundling
    // This is mainly for error recovery
    log::info!("Database initialization requested");
    
    if check_database_exists().await? {
        log::info!("Database already exists");
        return Ok(true);
    }
    
    // If database doesn't exist, this is likely an error in production
    // In development, we could try to create it, but for now just return false
    log::error!("Database not found - this should not happen in production builds");
    Ok(false)
}

#[command]
async fn test_pipeline_connectivity() -> Result<bool, String> {
    // Get the current port
    let port = match get_api_port().await {
        Ok(port) => port,
        Err(_) => return Ok(false), // Server not started
    };
    
    // Test basic connectivity to the pipeline API
    let client = reqwest::Client::new();
    
    let response = client
        .get(format!("http://localhost:{}/api/pipeline-info", port))
        .timeout(std::time::Duration::from_secs(5))
        .send()
        .await;
    
    match response {
        Ok(resp) => Ok(resp.status().is_success()),
        Err(_) => Ok(false),
    }
}

#[command]
async fn get_job_results(job_id: String) -> Result<serde_json::Value, String> {
    // Get the current port
    let port = get_api_port().await?;
    
    let client = reqwest::Client::new();
    
    let response = client
        .get(format!("http://localhost:{}/api/results/{}", port, job_id))
        .send()
        .await
        .map_err(|e| format!("Request failed: {}", e))?;
    
    if response.status().is_success() {
        response
            .json()
            .await
            .map_err(|e| format!("Failed to parse response: {}", e))
    } else {
        Err("Failed to get job results".to_string())
    }
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
  tauri::Builder::default()
    .setup(|app| {
      // Enable logging for both debug and release builds
      app.handle().plugin(
        tauri_plugin_log::Builder::default()
          .level(log::LevelFilter::Info)
          .targets([
            tauri_plugin_log::Target::new(tauri_plugin_log::TargetKind::Stdout),
            tauri_plugin_log::Target::new(tauri_plugin_log::TargetKind::LogDir { 
              file_name: Some("GenePredict.log".to_string()) 
            }),
          ])
          .build(),
      )?;
      
      log::info!("=== GenePredict Starting ===");
      log::info!("Build mode: {}", if cfg!(debug_assertions) { "DEBUG" } else { "RELEASE" });
      log::info!("Executable path: {:?}", std::env::current_exe());
      
      // Initialize plugin registry
      log::info!("Initializing plugin registry during app setup...");
      if let Err(e) = PluginRegistryApi::initialize() {
        log::error!("Failed to initialize plugin registry: {}", e);
        // Don't fail the app startup, just log the error
        log::warn!("Continuing without plugin system initialization");
      } else {
        log::info!("Plugin registry initialized successfully");
      }
      
      // Auto-start the GeneKnow API server in production
      #[cfg(not(debug_assertions))]
      {
        log::info!("Production mode detected - scheduling API server startup");
        let app_handle = app.handle().clone();
        
        // Start the server after a short delay to ensure app is fully initialized
        tauri::async_runtime::spawn(async move {
            // Wait a bit for the app to fully initialize
            log::info!("Waiting 1 second before starting API server...");
            tokio::time::sleep(tokio::time::Duration::from_millis(1000)).await;
            
            log::info!("Now attempting to start GeneKnow API server...");
            let mut attempts = 0;
            const MAX_ATTEMPTS: u32 = 3;
            
            while attempts < MAX_ATTEMPTS {
                log::info!("API server start attempt {} of {}", attempts + 1, MAX_ATTEMPTS);
                match start_api_server().await {
                    Ok(_) => {
                        log::info!("✅ GeneKnow API server started successfully on attempt {}", attempts + 1);
                        
                        // Notify the frontend that the server is ready
                        if let Err(e) = app_handle.emit("api-server-ready", true) {
                            log::warn!("Failed to emit server ready event: {}", e);
                        }
                        break;
                    }
                    Err(e) => {
                        attempts += 1;
                        log::error!("❌ Failed to start API server (attempt {}/{}): {}", attempts, MAX_ATTEMPTS, e);
                        
                        if attempts < MAX_ATTEMPTS {
                            log::info!("Retrying in 2 seconds...");
                            tokio::time::sleep(tokio::time::Duration::from_secs(2)).await;
                        } else {
                            log::error!("Failed to start API server after {} attempts", MAX_ATTEMPTS);
                            // Notify the frontend about the failure
                            if let Err(emit_err) = app_handle.emit("api-server-error", e.clone()) {
                                log::warn!("Failed to emit server error event: {}", emit_err);
                            }
                        }
                    }
                }
            }
            
            log::info!("API server startup sequence completed");
        });
      }
      
      // Also start the server in development mode for easier testing
      #[cfg(debug_assertions)]
      {
        log::info!("Development mode: API server will be started on demand");
        // Optionally auto-start in dev mode too
        let app_handle = app.handle().clone();
        tauri::async_runtime::spawn(async move {
            tokio::time::sleep(tokio::time::Duration::from_millis(2000)).await;
            log::info!("Attempting to start API server in development mode...");
            match start_api_server().await {
                Ok(_) => {
                    log::info!("✅ API server started in development mode");
                    let _ = app_handle.emit("api-server-ready", true);
                }
                Err(e) => {
                    log::error!("❌ Failed to start API server in dev mode: {}", e);
                }
            }
        });
      }
      
      Ok(())
    })
    .invoke_handler(tauri::generate_handler![
        // Legacy commands (backwards compatibility)
        convert_fastq_to_vcf,
        extract_genomic_regions,
        check_genomic_dependencies,
        get_available_vcf_files,
        generate_test_fastq,
        // New plugin management commands
        list_plugins,
        run_plugin,
        has_plugin,
        get_plugin_metadata,
        reload_plugins,
        get_plugin_registry_stats,
        // GeneKnow Pipeline API commands
        check_api_health,
        start_api_server,
        stop_api_server,
        get_api_server_status,
        get_api_port,
        process_genomic_file,
        get_job_status,
        get_job_results,
        // First-run setup commands
        check_environment,
        check_database_exists,
        initialize_database,
        test_pipeline_connectivity,
        // File processing helpers
        save_temp_file,
        delete_temp_file
    ])
    .on_window_event(|_window, event| {
        // Stop API server when app closes
        if let tauri::WindowEvent::Destroyed = event {
            let mut process_guard = API_SERVER_PROCESS.lock().unwrap();
            if let Some(mut child) = process_guard.take() {
                let _ = child.kill();
                log::info!("Stopped GeneKnow API server on app exit");
            }
        }
    })
    .run(tauri::generate_context!())
    .expect("error while running tauri application");
}
