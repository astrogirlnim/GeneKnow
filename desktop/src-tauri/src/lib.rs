use std::collections::HashMap;
use std::process::Command;
use serde::{Deserialize, Serialize};
use tauri::command;

// Import modules
pub mod utils;
pub mod plugin;
pub mod python_script_plugin;
pub mod plugin_registry;

use utils::execute_python;
use plugin::{PluginSummary, PluginError};
use plugin_registry::PluginRegistryApi;

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

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
  tauri::Builder::default()
    .setup(|app| {
      if cfg!(debug_assertions) {
        app.handle().plugin(
          tauri_plugin_log::Builder::default()
            .level(log::LevelFilter::Info)
            .build(),
        )?;
      }
      
      // Initialize plugin registry
      log::info!("Initializing plugin registry during app setup...");
      if let Err(e) = PluginRegistryApi::initialize() {
        log::error!("Failed to initialize plugin registry: {}", e);
        // Don't fail the app startup, just log the error
        log::warn!("Continuing without plugin system initialization");
      } else {
        log::info!("Plugin registry initialized successfully");
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
        get_plugin_registry_stats
    ])
    .run(tauri::generate_context!())
    .expect("error while running tauri application");
}
