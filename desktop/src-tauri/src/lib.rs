use std::collections::HashMap;
use std::process::Command;
use serde::{Deserialize, Serialize};
use tauri::command;

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
    // Build the Python command
    let mut cmd = Command::new("python3");
    cmd.arg("src/api/fastq_to_vcf_pipeline.py");
    cmd.arg("-r").arg(&options.reference);
    cmd.arg("-1").arg(&options.fastq1);
    
    if let Some(fastq2) = &options.fastq2 {
        cmd.arg("-2").arg(fastq2);
    }
    
    cmd.arg("-o").arg(&options.output_prefix);
    
    if let Some(threads) = options.threads {
        cmd.arg("-t").arg(threads.to_string());
    }
    
    if let Some(aligner) = &options.aligner {
        cmd.arg("-a").arg(aligner);
    }
    
    // Execute the command
    let start_time = std::time::Instant::now();
    let output = cmd.output().map_err(|e| e.to_string())?;
    let execution_time = start_time.elapsed().as_secs_f64();
    
    if output.status.success() {
        // Parse the output to extract results
        let stdout = String::from_utf8_lossy(&output.stdout);
        
        // Extract variant count from output
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
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        Err(format!("Pipeline failed: {}", stderr))
    }
}

#[command]
async fn extract_genomic_regions(options: ExtractRegionOptions) -> Result<ExtractRegionResult, String> {
    // First, update the Python script to use the provided BED file
    let mut cmd = Command::new("python3");
    cmd.arg("src/api/extract_by_region.py");
    
    // Set environment variables for the script
    cmd.env("BED_PATH", &options.bed_file);
    cmd.env("OUTPUT_DIR", &options.output_dir);
    
    if let Some(max_processes) = options.max_processes {
        cmd.env("MAX_PROCESSES", max_processes.to_string());
    }
    
    // Execute the command
    let start_time = std::time::Instant::now();
    let output = cmd.output().map_err(|e| e.to_string())?;
    let execution_time = start_time.elapsed().as_secs_f64();
    
    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        
        // Parse the output to extract results
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
    // This would call the Python config module to get available VCF files
    let output = Command::new("python3")
        .arg("-c")
        .arg("import sys; sys.path.append('src/api'); from config_data_source import get_vcf_files; import json; print(json.dumps(get_vcf_files()))")
        .output()
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
    let output = Command::new("python3")
        .arg("-c")
        .arg(format!(r#"
import random
import gzip
import os
import json

output_dir = '{}'
num_reads = {}
os.makedirs(output_dir, exist_ok=True)

def generate_fastq(filename, num_reads=1000, read_length=150):
    bases = ['A', 'T', 'G', 'C']
    filepath = os.path.join(output_dir, filename)
    with gzip.open(filepath, 'wt') as f:
        for i in range(num_reads):
            seq = ''.join(random.choice(bases) for _ in range(read_length))
            qual = 'I' * read_length
            f.write(f'@read_{{i}}\n{{seq}}\n+\n{{qual}}\n')
    return filepath

file1 = generate_fastq('test_R1.fastq.gz', num_reads)
file2 = generate_fastq('test_R2.fastq.gz', num_reads)

result = {{
    'success': True,
    'file1': file1,
    'file2': file2
}}
print(json.dumps(result))
"#, output_dir, num_reads))
        .output()
        .map_err(|e| e.to_string())?;
    
    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        serde_json::from_str(&stdout).map_err(|e| e.to_string())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        Err(format!("Failed to generate test FASTQ: {}", stderr))
    }
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
      Ok(())
    })
    .invoke_handler(tauri::generate_handler![
        convert_fastq_to_vcf,
        extract_genomic_regions,
        check_genomic_dependencies,
        get_available_vcf_files,
        generate_test_fastq
    ])
    .run(tauri::generate_context!())
    .expect("error while running tauri application");
}
