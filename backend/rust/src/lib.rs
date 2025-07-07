// Learn more about Tauri commands at https://tauri.app/develop/calling-rust/
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use tokio::sync::RwLock;
use serde::{Deserialize, Serialize};
use tauri::{Manager, State};
use tracing::{info, warn, error, debug};
use uuid::Uuid;

pub mod plugins;
pub mod file_processing;
pub mod ml_engine;
pub mod config;
pub mod error;

pub use error::{GenePredicateError, Result};

/// Represents a genomic file uploaded by the user
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomicFile {
    pub id: String,
    pub name: String,
    pub file_type: String,
    pub size: u64,
    pub path: PathBuf,
    pub status: ProcessingStatus,
    pub created_at: chrono::DateTime<chrono::Utc>,
}

/// Processing status of a genomic file
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ProcessingStatus {
    Uploaded,
    Validating,
    Processing,
    Completed,
    Failed(String),
}

/// Risk assessment result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RiskAssessment {
    pub file_id: String,
    pub overall_risk: String,
    pub risk_score: f64,
    pub confidence: f64,
    pub variants_count: u32,
    pub genes: Vec<String>,
    pub created_at: chrono::DateTime<chrono::Utc>,
}

/// Application state
#[derive(Debug)]
pub struct AppState {
    pub files: Arc<RwLock<HashMap<String, GenomicFile>>>,
    pub assessments: Arc<RwLock<HashMap<String, RiskAssessment>>>,
    pub config: config::AppConfig,
    pub plugins: plugins::PluginManager,
}

impl AppState {
    pub fn new() -> Result<Self> {
        let config = config::AppConfig::load()?;
        let plugins = plugins::PluginManager::new();
        
        Ok(AppState {
            files: Arc::new(RwLock::new(HashMap::new())),
            assessments: Arc::new(RwLock::new(HashMap::new())),
            config,
            plugins,
        })
    }
}

/// Upload a genomic file
#[tauri::command]
async fn upload_genomic_file(
    file_path: String,
    file_name: String,
    file_size: u64,
    state: State<'_, AppState>,
) -> Result<GenomicFile, String> {
    info!("📁 Uploading genomic file: {} ({}KB)", file_name, file_size / 1024);
    
    let file_id = Uuid::new_v4().to_string();
    let file_type = determine_file_type(&file_name);
    
    let genomic_file = GenomicFile {
        id: file_id.clone(),
        name: file_name.clone(),
        file_type: file_type.clone(),
        size: file_size,
        path: PathBuf::from(file_path),
        status: ProcessingStatus::Uploaded,
        created_at: chrono::Utc::now(),
    };
    
    // Store the file in state
    state.files.write().await.insert(file_id.clone(), genomic_file.clone());
    
    info!("✅ File uploaded successfully: {}", file_name);
    debug!("📊 File details: {:?}", genomic_file);
    
    Ok(genomic_file)
}

/// Process a genomic file for risk assessment
#[tauri::command]
async fn process_genomic_file(
    file_id: String,
    state: State<'_, AppState>,
) -> Result<RiskAssessment, String> {
    info!("🔄 Processing genomic file: {}", file_id);
    
    // Get file from state
    let file = {
        let files = state.files.read().await;
        files.get(&file_id).cloned()
    };
    
    let file = match file {
        Some(file) => file,
        None => {
            error!("❌ File not found: {}", file_id);
            return Err(format!("File not found: {}", file_id));
        }
    };
    
    // Update file status to processing
    {
        let mut files = state.files.write().await;
        if let Some(file) = files.get_mut(&file_id) {
            file.status = ProcessingStatus::Processing;
        }
    }
    
    info!("🔬 Starting file processing pipeline for: {}", file.name);
    
    // Process the file using the plugin system
    let result = match file.file_type.as_str() {
        "vcf" => process_vcf_file(&file, &state).await,
        "bam" => process_bam_file(&file, &state).await,
        "fastq" => process_fastq_file(&file, &state).await,
        _ => {
            error!("❌ Unsupported file type: {}", file.file_type);
            Err("Unsupported file type".to_string())
        }
    };
    
    match result {
        Ok(assessment) => {
            // Update file status to completed
            {
                let mut files = state.files.write().await;
                if let Some(file) = files.get_mut(&file_id) {
                    file.status = ProcessingStatus::Completed;
                }
            }
            
            // Store assessment
            state.assessments.write().await.insert(file_id.clone(), assessment.clone());
            
            info!("✅ File processed successfully: {}", file.name);
            debug!("🎯 Assessment results: {:?}", assessment);
            
            Ok(assessment)
        }
        Err(e) => {
            error!("❌ Processing failed for file {}: {}", file.name, e);
            
            // Update file status to failed
            {
                let mut files = state.files.write().await;
                if let Some(file) = files.get_mut(&file_id) {
                    file.status = ProcessingStatus::Failed(e.clone());
                }
            }
            
            Err(e)
        }
    }
}

/// Get application status and info
#[tauri::command]
async fn get_app_info(state: State<'_, AppState>) -> Result<serde_json::Value, String> {
    let files_count = state.files.read().await.len();
    let assessments_count = state.assessments.read().await.len();
    
    let info = serde_json::json!({
        "app_name": "GenePredict",
        "version": "0.1.0",
        "description": "AI for Genomic Risk Assessment",
        "files_processed": files_count,
        "assessments_generated": assessments_count,
        "privacy_mode": "local_only",
        "status": "ready"
    });
    
    debug!("📊 App info requested: {:?}", info);
    Ok(info)
}

/// List all uploaded files
#[tauri::command]
async fn list_files(state: State<'_, AppState>) -> Result<Vec<GenomicFile>, String> {
    let files = state.files.read().await;
    let files_vec: Vec<GenomicFile> = files.values().cloned().collect();
    
    debug!("📋 Files list requested: {} files", files_vec.len());
    Ok(files_vec)
}

/// Get risk assessment for a file
#[tauri::command]
async fn get_risk_assessment(
    file_id: String,
    state: State<'_, AppState>,
) -> Result<RiskAssessment, String> {
    let assessments = state.assessments.read().await;
    
    match assessments.get(&file_id) {
        Some(assessment) => {
            debug!("🎯 Risk assessment retrieved for file: {}", file_id);
            Ok(assessment.clone())
        }
        None => {
            warn!("⚠️  No risk assessment found for file: {}", file_id);
            Err(format!("No risk assessment found for file: {}", file_id))
        }
    }
}

// Helper functions
fn determine_file_type(file_name: &str) -> String {
    let path = PathBuf::from(file_name);
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("vcf") => "vcf".to_string(),
        Some("bam") => "bam".to_string(),
        Some("fastq") | Some("fq") => "fastq".to_string(),
        _ => "unknown".to_string(),
    }
}

async fn process_vcf_file(file: &GenomicFile, state: &AppState) -> std::result::Result<RiskAssessment, String> {
    info!("🧬 Processing VCF file: {}", file.name);
    
    // Simulate VCF processing
    tokio::time::sleep(std::time::Duration::from_secs(2)).await;
    
    // Mock assessment for VCF
    let assessment = RiskAssessment {
        file_id: file.id.clone(),
        overall_risk: "medium".to_string(),
        risk_score: 6.3,
        confidence: 87.0,
        variants_count: 1247,
        genes: vec!["BRCA1".to_string(), "BRCA2".to_string(), "TP53".to_string()],
        created_at: chrono::Utc::now(),
    };
    
    info!("✅ VCF processing completed for: {}", file.name);
    Ok(assessment)
}

async fn process_bam_file(file: &GenomicFile, state: &AppState) -> std::result::Result<RiskAssessment, String> {
    info!("🧬 Processing BAM file: {}", file.name);
    
    // Simulate BAM processing
    tokio::time::sleep(std::time::Duration::from_secs(3)).await;
    
    // Mock assessment for BAM
    let assessment = RiskAssessment {
        file_id: file.id.clone(),
        overall_risk: "high".to_string(),
        risk_score: 7.8,
        confidence: 92.0,
        variants_count: 2156,
        genes: vec!["BRCA1".to_string(), "BRCA2".to_string(), "TP53".to_string(), "PALB2".to_string()],
        created_at: chrono::Utc::now(),
    };
    
    info!("✅ BAM processing completed for: {}", file.name);
    Ok(assessment)
}

async fn process_fastq_file(file: &GenomicFile, state: &AppState) -> std::result::Result<RiskAssessment, String> {
    info!("🧬 Processing FASTQ file: {}", file.name);
    
    // Simulate FASTQ processing
    tokio::time::sleep(std::time::Duration::from_secs(4)).await;
    
    // Mock assessment for FASTQ
    let assessment = RiskAssessment {
        file_id: file.id.clone(),
        overall_risk: "low".to_string(),
        risk_score: 3.2,
        confidence: 78.0,
        variants_count: 543,
        genes: vec!["ATM".to_string(), "CHEK2".to_string()],
        created_at: chrono::Utc::now(),
    };
    
    info!("✅ FASTQ processing completed for: {}", file.name);
    Ok(assessment)
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    info!("🚀 Initializing GenePredict application...");
    
    // Initialize application state
    let app_state = AppState::new().expect("Failed to initialize application state");
    
    info!("📋 Registered Tauri commands:");
    info!("  - upload_genomic_file");
    info!("  - process_genomic_file");
    info!("  - get_app_info");
    info!("  - list_files");
    info!("  - get_risk_assessment");
    
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .manage(app_state)
        .invoke_handler(tauri::generate_handler![
            upload_genomic_file,
            process_genomic_file,
            get_app_info,
            list_files,
            get_risk_assessment
        ])
        .setup(|app| {
            info!("🎯 GenePredict application setup completed");
            info!("🔒 Privacy mode: Local processing only");
            info!("🧬 Ready for genomic risk assessment");
            Ok(())
        })
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
