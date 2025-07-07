use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use crate::error::{GenePredicateError, Result};
use tracing::{info, debug, warn};

/// Application configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AppConfig {
    pub app: AppSettings,
    pub processing: ProcessingSettings,
    pub security: SecuritySettings,
    pub paths: PathSettings,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AppSettings {
    pub name: String,
    pub version: String,
    pub debug_mode: bool,
    pub max_file_size_mb: u64,
    pub supported_formats: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingSettings {
    pub max_concurrent_jobs: u32,
    pub timeout_seconds: u64,
    pub enable_gpu: bool,
    pub temp_dir_cleanup: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SecuritySettings {
    pub local_only: bool,
    pub encrypt_temp_files: bool,
    pub audit_logging: bool,
    pub data_retention_hours: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathSettings {
    pub temp_dir: PathBuf,
    pub data_dir: PathBuf,
    pub plugins_dir: PathBuf,
    pub reference_data_dir: PathBuf,
}

impl Default for AppConfig {
    fn default() -> Self {
        let home_dir = dirs::home_dir().unwrap_or_else(|| PathBuf::from("."));
        let app_dir = home_dir.join(".genepredict");
        
        AppConfig {
            app: AppSettings {
                name: "GenePredict".to_string(),
                version: "0.1.0".to_string(),
                debug_mode: cfg!(debug_assertions),
                max_file_size_mb: 500, // 500MB max file size
                supported_formats: vec![
                    "vcf".to_string(),
                    "bam".to_string(),
                    "fastq".to_string(),
                    "fq".to_string(),
                ],
            },
            processing: ProcessingSettings {
                max_concurrent_jobs: 4,
                timeout_seconds: 300, // 5 minutes
                enable_gpu: false, // Disable GPU by default for compatibility
                temp_dir_cleanup: true,
            },
            security: SecuritySettings {
                local_only: true,
                encrypt_temp_files: true,
                audit_logging: true,
                data_retention_hours: 24, // Delete processed data after 24 hours
            },
            paths: PathSettings {
                temp_dir: app_dir.join("temp"),
                data_dir: app_dir.join("data"),
                plugins_dir: app_dir.join("plugins"),
                reference_data_dir: app_dir.join("reference"),
            },
        }
    }
}

impl AppConfig {
    /// Load configuration from file or create default
    pub fn load() -> Result<Self> {
        info!("📋 Loading application configuration...");
        
        let home_dir = dirs::home_dir().unwrap_or_else(|| PathBuf::from("."));
        let config_path = home_dir.join(".genepredict").join("config.toml");
        
        if config_path.exists() {
            info!("📁 Loading config from: {:?}", config_path);
            let config_str = std::fs::read_to_string(&config_path)
                .map_err(|e| GenePredicateError::Config(format!("Failed to read config file: {}", e)))?;
            
            let config: AppConfig = toml::from_str(&config_str)
                .map_err(|e| GenePredicateError::Config(format!("Failed to parse config: {}", e)))?;
            
            debug!("✅ Configuration loaded from file");
            Ok(config)
        } else {
            warn!("⚠️  Config file not found, using default configuration");
            let config = AppConfig::default();
            
            // Create config directory
            if let Some(parent) = config_path.parent() {
                std::fs::create_dir_all(parent)
                    .map_err(|e| GenePredicateError::Config(format!("Failed to create config directory: {}", e)))?;
            }
            
            // Save default config
            config.save(&config_path)?;
            
            info!("📄 Default configuration created at: {:?}", config_path);
            Ok(config)
        }
    }
    
    /// Save configuration to file
    pub fn save(&self, path: &PathBuf) -> Result<()> {
        let config_str = toml::to_string_pretty(self)
            .map_err(|e| GenePredicateError::Config(format!("Failed to serialize config: {}", e)))?;
        
        std::fs::write(path, config_str)
            .map_err(|e| GenePredicateError::Config(format!("Failed to write config file: {}", e)))?;
        
        debug!("💾 Configuration saved to: {:?}", path);
        Ok(())
    }
    
    /// Ensure all required directories exist
    pub fn ensure_directories(&self) -> Result<()> {
        let dirs = [
            &self.paths.temp_dir,
            &self.paths.data_dir,
            &self.paths.plugins_dir,
            &self.paths.reference_data_dir,
        ];
        
        for dir in &dirs {
            if !dir.exists() {
                info!("📁 Creating directory: {:?}", dir);
                std::fs::create_dir_all(dir)
                    .map_err(|e| GenePredicateError::Config(format!("Failed to create directory {:?}: {}", dir, e)))?;
            } else {
                debug!("✅ Directory exists: {:?}", dir);
            }
        }
        
        Ok(())
    }
}

// Add chrono and dirs dependencies to Cargo.toml if not already there
#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;
    
    #[test]
    fn test_default_config() {
        let config = AppConfig::default();
        assert_eq!(config.app.name, "GenePredict");
        assert!(config.security.local_only);
        assert!(config.processing.temp_dir_cleanup);
    }
    
    #[test]
    fn test_config_save_load() {
        let temp_dir = TempDir::new().unwrap();
        let config_path = temp_dir.path().join("test_config.toml");
        
        let original_config = AppConfig::default();
        original_config.save(&config_path).unwrap();
        
        assert!(config_path.exists());
        
        // Note: We can't easily test loading without implementing a custom load method
        // that accepts a path parameter
    }
} 