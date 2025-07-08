use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use thiserror::Error;

/// Core trait that all genomic plugins must implement
/// 
/// This trait defines the contract for all genomic analysis plugins,
/// providing a standardized interface for plugin discovery, execution,
/// and metadata management.
pub trait GenomicPlugin: Send + Sync {
    /// Unique identifier for this plugin (e.g., "fastq_to_vcf")
    fn id(&self) -> &'static str;
    
    /// Human-readable name for this plugin
    fn name(&self) -> &'static str;
    
    /// Short description of what this plugin does
    fn description(&self) -> &'static str;
    
    /// Plugin version
    fn version(&self) -> &'static str;
    
    /// Execute the plugin with arbitrary arguments encoded as JSON
    /// 
    /// # Arguments
    /// * `args` - JSON value containing plugin-specific arguments
    /// 
    /// # Returns
    /// Returns the plugin execution result as JSON or an error
    fn run(&self, args: serde_json::Value) -> Result<serde_json::Value, PluginError>;
    
    /// Get plugin manifest information
    fn manifest(&self) -> &PluginManifest;
    
    /// Validate plugin arguments against schema
    fn validate_args(&self, args: &serde_json::Value) -> Result<(), PluginError> {
        // Default implementation - can be overridden by specific plugins
        log::debug!("Validating arguments for plugin {}: {:?}", self.id(), args);
        Ok(())
    }
}

/// Plugin manifest structure parsed from manifest.json
/// 
/// Contains metadata about the plugin including its configuration,
/// input/output schemas, and execution parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PluginManifest {
    /// Plugin unique identifier
    pub id: String,
    
    /// Plugin human-readable name
    pub name: String,
    
    /// Plugin description
    pub description: String,
    
    /// Plugin version
    pub version: String,
    
    /// Plugin author information
    pub author: Option<String>,
    
    /// Plugin script file path (relative to plugin directory)
    pub script_path: String,
    
    /// Plugin category (e.g., "genomic_processing", "data_analysis")
    pub category: String,
    
    /// Plugin configuration parameters
    pub config: HashMap<String, serde_json::Value>,
    
    /// Input argument schema (JSON Schema)
    pub input_schema: Option<serde_json::Value>,
    
    /// Output result schema (JSON Schema)
    pub output_schema: Option<serde_json::Value>,
    
    /// Plugin requirements (Python packages, system tools)
    pub requirements: Vec<String>,
    
    /// Plugin tags for organization
    pub tags: Vec<String>,
    
    /// Whether plugin is enabled
    pub enabled: bool,
    
    /// Cross-platform compatibility flags
    pub platform_support: PlatformSupport,
}

/// Platform support configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlatformSupport {
    pub windows: bool,
    pub macos: bool,
    pub linux: bool,
}

impl Default for PlatformSupport {
    fn default() -> Self {
        Self {
            windows: true,
            macos: true,
            linux: true,
        }
    }
}

/// Summary information about a plugin for frontend display
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PluginSummary {
    pub id: String,
    pub name: String,
    pub description: String,
    pub version: String,
    pub category: String,
    pub tags: Vec<String>,
    pub enabled: bool,
    pub platform_compatible: bool,
}

/// Plugin execution result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PluginExecutionResult {
    pub success: bool,
    pub plugin_id: String,
    pub execution_time: f64,
    pub result: Option<serde_json::Value>,
    pub error: Option<String>,
    pub logs: Vec<String>,
}

/// Plugin error types
#[derive(Error, Debug)]
pub enum PluginError {
    #[error("Plugin not found: {0}")]
    PluginNotFound(String),
    
    #[error("Plugin execution failed: {0}")]
    ExecutionFailed(String),
    
    #[error("Plugin manifest invalid: {0}")]
    InvalidManifest(String),
    
    #[error("Plugin argument validation failed: {0}")]
    ValidationFailed(String),
    
    #[error("Plugin initialization failed: {0}")]
    InitializationFailed(String),
    
    #[error("Plugin platform incompatible: {0}")]
    PlatformIncompatible(String),
    
    #[error("Plugin disabled: {0}")]
    PluginDisabled(String),
    
    #[error("JSON parsing error: {0}")]
    JsonError(#[from] serde_json::Error),
    
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

impl PluginError {
    /// Convert plugin error to string for Tauri command responses
    pub fn to_string(&self) -> String {
        format!("{}", self)
    }
}

/// Plugin configuration for runtime parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PluginConfig {
    /// Plugin-specific configuration values
    pub values: HashMap<String, serde_json::Value>,
    
    /// Remote configuration overrides (from Firebase Remote Config)
    pub remote_overrides: Option<HashMap<String, serde_json::Value>>,
    
    /// Whether to use remote configuration
    pub use_remote_config: bool,
}

impl PluginConfig {
    /// Get a configuration value with fallback to default
    pub fn get_value(&self, key: &str, default: serde_json::Value) -> serde_json::Value {
        // Check remote overrides first if enabled
        if self.use_remote_config {
            if let Some(remote) = &self.remote_overrides {
                if let Some(value) = remote.get(key) {
                    return value.clone();
                }
            }
        }
        
        // Fallback to local configuration
        self.values.get(key).cloned().unwrap_or(default)
    }
    
    /// Merge remote configuration with local configuration
    pub fn merge_remote_config(&mut self, remote_config: HashMap<String, serde_json::Value>) {
        self.remote_overrides = Some(remote_config);
    }
}

/// Plugin metadata for discovery and management
#[derive(Debug, Clone)]
pub struct PluginMetadata {
    pub manifest: PluginManifest,
    pub config: PluginConfig,
    pub directory_path: std::path::PathBuf,
    pub last_loaded: std::time::SystemTime,
}

impl PluginMetadata {
    /// Check if plugin is compatible with current platform
    pub fn is_platform_compatible(&self) -> bool {
        let platform = &self.manifest.platform_support;
        
        #[cfg(target_os = "windows")]
        return platform.windows;
        
        #[cfg(target_os = "macos")]
        return platform.macos;
        
        #[cfg(target_os = "linux")]
        return platform.linux;
        
        #[cfg(not(any(target_os = "windows", target_os = "macos", target_os = "linux")))]
        return false;
    }
    
    /// Convert to plugin summary for frontend
    pub fn to_summary(&self) -> PluginSummary {
        PluginSummary {
            id: self.manifest.id.clone(),
            name: self.manifest.name.clone(),
            description: self.manifest.description.clone(),
            version: self.manifest.version.clone(),
            category: self.manifest.category.clone(),
            tags: self.manifest.tags.clone(),
            enabled: self.manifest.enabled,
            platform_compatible: self.is_platform_compatible(),
        }
    }
}

/// Plugin loading and validation utilities
pub struct PluginLoader;

impl PluginLoader {
    /// Load plugin manifest from JSON file
    pub fn load_manifest(manifest_path: &std::path::Path) -> Result<PluginManifest, PluginError> {
        log::info!("Loading plugin manifest from: {:?}", manifest_path);
        
        let manifest_content = std::fs::read_to_string(manifest_path)?;
        let manifest: PluginManifest = serde_json::from_str(&manifest_content)?;
        
        // Validate manifest
        Self::validate_manifest(&manifest)?;
        
        log::info!("Successfully loaded plugin manifest: {} ({})", manifest.name, manifest.id);
        Ok(manifest)
    }
    
    /// Validate plugin manifest structure and required fields
    fn validate_manifest(manifest: &PluginManifest) -> Result<(), PluginError> {
        if manifest.id.is_empty() {
            return Err(PluginError::InvalidManifest("Plugin ID cannot be empty".to_string()));
        }
        
        if manifest.name.is_empty() {
            return Err(PluginError::InvalidManifest("Plugin name cannot be empty".to_string()));
        }
        
        if manifest.script_path.is_empty() {
            return Err(PluginError::InvalidManifest("Plugin script path cannot be empty".to_string()));
        }
        
        // Validate ID format (alphanumeric and underscores only)
        if !manifest.id.chars().all(|c| c.is_alphanumeric() || c == '_') {
            return Err(PluginError::InvalidManifest(
                "Plugin ID must contain only alphanumeric characters and underscores".to_string()
            ));
        }
        
        log::debug!("Plugin manifest validation passed for: {}", manifest.id);
        Ok(())
    }
    
    /// Create plugin configuration from manifest
    pub fn create_config(manifest: &PluginManifest) -> PluginConfig {
        PluginConfig {
            values: manifest.config.clone(),
            remote_overrides: None,
            use_remote_config: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    
    fn create_test_manifest() -> PluginManifest {
        let mut config = HashMap::new();
        config.insert("threads".to_string(), serde_json::json!(4));
        config.insert("memory_limit".to_string(), serde_json::json!("8GB"));
        
        PluginManifest {
            id: "test_plugin".to_string(),
            name: "Test Plugin".to_string(),
            description: "A test plugin for unit testing".to_string(),
            version: "1.0.0".to_string(),
            author: Some("Test Author".to_string()),
            script_path: "test_plugin.py".to_string(),
            category: "testing".to_string(),
            config,
            input_schema: None,
            output_schema: None,
            requirements: vec!["python>=3.8".to_string()],
            tags: vec!["test".to_string(), "validation".to_string()],
            enabled: true,
            platform_support: PlatformSupport::default(),
        }
    }
    
    #[test]
    fn test_plugin_manifest_validation() {
        let manifest = create_test_manifest();
        assert!(PluginLoader::validate_manifest(&manifest).is_ok());
    }
    
    #[test]
    fn test_plugin_config_value_retrieval() {
        let manifest = create_test_manifest();
        let config = PluginLoader::create_config(&manifest);
        
        let threads = config.get_value("threads", serde_json::json!(1));
        assert_eq!(threads, serde_json::json!(4));
        
        let nonexistent = config.get_value("nonexistent", serde_json::json!("default"));
        assert_eq!(nonexistent, serde_json::json!("default"));
    }
    
    #[test]
    fn test_plugin_metadata_platform_compatibility() {
        let manifest = create_test_manifest();
        let config = PluginLoader::create_config(&manifest);
        
        let metadata = PluginMetadata {
            manifest,
            config,
            directory_path: std::path::PathBuf::from("/test/path"),
            last_loaded: std::time::SystemTime::now(),
        };
        
        // Should be compatible on all platforms by default
        assert!(metadata.is_platform_compatible());
    }
    
    #[test]
    fn test_plugin_summary_conversion() {
        let manifest = create_test_manifest();
        let config = PluginLoader::create_config(&manifest);
        
        let metadata = PluginMetadata {
            manifest,
            config,
            directory_path: std::path::PathBuf::from("/test/path"),
            last_loaded: std::time::SystemTime::now(),
        };
        
        let summary = metadata.to_summary();
        assert_eq!(summary.id, "test_plugin");
        assert_eq!(summary.name, "Test Plugin");
        assert!(summary.enabled);
        assert!(summary.platform_compatible);
    }
} 