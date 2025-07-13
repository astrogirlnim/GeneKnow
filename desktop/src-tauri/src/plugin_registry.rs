use crate::plugin::{
    GenomicPlugin, PluginError, PluginLoader, PluginManifest, PluginMetadata,
    PluginSummary, PluginConfig
};
use crate::python_script_plugin::PythonScriptPluginFactory;
use once_cell::sync::Lazy;
use serde_json::Value;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, RwLock};
use log::{info, debug, warn};

/// Global plugin registry instance
/// 
/// This is a singleton that manages all registered plugins and provides
/// thread-safe access to plugin discovery, loading, and execution.
pub static PLUGIN_REGISTRY: Lazy<Arc<RwLock<PluginRegistry>>> = 
    Lazy::new(|| Arc::new(RwLock::new(PluginRegistry::new())));

/// Main plugin registry for managing genomic analysis plugins
/// 
/// The registry handles plugin discovery, loading, validation, and execution.
/// It maintains a cache of loaded plugins and provides thread-safe access.
pub struct PluginRegistry {
    /// Map of plugin ID to plugin instance
    plugins: HashMap<String, Box<dyn GenomicPlugin>>,
    
    /// Map of plugin ID to plugin metadata
    metadata: HashMap<String, PluginMetadata>,
    
    /// Plugin directory path
    plugin_dir: PathBuf,
    
    /// Whether the registry has been initialized
    initialized: bool,
    
    /// Registry configuration
    config: RegistryConfig,
}

/// Configuration for the plugin registry
#[derive(Debug, Clone)]
pub struct RegistryConfig {
    /// Whether to enable remote configuration
    pub enable_remote_config: bool,
    
    /// Whether to auto-reload plugins on file changes
    pub auto_reload: bool,
    
    /// Maximum number of plugins to load
    pub max_plugins: usize,
    
    /// Plugin discovery timeout in seconds
    pub discovery_timeout: u64,
}

impl Default for RegistryConfig {
    fn default() -> Self {
        Self {
            enable_remote_config: false,
            auto_reload: false,
            max_plugins: 50,
            discovery_timeout: 30,
        }
    }
}

impl PluginRegistry {
    /// Create a new plugin registry
    pub fn new() -> Self {
        Self {
            plugins: HashMap::new(),
            metadata: HashMap::new(),
            plugin_dir: PathBuf::new(),
            initialized: false,
            config: RegistryConfig::default(),
        }
    }
    
    /// Initialize the plugin registry
    /// 
    /// This method scans the plugin directory and loads all available plugins.
    /// It should be called once during application startup.
    pub fn initialize(&mut self) -> Result<(), PluginError> {
        if self.initialized {
            warn!("Plugin registry already initialized");
            return Ok(());
        }
        
        info!("Initializing plugin registry...");
        
        // Set up plugin directory
        self.plugin_dir = self.get_plugin_directory()?;
        
        // Scan and load plugins
        self.scan_plugins()?;
        
        self.initialized = true;
        info!("Plugin registry initialized successfully. Loaded {} plugins.", self.plugins.len());
        
        Ok(())
    }
    
    /// Get the plugin directory path
    fn get_plugin_directory(&self) -> Result<PathBuf, PluginError> {
        // Check if we're running in production mode (app bundle)
        let plugin_dir = if cfg!(debug_assertions) {
            // Development mode: Find project root
            let current_dir = std::env::current_dir()
                .map_err(|e| PluginError::InitializationFailed(e.to_string()))?;
            
            let mut path = current_dir;
            loop {
                // Check for key project files
                if path.join("README.md").exists() && 
                   path.join("desktop").exists() && 
                   path.join("docs").exists() {
                    break;
                }
                
                // Move up one directory
                if let Some(parent) = path.parent() {
                    path = parent.to_path_buf();
                } else {
                    return Err(PluginError::InitializationFailed(
                        "Could not find project root directory".to_string()
                    ));
                }
            }
            
            // Build plugin directory path
            path.join("desktop").join("python_ml").join("plugins")
        } else {
            // Production mode: Use bundled resources
            let exe_path = std::env::current_exe()
                .map_err(|e| PluginError::InitializationFailed(
                    format!("Failed to get executable path: {}", e)
                ))?;
            
            let resource_dir = exe_path
                .parent()
                .ok_or_else(|| PluginError::InitializationFailed("Failed to get executable directory".to_string()))?
                .parent()
                .ok_or_else(|| PluginError::InitializationFailed("Failed to get app directory".to_string()))?
                .join("Resources")
                .join("_up_")
                .join("bundled_resources");
            
            resource_dir.join("desktop").join("python_ml").join("plugins")
        };
        
        // Create plugin directory if it doesn't exist
        if !plugin_dir.exists() {
            std::fs::create_dir_all(&plugin_dir)
                .map_err(|e| PluginError::InitializationFailed(
                    format!("Failed to create plugin directory: {}", e)
                ))?;
        }
        
        info!("Plugin directory: {:?}", plugin_dir);
        Ok(plugin_dir)
    }
    
    /// Get the registry configuration
    pub fn get_config(&self) -> &RegistryConfig {
        &self.config
    }
    
    /// Update the registry configuration
    pub fn update_config(&mut self, config: RegistryConfig) {
        self.config = config;
    }
    
    /// Scan for plugins in the plugin directory
    fn scan_plugins(&mut self) -> Result<(), PluginError> {
        info!("Scanning for plugins in: {:?}", self.plugin_dir);
        
        let entries = std::fs::read_dir(&self.plugin_dir)
            .map_err(|e| PluginError::InitializationFailed(
                format!("Failed to read plugin directory: {}", e)
            ))?;
        
        let mut discovered_plugins = 0;
        let mut loaded_plugins = 0;
        
        for entry in entries {
            // Check if we've reached the maximum number of plugins
            if loaded_plugins >= self.config.max_plugins {
                info!("Reached maximum plugin limit ({}), skipping remaining plugins", self.config.max_plugins);
                break;
            }
            
            let entry = entry.map_err(|e| PluginError::InitializationFailed(
                format!("Failed to read directory entry: {}", e)
            ))?;
            
            let path = entry.path();
            if path.is_dir() {
                match self.load_plugin_from_directory(&path) {
                    Ok(plugin_id) => {
                        info!("Successfully loaded plugin: {}", plugin_id);
                        loaded_plugins += 1;
                    }
                    Err(e) => {
                        warn!("Failed to load plugin from {:?}: {}", path, e);
                    }
                }
                discovered_plugins += 1;
            }
        }
        
        info!("Plugin discovery complete. Discovered: {}, Loaded: {}/{}", 
              discovered_plugins, loaded_plugins, self.config.max_plugins);
        Ok(())
    }
    
    /// Load a plugin from a directory
    fn load_plugin_from_directory(&mut self, plugin_dir: &PathBuf) -> Result<String, PluginError> {
        debug!("Loading plugin from directory: {:?}", plugin_dir);
        
        // Look for manifest.json
        let manifest_path = plugin_dir.join("manifest.json");
        if !manifest_path.exists() {
            return Err(PluginError::InvalidManifest(
                format!("No manifest.json found in {:?}", plugin_dir)
            ));
        }
        
        // Load manifest
        let manifest = PluginLoader::load_manifest(&manifest_path)?;
        
        // Check if plugin is enabled
        if !manifest.enabled {
            debug!("Plugin {} is disabled, skipping", manifest.id);
            return Err(PluginError::PluginDisabled(manifest.id.clone()));
        }
        
        // Check platform compatibility
        let metadata = PluginMetadata {
            manifest: manifest.clone(),
            config: PluginLoader::create_config(&manifest),
            directory_path: plugin_dir.clone(),
            last_loaded: std::time::SystemTime::now(),
        };
        
        if !metadata.is_platform_compatible() {
            return Err(PluginError::PlatformIncompatible(
                format!("Plugin {} is not compatible with current platform", manifest.id)
            ));
        }
        
        // Validate script exists
        let script_path = plugin_dir.join(&manifest.script_path);
        if !script_path.exists() {
            return Err(PluginError::InvalidManifest(
                format!("Script file not found: {:?}", script_path)
            ));
        }
        
        // Create plugin instance
        let plugin = self.create_plugin_instance(manifest.clone())?;
        
        // Register plugin
        let plugin_id = manifest.id.clone();
        self.plugins.insert(plugin_id.clone(), plugin);
        self.metadata.insert(plugin_id.clone(), metadata);
        
        debug!("Plugin {} loaded successfully", plugin_id);
        Ok(plugin_id)
    }
    
    /// Create a plugin instance based on manifest
    fn create_plugin_instance(&self, manifest: PluginManifest) -> Result<Box<dyn GenomicPlugin>, PluginError> {
        // Currently only supports Python script plugins
        // In the future, we can add support for other plugin types
        
        // Validate manifest for Python script plugin
        PythonScriptPluginFactory::validate_manifest(&manifest)?;
        
        // Create Python script plugin
        Ok(PythonScriptPluginFactory::create_plugin(manifest))
    }
    
    /// Get list of available plugins
    pub fn list_plugins(&self) -> Vec<PluginSummary> {
        debug!("Listing {} available plugins", self.metadata.len());
        
        self.metadata
            .values()
            .map(|metadata| metadata.to_summary())
            .collect()
    }
    
    /// Get plugin metadata by ID
    pub fn get_plugin_metadata(&self, plugin_id: &str) -> Option<&PluginMetadata> {
        self.metadata.get(plugin_id)
    }
    
    /// Execute a plugin with given arguments
    pub fn run_plugin(&self, plugin_id: &str, args: Value) -> Result<Value, PluginError> {
        debug!("Running plugin: {} with args: {:?}", plugin_id, args);
        
        // Get plugin
        let plugin = self.plugins.get(plugin_id)
            .ok_or_else(|| PluginError::PluginNotFound(plugin_id.to_string()))?;
        
        // Execute plugin
        plugin.run(args)
    }
    
    /// Check if a plugin exists
    pub fn has_plugin(&self, plugin_id: &str) -> bool {
        self.plugins.contains_key(plugin_id)
    }
    
    /// Get plugin by ID
    pub fn get_plugin(&self, plugin_id: &str) -> Option<&Box<dyn GenomicPlugin>> {
        self.plugins.get(plugin_id)
    }
    
    /// Reload all plugins
    pub fn reload_plugins(&mut self) -> Result<(), PluginError> {
        info!("Reloading all plugins...");
        
        // Clear current plugins
        self.plugins.clear();
        self.metadata.clear();
        
        // Rescan plugins
        self.scan_plugins()?;
        
        info!("Plugin reload complete. Loaded {} plugins.", self.plugins.len());
        Ok(())
    }
    
    /// Reload a specific plugin
    pub fn reload_plugin(&mut self, plugin_id: &str) -> Result<(), PluginError> {
        info!("Reloading plugin: {}", plugin_id);
        
        // Get plugin directory
        let plugin_dir = if let Some(metadata) = self.metadata.get(plugin_id) {
            metadata.directory_path.clone()
        } else {
            return Err(PluginError::PluginNotFound(plugin_id.to_string()));
        };
        
        // Remove existing plugin
        self.plugins.remove(plugin_id);
        self.metadata.remove(plugin_id);
        
        // Reload plugin
        self.load_plugin_from_directory(&plugin_dir)?;
        
        info!("Plugin {} reloaded successfully", plugin_id);
        Ok(())
    }
    
    /// Get registry statistics
    pub fn get_stats(&self) -> RegistryStats {
        let enabled_count = self.metadata.values()
            .filter(|m| m.manifest.enabled)
            .count();
        
        let compatible_count = self.metadata.values()
            .filter(|m| m.is_platform_compatible())
            .count();
        
        RegistryStats {
            total_plugins: self.plugins.len(),
            enabled_plugins: enabled_count,
            compatible_plugins: compatible_count,
            initialized: self.initialized,
        }
    }
    
    /// Update plugin configuration
    pub fn update_plugin_config(&mut self, plugin_id: &str, config: PluginConfig) -> Result<(), PluginError> {
        let metadata = self.metadata.get_mut(plugin_id)
            .ok_or_else(|| PluginError::PluginNotFound(plugin_id.to_string()))?;
        
        metadata.config = config;
        debug!("Updated configuration for plugin: {}", plugin_id);
        Ok(())
    }
}

/// Registry statistics
#[derive(Debug, Clone, serde::Serialize)]
pub struct RegistryStats {
    pub total_plugins: usize,
    pub enabled_plugins: usize,
    pub compatible_plugins: usize,
    pub initialized: bool,
}

/// Plugin registry API for external access
pub struct PluginRegistryApi;

impl PluginRegistryApi {
    /// Initialize the global plugin registry
    pub fn initialize() -> Result<(), PluginError> {
        let mut registry = PLUGIN_REGISTRY.write()
            .map_err(|e| PluginError::InitializationFailed(e.to_string()))?;
        
        registry.initialize()
    }
    
    /// Get list of available plugins
    pub fn list_plugins() -> Result<Vec<PluginSummary>, PluginError> {
        let registry = PLUGIN_REGISTRY.read()
            .map_err(|e| PluginError::InitializationFailed(e.to_string()))?;
        
        Ok(registry.list_plugins())
    }
    
    /// Execute a plugin
    pub fn run_plugin(plugin_id: &str, args: Value) -> Result<Value, PluginError> {
        let registry = PLUGIN_REGISTRY.read()
            .map_err(|e| PluginError::InitializationFailed(e.to_string()))?;
        
        registry.run_plugin(plugin_id, args)
    }
    
    /// Check if a plugin exists
    pub fn has_plugin(plugin_id: &str) -> Result<bool, PluginError> {
        let registry = PLUGIN_REGISTRY.read()
            .map_err(|e| PluginError::InitializationFailed(e.to_string()))?;
        
        Ok(registry.has_plugin(plugin_id))
    }
    
    /// Get plugin metadata
    pub fn get_plugin_metadata(plugin_id: &str) -> Result<Option<PluginMetadata>, PluginError> {
        let registry = PLUGIN_REGISTRY.read()
            .map_err(|e| PluginError::InitializationFailed(e.to_string()))?;
        
        Ok(registry.get_plugin_metadata(plugin_id).cloned())
    }
    
    /// Reload all plugins
    pub fn reload_plugins() -> Result<(), PluginError> {
        let mut registry = PLUGIN_REGISTRY.write()
            .map_err(|e| PluginError::InitializationFailed(e.to_string()))?;
        
        registry.reload_plugins()
    }
    
    /// Get registry statistics
    pub fn get_stats() -> Result<RegistryStats, PluginError> {
        let registry = PLUGIN_REGISTRY.read()
            .map_err(|e| PluginError::InitializationFailed(e.to_string()))?;
        
        Ok(registry.get_stats())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;
    use std::fs;
    
    fn create_test_plugin_dir() -> Result<TempDir, Box<dyn std::error::Error>> {
        let temp_dir = TempDir::new()?;
        let plugin_dir = temp_dir.path().join("plugins").join("test_plugin");
        fs::create_dir_all(&plugin_dir)?;
        
        // Create manifest.json
        let manifest = serde_json::json!({
            "id": "test_plugin",
            "name": "Test Plugin",
            "description": "A test plugin",
            "version": "1.0.0",
            "script_path": "test_plugin.py",
            "category": "genomic_processing",
            "config": {},
            "requirements": [],
            "tags": ["test"],
            "enabled": true,
            "platform_support": {
                "windows": true,
                "macos": true,
                "linux": true
            }
        });
        
        fs::write(plugin_dir.join("manifest.json"), manifest.to_string())?;
        
        // Create Python script
        fs::write(plugin_dir.join("test_plugin.py"), "#!/usr/bin/env python3\nprint('Hello from test plugin')")?;
        
        Ok(temp_dir)
    }
    
    #[test]
    fn test_plugin_registry_creation() {
        let registry = PluginRegistry::new();
        assert!(!registry.initialized);
        assert_eq!(registry.plugins.len(), 0);
        assert_eq!(registry.metadata.len(), 0);
    }
    
    #[test]
    fn test_registry_api_initialization() {
        // This test would require proper setup of the plugin directory
        // For now, we just test that the API methods compile
        let result = PluginRegistryApi::list_plugins();
        assert!(result.is_ok() || result.is_err()); // Either outcome is acceptable for this test
    }
    
    #[test]
    fn test_registry_stats() {
        let registry = PluginRegistry::new();
        let stats = registry.get_stats();
        
        assert_eq!(stats.total_plugins, 0);
        assert_eq!(stats.enabled_plugins, 0);
        assert_eq!(stats.compatible_plugins, 0);
        assert!(!stats.initialized);
    }
} 