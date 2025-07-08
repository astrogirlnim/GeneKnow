use crate::plugin::{GenomicPlugin, PluginError, PluginManifest, PluginExecutionResult};
use crate::utils::execute_python;
use serde_json::Value;
use std::time::Instant;
use log::{info, debug, error};

/// Python script plugin implementation
/// 
/// This plugin adapter bridges the GenomicPlugin trait to the existing
/// execute_python utility, allowing Python scripts to be executed through
/// the plugin system.
pub struct PythonScriptPlugin {
    manifest: PluginManifest,
}

impl PythonScriptPlugin {
    /// Create a new Python script plugin instance
    /// 
    /// # Arguments
    /// * `manifest` - Plugin manifest containing metadata and configuration
    /// 
    /// # Returns
    /// Returns a new PythonScriptPlugin instance
    pub fn new(manifest: PluginManifest) -> Self {
        debug!("Creating PythonScriptPlugin for: {} ({})", manifest.name, manifest.id);
        Self { manifest }
    }
    
    /// Convert JSON arguments to command line arguments
    /// 
    /// This method converts the JSON plugin arguments into command line
    /// arguments suitable for passing to the Python script.
    /// 
    /// # Arguments
    /// * `args` - JSON arguments from plugin execution
    /// 
    /// # Returns
    /// Returns a vector of command line arguments
    fn json_to_args(&self, args: &Value) -> Vec<String> {
        let mut cmd_args = Vec::new();
        
        // Always add --json flag to request JSON output
        cmd_args.push("--json".to_string());
        
        // Convert JSON object to command line arguments
        if let Some(obj) = args.as_object() {
            for (key, value) in obj {
                // Skip internal plugin parameters
                if key.starts_with("_") {
                    continue;
                }
                
                // Add argument flag
                cmd_args.push(format!("--{}", key.replace("_", "-")));
                
                // Add argument value based on type
                match value {
                    Value::String(s) => cmd_args.push(s.clone()),
                    Value::Number(n) => cmd_args.push(n.to_string()),
                    Value::Bool(b) => {
                        if *b {
                            // For boolean flags, we don't add a value
                            continue;
                        } else {
                            // Remove the flag if boolean is false
                            cmd_args.pop();
                        }
                    }
                    Value::Array(arr) => {
                        // For arrays, add multiple values
                        for item in arr {
                            if let Some(s) = item.as_str() {
                                cmd_args.push(s.to_string());
                            }
                        }
                    }
                    _ => {
                        // For other types, convert to string
                        cmd_args.push(value.to_string());
                    }
                }
            }
        }
        
        debug!("Converted JSON args to command line: {:?}", cmd_args);
        cmd_args
    }
    
    /// Parse Python script output into plugin execution result
    /// 
    /// # Arguments
    /// * `output` - Command output from Python script execution
    /// * `execution_time` - Time taken for script execution
    /// 
    /// # Returns
    /// Returns a PluginExecutionResult
    fn parse_output(&self, output: &std::process::Output, execution_time: f64) -> PluginExecutionResult {
        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        
        // Split logs from actual output
        let mut logs = Vec::new();
        let mut json_output = None;
        
        // Parse stdout line by line
        for line in stdout.lines() {
            if line.trim().starts_with('{') || line.trim().starts_with('[') {
                // Potential JSON output - try to parse
                if let Ok(json_value) = serde_json::from_str::<Value>(line.trim()) {
                    json_output = Some(json_value);
                } else {
                    logs.push(line.to_string());
                }
            } else {
                logs.push(line.to_string());
            }
        }
        
        // Add stderr to logs
        if !stderr.is_empty() {
            logs.push(format!("STDERR: {}", stderr));
        }
        
        // Determine success based on exit code and JSON output
        let success = output.status.success() && json_output.is_some();
        
        // Create execution result
        PluginExecutionResult {
            success,
            plugin_id: self.manifest.id.clone(),
            execution_time,
            result: json_output,
            error: if !success {
                Some(if !stderr.is_empty() {
                    stderr.to_string()
                } else {
                    "Script execution failed".to_string()
                })
            } else {
                None
            },
            logs,
        }
    }
    
    /// Get script name without extension from manifest
    fn get_script_name(&self) -> String {
        let script_path = &self.manifest.script_path;
        if script_path.ends_with(".py") {
            script_path[..script_path.len() - 3].to_string()
        } else {
            script_path.clone()
        }
    }
}

impl GenomicPlugin for PythonScriptPlugin {
    fn id(&self) -> &'static str {
        // Note: This leaks memory, but it's acceptable for plugin IDs
        // which are typically created once and used for the lifetime of the app
        Box::leak(self.manifest.id.clone().into_boxed_str())
    }
    
    fn name(&self) -> &'static str {
        Box::leak(self.manifest.name.clone().into_boxed_str())
    }
    
    fn description(&self) -> &'static str {
        Box::leak(self.manifest.description.clone().into_boxed_str())
    }
    
    fn version(&self) -> &'static str {
        Box::leak(self.manifest.version.clone().into_boxed_str())
    }
    
    fn run(&self, args: Value) -> Result<Value, PluginError> {
        info!("Executing Python plugin: {} with args: {:?}", self.manifest.id, args);
        
        // Check if plugin is enabled
        if !self.manifest.enabled {
            return Err(PluginError::PluginDisabled(self.manifest.id.clone()));
        }
        
        // Validate arguments
        self.validate_args(&args)?;
        
        // Convert JSON args to command line args
        let cmd_args = self.json_to_args(&args);
        let cmd_args_str: Vec<&str> = cmd_args.iter().map(|s| s.as_str()).collect();
        
        // Execute Python script
        let start_time = Instant::now();
        let script_name = self.get_script_name();
        
        info!("Running Python script: {} with args: {:?}", script_name, cmd_args_str);
        
        let output = execute_python(&script_name, &cmd_args_str)
            .map_err(|e| PluginError::ExecutionFailed(e.to_string()))?;
        
        let execution_time = start_time.elapsed().as_secs_f64();
        
        // Parse output into plugin result
        let result = self.parse_output(&output, execution_time);
        
        if result.success {
            info!("Plugin {} executed successfully in {:.2}s", self.manifest.id, execution_time);
            Ok(serde_json::to_value(result)?)
        } else {
            error!("Plugin {} execution failed: {:?}", self.manifest.id, result.error);
            Err(PluginError::ExecutionFailed(
                result.error.unwrap_or_else(|| "Unknown error".to_string())
            ))
        }
    }
    
    fn manifest(&self) -> &PluginManifest {
        &self.manifest
    }
    
    fn validate_args(&self, args: &Value) -> Result<(), PluginError> {
        debug!("Validating arguments for plugin {}: {:?}", self.manifest.id, args);
        
        // Basic validation - check if args is an object
        if !args.is_object() && !args.is_null() {
            return Err(PluginError::ValidationFailed(
                "Plugin arguments must be a JSON object".to_string()
            ));
        }
        
        // If input schema is provided, validate against it
        if let Some(_input_schema) = &self.manifest.input_schema {
            // TODO: Implement JSON schema validation
            // For now, we'll just do basic validation
            debug!("Input schema validation not yet implemented for plugin {}", self.manifest.id);
        }
        
        // Plugin-specific validation can be added here
        // For now, we'll just log the validation
        debug!("Argument validation passed for plugin {}", self.manifest.id);
        
        Ok(())
    }
}

/// Factory for creating Python script plugins
pub struct PythonScriptPluginFactory;

impl PythonScriptPluginFactory {
    /// Create a new Python script plugin from manifest
    /// 
    /// # Arguments
    /// * `manifest` - Plugin manifest
    /// 
    /// # Returns
    /// Returns a boxed GenomicPlugin trait object
    pub fn create_plugin(manifest: PluginManifest) -> Box<dyn GenomicPlugin> {
        Box::new(PythonScriptPlugin::new(manifest))
    }
    
    /// Validate that a manifest is suitable for Python script plugin
    /// 
    /// # Arguments
    /// * `manifest` - Plugin manifest to validate
    /// 
    /// # Returns
    /// Returns Ok if manifest is valid for Python script plugin
    pub fn validate_manifest(manifest: &PluginManifest) -> Result<(), PluginError> {
        // Check that script path ends with .py
        if !manifest.script_path.ends_with(".py") {
            return Err(PluginError::InvalidManifest(
                "Python script plugin must have a .py script file".to_string()
            ));
        }
        
        // Check that category is appropriate
        let valid_categories = vec![
            "genomic_processing",
            "data_analysis", 
            "file_conversion",
            "quality_control",
            "visualization",
            "utilities"
        ];
        
        if !valid_categories.contains(&manifest.category.as_str()) {
            debug!("Plugin category '{}' is not in recommended categories: {:?}", 
                   manifest.category, valid_categories);
        }
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin::{PluginManifest, PlatformSupport};
    use std::collections::HashMap;
    
    fn create_test_manifest() -> PluginManifest {
        let mut config = HashMap::new();
        config.insert("threads".to_string(), serde_json::json!(4));
        
        PluginManifest {
            id: "test_python_plugin".to_string(),
            name: "Test Python Plugin".to_string(),
            description: "A test Python plugin".to_string(),
            version: "1.0.0".to_string(),
            author: Some("Test Author".to_string()),
            script_path: "test_script.py".to_string(),
            category: "genomic_processing".to_string(),
            config,
            input_schema: None,
            output_schema: None,
            requirements: vec!["python>=3.8".to_string()],
            tags: vec!["test".to_string()],
            enabled: true,
            platform_support: PlatformSupport::default(),
        }
    }
    
    #[test]
    fn test_python_plugin_creation() {
        let manifest = create_test_manifest();
        let plugin = PythonScriptPlugin::new(manifest.clone());
        
        assert_eq!(plugin.id(), "test_python_plugin");
        assert_eq!(plugin.name(), "Test Python Plugin");
        assert_eq!(plugin.description(), "A test Python plugin");
        assert_eq!(plugin.version(), "1.0.0");
    }
    
    #[test]
    fn test_json_to_args_conversion() {
        let manifest = create_test_manifest();
        let plugin = PythonScriptPlugin::new(manifest);
        
        let args = serde_json::json!({
            "input_file": "/path/to/input.txt",
            "output_dir": "/path/to/output",
            "threads": 4,
            "verbose": true,
            "skip_validation": false
        });
        
        let cmd_args = plugin.json_to_args(&args);
        
        // Should contain --json flag
        assert!(cmd_args.contains(&"--json".to_string()));
        
        // Should contain input file argument
        assert!(cmd_args.contains(&"--input-file".to_string()));
        assert!(cmd_args.contains(&"/path/to/input.txt".to_string()));
        
        // Should contain threads argument
        assert!(cmd_args.contains(&"--threads".to_string()));
        assert!(cmd_args.contains(&"4".to_string()));
        
        // Should contain verbose flag but not skip-validation
        assert!(cmd_args.contains(&"--verbose".to_string()));
        assert!(!cmd_args.contains(&"--skip-validation".to_string()));
    }
    
    #[test]
    fn test_script_name_extraction() {
        let manifest = create_test_manifest();
        let plugin = PythonScriptPlugin::new(manifest);
        
        assert_eq!(plugin.get_script_name(), "test_script");
    }
    
    #[test]
    fn test_manifest_validation() {
        let manifest = create_test_manifest();
        assert!(PythonScriptPluginFactory::validate_manifest(&manifest).is_ok());
        
        let mut invalid_manifest = manifest.clone();
        invalid_manifest.script_path = "invalid_script.txt".to_string();
        assert!(PythonScriptPluginFactory::validate_manifest(&invalid_manifest).is_err());
    }
    
    #[test]
    fn test_plugin_factory() {
        let manifest = create_test_manifest();
        let plugin = PythonScriptPluginFactory::create_plugin(manifest.clone());
        
        assert_eq!(plugin.id(), "test_python_plugin");
        assert_eq!(plugin.name(), "Test Python Plugin");
    }
} 