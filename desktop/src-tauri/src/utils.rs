use std::path::PathBuf;
use std::process::Command;
use std::io::Result;
use log::{info, error, debug};

/// Execute a Python script with the given arguments
/// 
/// This function provides a centralized way to execute Python scripts with proper
/// cross-platform path handling, logging, and error reporting.
/// 
/// # Arguments
/// * `script_name` - The name of the Python script (without .py extension)
/// * `args` - A slice of string arguments to pass to the script
/// 
/// # Returns
/// Returns the command output as a Result<std::process::Output, String>
/// 
/// # Examples
/// ```
/// let result = execute_python("fastq_to_vcf_pipeline", &["--help"]);
/// match result {
///     Ok(output) => println!("Success: {}", String::from_utf8_lossy(&output.stdout)),
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
pub fn execute_python(script_name: &str, args: &[&str]) -> Result<std::process::Output> {
    info!("Executing Python script: {} with args: {:?}", script_name, args);
    
    // Build the script path - scripts are in desktop/python_ml/
    let script_path = get_python_script_path(script_name)?;
    
    debug!("Script path resolved to: {:?}", script_path);
    
    // Build the command
    let mut cmd = Command::new("python3");
    cmd.arg(&script_path);
    
    // Add all arguments
    for arg in args {
        cmd.arg(arg);
    }
    
    // Set working directory to the project root for relative path access
    if let Some(project_root) = get_project_root() {
        debug!("Working directory set to: {:?}", project_root);
        cmd.current_dir(project_root);
    }
    
    // Log the full command being executed
    info!("Running command: python3 {:?} {:?}", script_path, args);
    
    // Execute the command
    let output = cmd.output();
    
    match &output {
        Ok(result) => {
            if result.status.success() {
                info!("Python script executed successfully in {:.2}s", 
                     get_execution_time_placeholder());
                debug!("Script stdout: {}", String::from_utf8_lossy(&result.stdout));
            } else {
                error!("Python script failed with exit code: {:?}", result.status.code());
                error!("Script stderr: {}", String::from_utf8_lossy(&result.stderr));
            }
        }
        Err(e) => {
            error!("Failed to execute Python script: {}", e);
        }
    }
    
    output
}

/// Get the full path to a Python script in the python_ml directory
/// 
/// # Arguments
/// * `script_name` - The name of the script without .py extension
/// 
/// # Returns
/// Returns the full path to the script as a PathBuf
fn get_python_script_path(script_name: &str) -> Result<PathBuf> {
    // Get the directory where the Tauri app is running from
    let mut script_path = get_project_root()
        .ok_or_else(|| std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "Could not determine project root"
        ))?;
    
    // Add the python_ml directory and script name
    script_path.push("desktop");
    script_path.push("python_ml");
    script_path.push(format!("{}.py", script_name));
    
    // Verify the script exists
    if !script_path.exists() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("Python script not found: {:?}", script_path)
        ));
    }
    
    Ok(script_path)
}

/// Get the project root directory
/// 
/// # Returns
/// Returns the project root directory as an Option<PathBuf>
fn get_project_root() -> Option<PathBuf> {
    // In development, we want to find the project root by looking for key files
    let current_dir = std::env::current_dir().ok()?;
    
    // Look for indicators that we're in the project root
    let mut path = current_dir;
    loop {
        // Check for key project files
        if path.join("README.md").exists() && 
           path.join("desktop").exists() && 
           path.join("docs").exists() {
            return Some(path);
        }
        
        // Move up one directory
        if let Some(parent) = path.parent() {
            path = parent.to_path_buf();
        } else {
            break;
        }
    }
    
    // Fallback to current directory
    std::env::current_dir().ok()
}

/// Placeholder for execution time calculation
/// TODO: Implement proper timing mechanism
fn get_execution_time_placeholder() -> f64 {
    0.0
}

/// Execute a Python script and parse the JSON output
/// 
/// This is a convenience function that executes a Python script and parses
/// the JSON output into a serde_json::Value.
/// 
/// # Arguments
/// * `script_name` - The name of the Python script (without .py extension)
/// * `args` - A slice of string arguments to pass to the script
/// 
/// # Returns
/// Returns the parsed JSON output as a Result<serde_json::Value, String>
pub fn execute_python_json(script_name: &str, args: &[&str]) -> std::result::Result<serde_json::Value, String> {
    let output = execute_python(script_name, args)
        .map_err(|e| format!("Failed to execute Python script: {}", e))?;
    
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("Python script failed: {}", stderr));
    }
    
    let stdout = String::from_utf8_lossy(&output.stdout);
    debug!("Parsing JSON output: {}", stdout);
    
    // Parse the JSON output
    serde_json::from_str(stdout.trim())
        .map_err(|e| format!("Failed to parse JSON output: {}", e))
}

/// Convert file paths to OS-specific format
/// 
/// This function takes a path string and converts it to the appropriate
/// format for the current operating system.
/// 
/// # Arguments
/// * `path` - The path string to convert
/// 
/// # Returns
/// Returns the converted path as a PathBuf
pub fn normalize_path(path: &str) -> PathBuf {
    PathBuf::from(path)
}

/// Escape spaces in file paths for Windows compatibility
/// 
/// # Arguments
/// * `path` - The path to escape
/// 
/// # Returns
/// Returns the escaped path as a String
pub fn escape_path_spaces(path: &str) -> String {
    if cfg!(target_os = "windows") && path.contains(' ') {
        format!("\"{}\"", path)
    } else {
        path.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_path() {
        let path = normalize_path("test/path/file.txt");
        assert!(path.to_string_lossy().contains("file.txt"));
    }

    #[test]
    fn test_escape_path_spaces() {
        let path_with_spaces = "path with spaces/file.txt";
        let escaped = escape_path_spaces(path_with_spaces);
        
        if cfg!(target_os = "windows") {
            assert!(escaped.starts_with('"') && escaped.ends_with('"'));
        } else {
            assert_eq!(escaped, path_with_spaces);
        }
    }
} 