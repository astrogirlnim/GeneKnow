use thiserror::Error;

/// Custom error type for GenePredict
#[derive(Error, Debug)]
pub enum GenePredicateError {
    #[error("Configuration error: {0}")]
    Config(String),

    #[error("File processing error: {0}")]
    FileProcessing(String),

    #[error("Plugin error: {0}")]
    Plugin(String),

    #[error("ML engine error: {0}")]
    MLEngine(String),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Serialization error: {0}")]
    Serialization(#[from] serde_json::Error),

    #[error("Configuration loading error: {0}")]
    ConfigError(#[from] config::ConfigError),

    #[error("Unknown error: {0}")]
    Unknown(String),
}

/// Result type alias for GenePredict operations
pub type Result<T> = std::result::Result<T, GenePredicateError>;

impl From<String> for GenePredicateError {
    fn from(s: String) -> Self {
        GenePredicateError::Unknown(s)
    }
}

impl From<&str> for GenePredicateError {
    fn from(s: &str) -> Self {
        GenePredicateError::Unknown(s.to_string())
    }
} 