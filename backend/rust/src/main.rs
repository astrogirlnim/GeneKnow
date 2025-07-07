// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

fn main() {
    // Initialize logging
    tracing_subscriber::fmt()
        .with_env_filter("genepredict=debug,info")
        .init();
    
    tracing::info!("🧬 Starting GenePredict - AI for Genomic Risk Assessment");
    tracing::info!("🚀 Initializing Tauri application...");
    
    genepredict_lib::run()
}
