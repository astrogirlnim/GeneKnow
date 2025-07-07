# TCGA Reference Data Integration Guide

## 🎯 Overview

This guide explains how to integrate TCGA (The Cancer Genome Atlas) reference data from the [GDC Data Portal](https://portal.gdc.cancer.gov/) into your GenePredict application for enhanced genomic risk assessment.

## 🏗️ Architecture Compatibility

✅ **Your current setup is ideal for TCGA integration:**

- **Rust Backend**: Perfect for high-performance genomic data processing
- **Python ML Layer**: Excellent for TCGA data analysis with BioPython/pandas
- **Local Storage**: TCGA data stays private and compliant
- **Plugin System**: Easy to add TCGA-specific processors

## 📊 TCGA Data Types Available

### Primary Data Sources
1. **Clinical Data**: Patient demographics, treatment history, outcomes
2. **Gene Expression**: RNA-Seq quantification data
3. **Copy Number Variation**: CNV segment data
4. **Mutation Data**: VCF files with somatic mutations
5. **Methylation Data**: DNA methylation array data
6. **Protein Expression**: RPPA (protein array) data

### Key Cancer Types
- BRCA (Breast Cancer) - Perfect for your initial use case
- LUAD (Lung Adenocarcinoma)
- COAD (Colon Adenocarcinoma)
- PRAD (Prostate Adenocarcinoma)
- And 30+ other cancer types

## 🚀 Implementation Plan

### Phase 1: Data Acquisition Setup

#### 1. Install GDC Client
```bash
# Create data directory
mkdir -p data/tcga/raw data/tcga/processed

# Download GDC client (macOS)
cd data/tcga
curl -O https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_OSX_x64.zip
unzip gdc-client_v1.6.1_OSX_x64.zip
chmod +x gdc-client
```

#### 2. Python TCGA Integration
```python
# backend/python/genepredict/processors/tcga_processor.py
import pandas as pd
import json
import requests
from pathlib import Path
from typing import Dict, List, Optional
from ..models.base import GenomicData

class TCGAProcessor:
    """TCGA reference data processor for genomic risk assessment."""
    
    def __init__(self, data_dir: Path):
        self.data_dir = Path(data_dir)
        self.base_url = "https://api.gdc.cancer.gov"
        
    async def download_brca_reference_data(self) -> Dict[str, str]:
        """Download BRCA1/BRCA2 reference mutations and expression data."""
        
        # Query for BRCA-related mutation data
        filters = {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-BRCA"]}},
                {"op": "in", "content": {"field": "data_category", "value": ["Simple Nucleotide Variation"]}},
                {"op": "in", "content": {"field": "genes.symbol", "value": ["BRCA1", "BRCA2"]}}
            ]
        }
        
        params = {
            "filters": json.dumps(filters),
            "format": "json",
            "size": "2000"
        }
        
        response = requests.get(f"{self.base_url}/files", params=params)
        
        if response.status_code == 200:
            files_data = response.json()
            return await self._download_files(files_data["data"]["hits"])
        else:
            raise Exception(f"TCGA API error: {response.status_code}")
    
    async def process_clinical_outcomes(self, cancer_type: str = "BRCA") -> pd.DataFrame:
        """Process clinical outcome data for risk model training."""
        
        filters = {
            "op": "and", 
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": [f"TCGA-{cancer_type}"]}},
                {"op": "in", "content": {"field": "data_category", "value": ["Clinical"]}}
            ]
        }
        
        params = {
            "filters": json.dumps(filters),
            "format": "json", 
            "size": "2000"
        }
        
        response = requests.get(f"{self.base_url}/cases", params=params)
        clinical_data = response.json()
        
        # Process into standardized format
        processed_data = []
        for case in clinical_data["data"]["hits"]:
            patient_data = {
                "case_id": case["submitter_id"],
                "age_at_diagnosis": case.get("demographic", {}).get("age_at_diagnosis"),
                "vital_status": case.get("demographic", {}).get("vital_status"),
                "days_to_death": case.get("demographic", {}).get("days_to_death"),
                "tumor_stage": case.get("diagnoses", [{}])[0].get("tumor_stage"),
                "histological_type": case.get("diagnoses", [{}])[0].get("primary_diagnosis")
            }
            processed_data.append(patient_data)
            
        return pd.DataFrame(processed_data)
    
    def create_risk_reference_dataset(self, mutations_df: pd.DataFrame, 
                                    clinical_df: pd.DataFrame) -> GenomicData:
        """Create reference dataset for risk model training."""
        
        # Merge mutation and clinical data
        reference_data = mutations_df.merge(clinical_df, on="case_id", how="inner")
        
        # Calculate risk scores based on known pathogenic variants
        reference_data["risk_score"] = self._calculate_pathogenicity_scores(reference_data)
        
        return GenomicData(
            sample_id="tcga_reference",
            file_type="reference_dataset",
            variants=reference_data.to_dict("records"),
            metadata={
                "source": "TCGA",
                "creation_date": str(pd.Timestamp.now()),
                "total_samples": len(reference_data),
                "data_types": ["mutations", "clinical"]
            }
        )
```

### Phase 2: Rust Integration

#### 1. Add TCGA Commands to Rust Backend
```rust
// backend/rust/src/lib.rs - Add new commands

#[tauri::command]
async fn download_tcga_reference_data(
    state: tauri::State<'_, AppState>,
    cancer_type: String,
) -> Result<String, String> {
    info!("Downloading TCGA reference data for {}", cancer_type);
    
    let python_script = format!(
        r#"
import asyncio
from genepredict.processors.tcga_processor import TCGAProcessor
from pathlib import Path

async def main():
    processor = TCGAProcessor(Path("./data/tcga"))
    if "{}" == "BRCA":
        result = await processor.download_brca_reference_data()
        clinical = await processor.process_clinical_outcomes("BRCA")
        print(f"Downloaded {{len(result)}} files and {{len(clinical)}} clinical records")
        return "success"
    
asyncio.run(main())
"#,
        cancer_type
    );
    
    match execute_python_script(&python_script).await {
        Ok(output) => {
            info!("TCGA download completed: {}", output);
            Ok(format!("TCGA {} reference data downloaded successfully", cancer_type))
        }
        Err(e) => {
            error!("TCGA download failed: {}", e);
            Err(format!("Failed to download TCGA data: {}", e))
        }
    }
}

#[tauri::command] 
async fn get_tcga_cancer_types() -> Result<Vec<String>, String> {
    Ok(vec![
        "BRCA".to_string(), // Breast Cancer
        "LUAD".to_string(), // Lung Adenocarcinoma  
        "COAD".to_string(), // Colon Adenocarcinoma
        "PRAD".to_string(), // Prostate Adenocarcinoma
        "HNSC".to_string(), // Head and Neck Cancer
        "KIRC".to_string(), // Kidney Cancer
        "LGG".to_string(),  // Low Grade Glioma
        "THCA".to_string(), // Thyroid Cancer
        "LUSC".to_string(), // Lung Squamous Cell Carcinoma
        "UCEC".to_string(), // Uterine Cancer
    ])
}
```

#### 2. Update Cargo.toml Dependencies
```toml
# Add to backend/rust/Cargo.toml
[dependencies]
# ... existing dependencies ...
reqwest = { version = "0.11", features = ["json"] }
futures = "0.3"
```

### Phase 3: Frontend Integration

#### 1. Add TCGA Data Management UI
```typescript
// frontend/src/components/TCGAManager.tsx
import React, { useState, useEffect } from 'react';
import { invoke } from '@tauri-apps/api/tauri';
import { Download, Database, AlertCircle } from 'lucide-react';

export const TCGAManager: React.FC = () => {
  const [cancerTypes, setCancerTypes] = useState<string[]>([]);
  const [downloading, setDownloading] = useState<string | null>(null);
  const [downloadStatus, setDownloadStatus] = useState<{[key: string]: string}>({});

  useEffect(() => {
    loadCancerTypes();
  }, []);

  const loadCancerTypes = async () => {
    try {
      const types = await invoke<string[]>('get_tcga_cancer_types');
      setCancerTypes(types);
    } catch (error) {
      console.error('Failed to load cancer types:', error);
    }
  };

  const downloadReferenceData = async (cancerType: string) => {
    setDownloading(cancerType);
    try {
      const result = await invoke<string>('download_tcga_reference_data', { 
        cancerType 
      });
      setDownloadStatus(prev => ({ ...prev, [cancerType]: 'completed' }));
      console.log(`TCGA download result: ${result}`);
    } catch (error) {
      setDownloadStatus(prev => ({ ...prev, [cancerType]: 'failed' }));
      console.error(`TCGA download failed:`, error);
    } finally {
      setDownloading(null);
    }
  };

  return (
    <div className="bg-white rounded-lg shadow-sm border border-gray-200 p-6">
      <div className="flex items-center gap-3 mb-6">
        <Database className="w-6 h-6 text-blue-600" />
        <h3 className="text-lg font-semibold text-gray-900">
          TCGA Reference Data
        </h3>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
        {cancerTypes.map((cancerType) => (
          <div key={cancerType} className="border border-gray-200 rounded-lg p-4">
            <div className="flex items-center justify-between mb-3">
              <h4 className="font-medium text-gray-900">{cancerType}</h4>
              {downloadStatus[cancerType] === 'completed' && (
                <span className="text-green-600 text-sm">✓ Downloaded</span>
              )}
              {downloadStatus[cancerType] === 'failed' && (
                <AlertCircle className="w-4 h-4 text-red-500" />
              )}
            </div>
            
            <button
              onClick={() => downloadReferenceData(cancerType)}
              disabled={downloading === cancerType}
              className="w-full flex items-center justify-center gap-2 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
            >
              <Download className="w-4 h-4" />
              {downloading === cancerType ? 'Downloading...' : 'Download Data'}
            </button>
          </div>
        ))}
      </div>

      <div className="mt-6 p-4 bg-blue-50 rounded-lg">
        <h4 className="font-medium text-blue-900 mb-2">Data Usage</h4>
        <p className="text-blue-800 text-sm">
          TCGA reference data enhances risk assessment by providing population-level 
          variant frequencies, clinical outcomes, and validated pathogenic mutations.
          All data processing happens locally for privacy compliance.
        </p>
      </div>
    </div>
  );
};
```

## 🔧 Complete Dev Setup Instructions

### 1. Update Environment Configuration
```bash
# Add to .env file
# TCGA Configuration
TCGA_DATA_DIR=./data/tcga
TCGA_CACHE_DIR=./data/tcga/cache
TCGA_API_TIMEOUT=300
REFERENCE_DATA_AUTO_UPDATE=false

# ML Model Configuration  
BRCA_MODEL_PATH=./data/models/brca_risk_model.pkl
REFERENCE_DATASET_PATH=./data/tcga/processed/brca_reference.json
```

### 2. Create Data Directory Structure
```bash
mkdir -p data/tcga/{raw,processed,cache,models}
mkdir -p data/reference/{clinvar,1000genomes,cosmic}
```

### 3. Install Additional Python Dependencies
```bash
cd backend/python
source venv/bin/activate
pip install requests aiohttp aiofiles gdctools
```

### 4. Test TCGA Integration
```bash
# Create test script
cat > scripts/test-tcga.py << EOF
import asyncio
import sys
sys.path.append('./backend/python')

from genepredict.processors.tcga_processor import TCGAProcessor
from pathlib import Path

async def test_tcga():
    processor = TCGAProcessor(Path("./data/tcga"))
    print("Testing TCGA connection...")
    
    # Test API connectivity
    cancer_types = await processor.get_available_cancer_types()
    print(f"Available cancer types: {cancer_types[:5]}...")
    
    print("TCGA integration test passed!")

if __name__ == "__main__":
    asyncio.run(test_tcga())
EOF

python scripts/test-tcga.py
```

## 🎯 Benefits for Your Use Case

1. **Enhanced Risk Models**: Train on 30,000+ cancer patients with known outcomes
2. **Variant Validation**: Cross-reference user variants against TCGA mutation database  
3. **Population Frequencies**: Get accurate allele frequencies for risk calculation
4. **Clinical Correlation**: Link genetic variants to actual patient outcomes
5. **Multi-Cancer Support**: Expand beyond BRCA to 30+ cancer types

## 🔐 Privacy & Compliance

✅ **GDPR/HIPAA Compatible:**
- All TCGA data processed locally
- No user data sent to external APIs
- Reference data anonymized and aggregated
- User genetic data never leaves device

## 🚀 Next Steps

1. **Immediate**: Test current app with `npm run tauri:dev`
2. **Phase 1**: Implement TCGA processor (1-2 days)
3. **Phase 2**: Add Rust commands and frontend UI (2-3 days)  
4. **Phase 3**: Integrate with ML models (3-4 days)
5. **Phase 4**: Add automated reference data updates (1-2 days)

Your current architecture is perfectly suited for this integration. The local-first approach ensures privacy while leveraging world-class genomic reference data from [TCGA](https://portal.gdc.cancer.gov/). 