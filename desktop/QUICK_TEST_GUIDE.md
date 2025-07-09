# Quick Test Guide - Frontend-Backend Integration

This guide helps you quickly verify that the frontend and backend are properly connected.

## Prerequisites

1. Install dependencies:
   ```bash
   # Python dependencies
   cd geneknow_pipeline
   pip install -r requirements.txt
   
   # Frontend dependencies
   cd ../desktop/ui
   pnpm install
   ```

2. Make sure you have a test file. You can use:
   - `test_data/tcga_downloads/*.maf.gz` - Real TCGA MAF files
   - Any VCF, FASTQ, or BAM file

## Quick Test Steps

### 1. Start the Application

```bash
cd desktop/ui
pnpm tauri-dev
```

This will:
- Start the React frontend
- Launch the Tauri desktop app
- The Python API server will auto-start

### 2. Test File Processing

1. **Navigate to Upload Page**
   - Click "Get Started" on the welcome page
   - Or go directly to the Upload page

2. **Select a File**
   - Click the upload area or drag & drop a file
   - Use your system's file browser to select a genomic file (MAF, VCF, FASTQ, BAM)
   - The file will be read and temporarily saved
   - The file name and temporary path will be displayed

3. **Start Analysis**
   - Click "Start Analysis"
   - You should see:
     - Progress bar showing processing steps
     - Real-time updates (e.g., "Processing variants", "Calculating risk scores")
     - Processing percentage

4. **View Results**
   - After completion, you'll be redirected to the dashboard
   - You should see:
     - Risk probability percentage
     - Hazard score
     - Analysis report sections
     - Variant details
     - Risk scores for different cancer types

### 3. Verify Backend Connection

Check the terminal running `pnpm tauri-dev` for:
```
[INFO] GeneKnow API server started successfully
[INFO] Processing genomic file: /path/to/your/file
```

### 4. Test Mock Data

1. On the upload page, click one of the mock patient cards
2. The system will process mock data and show results
3. This tests the pipeline even without a real file

## What to Expect

### Successful Processing
- File selection works with native dialog
- Progress updates show during processing
- Dashboard displays actual results:
  - Variant counts
  - Risk scores
  - Processing time
  - Report sections with findings

### Common Issues

**"Failed to start API server"**
- Check if port 5001 is already in use
- Ensure Python is in your PATH

**"Failed to process file"**
- Verify the file format is supported
- Check file permissions
- Look at console logs for specific errors

**No progress updates**
- The backend might be processing - check terminal logs
- Processing can take 10-60 seconds depending on file size

## Debugging

1. **Check API Health**
   ```bash
   curl http://localhost:5001/api/health
   ```

2. **View Tauri Console**
   - Press F12 in the desktop app
   - Check for JavaScript errors
   - Look at Network tab for API calls

3. **Python API Logs**
   - Check the terminal for Python errors
   - Look for LangGraph processing steps

## Next Steps

Once basic processing works:
1. Try different file formats
2. Test larger files
3. Check error handling with invalid files
4. Verify all dashboard metrics display correctly

The integration is working if you can:
✅ Select files using native dialog
✅ See real-time progress updates
✅ View actual genomic analysis results
✅ Navigate between pages with data persistence 