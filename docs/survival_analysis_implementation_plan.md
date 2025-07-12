# Survival Analysis Visualization Implementation Plan

## Executive Summary

This document outlines the implementation plan for adding survival analysis visualization to the GeneKnow desktop application. The feature will provide patients with personalized survival curves based on their genetic risk profile compared to population averages.

## Current Architecture Analysis

### Backend Components Status

#### ✅ **Existing Infrastructure**
- **`survival_analyzer.py`** - Complete survival analysis node with:
  - Hazard ratio calculations for 15+ key cancer genes
  - Pathway-based survival impact assessment
  - Survival curve generation using exponential models
  - Median survival calculations
  - Prognostic factor analysis
  - Clinical recommendations generation

#### ⚠️ **Integration Issues**
- **`formatter.py`** - Survival analysis data not being properly added to structured_json
- **`enhanced_api_server.py`** - Survival analyzer node not being invoked in pipeline
- **`graph.py`** - Survival analyzer not included in main pipeline flow

#### ✅ **Frontend Components**
- **`ClinicalViewPage.tsx`** - Already expects survival analysis data structure
- **Survival Analysis Tab** - UI components ready for data visualization
- **Data Structure Types** - TypeScript interfaces defined

## Data Flow Architecture

```
┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
│   Risk Assessment   │───▶│  Survival Analyzer  │───▶│   Data Formatter    │
│  (Risk Scores +     │    │  (Hazard Ratios +   │    │  (Structured JSON)  │
│   Pathway Burden)   │    │   Survival Curves)  │    │                     │
└─────────────────────┘    └─────────────────────┘    └─────────────────────┘
                                                                │
                                                                ▼
┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
│   Clinical View     │◀───│      API Server     │◀───│   Report Writer     │
│  (Survival Charts)  │    │  (Results Endpoint) │    │  (Final Output)     │
└─────────────────────┘    └─────────────────────┘    └─────────────────────┘
```

## Implementation Plan

### Phase 1: Backend Integration (Priority: High)

#### 1.1 Update Pipeline Graph
**File**: `geneknow_pipeline/graph.py`
**Changes**:
- Add survival_analyzer node to pipeline flow
- Insert after pathway_burden and before formatter
- Add conditional logic to run only when risk scores are available

```python
# Add survival analysis after pathway analysis
if state.get("risk_scores"):
    state = survival_analyzer.process(state)
```

#### 1.2 Fix Data Formatter
**File**: `geneknow_pipeline/nodes/formatter.py`
**Changes**:
- Ensure survival_analysis data is properly added to structured_json
- Add survival data validation
- Include survival curves in API response

```python
# Add to structured_json
"survival_analysis": state.get("survival_analysis", {})
```

#### 1.3 Update API Server
**File**: `geneknow_pipeline/enhanced_api_server.py`
**Changes**:
- Ensure survival analysis results are included in /api/results endpoint
- Add survival data to job status responses
- Include survival analysis in pipeline-info endpoint

### Phase 2: Frontend Visualization (Priority: High)

#### 2.1 Complete Survival Chart Component
**File**: `desktop/ui/src/pages/ClinicalViewPage.tsx`
**Changes**:
- Implement actual survival curve visualization using Chart.js or D3.js
- Add interactive features (hover, zoom, time range selection)
- Display median survival comparisons
- Show confidence intervals

#### 2.2 Add Chart Dependencies
**File**: `desktop/ui/package.json`
**Changes**:
- Add Chart.js or D3.js for survival curve visualization
- Add chart utilities and plugins

```json
{
  "dependencies": {
    "chart.js": "^4.4.0",
    "react-chartjs-2": "^5.2.0"
  }
}
```

#### 2.3 Survival Analysis Data Processing
**File**: `desktop/ui/src/pages/ClinicalViewPage.tsx`
**Changes**:
- Process survival analysis data from backend
- Transform data for chart visualization
- Handle missing data gracefully
- Add loading states and error handling

### Phase 3: Enhanced Features (Priority: Medium)

#### 3.1 Interactive Survival Curves
**Features**:
- Hover tooltips showing exact survival percentages
- Time range selection (1-year, 5-year, 10-year views)
- Confidence interval display
- Risk factor annotations

#### 3.2 Survival Statistics Display
**Components**:
- Median survival comparison table
- Life years gained/lost calculations
- Hazard ratio explanations
- Prognostic factor breakdown

#### 3.3 Clinical Recommendations
**Features**:
- Personalized screening recommendations
- Risk mitigation strategies
- Survival improvement options
- Treatment response predictions

### Phase 4: Data Validation & Testing (Priority: High)

#### 4.1 Backend Testing
**File**: `geneknow_pipeline/test_survival_analysis.py`
**Tests**:
- Survival curve generation accuracy
- Hazard ratio calculations
- Edge case handling (no mutations, multiple high-risk variants)
- Performance testing with large datasets

#### 4.2 Frontend Testing
**File**: `desktop/ui/src/components/SurvivalChart.test.tsx`
**Tests**:
- Chart rendering with real data
- Interactive features (hover, zoom)
- Responsive design
- Error state handling

#### 4.3 Integration Testing
**Tests**:
- End-to-end survival analysis pipeline
- API response validation
- Data transformation accuracy
- Cross-browser compatibility

## Technical Implementation Details

### Backend Data Structure

```python
{
  "survival_analysis": {
    "survival_curves": {
      "breast": {
        "time_points": [0, 0.2, 0.4, ...],  # Years
        "population_survival": [1.0, 0.95, 0.89, ...],  # Percentages
        "patient_survival": [1.0, 0.91, 0.83, ...],  # Percentages
        "confidence_interval": {
          "lower": [1.0, 0.88, 0.79, ...],
          "upper": [1.0, 0.94, 0.87, ...]
        },
        "median_survival": {
          "population": 18.5,  # Years
          "patient": 14.2,     # Years
          "difference": 4.3    # Years lost
        },
        "hazard_ratio": 1.7,
        "five_year_survival": {
          "population": 89.0,  # Percent
          "patient": 76.5      # Percent
        }
      }
    },
    "prognostic_factors": {
      "positive_prognostic_factors": [
        {
          "gene": "ERBB2",
          "hazard_ratio": 0.7,
          "impact": "positive",
          "confidence": 0.9
        }
      ],
      "negative_prognostic_factors": [
        {
          "gene": "TP53",
          "hazard_ratio": 2.1,
          "impact": "negative",
          "confidence": 0.95
        }
      ]
    },
    "clinical_interpretation": {
      "mutations_analyzed": 5,
      "pathways_analyzed": 3,
      "recommendation": "Enhanced surveillance recommended due to 1.7x increased risk"
    }
  }
}
```

### Frontend Chart Configuration

```typescript
interface SurvivalChartData {
  labels: number[];  // Time points
  datasets: [
    {
      label: 'Population Average';
      data: number[];
      borderColor: '#6B7280';
      backgroundColor: 'rgba(107, 114, 128, 0.1)';
    },
    {
      label: 'Your Profile';
      data: number[];
      borderColor: '#EF4444';
      backgroundColor: 'rgba(239, 68, 68, 0.1)';
    }
  ];
}
```

## Testing Strategy

### 1. Unit Tests
- **Backend**: Test survival curve calculations with known inputs
- **Frontend**: Test chart rendering with mock data
- **Integration**: Test API data flow

### 2. Integration Tests
- **End-to-end**: Complete pipeline with survival analysis
- **API**: Verify survival data in response
- **UI**: Test chart interactivity

### 3. User Acceptance Tests
- **Accuracy**: Verify survival predictions make clinical sense
- **Usability**: Test chart readability and interaction
- **Performance**: Ensure fast rendering with large datasets

## Deployment Strategy

### Phase 1 Deployment (Backend)
1. Deploy survival analyzer integration
2. Update API endpoints
3. Test with existing frontend (should show "data available")

### Phase 2 Deployment (Frontend)
1. Deploy chart visualization components
2. Update clinical view page
3. Test with real backend data

### Phase 3 Deployment (Polish)
1. Add interactive features
2. Deploy enhanced statistics
3. Add clinical recommendations

## Risk Assessment

### Technical Risks
- **Chart Performance**: Large datasets may cause rendering issues
- **Data Accuracy**: Survival predictions must be clinically sound
- **Browser Compatibility**: Chart libraries may have compatibility issues

### Mitigation Strategies
- **Performance**: Implement data sampling and virtualization
- **Accuracy**: Validate against published survival statistics
- **Compatibility**: Test across major browsers and devices

## Success Metrics

### Technical Metrics
- **Backend**: Survival analysis completes in <2 seconds
- **Frontend**: Charts render in <1 second
- **Accuracy**: Survival predictions within 5% of published data

### User Experience Metrics
- **Engagement**: Users spend >30 seconds on survival analysis tab
- **Comprehension**: Users understand survival impact (measured via feedback)
- **Clinical Value**: Healthcare providers find recommendations actionable

## Timeline

### Week 1-2: Backend Integration
- Update pipeline graph
- Fix formatter
- Update API server
- Basic testing

### Week 3-4: Frontend Visualization
- Implement chart components
- Add data processing
- Basic interactivity
- Integration testing

### Week 5-6: Polish & Testing
- Enhanced features
- Comprehensive testing
- Performance optimization
- Documentation

### Week 7: Deployment & Monitoring
- Production deployment
- User feedback collection
- Performance monitoring
- Bug fixes

## Conclusion

The survival analysis visualization feature will significantly enhance GeneKnow's clinical value by providing patients with personalized survival insights based on their genetic profile. The implementation leverages existing backend infrastructure while adding comprehensive frontend visualization capabilities.

The modular approach ensures that the feature can be developed incrementally, with each phase providing increasing value to users. The robust testing strategy ensures clinical accuracy and user experience quality.

This feature positions GeneKnow as a comprehensive genomic risk assessment platform that not only identifies risks but also provides actionable survival insights for patients and healthcare providers. 