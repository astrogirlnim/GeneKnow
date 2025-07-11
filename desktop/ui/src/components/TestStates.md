# Quick Test States Reference

To test different states, change the `status` in `MOCK_SHAP_VALIDATION` in both files:

## FLAG_FOR_REVIEW (Default)
```typescript
status: 'FLAG_FOR_REVIEW' as const,
```

## PASS State
```typescript
status: 'PASS' as const,
reasons: [], // Empty for PASS
```

## ERROR State  
```typescript
status: 'ERROR' as const,
reasons: ['SHAP validation error: Model structure incompatible'],
top_contributors: [], // Empty for ERROR
```

## SKIPPED State
```typescript
status: 'SKIPPED' as const,
reasons: ['ML fusion model or features not available'],
top_contributors: [], // Empty for SKIPPED
```

## Files to Update:
- `desktop/ui/src/pages/DashboardPage.tsx` 
- `desktop/ui/src/pages/ClinicalViewPage.tsx`

## Location in Code:
Look for `MOCK_SHAP_VALIDATION` object around line 190-220 