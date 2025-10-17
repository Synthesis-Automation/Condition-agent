# ML-Based Recommendation API - Complete Guide

## Overview

The ML-based recommendation endpoint (`/api/v1/recommend/conditions`) uses machine learning to find structurally similar precedent reactions and recommend optimal reaction conditions based on:

- **DRFP Fingerprints**: Structural similarity between reactions
- **Precedent Database**: Historical reaction data with known outcomes
- **Catalyst Filtering**: Optional filtering by metal type (Pd, Cu, Ni, etc.)
- **Reranking**: Boost results by chemical rules or reagent popularity

## Endpoint Details

**URL:** `POST http://localhost:8000/api/v1/recommend/conditions`  
**Content-Type:** `application/json`

## Request Format

### Basic Request
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "reaction_type": "C_N_Coupling",
  "k": 50,
  "limit": 5
}
```

### Request with Catalyst Filter
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "reaction_type": "C_N_Coupling",
  "k": 50,
  "limit": 5,
  "relax": {
    "catalyst_class": "Pd"
  }
}
```

### Request with Auto-Detection
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "k": 50,
  "limit": 5,
  "relax": {
    "catalyst_class": "Cu"
  }
}
```

## Request Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `reaction` | string | ✅ Yes | - | Reaction SMILES (reactants>>products) |
| `reaction_type` | string | ❌ No | null | Reaction family ("Suzuki", "C_N_Coupling", etc.) - omit for auto-detection |
| `k` | integer | ❌ No | 50 | Number of similar precedents to retrieve |
| `limit` | integer | ❌ No | 5 | Number of recommendations to return |
| `relax` | object | ❌ No | {} | Filtering options |
| `relax.catalyst_class` | string | ❌ No | null | Catalyst filter ("Pd", "Cu", "Ni", "Ru", "Co", etc.) |
| `rerank_strategy` | string | ❌ No | "rule" | Reranking method ("rule", "analytics", "none") |
| `filter_unknown_reagents` | boolean | ❌ No | false | Exclude precedents with unknown reagents |

## Response Format

### Successful Response Structure

```json
{
  "status": "success",
  "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "recommendations": [
    {
      "rank": 1,
      "confidence": 0.92,
      "reagents": [
        {
          "id": "14221-01-3",
          "role": "metal_precursor",
          "name": "Pd(OAc)2",
          "abbreviation": null,
          "cas": "14221-01-3",
          "smiles": "[Pd+2]",
          "equivalents": null
        },
        {
          "id": "11099-03-9",
          "role": "ligand",
          "name": "XPhos",
          "abbreviation": "XPhos",
          "cas": "11099-03-9",
          "smiles": "CC(C)c1cc(C(C)C)c(-c2ccccc2P(C3CCCCC3)C4CCCCC4)c(C(C)C)c1",
          "equivalents": null
        },
        {
          "id": "917-69-1",
          "role": "base",
          "name": "Potassium tert-butoxide",
          "cas": "917-69-1",
          "smiles": "CC(C)(C)[O-].[K+]",
          "equivalents": null
        },
        {
          "id": "1310-58-3",
          "role": "solvent",
          "name": "Toluene",
          "cas": "1310-58-3",
          "smiles": "Cc1ccccc1",
          "equivalents": null
        }
      ],
      "conditions": {
        "temperature": {
          "value": 100,
          "unit": "°C"
        },
        "time": {
          "value": 16,
          "unit": "hours"
        }
      },
      "precedent_count": 23
    }
  ],
  "meta": {
    "model": "ML-precedent-knn",
    "reaction_type": "C_N_Coupling",
    "detection_source": "user_supplied",
    "detection_confidence": 1.0,
    "k_precedents_retrieved": 50,
    "total_recommendations": 5,
    "catalyst_filter": "Pd",
    "processing_time_ms": 312.45
  }
}
```

## Response Field Descriptions

### Top Level
- `status`: "success" or "error"
- `reaction_smiles`: The input reaction (normalized)
- `recommendations`: Array of recommended conditions (sorted by rank/confidence)
- `meta`: Metadata about the recommendation process

### Recommendation Object
- `rank`: Position in recommendations (1 = best)
- `confidence`: Confidence score (0.0-1.0)
  - > 0.9: Very high confidence
  - 0.7-0.9: High confidence
  - 0.5-0.7: Medium confidence
  - < 0.5: Low confidence
- `reagents`: List of chemicals needed
- `conditions`: Reaction conditions
- `precedent_count`: Number of similar precedents supporting this recommendation

### Reagent Object
- `id`: Unique identifier (usually CAS number)
- `role`: Chemical role in reaction
  - `"metal_precursor"`: Metal catalyst source
  - `"ligand"`: Catalyst ligand
  - `"base"`: Base reagent
  - `"solvent"`: Reaction solvent
  - `"starting_material"`: Reactant
  - `"additive"`: Additional reagent
- `name`: Chemical name
- `abbreviation`: Short name (if available)
- `cas`: CAS registry number
- `smiles`: Chemical structure (SMILES)
- `equivalents`: Molar equivalents (may be null)

### Conditions Object
- `temperature`: Temperature value and unit
- `time`: Duration value and unit
- `atmosphere`: Atmospheric conditions (may be absent)

### Meta Object
- `model`: Recommendation model used ("ML-precedent-knn")
- `reaction_type`: Detected or specified reaction type
- `detection_source`: How the type was determined
  - `"user_supplied"`: User specified reaction_type
  - `"auto"`: ML-based auto-detection (rxn-insight)
  - `"router_fallback"`: Rule-based detection
- `detection_confidence`: Confidence in type detection (0.0-1.0)
- `k_precedents_retrieved`: Number of similar reactions found
- `total_recommendations`: Number of recommendations returned
- `catalyst_filter`: Active catalyst filter (if any)
- `processing_time_ms`: Time taken to generate recommendations
- `warning`: Warning message (if any)
- `info`: Informational message (if any)

## PowerShell Testing Examples

### Example 1: Basic ML Recommendation
```powershell
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "C_N_Coupling"
    k = 50
    limit = 5
} | ConvertTo-Json

$result = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $body -ContentType "application/json"

# Display results
Write-Host "Status: $($result.status)" -ForegroundColor Green
Write-Host "Reaction Type: $($result.meta.reaction_type)" -ForegroundColor Cyan
Write-Host "Recommendations: $($result.recommendations.Count)" -ForegroundColor Yellow

foreach ($rec in $result.recommendations) {
    Write-Host "`nRank $($rec.rank): Confidence $($rec.confidence)" -ForegroundColor Magenta
    Write-Host "  Catalyst: $($rec.reagents | Where-Object {$_.role -eq 'metal_precursor'} | Select-Object -ExpandProperty name)"
    Write-Host "  Base: $($rec.reagents | Where-Object {$_.role -eq 'base'} | Select-Object -ExpandProperty name)"
    Write-Host "  Solvent: $($rec.reagents | Where-Object {$_.role -eq 'solvent'} | Select-Object -ExpandProperty name)"
    Write-Host "  Temp: $($rec.conditions.temperature.value)$($rec.conditions.temperature.unit)"
    Write-Host "  Time: $($rec.conditions.time.value) $($rec.conditions.time.unit)"
}
```

### Example 2: Catalyst Filtering
```powershell
# Test multiple catalysts
$catalysts = @("Pd", "Cu", "Ni")
$reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

foreach ($cat in $catalysts) {
    Write-Host "`n=== Testing $cat Catalysts ===" -ForegroundColor Cyan
    
    $body = @{
        reaction = $reaction
        reaction_type = "C_N_Coupling"
        k = 50
        limit = 3
        relax = @{ catalyst_class = $cat }
    } | ConvertTo-Json -Depth 5
    
    $result = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
        -Method POST -Body $body -ContentType "application/json"
    
    if ($result.recommendations.Count -eq 0) {
        Write-Host "No $cat catalysts found" -ForegroundColor Red
    } else {
        Write-Host "Found $($result.recommendations.Count) recommendations" -ForegroundColor Green
        Write-Host "Top recommendation confidence: $($result.recommendations[0].confidence)" -ForegroundColor Yellow
        Write-Host "Precedent count: $($result.recommendations[0].precedent_count)" -ForegroundColor Gray
    }
}
```

### Example 3: Auto-Detection
```powershell
# Auto-detect reaction type
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    # No reaction_type specified - will auto-detect
    k = 50
    limit = 5
} | ConvertTo-Json

$result = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $body -ContentType "application/json"

Write-Host "Detected Type: $($result.meta.reaction_type)" -ForegroundColor Cyan
Write-Host "Detection Source: $($result.meta.detection_source)" -ForegroundColor Yellow
Write-Host "Confidence: $($result.meta.detection_confidence)" -ForegroundColor Green

if ($result.meta.detection_confidence -lt 0.5) {
    Write-Host "WARNING: Low detection confidence!" -ForegroundColor Red
}
```

### Example 4: Auto-Detection + Catalyst Filter
```powershell
# Auto-detect type but prefer specific catalyst
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    # reaction_type omitted for auto-detection
    k = 50
    limit = 5
    relax = @{ catalyst_class = "Pd" }
    rerank_strategy = "rule"
} | ConvertTo-Json -Depth 5

$result = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $body -ContentType "application/json"

Write-Host "Auto-Detected: $($result.meta.reaction_type)" -ForegroundColor Cyan
Write-Host "Catalyst Filter: $($result.meta.catalyst_filter)" -ForegroundColor Magenta
Write-Host "Results: $($result.recommendations.Count)" -ForegroundColor Green
```

## Error Scenarios

### Error: Invalid Reaction SMILES
**Request:**
```json
{
  "reaction": "invalid_smiles>>product"
}
```

**Response:**
```json
{
  "status": "error",
  "error": {
    "code": "INVALID_SMILES",
    "message": "Could not parse reaction SMILES",
    "details": "Invalid SMILES syntax in input"
  },
  "meta": {
    "processing_time_ms": 12.34
  }
}
```

### Warning: No Precedents Found
**Request:**
```json
{
  "reaction": "CC>>CCC",
  "reaction_type": "Suzuki"
}
```

**Response:**
```json
{
  "status": "success",
  "reaction_smiles": "CC>>CCC",
  "recommendations": [],
  "meta": {
    "model": "ML-precedent-knn",
    "reaction_type": "Suzuki",
    "detection_source": "user_supplied",
    "detection_confidence": 1.0,
    "k_precedents_retrieved": 0,
    "total_recommendations": 0,
    "warning": "No similar precedents found for this reaction",
    "processing_time_ms": 123.45
  }
}
```

### Warning: Low Confidence Detection
**Response:**
```json
{
  "status": "success",
  "recommendations": [ /* ... */ ],
  "meta": {
    "reaction_type": "C_N_Coupling",
    "detection_source": "auto",
    "detection_confidence": 0.45,
    "warning": "Low confidence auto-detection. Consider manually specifying reaction_type."
  }
}
```

### Info: No Results with Catalyst Filter
**Response:**
```json
{
  "status": "success",
  "recommendations": [],
  "meta": {
    "catalyst_filter": "Ni",
    "k_precedents_retrieved": 0,
    "total_recommendations": 0,
    "info": "No Ni-catalyzed precedents found. Try 'Pd' or 'Cu' instead."
  }
}
```

## Comparison: Pd vs Cu vs Ni

### Example Reaction: C-N Coupling
`Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

| Aspect | Palladium (Pd) | Copper (Cu) | Nickel (Ni) |
|--------|----------------|-------------|-------------|
| **Typical Confidence** | 0.85-0.95 | 0.75-0.90 | 0.70-0.85 |
| **Temperature** | 80-110°C | 80-120°C | 60-100°C |
| **Reaction Time** | 8-24h | 12-48h | 12-36h |
| **Common Ligands** | XPhos, BINAP, dppf | L-Proline, DMEDA | dtbbpy, bipyridine |
| **Cost** | High | Low | Medium |
| **Selectivity** | Very high | Medium-High | Medium |
| **Substrate Scope** | Broad | Narrower | Growing |
| **Best For** | Complex molecules | Cost-sensitive | Emerging area |

## Best Practices

### 1. Always Check Detection Confidence
```powershell
if ($result.meta.detection_confidence -lt 0.7) {
    Write-Host "WARNING: Consider specifying reaction_type manually" -ForegroundColor Yellow
}
```

### 2. Compare Multiple Catalysts
```powershell
# Get recommendations for all catalysts, pick the best
$results = @{}
foreach ($cat in @("Pd", "Cu", "Ni")) {
    $results[$cat] = Invoke-RestMethod ... # with catalyst_class = $cat
}

# Pick highest confidence
$best = $results.GetEnumerator() | Sort-Object {$_.Value.recommendations[0].confidence} -Descending | Select-Object -First 1
Write-Host "Best catalyst: $($best.Key) with confidence $($best.Value.recommendations[0].confidence)"
```

### 3. Use Auto-Detection for Exploration
```powershell
# Let API detect type, then refine with catalyst
$auto_result = Invoke-RestMethod ... # without reaction_type
$detected_type = $auto_result.meta.reaction_type

# Use detected type for refined search
$refined_body = @{
    reaction = $reaction
    reaction_type = $detected_type  # Use auto-detected type
    relax = @{ catalyst_class = "Pd" }
} | ConvertTo-Json -Depth 5
```

### 4. Handle Empty Results Gracefully
```powershell
if ($result.recommendations.Count -eq 0) {
    Write-Host "No recommendations found" -ForegroundColor Red
    
    # Try without catalyst filter
    $body_no_filter = @{
        reaction = $reaction
        reaction_type = $reaction_type
        # No relax parameter
    } | ConvertTo-Json
    
    $result2 = Invoke-RestMethod ...
    if ($result2.recommendations.Count -gt 0) {
        Write-Host "Found $($result2.recommendations.Count) without catalyst filter" -ForegroundColor Green
    }
}
```

## Troubleshooting

### Issue: Empty Recommendations
**Possible Causes:**
1. Catalyst filter too restrictive → Try without filter or different catalyst
2. Novel reaction type → Check with auto-detection
3. Invalid SMILES → Verify reaction syntax

### Issue: Low Confidence Scores
**Possible Causes:**
1. Unclear reaction type → Specify `reaction_type` manually
2. Novel substrate → Expected for new chemistry
3. Need more precedents → Increase `k` parameter

### Issue: Unexpected Reaction Type
**Solutions:**
1. Manually specify `reaction_type`
2. Check reaction SMILES for errors
3. Use `/api/v1/reaction/detect-type` to debug detection

## Related Endpoints

- `POST /match` - Rule-based matching (also supports catalyst filtering)
- `POST /api/v1/reaction/detect-type` - Standalone type detection
- `POST /api/v1/recommend` - Full recommendation workflow
- `GET /docs` - Interactive API documentation

## Performance Notes

- **Average Response Time:** 200-500ms
- **With Catalyst Filter:** Similar (filtering is fast)
- **Large k values:** May increase time (500+ ms)
- **Auto-Detection:** Adds ~50-100ms overhead

## Further Reading

- `CATALYST_ENDPOINT_SUPPORT.md` - Technical details on catalyst filtering
- `API_TESTING_GUIDE.md` - Full API testing guide
- `WEB_CLI_SYNC_COMPLETE.md` - CLI usage examples

---

**Last Updated:** October 17, 2025  
**API Version:** v1  
**Status:** Production Ready ✅
