# API Testing Quick Reference - PowerShell

## Prerequisites

Start the FastAPI server in a separate terminal:

```powershell
python -m uvicorn app.main:app --reload --port 8000
```

## Test Commands (Copy & Paste)

### 1. Health Check

```powershell
Invoke-WebRequest -Uri "http://localhost:8000/health" | Select-Object StatusCode, Content
```

### 2. Rule-Based Recommendation (/match)

**Correct field name: `reaction` (not `reaction_smiles`)**

**âš ï¸ Important: You must specify the `db` path to an existing database file.**

#### Basic Usage (No Catalyst Filter)

```powershell
$body = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    db = "data/rule_db/suzuki_db.json"
    include_trace = $true
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/match" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

#### With Catalyst Filter (NEW!)

Filter recommendations by catalyst type (Pd, Cu, Ni, etc.):

```powershell
# Filter for Copper catalysts only
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    db = "data/rule_db/C_N_Coupling_Cu_db.json"
    include_trace = $true
    relax = @{
        catalyst_class = "Cu"
    }
} | ConvertTo-Json -Depth 5

Invoke-RestMethod -Uri "http://localhost:8000/match" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

```powershell
# Filter for Palladium catalysts only
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    db = "data/rule_db/C_N_Coupling_Pd_db.json"
    relax = @{
        catalyst_class = "Pd"
    }
} | ConvertTo-Json -Depth 5

Invoke-RestMethod -Uri "http://localhost:8000/match" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

**Available databases:**

- `data/rule_db/suzuki_db.json` - Suzuki coupling
- `data/rule_db/C_N_Coupling_Pd_db.json` - C-N coupling (Pd)
- `data/rule_db/C_N_Coupling_Ni_db.json` - C-N coupling (Ni)  
- `data/rule_db/C_N_Coupling_Cu_db.json` - C-N coupling (Cu)
- `data/rule_db/amide_formation_db.json` - Amide formation

**Parameters:**

- `reaction`: Reaction SMILES (required)
- `db`: Database file path (default: "cn_coupling_pd_db.json" - may not exist)
- `include_trace`: Include matching trace (default: true)
- `relax`: Optional filters (NEW!)
  - `catalyst_class`: Filter by catalyst type ("Pd", "Cu", "Ni", "Ru", "Co", etc.)

**Supported Catalyst Values:**

- `"Pd"` - Palladium catalysts
- `"Cu"` - Copper catalysts
- `"Ni"` - Nickel catalysts
- `"Ru"` - Ruthenium catalysts
- `"Co"` - Cobalt catalysts
- `"other"` - Other catalyst types
- Omit or set to `null` for all catalysts

### 3. ML-Based Recommendation (/api/v1/recommend/conditions)

The ML-based recommendation endpoint uses machine learning to find similar precedent reactions and recommend conditions based on structural similarity (DRFP fingerprints).

#### Basic Usage (No Catalyst Filter)

```powershell
$body = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    reaction_type = "Suzuki"
    k = 50
    limit = 3
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

**Expected Response Structure:**

```json
{
  "status": "success",
  "reaction_smiles": "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1",
  "recommendations": [
    {
      "rank": 1,
      "confidence": 0.85,
      "reagents": [
        {
          "id": "7440-50-8",
          "role": "metal_precursor",
          "name": "Pd(PPh3)4",
          "abbreviation": null,
          "cas": "7440-50-8",
          "smiles": "[Pd]",
          "equivalents": null
        },
        {
          "id": "1310-73-2",
          "role": "base",
          "name": "Sodium hydroxide",
          "cas": "1310-73-2",
          "smiles": "[OH-].[Na+]",
          "equivalents": null
        },
        {
          "id": "67-56-1",
          "role": "solvent",
          "name": "Methanol",
          "cas": "67-56-1",
          "smiles": "CO",
          "equivalents": null
        }
      ],
      "conditions": {
        "temperature": {
          "value": 80,
          "unit": "Â°C"
        },
        "time": {
          "value": 12,
          "unit": "hours"
        }
      },
      "precedent_count": 15
    },
    {
      "rank": 2,
      "confidence": 0.78,
      "reagents": [ /* ... */ ],
      "conditions": { /* ... */ },
      "precedent_count": 12
    }
  ],
  "meta": {
    "model": "ML-precedent-knn",
    "reaction_type": "Suzuki",
    "detection_source": "user_supplied",
    "detection_confidence": 1.0,
    "k_precedents_retrieved": 50,
    "total_recommendations": 3,
    "processing_time_ms": 245.67
  }
}
```

**Response Field Descriptions:**

- `status`: "success" or "error"
- `reaction_smiles`: The input reaction (normalized)
- `recommendations`: Array of recommended conditions (sorted by rank)
  - `rank`: Position in recommendations (1 = best)
  - `confidence`: Confidence score (0.0-1.0)
  - `reagents`: List of chemicals needed
    - `role`: "metal_precursor", "ligand", "base", "solvent", "starting_material"
    - `name`: Chemical name
    - `cas`: CAS registry number
    - `smiles`: Chemical structure
  - `conditions`: Reaction conditions
    - `temperature`: Value and unit
    - `time`: Duration and unit
  - `precedent_count`: Number of similar precedents supporting this recommendation
- `meta`: Metadata about the recommendation
  - `reaction_type`: Detected or specified reaction type
  - `detection_source`: "user_supplied", "auto", or "router_fallback"
  - `detection_confidence`: How confident the detection is (0.0-1.0)
  - `k_precedents_retrieved`: Number of similar reactions found
  - `processing_time_ms`: Time taken to generate recommendations

#### With Catalyst Filter (NEW!)

Filter recommendations by catalyst type (Pd, Cu, Ni, etc.):

```powershell
# Filter for Palladium catalysts only
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "C_N_Coupling"
    k = 50
    limit = 5
    relax = @{
        catalyst_class = "Pd"
    }
} | ConvertTo-Json -Depth 5

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

**Expected Response with Catalyst Filter:**

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
          "unit": "Â°C"
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

**Note:** With catalyst filtering:

- Only catalysts matching the specified type are returned
- `meta.catalyst_filter` shows the active filter
- Precedent count may be lower due to filtering
- All recommended conditions will use Pd-based catalysts

```powershell
# Filter for Copper catalysts only
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "C_N_Coupling"
    k = 50
    limit = 5
    relax = @{
        catalyst_class = "Cu"
    }
} | ConvertTo-Json -Depth 5

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

**Expected Response with Cu Filter:**

```json
{
  "status": "success",
  "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "recommendations": [
    {
      "rank": 1,
      "confidence": 0.88,
      "reagents": [
        {
          "id": "7758-89-6",
          "role": "metal_precursor",
          "name": "CuI",
          "abbreviation": null,
          "cas": "7758-89-6",
          "smiles": "[Cu+]",
          "equivalents": null
        },
        {
          "id": "2404-44-6",
          "role": "ligand",
          "name": "L-Proline",
          "abbreviation": "Pro",
          "cas": "2404-44-6",
          "smiles": "O=C(O)C1CCCN1",
          "equivalents": null
        },
        {
          "id": "1310-73-2",
          "role": "base",
          "name": "K3PO4",
          "cas": "1310-73-2",
          "smiles": "[K+].[K+].[K+].[O-]P([O-])([O-])=O",
          "equivalents": null
        },
        {
          "id": "67-68-5",
          "role": "solvent",
          "name": "DMSO",
          "cas": "67-68-5",
          "smiles": "CS(C)=O",
          "equivalents": null
        }
      ],
      "conditions": {
        "temperature": {
          "value": 90,
          "unit": "Â°C"
        },
        "time": {
          "value": 24,
          "unit": "hours"
        }
      },
      "precedent_count": 18
    }
  ],
  "meta": {
    "model": "ML-precedent-knn",
    "reaction_type": "C_N_Coupling",
    "detection_source": "user_supplied",
    "detection_confidence": 1.0,
    "k_precedents_retrieved": 50,
    "total_recommendations": 5,
    "catalyst_filter": "Cu",
    "processing_time_ms": 298.12
  }
}
```

**Comparing Pd vs Cu Results:**

- **Pd**: Higher confidence (0.92), shorter time (16h), moderate temp (100Â°C), ligand: XPhos
- **Cu**: Lower confidence (0.88), longer time (24h), lower temp (90Â°C), ligand: L-Proline
- **Pd**: Typically more selective, faster reactions
- **Cu**: More economical, different selectivity profile

#### With Multiple Filters

```powershell
# Combine catalyst filter with reranking strategy
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "C_N_Coupling"
    k = 50
    limit = 5
    relax = @{
        catalyst_class = "Ni"
    }
    rerank_strategy = "rule"  # Options: "rule", "analytics", "none"
    filter_unknown_reagents = $false
} | ConvertTo-Json -Depth 5

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

**Parameters:**

- `reaction`: Reaction SMILES (required)
- `reaction_type`: Reaction family (optional, auto-detects if omitted)
- `k`: Number of precedents to retrieve (default: 50)
- `limit`: Number of recommendations to return (default: 5)
- `relax`: Optional filters (NEW!)
  - `catalyst_class`: Filter by catalyst type ("Pd", "Cu", "Ni", etc.)
- `rerank_strategy`: Reranking method (default: "rule")
  - `"rule"`: Boost by chemical rule matching
  - `"analytics"`: Boost by reagent popularity
  - `"none"`: DRFP similarity only
- `filter_unknown_reagents`: Exclude precedents with unknown reagents (default: false)

**Supported Reaction Types:**

- `"Suzuki"` - Suzuki-Miyaura coupling
- `"C_N_Coupling"` - C-N coupling (auto-selects Pd/Cu/Ni database)
- `"Amide_formation"` - Amide formation
- Omit or set to `null` for auto-detection

#### Error Responses and Edge Cases

**Error: Invalid Reaction SMILES**

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

**Error: No Precedents Found**

```json
{
  "status": "success",
  "reaction_smiles": "...",
  "recommendations": [],
  "meta": {
    "model": "ML-precedent-knn",
    "reaction_type": "Unknown",
    "detection_source": "router_fallback",
    "detection_confidence": 0.0,
    "k_precedents_retrieved": 0,
    "total_recommendations": 0,
    "warning": "No similar precedents found for this reaction",
    "processing_time_ms": 123.45
  }
}
```

**Warning: Low Confidence Detection**

```json
{
  "status": "success",
  "reaction_smiles": "...",
  "recommendations": [ /* ... */ ],
  "meta": {
    "model": "ML-precedent-knn",
    "reaction_type": "C_N_Coupling",
    "detection_source": "auto",
    "detection_confidence": 0.45,
    "warning": "Low confidence auto-detection. Consider manually specifying reaction_type.",
    "k_precedents_retrieved": 50,
    "total_recommendations": 3,
    "processing_time_ms": 234.56
  }
}
```

**Success: Empty Results with Catalyst Filter**

```json
{
  "status": "success",
  "reaction_smiles": "...",
  "recommendations": [],
  "meta": {
    "model": "ML-precedent-knn",
    "reaction_type": "C_N_Coupling",
    "detection_source": "user_supplied",
    "detection_confidence": 1.0,
    "catalyst_filter": "Ni",
    "k_precedents_retrieved": 0,
    "total_recommendations": 0,
    "info": "No Ni-catalyzed precedents found. Try 'Pd' or 'Cu' instead.",
    "processing_time_ms": 156.78
  }
}
```

### 4. Fusion Recommendation (/api/v1/recommend/fusion) - Deprecated

```powershell
$body = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    k = 50
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/fusion" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

### 5. Reaction Type Detection (/api/v1/reaction/detect-type)

Test the auto-detection capability separately:

```powershell
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/reaction/detect-type" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

**Response shows:**

- `selected_family`: The detected reaction type
- `rxn_insight`: Results from ML-based detection (if available)
- `router_fallback`: Results from rule-based detection
- `rxn_insight_available`: Whether ML detector is installed

### 6. Auto-Detection with Recommendations

#### ML Recommendations with Auto-Detection

Omit `reaction_type` to let the API auto-detect:

```powershell
# Auto-detect reaction type, no catalyst filter
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    # reaction_type is omitted - will auto-detect!
    k = 50
    limit = 5
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

#### Auto-Detection + Catalyst Filter

Combine auto-detection with catalyst preference:

```powershell
# Auto-detect type, but prefer Palladium catalysts
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    # reaction_type omitted for auto-detection
    k = 50
    limit = 5
    relax = @{
        catalyst_class = "Pd"
    }
    rerank_strategy = "rule"
} | ConvertTo-Json -Depth 5

$result = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $body -ContentType "application/json"

# Show what was detected
Write-Host "Detected Type: $($result.meta.reaction_type)" -ForegroundColor Cyan
Write-Host "Detection Source: $($result.meta.detection_source)" -ForegroundColor Yellow
Write-Host "Confidence: $($result.meta.detection_confidence)" -ForegroundColor Green

$result | ConvertTo-Json -Depth 10
```

#### Test Different Reaction Types with Auto-Detection

```powershell
# Test 1: C-N Coupling (should detect automatically)
Write-Host "`n=== Testing C-N Coupling Auto-Detection ===" -ForegroundColor Cyan
$cn_body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    k = 50
    limit = 3
} | ConvertTo-Json

$cn_result = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $cn_body -ContentType "application/json"
Write-Host "Detected as: $($cn_result.meta.reaction_type)" -ForegroundColor Green

# Test 2: Suzuki Coupling (should detect automatically)
Write-Host "`n=== Testing Suzuki Coupling Auto-Detection ===" -ForegroundColor Cyan
$suzuki_body = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    k = 50
    limit = 3
} | ConvertTo-Json

$suzuki_result = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $suzuki_body -ContentType "application/json"
Write-Host "Detected as: $($suzuki_result.meta.reaction_type)" -ForegroundColor Green

# Test 3: Amide Formation (should detect automatically)
Write-Host "`n=== Testing Amide Formation Auto-Detection ===" -ForegroundColor Cyan
$amide_body = @{
    reaction = "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1"
    k = 50
    limit = 3
} | ConvertTo-Json

$amide_result = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $amide_body -ContentType "application/json"
Write-Host "Detected as: $($amide_result.meta.reaction_type)" -ForegroundColor Green
```

### 7. Protocol Recommendation (No API endpoint - uses CLI)

**Protocol doesn't have an API endpoint yet. Use the CLI tool:**

```powershell
python scripts/web_recommendation_cli.py `
  --rxn "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1" `
  --strategy protocol `
  --k 3
```

## API Endpoint Summary

| Endpoint | Method | Key Fields | Catalyst Filter | Description |
|----------|--------|------------|-----------------|-------------|
| `/health` | GET | - | N/A | Health check |
| `/match` | POST | `reaction`, `db` | âœ?`relax.catalyst_class` | Rule-based (SCDB) |
| `/api/v1/recommend/conditions` | POST | `reaction`, `reaction_type` | âœ?`relax.catalyst_class` | ML-based (precedents) |
| `/api/v1/recommend` | POST | `reaction` | âœ?`relax.catalyst_class` | Full recommendation |
| `/api/v1/recommend/fusion` | POST | `reaction` | â?Deprecated | Fusion (deprecated) |
| `/api/v1/reaction/detect-type` | POST | `reaction` | N/A | Type detection |

## Catalyst Filtering Examples

### Example 1: Compare Different Catalysts for C-N Coupling

```powershell
# Test with Palladium
$bodyPd = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "C_N_Coupling"
    k = 50
    limit = 3
    relax = @{ catalyst_class = "Pd" }
} | ConvertTo-Json -Depth 5

$resultPd = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $bodyPd -ContentType "application/json"

Write-Host "=== Palladium Results ===" -ForegroundColor Cyan
$resultPd | ConvertTo-Json -Depth 10

# Test with Copper
$bodyCu = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "C_N_Coupling"
    k = 50
    limit = 3
    relax = @{ catalyst_class = "Cu" }
} | ConvertTo-Json -Depth 5

$resultCu = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $bodyCu -ContentType "application/json"

Write-Host "=== Copper Results ===" -ForegroundColor Green
$resultCu | ConvertTo-Json -Depth 10

# Test with Nickel
$bodyNi = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "C_N_Coupling"
    k = 50
    limit = 3
    relax = @{ catalyst_class = "Ni" }
} | ConvertTo-Json -Depth 5

$resultNi = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $bodyNi -ContentType "application/json"

Write-Host "=== Nickel Results ===" -ForegroundColor Yellow
$resultNi | ConvertTo-Json -Depth 10
```

### Example 2: Rule-Based with Catalyst Filter

```powershell
# Get Cu-specific rule-based recommendations
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    db = "data/rule_db/C_N_Coupling_Cu_db.json"
    include_trace = $true
    relax = @{ catalyst_class = "Cu" }
} | ConvertTo-Json -Depth 5

$result = Invoke-RestMethod -Uri "http://localhost:8000/match" `
    -Method POST -Body $body -ContentType "application/json"

# Check filter metadata
Write-Host "Original count: $($result._catalyst_filtered.original_count)" -ForegroundColor Yellow
Write-Host "Filtered count: $($result._catalyst_filtered.filtered_count)" -ForegroundColor Green
$result | ConvertTo-Json -Depth 10
```

### Example 3: Auto-Detection + Catalyst Filter

```powershell
# Let the API auto-detect reaction type, but specify catalyst preference
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    # reaction_type omitted - will auto-detect
    k = 50
    limit = 5
    relax = @{ catalyst_class = "Pd" }
    rerank_strategy = "rule"
} | ConvertTo-Json -Depth 5

$result = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
    -Method POST -Body $body -ContentType "application/json"

Write-Host "Detected type: $($result.meta.reaction_type)" -ForegroundColor Cyan
Write-Host "Detection source: $($result.meta.detection_source)" -ForegroundColor Gray
$result | ConvertTo-Json -Depth 10
```

## Common Mistakes

### â?Wrong Field for Catalyst Filtering

**CRITICAL: Use `relax.catalyst_class`, NOT `constraints.metal_preference`**

```powershell
# â?WRONG - Using constraints (will be IGNORED!)
$body = @{
    reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    reaction_type = "Suzuki"
    k = 10
    limit = 5
    constraints = @{
        metal_preference = "Cu"  # â?WRONG FIELD - DOES NOT WORK!
    }
    relax = @{}
} | ConvertTo-Json -Depth 5
```

```powershell
# âœ?CORRECT - Using relax.catalyst_class
$body = @{
    reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    reaction_type = "Suzuki"
    k = 10
    limit = 5
    relax = @{
        catalyst_class = "Cu"  # âœ?CORRECT - THIS WORKS!
    }
    constraints = @{}
} | ConvertTo-Json -Depth 5
```

**Why this matters:**

- `relax` parameters control the **precedent search** (filtering BEFORE recommendation)
- `constraints` parameters filter the **final recommendations** (filtering AFTER recommendation)
- Catalyst filtering MUST happen during search, so use `relax.catalyst_class`

### â?Wrong Field Names

```powershell
# WRONG for /match endpoint
$body = @{
    reaction_smiles = "..."  # Should be "reaction"
    reaction_type = "Suzuki"  # Not used in /match
} | ConvertTo-Json
```

```powershell
# CORRECT for /match endpoint
$body = @{
    reaction = "..."  # Correct!
    db = "data/rule_db/suzuki_db.json"
} | ConvertTo-Json
```

### â?Forgetting -Depth for Nested Objects

```powershell
# WRONG - catalyst_class won't be included
$body = @{
    reaction = "..."
    relax = @{ catalyst_class = "Pd" }
} | ConvertTo-Json  # Missing -Depth!

# CORRECT - nested objects properly serialized
$body = @{
    reaction = "..."
    relax = @{ catalyst_class = "Pd" }
} | ConvertTo-Json -Depth 5
```

### â?Using Wrong Database with Catalyst Filter

```powershell
# WRONG - Asking for Cu catalysts from Pd database
$body = @{
    reaction = "..."
    db = "data/rule_db/C_N_Coupling_Pd_db.json"  # Pd database
    relax = @{ catalyst_class = "Cu" }  # Cu filter - mismatch!
} | ConvertTo-Json -Depth 5

# CORRECT - Match database to catalyst
$body = @{
    reaction = "..."
    db = "data/rule_db/C_N_Coupling_Cu_db.json"  # Cu database
    relax = @{ catalyst_class = "Cu" }  # Cu filter - match!
} | ConvertTo-Json -Depth 5
```

## Quick Reference Card

### Auto-Detection Usage

```powershell
# Option 1: Omit reaction_type completely
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    # No reaction_type field - will auto-detect
    k = 50
    limit = 5
} | ConvertTo-Json

# Option 2: Set reaction_type to null
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = $null  # Explicit auto-detect
    k = 50
    limit = 5
} | ConvertTo-Json

# Option 3: Auto-detect + catalyst filter
$body = @{
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    # reaction_type omitted
    k = 50
    limit = 5
    relax = @{ catalyst_class = "Pd" }
} | ConvertTo-Json -Depth 5
```

### Catalyst Filter Syntax

```powershell
# Add this to any request that supports filtering:
relax = @{
    catalyst_class = "Pd"  # or "Cu", "Ni", "Ru", "Co", etc.
}
```

### Common Catalyst Values

| Value | Description | Example Reactions |
|-------|-------------|-------------------|
| `"Pd"` | Palladium | Suzuki, Buchwald-Hartwig |
| `"Cu"` | Copper | Ullmann, Chan-Lam |
| `"Ni"` | Nickel | Negishi, Kumada |
| `"Ru"` | Ruthenium | Metathesis, Hydrogenation |
| `"Co"` | Cobalt | C-H activation |

### Reaction Type Values

| Value | Description | Databases |
|-------|-------------|-----------|
| `"Suzuki"` | Suzuki-Miyaura C-C coupling | `suzuki_db.json` |
| `"C_N_Coupling"` | C-N coupling (all catalysts) | `C_N_Coupling_Pd_db.json`, `C_N_Coupling_Cu_db.json`, `C_N_Coupling_Ni_db.json` |
| `"Amide_formation"` | Amide bond formation | `amide_formation_db.json` |
| `null` / omitted | **Auto-detect** | - |

### Auto-Detection Behavior

| Scenario | Detection Method | Notes |
|----------|------------------|-------|
| `reaction_type` omitted | ML-based (rxn-insight) â†?Rule-based fallback | Recommended approach |
| `reaction_type = null` | Same as omitted | Explicit auto-detect |
| `reaction_type = "Suzuki"` | User-specified, no detection | Manual override |
| Auto-detect + catalyst filter | Detects type, filters by catalyst | Best for exploration |

**Detection Priority:**

1. User-specified `reaction_type` (if provided)
2. ML-based detection via `rxn-insight` (if available)
3. Rule-based detection via `chemtools.router`
4. Fallback to "Unknown"

### Rerank Strategies

| Value | Description | Best For |
|-------|-------------|----------|
| `"rule"` | Boost by chemical rule matching | Default - balanced quality |
| `"analytics"` | Boost by reagent popularity | Popular, well-tested conditions |
| `"none"` | DRFP similarity only | Novel chemistry exploration |

## View OpenAPI Documentation

Open in browser: <http://localhost:8000/docs>

This interactive documentation shows all endpoints, request/response schemas, and allows testing directly in the browser with the new catalyst filtering options.
