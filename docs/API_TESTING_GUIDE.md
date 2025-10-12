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

**⚠️ Important: You must specify the `db` path to an existing database file.**

```powershell
$body = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    db = "data/conditionDB/suzuki_db.json"
    include_trace = $true
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/match" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

**Available databases:**
- `data/conditionDB/suzuki_db.json` - Suzuki coupling
- `data/conditionDB/C_N_Coupling_Pd_db.json` - C-N coupling (Pd)
- `data/conditionDB/C_N_Coupling_Ni_db.json` - C-N coupling (Ni)  
- `data/conditionDB/C_N_Coupling_Cu_db.json` - C-N coupling (Cu)
- `data/conditionDB/amide_formation_db.json` - Amide formation

**Optional parameters:**
- `db`: Database file path (default: "cn_coupling_pd_db.json" - may not exist)
- `include_trace`: Include matching trace (default: true)

### 3. ML-Based Recommendation (/api/v1/recommend/conditions)
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

**Parameters:**
- `reaction`: Reaction SMILES (required)
- `reaction_type`: Reaction family (optional)
- `k`: Number of precedents (default: 50)
- `limit`: Number of recommendations (default: 5)

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
```powershell
$body = @{
    reaction_smiles = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/reaction/detect-type" -Method POST `
    -Body $body -ContentType "application/json" | ConvertTo-Json -Depth 10
```

### 6. Protocol Recommendation (No API endpoint - uses CLI)
**Protocol doesn't have an API endpoint yet. Use the CLI tool:**

```powershell
python scripts/web_recommendation_cli.py `
  --rxn "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1" `
  --strategy protocol `
  --k 3
```

## API Endpoint Summary

| Endpoint | Method | Request Field | Description |
|----------|--------|---------------|-------------|
| `/health` | GET | - | Health check |
| `/match` | POST | `reaction` | Rule-based (SCDB) |
| `/api/v1/recommend/conditions` | POST | `reaction_smiles`, `reaction_type` | ML-based (precedents) |
| `/api/v1/recommend/fusion` | POST | `reaction_smiles`, `reaction_type` | Fusion (deprecated) |
| `/api/v1/reaction/detect-type` | POST | `reaction_smiles` | Type detection |

## Common Mistakes

❌ **WRONG** (for `/match`):
```powershell
$body = @{
    reaction_smiles = "..."  # Wrong field name!
    reaction_type = "Suzuki"
} | ConvertTo-Json
```

✅ **CORRECT** (for `/match`):
```powershell
$body = @{
    reaction = "..."  # Correct field name
    include_trace = $true
} | ConvertTo-Json
```

## View OpenAPI Documentation
Open in browser: http://localhost:8000/docs

This interactive documentation shows all endpoints, request/response schemas, and allows testing directly in the browser.
