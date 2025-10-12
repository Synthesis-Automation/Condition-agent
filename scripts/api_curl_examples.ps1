# Quick API Test Examples for PowerShell
# Make sure the server is running: python -m uvicorn app.main:app --reload --port 8000

# =============================================================================
# HEALTH CHECK
# =============================================================================
Invoke-WebRequest -Uri "http://localhost:8000/health" | Select-Object StatusCode, Content

# =============================================================================
# RULE-BASED RECOMMENDATION
# =============================================================================
$ruleBody = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    db = "data/conditionDB/suzuki_db.json"
    include_trace = $true
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/match" -Method POST `
    -Body $ruleBody -ContentType "application/json" | ConvertTo-Json -Depth 10

# =============================================================================
# ML-BASED RECOMMENDATION
# =============================================================================
$mlBody = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    reaction_type = "Suzuki"
    k = 50
    limit = 3
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" -Method POST `
    -Body $mlBody -ContentType "application/json" | ConvertTo-Json -Depth 10

# =============================================================================
# FUSION RECOMMENDATION
# =============================================================================
$fusionBody = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    k = 50
    max_variants = 3
} | ConvertTo-Json

Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/fusion" -Method POST `
    -Body $fusionBody -ContentType "application/json" | ConvertTo-Json -Depth 10

# =============================================================================
# PROTOCOL (Currently uses local module - no API endpoint yet)
# =============================================================================
# Use the CLI tool instead:
# python scripts/web_recommendation_cli.py --rxn "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1" --strategy protocol --k 3
