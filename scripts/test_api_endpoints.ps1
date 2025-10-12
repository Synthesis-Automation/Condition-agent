# PowerShell Script to Test Condition-Agent API Endpoints
# Usage: .\scripts\test_api_endpoints.ps1

$baseUrl = "http://localhost:8000"

Write-Host "Testing Condition-Agent API Endpoints" -ForegroundColor Cyan
Write-Host "======================================" -ForegroundColor Cyan
Write-Host ""

# Test 1: Health Check
Write-Host "1. Health Check" -ForegroundColor Yellow
$response = Invoke-WebRequest -Uri "$baseUrl/health" -Method GET
Write-Host "Status: $($response.StatusCode)" -ForegroundColor Green
$response.Content | ConvertFrom-Json | ConvertTo-Json -Depth 10
Write-Host ""

# Test 2: Rule-Based Recommendation
Write-Host "2. Rule-Based Recommendation (/match)" -ForegroundColor Yellow
$ruleBody = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    db = "data/conditionDB/suzuki_db.json"
    include_trace = $true
} | ConvertTo-Json

$response = Invoke-RestMethod -Uri "$baseUrl/match" -Method POST `
    -Body $ruleBody -ContentType "application/json"
Write-Host "Response:"
$response | ConvertTo-Json -Depth 10
Write-Host ""

# Test 3: ML-Based Recommendation
Write-Host "3. ML-Based Recommendation (/api/v1/recommend/conditions)" -ForegroundColor Yellow
$mlBody = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    reaction_type = "Suzuki"
    k = 50
    limit = 3
} | ConvertTo-Json

$response = Invoke-RestMethod -Uri "$baseUrl/api/v1/recommend/conditions" -Method POST `
    -Body $mlBody -ContentType "application/json"
Write-Host "Response:"
$response | ConvertTo-Json -Depth 10
Write-Host ""

# Test 4: Fusion Recommendation (Deprecated)
Write-Host "4. Fusion Recommendation (/api/v1/recommend/fusion)" -ForegroundColor Yellow
$fusionBody = @{
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    k = 50
    max_variants = 3
} | ConvertTo-Json

$response = Invoke-RestMethod -Uri "$baseUrl/api/v1/recommend/fusion" -Method POST `
    -Body $fusionBody -ContentType "application/json"
Write-Host "Response:"
$response | ConvertTo-Json -Depth 10
Write-Host ""

Write-Host "All tests completed!" -ForegroundColor Green
Write-Host ""
Write-Host "Note: Protocol endpoint not yet implemented." -ForegroundColor Yellow
Write-Host "Use CLI tools for protocol recommendations:" -ForegroundColor Yellow
Write-Host "  python scripts/web_recommendation_cli.py --strategy protocol --rxn '...'" -ForegroundColor Cyan
