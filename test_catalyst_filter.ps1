# Test script to verify catalyst filtering works correctly
# This script tests the CORRECT usage of catalyst_class filtering

Write-Host "`n=== Testing Catalyst Filtering ===" -ForegroundColor Cyan

# Test 1: WRONG usage (user's original request)
Write-Host "`n--- Test 1: WRONG - Using constraints.metal_preference ---" -ForegroundColor Red
$wrongBody = @{
    reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    reaction_type = "Suzuki"
    k = 10
    limit = 5
    constraints = @{
        required_reagents = @("copper catalyst", "Cu catalyst")
        exclude_reagents = @()
        exclude_roles = @()
        metal_preference = "Cu"
    }
    relax = @{}
} | ConvertTo-Json -Depth 5

try {
    $wrongResult = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
        -Method POST -Body $wrongBody -ContentType "application/json"
    
    Write-Host "Result:" -ForegroundColor Yellow
    $wrongResult | ConvertTo-Json -Depth 10
} catch {
    Write-Host "Error: $_" -ForegroundColor Red
}

# Test 2: CORRECT usage (relax.catalyst_class)
Write-Host "`n--- Test 2: CORRECT - Using relax.catalyst_class ---" -ForegroundColor Green
$correctBody = @{
    reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    reaction_type = "Suzuki"
    k = 10
    limit = 5
    relax = @{
        catalyst_class = "Cu"
    }
    constraints = @{}
} | ConvertTo-Json -Depth 5

try {
    $correctResult = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
        -Method POST -Body $correctBody -ContentType "application/json"
    
    Write-Host "Result:" -ForegroundColor Yellow
    $correctResult | ConvertTo-Json -Depth 10
} catch {
    Write-Host "Error: $_" -ForegroundColor Red
}

# Test 3: Compare Pd vs Cu catalysts for Suzuki
Write-Host "`n--- Test 3: Compare Pd vs Cu for Suzuki ---" -ForegroundColor Cyan

Write-Host "`nWith Pd catalyst:" -ForegroundColor Magenta
$pdBody = @{
    reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    reaction_type = "Suzuki"
    k = 10
    limit = 3
    relax = @{
        catalyst_class = "Pd"
    }
} | ConvertTo-Json -Depth 5

try {
    $pdResult = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
        -Method POST -Body $pdBody -ContentType "application/json"
    
    if ($pdResult.recommendations -and $pdResult.recommendations.Count -gt 0) {
        $firstCatalyst = $pdResult.recommendations[0].reagents | Where-Object { $_.role -eq "metal_precursor" } | Select-Object -First 1
        Write-Host "First recommendation catalyst: $($firstCatalyst.name)" -ForegroundColor Green
    } else {
        Write-Host "No recommendations found for Pd" -ForegroundColor Yellow
    }
} catch {
    Write-Host "Error: $_" -ForegroundColor Red
}

Write-Host "`nWith Cu catalyst:" -ForegroundColor Magenta
$cuBody = @{
    reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    reaction_type = "Suzuki"
    k = 10
    limit = 3
    relax = @{
        catalyst_class = "Cu"
    }
} | ConvertTo-Json -Depth 5

try {
    $cuResult = Invoke-RestMethod -Uri "http://localhost:8000/api/v1/recommend/conditions" `
        -Method POST -Body $cuBody -ContentType "application/json"
    
    if ($cuResult.recommendations -and $cuResult.recommendations.Count -gt 0) {
        $firstCatalyst = $cuResult.recommendations[0].reagents | Where-Object { $_.role -eq "metal_precursor" } | Select-Object -First 1
        Write-Host "First recommendation catalyst: $($firstCatalyst.name)" -ForegroundColor Green
    } else {
        Write-Host "No recommendations found for Cu" -ForegroundColor Yellow
    }
} catch {
    Write-Host "Error: $_" -ForegroundColor Red
}

Write-Host "`n=== Test Complete ===" -ForegroundColor Cyan
