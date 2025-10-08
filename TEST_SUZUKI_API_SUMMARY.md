# Test Suzuki API - Quick Summary

## Files Created

1. **`test_suzuki_api.py`** - Main test script for FastAPI endpoints
2. **`TEST_SUZUKI_API_README.md`** - Comprehensive documentation

## What Was Created

A complete test suite for validating Suzuki reaction recommendations via three different FastAPI endpoints:

### Endpoints Tested

1. **`POST /match`** - Rule-based (SMARTS pattern matching)
2. **`POST /api/v1/recommend/conditions`** - ML-based (DRFP k-NN)
3. **`POST /api/v1/recommend/fusion`** - Fusion (multi-source evidence)

### Test Coverage

- 5 diverse Suzuki reactions (simple → challenging)
- Automatic server health check
- Formatted output with color-coded results
- Interactive progression between tests
- Comprehensive result comparison

## How to Use

### 1. Start Server

```powershell
# Terminal 1
.\.venv\Scripts\Activate.ps1
uvicorn app.main:app --reload --port 8000
```

### 2. Run Tests

```powershell
# Terminal 2
.\.venv\Scripts\Activate.ps1
python test_suzuki_api.py
```

## Key Features

✅ **Health Check**: Verifies server is running before tests
✅ **Three Approaches**: Tests rule, ML, and fusion simultaneously
✅ **Diverse Substrates**: From simple to challenging reactions
✅ **Detailed Output**: Shows conditions, scores, and metadata
✅ **Interactive**: Pause between reactions for review
✅ **Error Handling**: Graceful handling of API errors

## Sample Output

```
================================================================================
  APPROACH 1: RULE-BASED (SMARTS Pattern Matching)
================================================================================
✓ Match Type: exact
  Entry: Ph-Br_standard_conditions
  Recommended Conditions:
    Catalyst: Pd(OAc)2 + XPhos
    Base: K3PO4
    Solvent: Dioxane

================================================================================
  APPROACH 3: FUSION (Multi-Source Evidence)
================================================================================
✓ Fusion System Active
  Adaptive Weights:
    α (precedents) = 0.650
    β (analytics)  = 0.250
    γ (rules)      = 0.080
    δ (ML)         = 0.020
```

## Differences from Previous Script

The previous `demo_fusion_suzuki.py` was **module-based** (direct Python imports):
- Called `recommend_from_reaction()` directly
- No HTTP/API layer
- Used for testing core logic

This new `test_suzuki_api.py` is **API-based** (HTTP requests):
- Tests actual FastAPI endpoints
- Validates REST API contracts
- Used for integration testing
- Robot/client-friendly testing

## When to Use Which

**Use `demo_fusion_suzuki.py`**:
- Testing recommendation logic directly
- Debugging fusion algorithm
- Local development without server

**Use `test_suzuki_api.py`**:
- Testing deployed API
- Integration testing
- Validating HTTP contracts
- Robot/automation testing

## Next Steps

1. ✅ Server running? → Run `python test_suzuki_api.py`
2. 📊 Review outputs for each approach
3. 🔍 Compare rule vs ML vs fusion recommendations
4. 📝 Check fusion adaptive weights and reasoning

## Related Documentation

- Full guide: `TEST_SUZUKI_API_README.md`
- API docs: http://localhost:8000/docs
- Repository guidelines: `AGENTS.md`
