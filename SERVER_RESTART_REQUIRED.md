# IMPORTANT: Server Restart Required!

## Current Status

The test script has been fixed with the following changes:

### ✅ Fixed Issues in `test_suzuki_api.py`:

1. **Parameter Naming Conflict** (Line 341-345)
   - Changed `test_rule`, `test_ml`, `test_fusion` parameters
   - To `run_rule`, `run_ml`, `run_fusion`
   - This fixes the `TypeError: 'bool' object is not callable` error

2. **Rule-Based Conditions Handling** (Lines 88-136)
   - Added support for both list and dict structures
   - Handles cases where `conditions` might be a list
   - More robust error handling for various data formats

3. **Better Error Tracebacks** (Lines 139-144)
   - Added `traceback.print_exc()` to exception handlers
   - This will help diagnose issues better

## ⚠️ ML Endpoint Still Returns 500

The `/api/v1/recommend/conditions` endpoint is returning a 500 Internal Server Error.

### Possible Causes:

1. **Server not restarted** after the middleware fix
2. **Missing DRFP data** or other prerequisite data
3. **Error in `recommend_conditions_structured()` function**

### How to Fix:

#### Step 1: RESTART THE SERVER (CRITICAL!)

```powershell
# In the uvicorn terminal:
# Press Ctrl+C to stop

# Then restart:
uvicorn app.main:app --reload --port 8000
```

The `--reload` flag should auto-reload, but sometimes it doesn't catch all changes, especially:
- Import errors
- Module-level code
- Middleware changes

**Manual restart is recommended!**

#### Step 2: Check Server Logs

After restarting, watch the server terminal for:
- Import errors during startup
- Stack traces when the endpoint is called
- Missing data file warnings

#### Step 3: Run Tests Again

```powershell
# Quick test first:
python quick_api_test.py

# If that works, run full test:
python test_suzuki_api.py
```

## Expected Behavior After Restart

### Rule-Based (`/match`)
```
✓ Match Type: scheme
  Entry: Aryl iodides/bromides + aryl boron (SPhos set)
  
  Recommended Conditions:
    Catalyst: Pd(OAc)2 + SPhos
    Base: K3PO4
    Solvent: Dioxane
```

### ML-Based (`/api/v1/recommend/conditions`)
```
✓ Detected Type: Suzuki (confidence: 98.50%)
  Found 5 recommendations
  
  1. Confidence: 95.20% | Support: 47 precedents
     Catalyst: Pd(OAc)2/XPhos
     Base: K3PO4
     Solvent: 1,4-Dioxane
```

### Fusion (`/api/v1/recommend/fusion`)
```
✓ Fusion System Active
  Adaptive Weights:
    α (precedents) = 0.650
    β (analytics)  = 0.250
    γ (rules)      = 0.080
    δ (ML)         = 0.020
```

## Debugging the 500 Error

If the error persists after restart, check the server terminal output for the actual Python exception. Look for:

```
ERROR:    Exception in ASGI application
Traceback (most recent call last):
  ...
  [The actual error will be here]
```

Common issues:
- `FileNotFoundError`: Missing DRFP .npz file or dataset
- `KeyError`: Missing required field in data
- `AttributeError`: Method not found (import issue)
- `ModuleNotFoundError`: Missing dependency

## Quick Diagnostic Commands

```powershell
# Check if DRFP file exists:
Test-Path "data\reaction_dataset\Suzuki_drfp.npz"

# Check if recommend module imports:
python -c "from chemtools import recommend; print('OK')"

# Check if chem.recommend works:
python -c "from chemtools import chem; print(chem.recommend)"

# Test recommend_conditions_structured directly:
python -c "from chemtools.recommend.core import recommend_conditions_structured; print('OK')"
```

---

**Action Required**: 
1. ✅ Stop the server (Ctrl+C)
2. ✅ Start it again: `uvicorn app.main:app --reload --port 8000`
3. ✅ Run `python quick_api_test.py`
4. ✅ Check results and server logs
