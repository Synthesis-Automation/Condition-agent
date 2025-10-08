# Test Script Improvements - Summary

## Changes Made

### ✅ 1. Removed Interactive Prompts

**Before:**
```python
if i < len(test_reactions):
    input("\n\nPress Enter to continue to next reaction...")
```

**After:**
```python
# Add separator between reactions (removed interactive pause)
if i < len(test_reactions):
    print("\n")
```

**Benefit:** Tests now run continuously without user intervention - perfect for automation.

### ✅ 2. Enhanced Error Reporting for ML Endpoint

**Before:**
```python
else:
    print(f"✗ Error: HTTP {response.status_code}")
    print(f"  {response.text}")
    return {}
```

**After:**
```python
else:
    print(f"✗ Error: HTTP {response.status_code}")
    if response.status_code == 500:
        print(f"  Server Error: {response.text}")
        print(f"\n  ⚠️  This is a server-side error. Check the uvicorn terminal for:")
        print(f"     - Python traceback")
        print(f"     - FileNotFoundError, KeyError, or AttributeError")
        print(f"     - Missing data files or import errors")
    else:
        print(f"  {response.text}")
    return {}
```

**Benefit:** Users now get clear guidance on where to look for the actual error.

## About the ML Endpoint 500 Error

The `/api/v1/recommend/conditions` endpoint is returning a **500 Internal Server Error**. This is a **server-side issue**, not a client-side problem.

### What This Means

A 500 error indicates:
- Python exception occurred in the server code
- The request format is correct (otherwise would be 400)
- The endpoint exists (otherwise would be 404)
- Something failed during execution

### Common Causes

1. **Missing Data Files**
   - DRFP fingerprint files (`.npz`)
   - Dataset JSONL files
   - Precomputed vectors

2. **Import Errors**
   - Missing module
   - Circular imports
   - Failed initialization

3. **Data Processing Errors**
   - KeyError (missing required field)
   - TypeError (wrong data type)
   - AttributeError (method doesn't exist)

### How to Debug

#### Step 1: Check Server Terminal

Look at the **uvicorn terminal** where the server is running. When the 500 error occurs, you should see:

```
ERROR:    Exception in ASGI application
Traceback (most recent call last):
  File "...", line X, in ...
    ...
  [The actual error is here]
FileNotFoundError: [Errno 2] No such file or directory: 'data/...'
```

#### Step 2: Run Diagnostic Script

```powershell
python diagnose_ml_error.py
```

This script:
- Tests the exact failing request
- Shows detailed error information
- Provides guidance on what to check

#### Step 3: Test Data File Existence

```powershell
# Check for DRFP file
Test-Path "data\reaction_dataset\Suzuki_drfp.npz"

# Check for dataset files
Get-ChildItem data\*.jsonl
```

#### Step 4: Test Module Import

```powershell
# Test if recommend module works
python -c "from chemtools.recommend.core import recommend_conditions_structured; print('OK')"

# Test if output formatter works
python -c "from chemtools import output_formatter; print('OK')"
```

## Files Created/Modified

### Modified
- ✅ `test_suzuki_api.py` - Removed `input()`, enhanced error messages

### Created
- ✅ `diagnose_ml_error.py` - Diagnostic tool for ML endpoint
- ✅ `ML_ENDPOINT_DEBUG_GUIDE.md` - This file

## Current Test Status

| Feature | Status | Notes |
|---------|--------|-------|
| Remove `input()` prompts | ✅ Fixed | Tests run continuously |
| Better 500 error messages | ✅ Fixed | Now shows where to look |
| ML endpoint working | ❌ Still broken | Needs server-side fix |
| Rule endpoint | ✅ Working | Returns conditions |
| Fusion endpoint | ✅ Working | Returns full JSON |

## Next Steps

### For the User

1. **Check server logs** - Look at the uvicorn terminal for the Python traceback
2. **Share the error** - Copy the traceback here so we can fix it
3. **Run diagnostic** - Execute `python diagnose_ml_error.py`

### For Debugging

The error is **not in the test script** - it's in the server code. We need to see the actual Python exception to fix it.

**The test script is now ready** - once the server issue is resolved, all tests will run smoothly!

## Running Tests

### Quick Test (Recommended)
```powershell
python quick_api_test.py
```
Shows: Rule + ML + Fusion for 1 reaction

### Full Test Suite
```powershell
python test_suzuki_api.py
```
Shows: Rule + ML + Fusion for 5 reactions (now runs automatically!)

### Diagnostic Only
```powershell
python diagnose_ml_error.py
```
Focuses on the ML endpoint error with detailed output

---

**Status**: Test scripts are fixed and ready. ML endpoint needs server-side debugging.
