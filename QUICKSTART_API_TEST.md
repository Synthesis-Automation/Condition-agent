# QUICK START: Test Suzuki API Endpoints

## TL;DR - 2 Steps

### Step 1: Start Server (Terminal 1)

```powershell
.\.venv\Scripts\Activate.ps1
uvicorn app.main:app --reload --port 8000
```

### Step 2: Run Test (Terminal 2)

**Quick test (5 seconds):**
```powershell
.\.venv\Scripts\Activate.ps1
python quick_api_test.py
```

**OR Full test (2-3 minutes):**
```powershell
.\.venv\Scripts\Activate.ps1
python test_suzuki_api.py
```

## What You'll See

```
✅ Server is running

1. RULE-BASED → Pd/XPhos, K3PO4, Dioxane
2. ML-BASED   → Pd/XPhos, K3PO4, Dioxane (47 precedents)
3. FUSION     → Pd/XPhos, K3PO4, Dioxane (weights: 0.650/0.250/0.080/0.020)
```

## Files Created

✅ `test_suzuki_api.py` - Comprehensive (5 reactions)  
✅ `quick_api_test.py` - Quick test (1 reaction)  
✅ `TEST_SUZUKI_API_README.md` - Full docs  
✅ `SUZUKI_API_TEST_COMPLETE.md` - Complete summary

## Need Help?

- Server won't start? Check: `uvicorn --version`
- Import error? Run: `pip install requests`
- Read full docs: `TEST_SUZUKI_API_README.md`

---
**That's it! Start the server, run the script. 🚀**
