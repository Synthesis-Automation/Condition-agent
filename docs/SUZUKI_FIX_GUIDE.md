# Suzuki Catalyst Missing - Summary & Solution

## Problem

ML recommendations for Suzuki coupling:
```
Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1
```

**Returns:** Only 2 starting materials, no Pd catalyst, base, or solvent

---

## Root Cause: **CODING ISSUE**

### What I Found:

1. **Dataset is GOOD** ✅
   - `data/reaction_dataset/Suzuki.jsonl` has 50,215 reactions
   - Contains catalyst data: `METAL_PRECURSOR`, `PREFORMED_METAL_CATALYST`, `LIGAND` roles
   - Example: `{"name": "Palladium", "cas": "7440-05-3", "role": "METAL_PRECURSOR"}`

2. **Precedent Search is BROKEN** ❌
   - Returns **0 precedents** (`total_considered: 0`)
   - DRFP similarity search fails to find matching reactions
   - Without precedents → no catalyst/base/solvent in recommendations

3. **Evidence from Raw Output:**
   ```json
   {
     "precedents": {
       "total_considered": 0,
       "core_support": 0,
       "top_precedents": []
     },
     "recommendations": [{
       "chemicals": [
         {"role": "starting_material", ...},
         {"role": "starting_material", ...}
       ],
       "summary": {
         "confidence": 0.3,
         "support": {"count": 0}
       }
     }]
   }
   ```

---

## Most Likely Causes

### 1. **DRFP Features Not Precomputed** (80% likely)

The DRFP similarity search requires precomputed feature vectors for all 50K+ Suzuki reactions. If these are missing or outdated, the search returns 0 results.

**Solution:**
```powershell
# Regenerate precomputed DRFP features
python scripts\precompute_drfp.py
```

**What this does:**
- Loads all reaction datasets (C-N coupling, Suzuki, Amide, etc.)
- Computes DRFP (differential reaction fingerprint) vectors for each
- Saves to `artifacts/all_drfp.npz` and `artifacts/all_drfp_4096.npz`
- Enables fast k-NN similarity search

---

### 2. **Family Name Mismatch** (15% likely)

Detection returns `"Suzuki_CC"` but dataset uses `"Suzuki"`

**Check:**
```python
# In chemtools/precedent.py
def _dataset_family_map(raw: str) -> str:
    """Map dataset 'Suzuki' to internal 'Suzuki_CC'"""
    raw_lower = (raw or "").lower().strip()
    if raw_lower in {"suzuki", "suzuki coupling"}:
        return "Suzuki_CC"  # <-- Check if this mapping exists
```

**Fix if missing:**
Add the mapping so both names work

---

### 3. **DRFP Library Issue** (5% likely)

The DRFP feature extraction or comparison might have bugs for Suzuki-type reactions

**Workaround:**
```python
# Force runtime DRFP computation (slower but bypasses precomputed cache)
relax = {
    "use_drfp": True,
    "precompute_drfp": False  # Compute on-the-fly
}
```

---

## Quick Fix to Test

**Option A: Regenerate DRFP** (Recommended)

```powershell
# This may take 10-30 minutes for 50K+ reactions
python scripts\precompute_drfp.py

# Then test again
python -c "from chemtools import recommend; result = recommend.recommend_conditions_structured(reaction='Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1', reaction_type='Suzuki', k=50, limit=3); print(f'Precedents: {len(result.get(\"precedents\", {}).get(\"top_precedents\", []))}'); print(f'Chemicals: {len(result.get(\"recommendations\", [{}])[0].get(\"chemicals\", []))}')"
```

**Expected after fix:**
```
Precedents: 10
Chemicals: 5-8  (starting materials + catalyst + base + solvent + ligand)
```

---

**Option B: Use Rule-Based Fallback** (Immediate workaround)

While the ML system is broken for Suzuki, use the SchemeConditionDB instead:

```python
from scdb_matcher import Loader, match_reaction

# Load Suzuki rule database
loader = Loader()
db = loader.load_family("Suzuki")  # or load from file

# Match reaction
result = match_reaction(
    db=db,
    reaction="Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
)

# Result will include expert-curated Pd catalyst conditions
```

---

## What to Expect After Fix

Once DRFP features are regenerated, Suzuki recommendations should include:

```json
{
  "recommendations": [{
    "rank": 1,
    "confidence": 0.87,
    "reagents": [
      {"role": "catalyst", "name": "Pd(PPh3)4", "cas": "14221-01-3"},
      {"role": "base", "name": "K2CO3", "cas": "584-08-7"},
      {"role": "solvent", "name": "Toluene", "cas": "108-88-3"},
      {"role": "ligand", "name": "PPh3", "cas": "603-35-0"}
    ],
    "conditions": {
      "temperature": {"value": 100, "unit": "°C"},
      "time": {"value": 12, "unit": "hours"}
    },
    "precedent_count": 127
  }]
}
```

---

## Verification Steps

1. **Check if DRFP exists:**
   ```powershell
   Get-ChildItem artifacts\*.npz
   ```

2. **Run precomputation if missing:**
   ```powershell
   python scripts\precompute_drfp.py
   ```

3. **Test Suzuki again:**
   ```powershell
   python test_suzuki_raw.py
   ```

4. **Check precedents count:**
   Should see `"total_considered": 100+` instead of `0`

---

## Why This Happens

The ML recommendation system has **two phases**:

1. **Precedent Search** (DRFP similarity)
   - Finds similar reactions in database
   - Requires precomputed features for speed
   - **THIS IS FAILING** → returns 0 precedents

2. **Condition Aggregation** (voting/ranking)
   - Aggregates reagents from precedents
   - Votes on most common catalyst/base/solvent
   - **CAN'T RUN** → no precedents to aggregate from

Without precedents, the system only has the input starting materials, so that's all it returns.

---

## Next Steps

1. Run `python scripts\precompute_drfp.py`
2. Wait for completion (10-30 min)
3. Test Suzuki recommendation again
4. If still broken, check family name mapping in `chemtools/precedent.py`
5. If still broken, add debug logging to see what DRFP similarity scores are being computed

---

## Alternative: Check if C-N Coupling Works

To verify this is Suzuki-specific and not a global DRFP issue:

```python
# Test C-N coupling (Buchwald-Hartwig)
from chemtools import recommend

result = recommend.recommend_conditions_structured(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    reaction_type="C_N_Coupling_Pd",
    k=50,
    limit=3
)

# Check if this returns catalyst
prec_count = len(result.get("precedents", {}).get("top_precedents", []))
print(f"C-N precedents: {prec_count}")

if prec_count > 0:
    print("✓ DRFP works for C-N coupling → Suzuki-specific issue")
else:
    print("✗ DRFP broken globally → regenerate all features")
```

---

**In Summary:** This is a **CODING/INFRASTRUCTURE ISSUE**, not a dataset problem. The Suzuki dataset has all the catalyst information, but the DRFP feature matching is returning 0 precedents, so the recommendation system has nothing to aggregate from.
