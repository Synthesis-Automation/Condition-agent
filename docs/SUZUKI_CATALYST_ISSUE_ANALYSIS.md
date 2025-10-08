# Suzuki Catalyst Missing - Root Cause Analysis

## Problem Statement

When requesting ML recommendations for a Suzuki coupling reaction:
```
Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1
```

The JSON response **does NOT include Pd catalyst** information. Only starting materials are returned.

---

## Root Cause: **CODING ISSUE** (Not Dataset Issue)

### Investigation Results

1. **Dataset Status: ✅ GOOD**
   - `data/reaction_dataset/Suzuki.jsonl` contains **50,215 reactions**
   - Catalyst information **IS present** in the dataset
   - Reagent roles include: `METAL_PRECURSOR`, `PREFORMED_METAL_CATALYST`, `LIGAND`
   - Sample reagent: `{"name": "Palladium", "cas": "7440-05-3", "role": "METAL_PRECURSOR"}`

2. **Precedent Search: ❌ FAILING**
   - DRFP similarity search returns **0 matching reactions**
   - `total_considered: 0`
   - `top_precedents: []`
   - This explains why recommendations only contain starting materials

3. **Recommendation Output: ❌ INCOMPLETE**
   - Only 2 chemicals in recommendation (both starting materials)
   - No catalyst, base, solvent, or ligand
   - Confidence: 0.3 (very low)
   - Support: 0 precedents

### Why Precedent Search Returns 0 Results

The issue is in the **DRFP feature matching** or **family mapping**:

| Component | Status | Notes |
|-----------|--------|-------|
| Dataset exists | ✅ | 50K+ Suzuki reactions available |
| Reaction type detection | ✅ | Correctly detected as "Suzuki_CC" |
| Feature extraction | ⚠️ | May be extracting features incorrectly |
| DRFP similarity search | ❌ | Not finding any matches |
| Precomputed DRFP | ⚠️ | May be missing or outdated |

---

## Possible Causes & Solutions

### Cause 1: DRFP Features Not Precomputed for Suzuki

**Symptoms:**
- Precedent search returns 0 results
- System falls back to generic features

**Solution:**
```bash
# Regenerate precomputed DRFP features
python scripts/precompute_drfp.py
```

**Files to check:**
- `artifacts/all_drfp.npz` - Should contain Suzuki reactions
- `artifacts/all_drfp_4096.npz` - High-dimensional version

---

### Cause 2: Family Name Mismatch

**The Code:**
- Detection returns: `"Suzuki_CC"`
- Dataset uses: `"Suzuki"` (in `reaction_type` field)

**Possible Fix:**

```python
# In chemtools/precedent.py - _dataset_family_map()
def _dataset_family_map(raw: str) -> str:
    """Map dataset family names to internal family names."""
    raw_lower = (raw or "").lower().strip()
    if raw_lower in {"suzuki", "suzuki coupling", "suzuki_cc"}:
        return "Suzuki_CC"
    # ... other mappings
```

---

### Cause 3: DRFP Similarity Threshold Too Strict

**The Issue:**
- DRFP similarity might be too low for generic aryl bromide + boronic acid
- System might be filtering out valid matches

**Possible Fix:**

```python
# In chemtools/precedent.py - _candidate_pool()
# Lower the similarity threshold for Suzuki reactions
if family_txt in {"Suzuki_CC", "Suzuki"}:
    min_similarity = 0.3  # Lower threshold for broader matching
else:
    min_similarity = 0.5  # Default threshold
```

---

### Cause 4: Reagent Role Aggregation Filters Out Catalyst

**Even if precedents are found**, the recommendation building logic might be dropping non-starting-material reagents.

**Check:**
```python
# In chemtools/recommend.py - look for chemical aggregation logic
# Ensure METAL_PRECURSOR, PREFORMED_METAL_CATALYST, LIGAND roles are preserved
```

---

## Recommended Action Plan

### Step 1: Verify DRFP Precomputation ⚠️ **DO THIS FIRST**

```powershell
# Check if DRFP artifacts exist and are recent
Get-ChildItem artifacts\*.npz | Select-Object Name, Length, LastWriteTime

# Regenerate DRFP features for all datasets
python scripts\precompute_drfp.py
```

### Step 2: Add Debug Logging

Add logging to see why precedent search fails:

```python
# In chemtools/precedent.py - knn() function
import logging
logger = logging.getLogger(__name__)

def knn(family: str, features: Dict[str, Any], k: int = 50, relax: Dict[str, Any] | None = None) -> Dict[str, Any]:
    logger.info(f"KNN search: family={family}, k={k}")
    result = _knn_impl(family, features, k, relax)
    logger.info(f"KNN results: {len(result.get('precedents', []))} precedents found")
    return result
```

### Step 3: Family Name Consistency Check

```python
# Quick test to see what family names exist in dataset
import json
families = set()
with open("data/reaction_dataset/Suzuki.jsonl", encoding="utf-8") as f:
    for line in f:
        rxn = json.loads(line)
        families.add(rxn.get("reaction_type"))

print(f"Unique family names in Suzuki dataset: {families}")
# Should show: {"Suzuki"}

# Now check what precedent.py expects:
from chemtools.precedent import _dataset_family_map
print(f"Mapped family: {_dataset_family_map('Suzuki')}")
# Should return: "Suzuki_CC" or similar
```

### Step 4: Force DRFP Computation at Runtime

If precomputed DRFP is the issue, force runtime computation:

```python
# In the API call or test
relax = {
    "use_drfp": True,
    "precompute_drfp": False,  # Force runtime computation
}

raw_data = recommend.recommend_conditions_structured(
    reaction=reaction,
    reaction_type="Suzuki",
    k=50,
    limit=3,
    relax=relax
)
```

---

## Quick Test Script

```python
"""Test if precedent search works for Suzuki."""

from chemtools import precedent
from chemtools.features import molecular as feat_molecular

# Suzuki reaction
elec = "Brc1ccccc1"  # Aryl bromide
nuc = "OB(O)c1ccccc1"  # Boronic acid

# Extract features
features = feat_molecular.featurize(elec, nuc)

# Search for precedents (try both family names)
for family_name in ["Suzuki", "Suzuki_CC"]:
    print(f"\n=== Testing family: {family_name} ===")
    
    pack = precedent.knn(
        family=family_name,
        features=features,
        k=10,
        relax={"use_drfp": True, "precompute_drfp": True}
    )
    
    precs = pack.get("precedents", [])
    print(f"Found {len(precs)} precedents")
    
    if precs:
        print("✅ Precedent search WORKS!")
        print(f"First precedent core: {precs[0].get('core')}")
        break
else:
    print("\n❌ NO PRECEDENTS FOUND with either family name!")
    print("   This confirms the DRFP feature matching is broken.")
```

---

## Expected Output After Fix

```json
{
  "recommendations": [
    {
      "rank": 1,
      "confidence": 0.85,
      "reagents": [
        {"role": "catalyst", "name": "Pd(PPh3)4", "cas": "14221-01-3"},
        {"role": "base", "name": "Potassium carbonate", "cas": "584-08-7"},
        {"role": "solvent", "name": "Toluene", "cas": "108-88-3"}
      ],
      "conditions": {
        "temperature": {"value": 100, "unit": "°C"},
        "time": {"value": 12, "unit": "hours"}
      },
      "precedent_count": 127
    }
  ]
}
```

---

## Summary

**Root Cause:** DRFP precedent search returns 0 results for Suzuki reactions

**Most Likely Issue:** Precomputed DRFP features missing or family name mismatch

**Immediate Action:** Run `python scripts/precompute_drfp.py` to regenerate features

**Verification:** Run test script above to confirm precedent search works

**Impact:** Once fixed, Suzuki recommendations will include full catalyst/ligand/base/solvent information from 50K+ precedent reactions
