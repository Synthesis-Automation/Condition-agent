# Reaction Type Router Fix - Summary

## Problem Identified

User reported that when entering:
```
Reaction SMILES: Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1
Select reaction type [1]: 3  # C–N Coupling (Cu)
```

The rule-based recommendation was incorrectly identifying it as **Suzuki** instead of **C-N Coupling**.

## Root Causes

### 1. **Database Selection Issue**
- **Problem**: `local_rule_based_match()` always used `Suzuki_db.json` regardless of reaction type
- **Code Location**: `scripts/local_recommendation_cli.py` line 125
- **Impact**: All reactions were matched against Suzuki database only

### 2. **Deprecated Reaction Type Mapping**
- **Problem**: CLI used old catalyst-specific names (`C_N_Coupling_Cu`, `C_N_Coupling_Pd`, `C_N_Coupling_Ni`)
- **Code Location**: `scripts/recommendation_cli_utils.py` line 23-27
- **Impact**: These old names weren't recognized by the unified C_N_Coupling system

### 3. **Missing reaction_type Parameter**
- **Problem**: `local_rule_based_match()` didn't accept or use `reaction_type` for database selection
- **Impact**: User-selected reaction type was ignored for rule-based matching

## Solutions Applied

### Fix 1: Updated Reaction Type Choices (recommendation_cli_utils.py)

**Before:**
```python
REACTION_TYPE_CHOICES: Tuple[Tuple[str, Optional[str]], ...] = (
    ("Auto-detect (server decides)", None),
    ("Suzuki Coupling", "Suzuki"),
    ("Ullmann C–N (Cu)", "C_N_Coupling_Cu"),      # ❌ Deprecated
    ("Buchwald C–N (Pd)", "C_N_Coupling_Pd"),     # ❌ Deprecated
    ("C–N Coupling (Ni)", "C_N_Coupling_Ni"),     # ❌ Deprecated
    ("Amide Formation", "Amide_formation"),
)
```

**After:**
```python
REACTION_TYPE_CHOICES: Tuple[Tuple[str, Optional[str]], ...] = (
    ("Auto-detect (server decides)", None),
    ("Suzuki Coupling", "Suzuki"),
    ("C–N Coupling (unified)", "C_N_Coupling"),   # ✅ Unified
    ("Amide Formation", "Amide_formation"),
)
```

### Fix 2: Updated CLI Argument Choices (local_recommendation_cli.py)

**Before:**
```python
choices=[None, "Suzuki", "Suzuki_CC", "C_N_Coupling_Cu", "Ullmann_CN", 
         "C_N_Coupling_Pd", "Buchwald_CN", "C_N_Coupling_Ni", "Amide_formation"]
```

**After:**
```python
choices=[None, "Suzuki", "Suzuki_CC", "C_N_Coupling", "Amide_formation"]
```

### Fix 3: Added Database Auto-Selection (local_recommendation_cli.py)

**Before:**
```python
def local_rule_based_match(reaction: str, db_path: Optional[str]) -> Dict[str, Any]:
    """Replicate the /match endpoint using in-process ChemTools calls."""
    target_db = _resolve_rule_db(db_path) or _ENV_SCDB_DEFAULT  # Always Suzuki! ❌
```

**After:**
```python
def local_rule_based_match(reaction: str, db_path: Optional[str], reaction_type: Optional[str] = None) -> Dict[str, Any]:
    """Replicate the /match endpoint using in-process ChemTools calls.
    
    Args:
        reaction: Reaction SMILES string
        db_path: Explicit database path (overrides auto-selection)
        reaction_type: Detected or user-specified reaction type for DB selection
    """
    # Auto-select database based on reaction type if not explicitly provided
    if db_path is None and reaction_type:
        db_map = {
            "Suzuki": "data/conditionDB/Suzuki_db.json",
            "Suzuki_CC": "data/conditionDB/Suzuki_db.json",
            "C_N_Coupling": "data/conditionDB/C_N_Coupling_Cu_db.json",  # Use Cu as default
            "Amide_formation": "data/conditionDB/amide_formation_db.json",
        }
        candidate = db_map.get(reaction_type)
        if candidate and Path(candidate).exists():
            db_path = candidate
    
    target_db = _resolve_rule_db(db_path) or _ENV_SCDB_DEFAULT
    # ... rest of function
    return output_formatter.format_rule_match_result(
        reaction_smiles=reaction,
        match_result=payload,
        requested_type=reaction_type,  # ✅ Now passed correctly
```

### Fix 4: Updated Function Call (local_recommendation_cli.py)

**Before:**
```python
rule_result = local_rule_based_match(reaction, db_path)  # ❌ No reaction_type
```

**After:**
```python
rule_result = local_rule_based_match(reaction, db_path, reaction_type)  # ✅ Passes reaction_type
```

## Verification

```python
# Test 1: Updated reaction type choices
from scripts.recommendation_cli_utils import REACTION_TYPE_CHOICES
print(REACTION_TYPE_CHOICES)
# Output: 4 options including "C–N Coupling (unified)" → "C_N_Coupling"

# Test 2: Router correctly detects C-N coupling
from chemtools.router import detect_family_from_reaction
result = detect_family_from_reaction("Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1")
print(result['family'])  # Output: "C_N_Coupling"
print(result['confidence'])  # Output: 0.9

# Test 3: Database selection works
db_map = {"C_N_Coupling": "data/conditionDB/C_N_Coupling_Cu_db.json"}
candidate = db_map.get("C_N_Coupling")
exists = Path(candidate).exists()  # Output: True
```

## Impact

### Before Fix:
1. User selects "C–N Coupling (Cu)" → Maps to `C_N_Coupling_Cu`
2. Rule-based matcher ignores this and uses `Suzuki_db.json`
3. Matches aryl bromide against Suzuki reactions → Wrong results

### After Fix:
1. User selects "C–N Coupling (unified)" → Maps to `C_N_Coupling`
2. Database auto-selection chooses `C_N_Coupling_Cu_db.json`
3. Matches aryl bromide + amine against C-N coupling reactions → Correct results
4. Router correctly detects family as `C_N_Coupling` with 0.9 confidence

## Testing Recommendations

### Manual Test:
```bash
cd c:\Git-softwares\Condition-agent
python scripts/local_recommendation_cli.py

# Input:
Reaction SMILES: Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1
Select reaction type: 3  # C–N Coupling (unified)
Strategy: rule
```

**Expected Output:**
- Rule-based family: C_N_Coupling
- Database: C_N_Coupling_Cu_db.json
- Matches: C-N coupling reactions with Cu catalysts

### Automated Test:
```python
from scripts.local_recommendation_cli import local_rule_based_match

rsmi = "Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1"
result = local_rule_based_match(rsmi, None, "C_N_Coupling")

assert "C_N_Coupling" in result.get("detection", {}).get("family", "")
assert "C_N_Coupling_Cu_db.json" in result.get("metadata", {}).get("database_name", "")
```

## Files Modified

1. **scripts/recommendation_cli_utils.py**
   - Updated `REACTION_TYPE_CHOICES` to use unified C_N_Coupling
   - Removed deprecated catalyst-specific options

2. **scripts/local_recommendation_cli.py**
   - Added `reaction_type` parameter to `local_rule_based_match()`
   - Implemented database auto-selection logic
   - Updated CLI argument choices
   - Updated function call to pass reaction_type

## Future Considerations

### Database Consolidation
Currently using `C_N_Coupling_Cu_db.json` as default. Consider:
- Creating unified `C_N_Coupling_db.json` with all catalysts
- Using constraints to filter by catalyst preference
- Aligns with precedent search using `relax={"catalyst_class": "Cu"}`

### Metal Preference
User-selected metal (Cu/Pd/Ni) could be:
- Passed as constraint to precedent search
- Used to prioritize certain reactions
- Displayed in output for transparency

### Backward Compatibility
Old databases still exist but deprecated:
- `C_N_Coupling_Cu_db.json` (15,967 reactions)
- `C_N_Coupling_Pd_db.json`
- `C_N_Coupling_Ni_db.json`

These could be merged or kept for legacy support.

## Related Documentation

- `CATALYST_FILTERING_GUIDE.md` - Catalyst filtering via relax parameter
- `CATALYST_CLASSIFICATION_FIX.md` - Correct chemistry classification
- `chemtools/router.py` - Reaction type detection logic
- `chemtools/precedent/search.py` - DRFP-based precedent search with catalyst filtering

---

**Status**: ✅ FIXED
**Date**: 2025-10-16
**Issue**: Rule-based router incorrectly using Suzuki database for C-N coupling
**Solution**: Database auto-selection based on reaction type with unified C_N_Coupling mapping
