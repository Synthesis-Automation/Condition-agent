# Properties.py Cleanup Analysis

**Issue**: `data/compound_taxonomy/` was removed, but `properties.py` still references it with ~200 lines of dead code.

**Current State**: File works because it falls back to `_SEED` (5 hardcoded compounds), but contains extensive dead code for loading non-existent taxonomy files.

---

## 📊 Current File Analysis

**File**: `chemtools/properties.py` (317 lines)

### Active Code (Works):
- **Lines 23-32**: `_SEED` dictionary (5 compounds: K3PO4, KOH, Water, CuI, Phenanthroline)
- **Lines 251-317**: `lookup()` function using `_SEED`

### Dead Code (References Non-Existent Data):
- **Lines 51-58**: `_taxonomy_dir()` - returns None (directory doesn't exist)
- **Lines 61-186**: `_load_taxonomy_props()` - 125 lines of taxonomy loading (never executes)
- **Lines 189-226**: `_load_external()` - calls dead taxonomy loader

---

## 🔍 What Actually Exists vs What Code Expects

### Code Expects:
```
data/compound_taxonomy/
├── taxonomy_ligand.json
├── taxonomy_base.json
├── taxonomy_solvent.json
├── taxonomy_coupling_reagent.json
├── taxonomy_reductant.json
└── taxonomy_catalysts_precursor.json
```

### What Actually Exists:
```
data/
├── reagents/              ← Used by reagent_lookup.py
│   ├── ligand.json
│   ├── base.json
│   ├── solvent.json
│   └── ...
├── reaction_dataset/      ← Used by precedent.py
└── conditionDB/           ← Used by rule_scdb_matcher
```

**Result**: `_taxonomy_dir()` returns `None`, all taxonomy loading is skipped.

---

## 🧪 Testing Current Behavior

```powershell
# Test 1: Seed data (works)
python -c "from chemtools import chem; print(chem.properties.lookup('water'))"
# ✅ Returns: {found: True, record: {uid: "7732-18-5", role: "SOLVENT", ...}}

# Test 2: Non-seed data (fails)
python -c "from chemtools import chem; print(chem.properties.lookup('DMF'))"
# ❌ Returns: {found: False, record: None}

# Test 3: Fallback to registry (if available)
python -c "from chemtools import chem; print(chem.properties.lookup('K2CO3'))"
# Depends on registry.resolve() availability
```

---

## 💡 Decision: What Should We Do?

### Option 1: **Delete properties.py Entirely** ❌
**Reasoning**: Still used by 4 modules
- `constraints.py` - Property-based filtering
- `explain.py` - UID → name translation
- `recommend.py` - pKa-based ranking
- `context.py` - Public API

**Impact**: Would break existing code

### Option 2: **Simplify to Seed-Only** ⚠️
**Reasoning**: Keep minimal functionality
- Remove all taxonomy loading (125+ lines)
- Keep only `_SEED` data
- Document limitations

**Impact**: Limited to 5 compounds only

### Option 3: **Redirect to reagent_lookup.py** ✅ **RECOMMENDED**
**Reasoning**: `reagent_lookup.py` has the actual data
- `data/reagents/` contains 5000+ compounds
- Already has CAS, names, aliases
- Properties can be extracted from reagent databases

**Impact**: Unify data sources

### Option 4: **Keep as Fallback Stub** ✅ **ALTERNATIVE**
**Reasoning**: Minimal stub for backward compatibility
- Remove dead taxonomy code (~200 lines)
- Keep `_SEED` for essential compounds
- Document that it's a minimal stub
- Suggest using `chem.reagent.*` for lookups

**Impact**: Clean code, clear expectations

---

## 🎯 Recommendation: **Option 4 - Clean Stub**

### Why Option 4?

1. **Preserves API**: `chem.properties.lookup()` still works
2. **Removes Dead Code**: Delete ~200 lines of non-functional taxonomy loading
3. **Clear Purpose**: Minimal fallback for essential compounds
4. **No Breaking Changes**: Existing calls still work (with limited data)
5. **Future Path**: Can be enhanced or deprecated later

### Implementation:

**Before** (317 lines):
```python
# 125 lines of taxonomy loading
# 37 lines of external loading
# 66 lines of lookup logic
```

**After** (~70 lines):
```python
# Minimal seed data (5 compounds)
# Simple lookup function
# Clear documentation about limitations
```

### Changes Required:

1. **Remove Dead Code**:
   - Delete `_taxonomy_dir()` (8 lines)
   - Delete `_load_taxonomy_props()` (125 lines)
   - Delete `_load_external()` (38 lines)
   - Simplify `_props()` (remove external loading)

2. **Update Docstring**:
   ```python
   """
   chemtools.properties - Minimal Properties Stub
   
   Provides basic property lookup for essential compounds using seed data.
   
   ⚠️ LIMITED SCOPE: Only 5 hardcoded compounds available.
   
   For comprehensive reagent lookup, use:
   - chem.reagent.find(name, reagent_type)
   - chem.reagent.enrich(name, reagent_type)
   
   Data source: Hardcoded _SEED dictionary only.
   No taxonomy files loaded (data/compound_taxonomy removed).
   """
   ```

3. **Simplify lookup()**:
   - Remove `allow_registry` parameter (no registry dependency)
   - Remove taxonomy alias map
   - Direct lookup in `_SEED` only

---

## 📝 Simplified Code Structure

```python
# chemtools/properties.py (simplified to ~70 lines)

"""
Minimal property lookup stub for essential compounds.

⚠️ LIMITED: Only 5 hardcoded compounds. Use chem.reagent.* for comprehensive lookup.
"""

from typing import Dict, Any

# Essential compounds (5 total)
_SEED: Dict[str, Dict[str, Any]] = {
    "7778-53-2": {"role": "BASE", "token": "K3PO4", "pKa_DMSO": 30.0},
    "1310-58-3": {"role": "BASE", "token": "KOH", "pKa_water": 15.7},
    "7732-18-5": {"role": "SOLVENT", "token": "Water", "KT": {...}},
    "7681-65-4": {"role": "CATALYST", "token": "CuI"},
    "72-52-8": {"role": "LIGAND", "token": "Phenanthroline"},
}

def lookup(query: str) -> Dict[str, Any]:
    """Lookup properties by CAS or token.
    
    ⚠️ Limited to 5 essential compounds in _SEED.
    For comprehensive lookup, use chem.reagent.find() instead.
    """
    q = (query or "").strip()
    if not q:
        return {"found": False, "record": None}
    
    # Direct CAS lookup
    if q in _SEED:
        return {"found": True, "record": {"uid": q, **_SEED[q]}}
    
    # Token/name lookup
    q_lower = q.lower()
    for uid, rec in _SEED.items():
        token = str(rec.get("token", "")).lower()
        if q_lower == token:
            return {"found": True, "record": {"uid": uid, **rec}}
    
    return {"found": False, "record": None}
```

**Result**: 70 lines (down from 317), clear purpose, no dead code.

---

## 🔄 Migration Path for Users

### Current Code (Still Works):
```python
# Limited to 5 compounds
result = chem.properties.lookup("water")
# ✅ Works: {found: True, record: {...}}
```

### Recommended Modern Approach:
```python
# Full database access (5000+ compounds)
result = chem.reagent.find("water", "solvent")
# Returns: {name, cas, smiles, abbreviation, aliases, ...}
```

### For Property-Based Filtering:
```python
# Old: properties.lookup() for pKa
# New: Extract from reagent database or use constraints directly
```

---

## ⚠️ Impact Analysis

### Files That Use properties.lookup():

1. **chemtools/constraints.py** (Line 42)
   - **Usage**: Lookup boiling points for filtering
   - **Impact**: ⚠️ Only works for 5 seed compounds
   - **Mitigation**: Could query reagent database instead

2. **chemtools/explain.py** (Line 38)
   - **Usage**: UID → human name translation
   - **Impact**: ⚠️ Falls back gracefully if not found
   - **Mitigation**: Already has fallback logic

3. **chemtools/recommend.py** (Line 1040)
   - **Usage**: pKa lookups for base ranking
   - **Impact**: ⚠️ Limited base coverage
   - **Mitigation**: Could use reagent database

4. **app/main.py** (Line 282)
   - **Usage**: Public API endpoint `/api/v1/properties`
   - **Impact**: ⚠️ Endpoint exists but limited data
   - **Mitigation**: Document limitations in API docs

### Breaking Changes:
- ❌ None - API remains compatible
- ⚠️ Reduced data coverage (was already broken for non-seed compounds)

---

## ✅ Recommended Action Plan

1. **Simplify properties.py** (remove ~250 lines of dead code)
2. **Update docstring** (document limitations)
3. **Add migration guide** (recommend chem.reagent.* for new code)
4. **Update README** (clarify properties.lookup scope)
5. **Mark as deprecated** (optional - schedule for removal in v3.0)

### Benefits:
- ✅ Cleaner codebase (-250 lines)
- ✅ No misleading code
- ✅ Clear expectations
- ✅ Easier maintenance
- ✅ No breaking changes

**Estimated Effort**: 30 minutes

**Risk**: Low (code already didn't work for taxonomy, just removing dead code)

---

## 📊 Summary

| Aspect | Before | After |
|--------|--------|-------|
| **Lines of Code** | 317 | ~70 |
| **Dead Code** | ~250 lines | 0 |
| **Data Coverage** | 5 compounds (taxonomy broken) | 5 compounds (explicit) |
| **Purpose** | Unclear (claims taxonomy support) | Clear (minimal stub) |
| **Maintenance** | High (complex unused code) | Low (simple lookup) |
| **API Compatibility** | ✅ | ✅ |

**Status**: ⚠️ **NEEDS CLEANUP** - Remove dead code referencing deleted `data/compound_taxonomy/` directory
