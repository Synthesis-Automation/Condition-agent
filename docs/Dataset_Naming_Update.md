# Dataset Naming Update Summary

**Date**: October 6, 2025  
**Update Type**: Standardization to PascalCase convention

---

## Overview

Updated all dataset filenames in `data/reaction_dataset/` to follow consistent PascalCase naming convention, matching the ML model family names used throughout the codebase.

---

## File Renaming Summary

### ✅ Completed Renames

| Old Filename | New Filename | Reaction Type |
|-------------|--------------|---------------|
| `C_N_coupling_Cu_Ullmann.jsonl` | `C_N_Coupling_Cu.jsonl` | Ullmann C-N coupling (Cu) |
| `C_N_coupling_Pd_Buchwald.jsonl` | `C_N_Coupling_Pd.jsonl` | Buchwald-Hartwig C-N coupling (Pd) |
| `C_N_coupling_Ni.jsonl` | `C_N_Coupling_Ni.jsonl` | Nickel C-N coupling |

### ✅ Already Correct

| Filename | Reaction Type |
|----------|---------------|
| `Amide_formation.jsonl` | Amide formation reactions |
| `Suzuki.jsonl` | Suzuki cross-coupling reactions |

---

## Final Dataset Structure

```
data/reaction_dataset/
├── Amide_formation.jsonl          # Amide formation reactions
├── C_N_Coupling_Cu.jsonl          # Cu-catalyzed C-N coupling (Ullmann)
├── C_N_Coupling_Ni.jsonl          # Ni-catalyzed C-N coupling
├── C_N_Coupling_Pd.jsonl          # Pd-catalyzed C-N coupling (Buchwald-Hartwig)
└── Suzuki.jsonl                   # Suzuki cross-coupling reactions
```

---

## Naming Convention

### Standard Format: `{ReactionType}_{Metal}.jsonl`

- **PascalCase with underscores**: `C_N_Coupling_Cu.jsonl`
- **Metal suffix**: `_Cu`, `_Pd`, `_Ni` (when multiple metals exist for same reaction type)
- **Simple names**: `Amide_formation.jsonl`, `Suzuki.jsonl` (when no metal variants)

### Rationale

1. **Consistency with ML Models**: Matches family names in `chemtools/recommend.py`
2. **Code Clarity**: Family aliases map directly to dataset filenames
3. **Maintainability**: Clear, predictable naming pattern
4. **Backward Compatibility**: Old names still work via family aliases in `precedent.py`

---

## Impact on Codebase

### ✅ Files Updated to Use New Names

1. **chemtools/precedent.py**
   - `_dataset_family_map()` updated to return new filenames
   - Family aliases maintain backward compatibility
   - Maps `"Ullmann_CN"` → `"C_N_Coupling_Cu"`
   - Maps `"Buchwald_CN"` → `"C_N_Coupling_Pd"`

2. **chemtools/recommend.py**
   - Family aliases defined: `C_N_Coupling_Cu`, `C_N_Coupling_Pd`, `C_N_Coupling_Ni`
   - `_nice_family_text()` displays friendly names
   - `_FAMILY_ROLE_EXPECTATIONS` uses new names

3. **Test Files**
   - `tests/chemtools/test_rule_feature_builders.py` - Uses `C_N_Coupling_Cu`
   - `tests/chemtools/integrations/test_mcp_tools.py` - Expects new family names

4. **Documentation**
   - `docs/Rule_Based_System_Naming_Update.md` - Complete migration guide
   - `docs/Ni_CN_Coupling_Test_Results.md` - Ni rule database test results

---

## Backward Compatibility

### Legacy Names Still Supported

The codebase maintains **full backward compatibility** via family aliases:

```python
# In chemtools/recommend.py
_FAMILY_ALIASES = {
    "Ullmann_CN": "C_N_Coupling_Cu",
    "Buchwald_CN": "C_N_Coupling_Pd",
    # ... other aliases
}
```

**This means:**
- Old API calls with `"Ullmann_CN"` still work
- Old code referencing legacy names continues to function
- Gradual migration path for external integrations
- No breaking changes for existing users

---

## Verification

### Dataset File Integrity

All dataset files maintain their original content:
- ✅ File contents unchanged (only filenames updated)
- ✅ JSONL structure preserved
- ✅ Reaction data intact
- ✅ All entries validated

### Testing Status

- ✅ Unit tests pass with new naming
- ✅ Integration tests updated
- ✅ ML model loading verified
- ✅ Rule-based system validated
- ✅ API endpoints functional

---

## Related Updates

### Complete Naming Standardization (October 6, 2025)

1. **ML Models**: `C_N_Coupling_Cu`, `C_N_Coupling_Pd`, `C_N_Coupling_Ni`
2. **Rule Databases**: `cn_coupling_cu_db.json`, `cn_coupling_pd_db.json`, `cn_coupling_ni.json`
3. **Datasets**: `C_N_Coupling_Cu.jsonl`, `C_N_Coupling_Pd.jsonl`, `C_N_Coupling_Ni.jsonl`
4. **Code References**: All module functions and tests updated

### Documentation

- `docs/Rule_Based_System_Naming_Update.md` - Rule-based system naming
- `docs/Ni_CN_Coupling_Test_Results.md` - Ni rule database validation
- `docs/Dataset_Naming_Update.md` - This document (dataset renaming)

---

## Migration Guide for External Users

### If You Reference Datasets Directly

**Old Code:**
```python
# Loading datasets by filename
cu_data = load_jsonl("data/reaction_dataset/C_N_coupling_Cu_Ullmann.jsonl")
pd_data = load_jsonl("data/reaction_dataset/C_N_coupling_Pd_Buchwald.jsonl")
ni_data = load_jsonl("data/reaction_dataset/C_N_coupling_Ni.jsonl")
```

**New Code:**
```python
# Use new filenames
cu_data = load_jsonl("data/reaction_dataset/C_N_Coupling_Cu.jsonl")
pd_data = load_jsonl("data/reaction_dataset/C_N_Coupling_Pd.jsonl")
ni_data = load_jsonl("data/reaction_dataset/C_N_Coupling_Ni.jsonl")
```

### If You Use API/Function Calls

**No changes needed!** The API layer handles family aliases automatically:

```python
# These all still work (backward compatible)
recommend("Ullmann_CN", reaction_smiles)  # Old name
recommend("C_N_Coupling_Cu", reaction_smiles)  # New name

# Family aliases map old → new internally
```

---

## Checklist

- [x] Rename `C_N_coupling_Cu_Ullmann.jsonl` → `C_N_Coupling_Cu.jsonl`
- [x] Rename `C_N_coupling_Pd_Buchwald.jsonl` → `C_N_Coupling_Pd.jsonl`
- [x] Rename `C_N_coupling_Ni.jsonl` → `C_N_Coupling_Ni.jsonl`
- [x] Update `chemtools/precedent.py` dataset mapping
- [x] Update `chemtools/recommend.py` family aliases
- [x] Update test files
- [x] Verify all files load correctly
- [x] Run test suite
- [x] Create documentation

---

## Summary

✅ **All dataset files now follow consistent PascalCase naming convention**

- Clear alignment with ML model family names
- Improved code readability and maintainability
- Full backward compatibility maintained
- Comprehensive documentation provided

The dataset naming update is **complete and production-ready**! 🎉

---

**Related Documentation**:
- Rule-based system naming: `docs/Rule_Based_System_Naming_Update.md`
- Ni database validation: `docs/Ni_CN_Coupling_Test_Results.md`
- Repository guidelines: `AGENTS.md`

