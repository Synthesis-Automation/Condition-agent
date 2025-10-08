# ChemTools Cleanup Plan

**Generated**: October 7, 2025  
**Purpose**: Identify and remove obsolete/unneeded files from `/chemtools`

---

## 🗑️ Files to Remove

### 1. **`chemtools/scdb_matcher/`** - Backward Compatibility Shim
**Status**: ❌ **REMOVE** (per user request - backward compatibility not needed)

**Files**:
- `chemtools/scdb_matcher/__init__.py`
- `chemtools/scdb_matcher/loader.py`

**Reason**: 
- Created solely for backward compatibility
- User confirmed backward compatibility not needed
- All imports should use `chemtools.rule_scdb_matcher` directly

**Impact**:
- Will break old imports: `from chemtools.scdb_matcher import ...`
- Need to update any remaining code to use `chemtools.rule_scdb_matcher`

---

### 2. **`chemtools/condition_core.py`** - Legacy Condition Parsing
**Status**: ⚠️ **CANDIDATE FOR REMOVAL** (still has API endpoint)

**Size**: 302 lines

**Current Usage**:
- ✅ **API endpoint**: `/api/v1/condition-core/parse` in `app/main.py`
- ✅ **API endpoint**: `/api/v1/condition-core/validate` in `app/main.py`
- ❌ Used by legacy scripts:
  - `scripts/condition_core_tester.py`
  - `scripts/validate_condition_core_ullmann.py`

**Functionality**:
- Metal/ligand pairing from reagent lists
- Dataset alias resolution
- Uses `reagent_lookup` internally (already available)

**Recommendation**: 
- ⚠️ **Keep for now** - Has active API endpoints
- 🔄 **Deprecate** - Mark as deprecated, plan removal in v3.0
- 📝 **Document** - Add migration guide to use `reagent_lookup` directly

---

### 3. **`chemtools/selector_payloads.py`** - Amide Selector Features
**Status**: ✅ **KEEP** (actively used)

**Size**: 302 lines

**Current Usage**:
- ✅ Used by `scripts/amide_selector_report.py`
- ✅ Referenced in test files
- ✅ Specific functionality for amide formation rules

**Recommendation**: 
- ✅ **Keep** - Active functionality for amide formation
- 📝 **Document** - Add better docstrings

---

### 4. **`chemtools/agent/`** - Deprecated Agent Config
**Status**: ⚠️ **CANDIDATE FOR REMOVAL** (deprecated stubs)

**Files**:
- `chemtools/agent/config.py` - Already marked as deprecated with warnings
- `chemtools/agent/features/mapping.py` - Unclear usage
- `chemtools/agent/__init__.py`

**Current Usage**:
- ❌ No active imports found in codebase
- ⚠️ `config.py` already emits `DeprecationWarning`

**Recommendation**:
- ❌ **Remove** - Already deprecated, emits warnings
- Legacy rule-based MCP integration has been removed

---

### 5. **Old Registry System** - Already Removed
**Status**: ✅ **Already removed** (confirmed `chemtools/registry.py` doesn't exist)

**Note**: References in scripts are to:
- `chemtools.cli.registry` - Does NOT exist (should be in `rule_scdb_matcher/cli.py`)
- `chemtools.features.role.registry` - Different file (feature registry, still used)

---

## 📋 Cleanup Actions Summary

### Immediate Removals ❌

1. **Remove `chemtools/scdb_matcher/`** (2 files)
   - `__init__.py`
   - `loader.py`
   
2. **Remove `chemtools/agent/`** (3+ files)
   - `config.py` (deprecated)
   - `features/mapping.py`
   - `__init__.py`

### Files to Keep ✅

1. **`condition_core.py`** - Keep (has active API endpoints)
2. **`selector_payloads.py`** - Keep (actively used for amide formation)
3. **`features/role/`** - Keep (used by recommend.py and app/main.py)
4. **`integrations/`** - Keep (optional integrations)
5. **`ml/`** - Keep (ML models)
6. **`featurizers/`** - Keep (active featurization)
7. **`util/`** - Keep (utility functions)

---

## 🔄 Required Code Updates

### Update 1: Remove `scdb_matcher` export from `chemtools/__init__.py`

**Before**:
```python
from . import rule_scdb_matcher as scdb_matcher  # backward compatibility
__all__ = ["chem", "ChemTools", "ResourceConfig", "scdb_matcher"]
```

**After**:
```python
# Remove scdb_matcher export
__all__ = ["chem", "ChemTools", "ResourceConfig"]
```

---

### Update 2: Update scripts to use `rule_scdb_matcher` directly

**Files to update**:
- `scripts/registry_resolver.py` - Uses `chemtools.cli.registry` (DOESN'T EXIST)
- `pyproject.toml` - References `chemtools.cli.registry:main` (DOESN'T EXIST)
- `Makefile` - Uses `chemtools.cli.registry` (DOESN'T EXIST)

**Note**: These references are broken - `cli.py` is in `rule_scdb_matcher/`, not `chemtools/cli/`

---

### Update 3: Fix broken `chemtools.cli.registry` references

**Current**: `chemtools.cli.registry` (doesn't exist)  
**Should be**: `chemtools.rule_scdb_matcher.cli`

**Files to fix**:
1. `scripts/registry_resolver.py`:
   ```python
   # Before
   from chemtools.cli.registry import main
   
   # After
   from chemtools.rule_scdb_matcher.cli import main
   ```

2. `pyproject.toml`:
   ```toml
   # Before
   chem-registry = "chemtools.cli.registry:main"
   
   # After
   chem-registry = "chemtools.rule_scdb_matcher.cli:main"
   ```

3. `Makefile`:
   ```makefile
   # Before
   python -m chemtools.cli.registry $(REG_ARGS)
   
   # After
   python -m chemtools.rule_scdb_matcher.cli $(REG_ARGS)
   ```

---

## 📊 Impact Analysis

### Disk Space Saved
- `scdb_matcher/`: ~50 lines total
- `agent/`: ~100+ lines total
- **Total**: ~150-200 lines removed

### Breaking Changes
- ❌ Old imports `from chemtools.scdb_matcher import ...` will break
- ✅ New imports `from chemtools.rule_scdb_matcher import ...` will work
- ✅ ChemTools v2.0 API (`chem.rules.*`) unaffected

### Risk Assessment
- 🟢 **Low Risk**: Removing `scdb_matcher/` (newly created, not in production)
- 🟢 **Low Risk**: Removing `agent/` (already deprecated)
- 🟡 **Medium Risk**: Need to fix broken `cli.registry` references

---

## ✅ Step-by-Step Execution Plan

### Phase 1: Fix Broken References (FIRST)
1. ✅ Update `scripts/registry_resolver.py`
2. ✅ Update `pyproject.toml`
3. ✅ Update `Makefile`
4. ✅ Update `AGENTS.md`

### Phase 2: Remove Backward Compatibility Shim
1. ❌ Delete `chemtools/scdb_matcher/__init__.py`
2. ❌ Delete `chemtools/scdb_matcher/loader.py`
3. ❌ Remove directory `chemtools/scdb_matcher/`
4. 🔄 Update `chemtools/__init__.py` (remove `scdb_matcher` export)

### Phase 3: Remove Deprecated Agent Module
1. ❌ Delete `chemtools/agent/config.py`
2. ❌ Delete `chemtools/agent/features/mapping.py`
3. ❌ Delete `chemtools/agent/features/__init__.py`
4. ❌ Delete `chemtools/agent/__init__.py`
5. ❌ Remove directory `chemtools/agent/`

### Phase 4: Update Documentation
1. 📝 Update `CHEMTOOLS_FILE_SUMMARY.md`
2. 📝 Update `RULE_BASED_INTEGRATION_COMPLETE.md`
3. 📝 Create `CLEANUP_COMPLETE.md` with summary

### Phase 5: Testing
1. ✅ Test imports: `from chemtools import chem`
2. ✅ Test API: `python -c "from app.main import app"`
3. ✅ Test CLI: `python -m chemtools.rule_scdb_matcher.cli --help`
4. ✅ Run pytest

---

## 📝 Files Analysis Summary

### Total Files in `/chemtools`: ~30-40 files

### Files to Remove: 5-7 files
- `scdb_matcher/__init__.py`
- `scdb_matcher/loader.py`
- `agent/config.py`
- `agent/features/mapping.py`
- `agent/features/__init__.py`
- `agent/__init__.py`

### Files to Keep: 25-30 files
- All core operations (smiles, router, properties, constraints)
- All data operations (precedent, recommend, recommend_ml, reagent_lookup, explain)
- All utilities (reaction_similarity, reaction_type_detector, output_formatter)
- `condition_core.py` (has active API endpoints)
- `selector_payloads.py` (actively used)
- All subdirectories: featurizers/, util/, ml/, features/, integrations/, rule_scdb_matcher/

---

## 🎯 Expected Outcome

**Before Cleanup**:
- Unnecessary backward compatibility layer
- Deprecated agent config emitting warnings
- Broken CLI references to non-existent paths

**After Cleanup**:
- ✅ Clean imports (no backward compatibility needed)
- ✅ No deprecated modules
- ✅ Fixed CLI references
- ✅ Smaller, cleaner codebase
- ✅ Clear migration path for users

**Estimated Time**: 15-20 minutes
**Risk Level**: 🟢 Low (mostly removing unused/deprecated code)
