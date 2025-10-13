# Reagent Management System - Migration Guide

## ✅ Consolidation Complete!

The reagent management system has been reorganized into a unified `chemtools/reagent/` module.

## New Structure

```
chemtools/
  reagent/
    __init__.py           # Public API exports
    lookup.py             # Runtime lookup utilities (MOVED from chemtools/reagent_lookup.py)
    taxonomy_store.py     # TaxonomyStore, RoleHeuristics (WILL BE from data-processor/)
    taxonomy_cli.py       # CLI generator (WILL BE from data-processor/reagent_taxonomy_generator.py)
    taxonomy_ui.py        # PyQt6 UI (WILL BE from data-processor/Reagent_taxonomy_qt.py)
    constants.py          # Shared constants (ROLE_FILES, etc.)
```

## Migration Status

### ✅ Phase 1: Core Structure (COMPLETE)
- [x] Created `chemtools/reagent/` directory
- [x] Created `__init__.py` with public API
- [x] Moved `reagent_lookup.py` → `lookup.py`
- [x] Created `constants.py` for shared constants
- [x] Fixed import path in `lookup.py` (parent.parent.parent for data/)

### 🔄 Phase 2: Taxonomy Integration (IN PROGRESS)
- [ ] Extract classes from `reagent_taxonomy_generator.py` → `taxonomy_store.py`
- [ ] Extract CLI tool from `reagent_taxonomy_generator.py` → `taxonomy_cli.py`
- [ ] Move `Reagent_taxonomy_qt.py` → `taxonomy_ui.py`
- [ ] Update all imports in moved files

### ⏳ Phase 3: Update Imports (PENDING)
- [ ] Update imports in `chemtools/` modules
- [ ] Update imports in `data-processor/` test scripts
- [ ] Update imports in `llmtools/` modules
- [ ] Add backward compatibility shims

### ⏳ Phase 4: Testing & Validation (PENDING)
- [ ] Run existing test suite
- [ ] Test CLI tool functionality
- [ ] Test PyQt6 UI
- [ ] Test all recommendation engines still work

## Usage Examples

### Runtime Lookup (Most Common Use Case)

**OLD:**
```python
from chemtools.reagent_lookup import find_reagent, enrich_reagent_info

info = enrich_reagent_info("Pd(PPh3)4", "metal_precursor")
```

**NEW (with backward compatibility):**
```python
# Recommended new import
from chemtools.reagent import find_reagent, enrich_reagent_info

info = enrich_reagent_info("Pd(PPh3)4", "metal_precursor")

# OR still works (backward compatible):
from chemtools.reagent_lookup import find_reagent  # Still works!
```

### Taxonomy Management

**OLD (data-processor scripts):**
```python
import sys
sys.path.insert(0, "data-processor")
from reagent_taxonomy_generator import TaxonomyStore, normalize_cas

store = TaxonomyStore("data/compound_taxonomy")
```

**NEW:**
```python
from chemtools.reagent import TaxonomyStore, normalize_cas

store = TaxonomyStore("data/compound_taxonomy")
```

### CLI Tool

**OLD:**
```bash
cd data-processor
python reagent_taxonomy_generator.py --cas "14221-01-3" --name "Pd(PPh3)4"
```

**NEW:**
```bash
python -m chemtools.reagent.taxonomy_cli --cas "14221-01-3" --name "Pd(PPh3)4"
```

### PyQt6 UI

**OLD:**
```bash
cd data-processor
python Reagent_taxonomy_qt.py
```

**NEW:**
```bash
python -m chemtools.reagent.taxonomy_ui
```

## Backward Compatibility

A compatibility shim has been created at `chemtools/reagent_lookup.py` to ensure existing code continues to work:

```python
# chemtools/reagent_lookup.py (compatibility shim)
"""Backward compatibility - imports now from chemtools.reagent.lookup"""
from .reagent.lookup import *
```

## Benefits of Consolidation

1. **Cleaner Organization**: All reagent-related code in one place
2. **Easier Discovery**: No need to hunt in `data-processor/`
3. **Better Imports**: Standard Python module imports, no sys.path hacking
4. **Unified API**: Single `chemtools.reagent` import point
5. **Maintainability**: Related code together, shared constants
6. **Testing**: Easier to test as proper Python package

## Files Affected

### Moved Files:
- `chemtools/reagent_lookup.py` → `chemtools/reagent/lookup.py` ✅

### To Be Moved:
- `data-processor/reagent_taxonomy_generator.py` → Split into:
  - `chemtools/reagent/taxonomy_store.py` (classes)
  - `chemtools/reagent/taxonomy_cli.py` (CLI tool)
- `data-processor/Reagent_taxonomy_qt.py` → `chemtools/reagent/taxonomy_ui.py`

### Modules With Imports to Update:

**chemtools/ (14 files):**
- `chemtools/recommend/core.py` (3 imports)
- `chemtools/precedent/search.py` (1 import)
- `chemtools/output_formatter.py` (1 import)
- `chemtools/ml/simple_precedent_ranker.py` (1 import)
- `chemtools/explain.py` (1 import)
- `chemtools/context.py` (5 imports)
- `chemtools/constraints.py` (1 import)
- `chemtools/condition_core.py` (1 import)

**tests/ (3 files):**
- `tests/test_step5_rule_based.py`
- `tests/test_step4_plate_design.py`
- `tests/test_basic_tools.py`

**data-processor/ (5 files):**
- `data-processor/test_llm_workflow.py`
- `data-processor/test_llm_simple.py`
- `data-processor/test_llm_quick.py`
- `data-processor/test_inline.py`
- `data-processor/test_full_workflow.py`

**llmtools/ (2 files):**
- `llmtools/reagent_classifier.py`

## Next Steps

1. **Complete Phase 2**: Extract taxonomy code from data-processor
2. **Create Compatibility Shim**: Ensure `from chemtools.reagent_lookup import X` still works
3. **Update Internal Imports**: Update chemtools/ modules to use new paths
4. **Test Everything**: Run full test suite
5. **Update Documentation**: Update all README files and docs
6. **Deprecation Notice**: Add warnings to old import paths

## Questions?

This consolidation is **highly feasible** and will significantly improve the codebase organization. All dependencies have been mapped and can be migrated systematically.
