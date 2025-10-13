# Reagent System Consolidation - Implementation Plan

## ✅ Status: Phase 1 Complete, Ready for Phase 2

This document tracks the complete consolidation of the reagent management system into `chemtools/reagent/`.

## Overview

**Goal**: Move all reagent-related code from `data-processor/` into `chemtools/reagent/` for better organization and easier imports.

**Benefits**:
- ✅ Single source of truth for reagent functionality
- ✅ No more `sys.path` hacking in scripts
- ✅ Standard Python package imports
- ✅ Better testability and maintainability
- ✅ Clearer API surface

## File Structure

```
chemtools/reagent/
├── __init__.py              ✅ DONE - Public API exports
├── constants.py             ✅ DONE - ROLE_FILES, priorities, patterns
├── lookup.py                ✅ DONE - Runtime database lookup (from reagent_lookup.py)
├── taxonomy_store.py        ⏳ TODO - TaxonomyStore, RoleHeuristics classes
├── taxonomy_utils.py        ⏳ TODO - Helper functions (CAS, tokenization, resolvers)
├── taxonomy_cli.py          ⏳ TODO - CLI tool (main() from generator)
└── taxonomy_ui.py           ⏳ TODO - PyQt6 GUI (from Reagent_taxonomy_qt.py)
```

## Implementation Phases

### Phase 1: Foundation ✅ COMPLETE

- [x] Create `chemtools/reagent/` directory
- [x] Create `__init__.py` with public exports
- [x] Create `constants.py` with all role/family constants
- [x] Move `chemtools/reagent_lookup.py` → `chemtools/reagent/lookup.py`
- [x] Fix import paths in `lookup.py` (data directory)
- [x] Create migration documentation

### Phase 2: Extract Taxonomy Core ⏳ IN PROGRESS

**Files to create:**

1. **taxonomy_store.py** - Core classes and data structures
   - `TaxonomyStore` class
   - `RoleHeuristics` class
   - `ROLE_KEYWORDS_RAW` patterns
   - `MANUAL_FAMILY_PATTERNS`
   - `EMBEDDING_FIELD_MAP`

2. **taxonomy_utils.py** - Utility functions
   - CAS handling: `normalize_cas()`, `_valid_cas_checksum()`
   - Tokenization: `tokenize_all()`, `_tokenize_text()`, `_sanitize_text()`
   - Resolvers: `resolve_identity_from_cas()`, `_resolve_via_pubchem()`, `_resolve_via_cactus()`
   - Entry building: `build_entry()`, `build_embedding_text()`, `dedupe_synonyms()`
   - HTTP helpers: `_http_get_json()`, `_http_get_text()`

3. **taxonomy_cli.py** - Command-line tool
   - `main()` function
   - `parse_args()` function  
   - Workflow logic

4. **taxonomy_ui.py** - PyQt6 GUI
   - Direct copy from `Reagent_taxonomy_qt.py`
   - Update imports to use `chemtools.reagent`

### Phase 3: Update Imports Throughout Codebase ⏳ PENDING

**chemtools/ modules (14 files to update):**

```python
# OLD
from . import reagent_lookup
from .reagent_lookup import find_reagent

# NEW
from .reagent import find_reagent, enrich_reagent_info
# OR
from . import reagent
```

Files affected:
- `chemtools/recommend/core.py` (3 imports)
- `chemtools/precedent/search.py` (1 import)
- `chemtools/output_formatter.py` (1 import)
- `chemtools/ml/simple_precedent_ranker.py` (1 import)
- `chemtools/explain.py` (1 import)
- `chemtools/context.py` (5 imports)
- `chemtools/constraints.py` (1 import)
- `chemtools/condition_core.py` (1 import)

**tests/ (3 files):**

```python
# OLD
from chemtools import reagent_lookup
from chemtools.reagent_lookup import find_reagent

# NEW
from chemtools.reagent import find_reagent, enrich_reagent_info
```

Files:
- `tests/test_step5_rule_based.py`
- `tests/test_step4_plate_design.py`
- `tests/test_basic_tools.py`

**data-processor/ scripts (5 files):**

```python
# OLD (with sys.path hacking)
import sys
sys.path.insert(0, "data-processor")
from reagent_taxonomy_generator import normalize_cas, TaxonomyStore

# NEW (clean imports)
from chemtools.reagent import normalize_cas, TaxonomyStore
```

Files:
- `data-processor/test_llm_workflow.py`
- `data-processor/test_llm_simple.py`
- `data-processor/test_llm_quick.py`
- `data-processor/test_inline.py`
- `data-processor/test_full_workflow.py`

**llmtools/ (2 files):**

```python
# Update ROLE references to use chemtools.reagent.constants
from chemtools.reagent.constants import ROLE_FILES
```

Files:
- `llmtools/reagent_classifier.py`

### Phase 4: Remove Old Files ⏳ PENDING

After verifying all imports work:

- [ ] Delete `chemtools/reagent_lookup.py` (moved to `chemtools/reagent/lookup.py`)
- [ ] Mark as deprecated: `data-processor/reagent_taxonomy_generator.py`
- [ ] Mark as deprecated: `data-processor/Reagent_taxonomy_qt.py`
- [ ] Update `data-processor/README.md` to point to new location

### Phase 5: Testing & Validation ⏳ PENDING

- [ ] Run full test suite: `pytest -q`
- [ ] Test CLI tool: `python -m chemtools.reagent.taxonomy_cli --list-families`
- [ ] Test PyQt6 UI: `python -m chemtools.reagent.taxonomy_ui`
- [ ] Test recommendation engines still work
- [ ] Test precedent search still works
- [ ] Verify data-processor scripts work with new imports

### Phase 6: Documentation ⏳ PENDING

- [ ] Update main README.md
- [ ] Update AGENTS.md (repository guidelines)
- [ ] Create `chemtools/reagent/README.md`
- [ ] Add examples to docstrings
- [ ] Update API_DOCUMENTATION.md

## API Design

### Public API (`from chemtools.reagent import ...`)

**Runtime Lookup (most common)**:
```python
from chemtools.reagent import find_reagent, enrich_reagent_info, enrich_conditions

# Find a reagent
info = find_reagent("Pd(PPh3)4", "metal_precursor")

# Enrich with details
details = enrich_reagent_info("PPh3", "ligand")

# Enrich full conditions dict
enriched = enrich_conditions({"catalyst": "Pd(PPh3)4", "base": "K2CO3"})
```

**Taxonomy Management**:
```python
from chemtools.reagent import TaxonomyStore, RoleHeuristics, normalize_cas

# Load taxonomy
store = TaxonomyStore("data/compound_taxonomy")

# Use heuristics
heuristics = RoleHeuristics(store)
role, family, matches = heuristics.infer_family("Pd(PPh3)4", ["tetrakis"])

# Validate CAS
cas = normalize_cas("14221-01-3")
```

**CLI Tool**:
```bash
# Add reagent to taxonomy
python -m chemtools.reagent.taxonomy_cli \
    --cas "14221-01-3" \
    --name "Tetrakis(triphenylphosphine)palladium(0)" \
    --abbr "Pd(PPh3)4" \
    --role metal_precursor \
    --family pd_0_complexes

# List all families
python -m chemtools.reagent.taxonomy_cli --list-families

# Auto-resolve from CAS (hits PubChem)
python -m chemtools.reagent.taxonomy_cli --cas "14221-01-3"
```

**GUI Tool**:
```bash
python -m chemtools.reagent.taxonomy_ui
```

## Migration Commands

### Step 1: Create taxonomy_store.py
```bash
# Extract classes and patterns from reagent_taxonomy_generator.py
# Lines to extract:
# - ROLE_KEYWORDS_RAW (61-168)
# - MANUAL_FAMILY_PATTERNS (170-197)
# - RoleHeuristics class (442-495)
# - TaxonomyStore class (497-620)
# - EMBEDDING_FIELD_MAP (624-730)
# - _format_family_field() (732-755)
# - build_embedding_text() (757-775)
```

### Step 2: Create taxonomy_utils.py
```bash
# Extract utility functions
# Lines to extract:
# - CAS functions (208-232)
# - Tokenization (198-206, 234-242)
# - HTTP helpers (234-273)
# - Resolvers (275-441)
# - Entry building (777-821)
```

### Step 3: Create taxonomy_cli.py
```bash
# Extract CLI tool
# Lines to extract:
# - parse_args() (823-841)
# - main() (843-978)
```

### Step 4: Update all imports
```bash
# Use find-and-replace or run update script
```

## Testing Checklist

After each phase:

- [ ] `pytest tests/test_basic_tools.py::test_reagent_lookup` passes
- [ ] `pytest tests/test_step4_plate_design.py` passes
- [ ] `pytest tests/test_step5_rule_based.py` passes
- [ ] Recommendation engines work: `python scripts/local_recommendation_cli.py --strategy rule`
- [ ] CLI tool works: `python -m chemtools.reagent.taxonomy_cli --list-families`
- [ ] No import errors in any module

## Rollback Plan

If issues arise:

1. Revert to commit before changes
2. OR: Copy old `reagent_lookup.py` back
3. OR: Restore `data-processor/` files

All changes are in new files, so rollback is safe.

## Success Criteria

✅ All tests pass
✅ No import errors
✅ CLI tool works
✅ GUI tool works
✅ Cleaner `chemtools/` structure
✅ No `sys.path` hacking needed
✅ Better discoverability

## Next Steps

1. **Create taxonomy_store.py** with all classes
2. **Create taxonomy_utils.py** with helper functions
3. **Create taxonomy_cli.py** as module with `__main__`
4. **Update imports** in chemtools/
5. **Update imports** in tests/
6. **Update imports** in data-processor/
7. **Test everything**
8. **Remove old files**
9. **Update documentation**

---

**Current Status**: Phase 1 complete ✅  
**Next Action**: Create taxonomy_store.py  
**Estimated Time**: 2-3 hours total for full migration
