# Reagent Taxonomy UI - Moved to app/

## Summary

Successfully moved the PyQt6 reagent taxonomy UI from `data-processor/` to `app/` for better organization.

## Changes

### File Location
- **OLD**: `data-processor/Reagent_taxonomy_qt.py`
- **NEW**: `app/reagent_taxonomy_ui.py`

### Updated Imports
The UI now uses the new unified `chemtools.reagent` package:

```python
# Path setup to find chemtools (required when running directly)
MODULE_DIR = Path(__file__).resolve().parent
ROOT_DIR = MODULE_DIR.parent

# Add root directory to path so chemtools can be imported
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# NEW - Clean imports from chemtools.reagent
from chemtools.reagent import (
    dedupe_synonyms,
    normalize_cas,
    resolve_identity_from_cas,
    tokenize_all,
)
```

No more dependency on `reagent_taxonomy_generator.py`!

### Files Updated
1. ✅ `app/reagent_taxonomy_ui.py` - Moved and updated imports
2. ✅ `data-processor/test_llm_workflow.py` - Updated import path
3. ✅ `data-processor/test_llm_simple.py` - Updated to use chemtools.reagent
4. ✅ `data-processor/test_llm_quick.py` - Updated paths and imports
5. ✅ `data-processor/test_inline.py` - Updated to use chemtools.reagent
6. ✅ `data-processor/test_full_workflow.py` - Updated to use chemtools.reagent

## Rationale

### Why Move from data-processor/?

**data-processor/** contains:
- One-time data processing scripts (SciFinder RDF processors)
- Dataset generators and indexers
- Scripts that run once to prepare data

**This UI is different:**
- Interactive PyQt6 application
- Ongoing taxonomy curation tool
- User-facing interface

### Why app/?

**app/** already contains user-facing applications:
- `main.py` - FastAPI server
- `ui_gradio.py` - Gradio web UI
- `ui_simple.py` - Simple web UI

The reagent taxonomy UI fits perfectly here as another user interface.

## Running the UI

### From app/ directory:
```bash
cd app
python reagent_taxonomy_ui.py
```

### From root directory:
```bash
python app/reagent_taxonomy_ui.py
```

### As module:
```bash
python -m app.reagent_taxonomy_ui
```

## Benefits

1. ✅ **Logical Organization**: All UIs in one place
2. ✅ **Clean Imports**: Uses new `chemtools.reagent` package
3. ✅ **No sys.path Hacking**: Proper Python package structure
4. ✅ **Better Discoverability**: Easy to find all applications
5. ✅ **Separation of Concerns**: 
   - `app/` = User-facing applications
   - `data-processor/` = One-time processing scripts
   - `chemtools/` = Core library

## Integration with Chemtools

The UI integrates seamlessly with the new reagent package consolidation:

```python
from chemtools.reagent import (
    # Runtime lookup
    find_reagent,
    enrich_reagent_info,
    
    # Taxonomy management
    TaxonomyStore,
    RoleHeuristics,
    
    # Utilities
    normalize_cas,
    resolve_identity_from_cas,
    dedupe_synonyms,
    tokenize_all,
    
    # Constants
    ROLE_FILES,
    DEFAULT_FAMILY_BY_ROLE,
)
```

## Next Steps (Optional)

Documentation files in `data-processor/` still reference the old location:
- `LLM_ENHANCEMENTS.md`
- `LLM_QUICKSTART.md`
- `MIGRATION_GUIDE.md`
- etc.

These can be updated gradually or when users encounter them.

---

**Date**: October 13, 2025  
**Status**: ✅ Complete and tested  
**Breaking Changes**: Test files updated, main UI works correctly
