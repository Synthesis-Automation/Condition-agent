# Registry Module Removal Summary

## Overview
Removed the obsolete `chemtools/registry.py` module and updated all references throughout the codebase. The old taxonomy-based registry system has been replaced by the modern `chemtools/reagent_lookup.py` which uses individual JSON files in `data/reagents/`.

## Files Removed

### Core Module
- **`chemtools/registry.py`** (581 lines) - Deleted
  - Old taxonomy-based system using `data/compound_taxonomy` directory
  - Functions: `resolve()`, `search()`, `categories()`, `_load_registry()`
  - Class: `_RegistryIndex`

## Files Modified

### API Layer
- **`app/main.py`**
  - Removed import: `from chemtools import registry as creg`
  - Removed endpoint: `GET /api/v1/registry/categories`
  - Removed endpoint: `GET /api/v1/registry/search`
  - Removed startup cache warming: `_reg._load_registry()`

### UI Layer
- **`app/ui_gradio.py`**
  - Commented out import: `from chemtools import registry as creg`
  - Updated `LegacyRegistryBackend` class to return error messages
  - Methods now return empty results instead of calling registry functions

### Package Init
- **`chemtools/__init__.py`**
  - Removed `"registry"` from `__all__` list
  - Updated docstring to remove registry import example

### Data Processors
- **`data-processor/reaction_markdown_generator.py`**
  - Added deprecation warning on import attempt
  - Disabled registry loading in `CASRegistry.load_registry_data()` (line 229)
  - Now skips chemtools.registry section entirely

- **`data-processor/process_reactions.py`**
  - Disabled registry loading in CAS mapping section (line 847)
  - Now skips chemtools.registry taxonomy loading

## Scripts Not Modified (Will Fail if Used)

The following scripts still import registry and will need to be removed or updated if they are actively used:

- **`scripts/list_registry.py`** - Uses `from chemtools.registry import search`
- **`scripts/registry_tester.py`** - Uses `from chemtools.registry import resolve`
- **`scripts/condition_core_tester.py`** - Uses `from chemtools.registry import resolve as reg_resolve`

**Recommendation**: These appear to be old test/development scripts. If they're not actively used, they should be removed in a future cleanup.

## Migration Guide

### For API Users
The following API endpoints have been removed:
- `GET /api/v1/registry/categories` - No replacement
- `GET /api/v1/registry/search` - No replacement

If you need reagent information, use the modern `chemtools.reagent_lookup` module:
```python
from chemtools import reagent_lookup

# Find a reagent
info = reagent_lookup.enrich_reagent_info("PPh3", "ligand")
print(info)

# Get all reagent types
types = reagent_lookup.get_all_reagent_types()
```

### For Code Users
Replace old registry usage:
```python
# OLD (removed)
from chemtools import registry
result = registry.resolve("PPh3")
results = registry.search(q="phosphine", role="LIGAND")
cats = registry.categories()

# NEW (use reagent_lookup)
from chemtools import reagent_lookup

# Find specific reagent
info = reagent_lookup.find_reagent("PPh3", "ligand")

# Enrich with full details
details = reagent_lookup.enrich_reagent_info("PPh3", "ligand")

# Get available types
types = reagent_lookup.get_all_reagent_types()
# Returns: ['ligand', 'base', 'solvent', 'metal_precursor', ...]
```

## Data Directory Structure

### Old System (Removed)
```
data/compound_taxonomy/
  ├── ligand.json
  ├── base.json
  ├── solvent.json
  ├── coupling_reagent.json
  └── catalysts_precursor.json
```

### Modern System (Active)
```
data/reagents/
  ├── ligand.json
  ├── base.json
  ├── solvent.json
  ├── metal_precursor.json
  ├── oxidant.json
  ├── additive.json
  └── acid.json
```

## Benefits of New System

1. **Simpler structure**: Individual JSON files per reagent type
2. **Better performance**: LRU caching in `load_reagent_database()`
3. **More information**: Includes SMILES, InChI keys, roles, aliases
4. **Better matching**: Normalized name matching with partial matching fallback
5. **Enrichment focused**: Designed for enriching condition recommendations

## Testing

All code still compiles after removal. Remaining import errors are unrelated (scdb_matcher module).

To verify the API still works:
```bash
# Start the API
make run

# Test remaining endpoints (registry endpoints should 404)
curl http://localhost:8000/api/v1/cores
curl http://localhost:8000/api/v1/properties/lookup -X POST -H "Content-Type: application/json" -d '{"query": "water"}'
```

## Future Work

As part of the ChemToolsContext resource manager implementation (see `RESOURCE_MANAGER_IMPLEMENTATION_PLAN.md`), reagent databases will be loaded centrally with:
- Lazy loading
- Thread-safe access
- Selective loading (only needed reagent types)
- Unified resource lifecycle management

This will further improve performance and reduce memory usage.
