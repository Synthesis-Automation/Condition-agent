# ChemTools v2.0 - Reagent Lookup Integration

## Overview

Successfully integrated the `reagent_lookup` module into ChemTools v2.0, replacing the old registry system with a cleaner, more maintainable reagent database system.

## What Changed

### 1. **New Reagent Namespace Added**
**File**: `chemtools/context.py`

Added `ReagentNamespace` class with methods:
```python
chem.reagent.find(name, reagent_type)          # Find reagent by name/CAS/abbr
chem.reagent.enrich(name, reagent_type)        # Get full reagent details
chem.reagent.enrich_conditions(conditions)     # Enrich complete condition set
chem.reagent.list_types()                      # List available reagent databases
```

### 2. **Fixed condition_core.py**
**File**: `chemtools/condition_core.py`

**Before** (using deleted registry):
```python
from .registry import resolve as registry_resolve
res = registry_resolve(uid or name)
```

**After** (using reagent_lookup):
```python
from . import reagent_lookup

def _lookup_reagent(name, uid, reagent_type="metal_precursor"):
    """Lookup reagent using reagent_lookup system."""
    query = uid if uid else name
    result = reagent_lookup.find_reagent(query, reagent_type)
    return result

# Usage:
res = _lookup_reagent(name, uid, "metal_precursor")
res = _lookup_reagent(name, uid, "ligand")
```

### 3. **Reagent Database Structure**
**Location**: `data/reagents/*.json`

Available databases:
- `acid.json` - Carboxylic acids, Lewis acids
- `additive.json` - Reaction additives
- `base.json` - Bases (organic, inorganic)
- `ligand.json` - Ligands for metal catalysis
- `metal_precursor.json` - Metal catalyst precursors
- `oxidant.json` - Oxidizing agents
- `reductant.json` - Reducing agents
- `solvent.json` - Reaction solvents

Each entry contains:
```json
{
  "name": "Triphenylphosphine",
  "abbreviation": ["PPh3"],
  "cas": "603-35-0",
  "smiles": "c1ccc(P(c2ccccc2)c2ccccc2)cc1",
  "inchi_key": "XQZTHPMLTXQOAI-UHFFFAOYSA-N",
  "aliases": ["triphenylphosphane"],
  "roles": {...}
}
```

## Usage Examples

### Find a Reagent
```python
from chemtools import chem

# Find by name
ligand = chem.reagent.find("PPh3", "ligand")
print(ligand["name"])  # → "Triphenylphosphine"
print(ligand["cas"])   # → "603-35-0"

# Find by CAS
base = chem.reagent.find("584-08-7", "base")
print(base["name"])    # → "Potassium carbonate"
```

### Enrich Reagent Info
```python
# Get full details
info = chem.reagent.enrich("K2CO3", "base")
print(info)
# {
#   "name": "Potassium carbonate",
#   "cas": "584-08-7",
#   "abbreviation": "K2CO3",
#   "smiles": "[K+].[K+].[O-]C([O-])=O",
#   "found": True,
#   ...
# }
```

### Enrich Full Conditions
```python
# Enrich a complete condition set
conditions = {
    "catalyst": "Pd(PPh3)4",
    "ligand": "PPh3",
    "base": "K2CO3",
    "solvent": "DMF"
}

enriched = chem.reagent.enrich_conditions(conditions)
# Adds *_details fields with full reagent information:
# {
#   "catalyst": "Pd(PPh3)4",
#   "catalyst_details": {...full info...},
#   "ligand": "PPh3",
#   "ligand_details": {...full info...},
#   ...
# }
```

### List Available Reagent Types
```python
types = chem.reagent.list_types()
print(types)
# ['acid', 'additive', 'base', 'ligand', 'metal_precursor', 
#  'oxidant', 'reductant', 'solvent']
```

## Benefits Over Old Registry

### ✅ **Structured Data**
- JSON databases with consistent schema
- Easy to add new reagents
- Version controllable

### ✅ **Multiple Lookup Methods**
- Name matching
- CAS number lookup
- Abbreviation search
- Alias support
- Partial name matching

### ✅ **Rich Information**
- SMILES structures
- InChI keys  
- CAS numbers
- Multiple names/aliases
- Role information

### ✅ **Performance**
- LRU caching with `@lru_cache`
- Fast in-memory lookups
- No database queries

### ✅ **Maintainability**
- Simple JSON files
- No complex taxonomy system
- Easy to update/extend

## Migration from Old Registry

The old `registry.py` system has been completely removed. All lookups now go through `reagent_lookup`:

| Old (Registry) | New (Reagent Lookup) |
|----------------|----------------------|
| `registry.resolve(name)` | `chem.reagent.find(name, type)` |
| `registry.search(query)` | `chem.reagent.enrich(name, type)` |
| N/A | `chem.reagent.enrich_conditions(conditions)` |

## Files Modified

1. ✅ `chemtools/context.py` - Added `ReagentNamespace` class
2. ✅ `chemtools/condition_core.py` - Replaced registry calls with reagent_lookup
3. ✅ `app/main.py` - Already updated to use ChemTools v2.0
4. ✅ `chemtools/__init__.py` - Clean exports (chem, ChemTools, ResourceConfig)

## Testing

```bash
# Test import
python -c "from chemtools import chem; print(chem.reagent.list_types())"

# Test lookup
python -c "from chemtools import chem; print(chem.reagent.find('PPh3', 'ligand'))"

# Test API
python -c "from app.main import app; print('API works!')"
```

## Next Steps

1. ✅ **Integration Complete** - Reagent namespace working
2. 📝 **Update UI Files** - Migrate `ui_simple.py` and `ui_gradio.py` to ChemTools v2.0
3. 🧪 **End-to-End Testing** - Test full application workflow
4. 📚 **Documentation** - Add reagent examples to guides

## Summary

The reagent lookup system is now fully integrated into ChemTools v2.0:

```python
from chemtools import chem

# All reagent operations under one namespace
chem.reagent.find()
chem.reagent.enrich()
chem.reagent.enrich_conditions()
chem.reagent.list_types()
```

This provides a clean, logical API for accessing reagent databases, replacing the old registry system with something more maintainable and powerful.
