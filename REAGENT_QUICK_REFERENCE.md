# Quick Reference: Reagent System After Cleanup

## Running the Applications

### PyQt6 Taxonomy Manager UI
```bash
# Run directly
python app/reagent_taxonomy_ui.py

# Or from app directory
cd app
python reagent_taxonomy_ui.py
```

### Command-line Taxonomy Tool
```bash
# List all reagent families
python -m chemtools.reagent.taxonomy_cli --list-families

# Add new reagent by CAS number
python -m chemtools.reagent.taxonomy_cli --cas "14221-01-3" --name "Pd(PPh3)4"

# Dry run mode
python -m chemtools.reagent.taxonomy_cli --cas "121-44-8" --dry-run
```

## Using in Your Code

### Runtime Reagent Lookup
```python
from chemtools.reagent import find_reagent, enrich_reagent_info

# Find reagent by name
result = find_reagent("Pd(PPh3)4")

# Enrich with full details
info = enrich_reagent_info("Pd(PPh3)4", "metal_precursor")
print(info["cas"])
print(info["smiles"])
```

### Taxonomy Management
```python
from chemtools.reagent import TaxonomyStore, RoleHeuristics

# Load taxonomy database
store = TaxonomyStore("data/compound_taxonomy")

# Find reagent by CAS
role, family, entry = store.find_by_cas("14221-01-3")

# Auto-detect role/family from name
role, family, tokens = RoleHeuristics.infer_family("Pd(PPh3)4", [])
```

### Utilities
```python
from chemtools.reagent import (
    normalize_cas,
    resolve_identity_from_cas,
    dedupe_synonyms,
    tokenize_all,
)

# Validate and format CAS number
cas = normalize_cas("14221013")  # Returns "14221-01-3"

# Resolve CAS to name/SMILES via PubChem/Cactus
info = resolve_identity_from_cas("121-44-8")
print(info["name"])          # "Triethylamine"
print(info["smiles"])        # "CCN(CC)CC"
print(info["synonyms"])      # ["N,N-Diethylethanamine", ...]

# Deduplicate synonym list
unique = dedupe_synonyms(["TEA", "tea", "Triethylamine", "TEA"])
# Returns: ["TEA", "Triethylamine"]
```

### Constants
```python
from chemtools.reagent import (
    ROLE_FILES,
    DEFAULT_FAMILY_BY_ROLE,
    ROLE_PRIORITY,
)

# Get all reagent roles
roles = list(ROLE_FILES.keys())
# ['acid', 'additive', 'base', 'condensation_agent', ...]

# Get default family for a role
default = DEFAULT_FAMILY_BY_ROLE["metal_precursor"]
# 'pd_ii_salts'

# Check role priority
priority = ROLE_PRIORITY["ligand"]  # 0 (highest)
```

## File Locations

### Applications
- **FastAPI Server**: `app/main.py`
- **Gradio UI**: `app/ui_gradio.py`
- **Simple UI**: `app/ui_simple.py`
- **Taxonomy UI**: `app/reagent_taxonomy_ui.py` ← PyQt6 GUI

### Libraries
- **Reagent Package**: `chemtools/reagent/`
  - `__init__.py` - Public API
  - `constants.py` - Shared constants
  - `lookup.py` - Runtime database lookup
  - `taxonomy_store.py` - TaxonomyStore & RoleHeuristics classes
  - `taxonomy_utils.py` - Utility functions
  - `taxonomy_cli.py` - Command-line tool

### Data
- **Reagent Database**: `data/reagents/`
- **Taxonomy Data**: `data/compound_taxonomy/` (create if needed)

## Import Patterns

### ❌ OLD (No longer works)
```python
# This will fail - file was removed
from chemtools import reagent_lookup
reagent_lookup.find_reagent("...")
```

### ✅ NEW (Current)
```python
# All from one unified package
from chemtools.reagent import find_reagent, TaxonomyStore, normalize_cas

find_reagent("...")
store = TaxonomyStore("...")
cas = normalize_cas("...")
```

## Common Tasks

### Add New Reagent via UI
1. Run: `python app/reagent_taxonomy_ui.py`
2. Enter CAS number
3. Click "Auto-resolve" (fetches name/SMILES from PubChem)
4. Select role and family
5. Fill in additional fields
6. Click "Save"

### Add New Reagent via CLI
```bash
# With auto-resolution
python -m chemtools.reagent.taxonomy_cli --cas "14221-01-3"

# Manual specification
python -m chemtools.reagent.taxonomy_cli \
    --cas "14221-01-3" \
    --name "Tetrakis(triphenylphosphine)palladium(0)" \
    --abbr "Pd(PPh3)4" \
    --role metal_precursor \
    --family pd_0_complexes
```

### Query Reagents in Code
```python
from chemtools.reagent import find_reagent

# Search by name
result = find_reagent("Pd(PPh3)4")
if result:
    print(f"Found: {result['name']}")
    print(f"CAS: {result['cas']}")
    print(f"Role: {result['role']}")
else:
    print("Not found in database")
```

## Troubleshooting

### ModuleNotFoundError: No module named 'chemtools'
**Solution**: Make sure you're running from the project root directory, or the code adds the root to `sys.path`:

```python
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent  # Adjust as needed
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
```

### ImportError: cannot import name 'find_reagent'
**Solution**: Update to use the new package:
```python
# OLD - won't work
from chemtools import reagent_lookup

# NEW - works
from chemtools.reagent import find_reagent
```

### Taxonomy directory not found
**Solution**: Create the directory if you're managing taxonomy data:
```bash
mkdir -p data/compound_taxonomy
```

---

**Last Updated**: October 13, 2025  
**Status**: Production-ready ✅
