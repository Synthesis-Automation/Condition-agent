# Validation Script Fix

## Issue

The `validate_reagent_db.py` script was using `DEFAULT_REGISTRY_DIR = "reagents"` which is a relative path without the "data/" prefix, causing it to fail when run from the project root.

## Fix Applied

Changed the script to construct the full default path:

```python
# Before
from chemtools.reagent import DEFAULT_REGISTRY_DIR  # Just "reagents"

# After
project_root = Path(__file__).parent.parent
DEFAULT_REGISTRY_PATH = project_root / "data" / "reagents"
```

## Now Works

```bash
# Run from anywhere in the project
python scripts/validate_reagent_db.py

# Output:
# Validating reagent database at: C:\Git-softwares\Condition-agent\data\reagents
# ✅ Works correctly
```

## Usage

```bash
# Default path (data/reagents)
python scripts/validate_reagent_db.py

# Custom path
python scripts/validate_reagent_db.py --registry-dir path/to/other/reagents

# With verbose output
python scripts/validate_reagent_db.py --verbose

# Specific role
python scripts/validate_reagent_db.py --role ligand --verbose
```

## Status

✅ **Fixed and tested**
- Default path now correctly points to `data/reagents`
- Script works when run from project root
- All validation features working properly
