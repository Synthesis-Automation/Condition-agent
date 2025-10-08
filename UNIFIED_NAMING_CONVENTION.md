# Unified Naming Convention

## Overview

All dataset files and rule database files now use **consistent PascalCase naming** with underscores for direct 1:1 correlation.

## File Naming Structure

### Pattern
```
Dataset:  {Family}.jsonl
Rule DB:  {Family}_db.json
DRFP:     {Family}_drfp.npz
```

### Current Files

| Family Name | Dataset File | Rule DB File | DRFP File |
|------------|--------------|--------------|-----------|
| **C_N_Coupling_Cu** | `C_N_Coupling_Cu.jsonl` | `C_N_Coupling_Cu_db.json` | `C_N_Coupling_Cu_drfp.npz` |
| **C_N_Coupling_Pd** | `C_N_Coupling_Pd.jsonl` | `C_N_Coupling_Pd_db.json` | `C_N_Coupling_Pd_drfp.npz` |
| **C_N_Coupling_Ni** | `C_N_Coupling_Ni.jsonl` | `C_N_Coupling_Ni_db.json` | `C_N_Coupling_Ni_drfp.npz` |
| **Suzuki** | `Suzuki.jsonl` | `Suzuki_db.json` | `Suzuki_drfp.npz` |
| **Amide_formation** | `Amide_formation.jsonl` | `Amide_formation_db.json` | `Amide_formation_drfp.npz` |

## Benefits

✅ **Direct Correlation**: Rule DB and dataset files match exactly (except suffix)  
✅ **Programmatic Lookup**: Easy to construct file paths: `f"{family}_db.json"`, `f"{family}.jsonl"`  
✅ **Clear Relationships**: Instantly see which files belong together  
✅ **Consistent Naming**: All files use same capitalization and separators  
✅ **Maintainability**: Reduces confusion and errors when adding new reaction types  

## Family Name Mapping

The system automatically maps between detection names and dataset file names:

### Detection Names → Dataset Names

| Detection Output | Dataset Family | Files Loaded |
|-----------------|----------------|--------------|
| `Ullmann_CN` | `C_N_Coupling_Cu` | `C_N_Coupling_Cu.jsonl`, `C_N_Coupling_Cu_db.json` |
| `Buchwald_CN` | `C_N_Coupling_Pd` | `C_N_Coupling_Pd.jsonl`, `C_N_Coupling_Pd_db.json` |
| `Suzuki_CC` | `Suzuki` | `Suzuki.jsonl`, `Suzuki_db.json` |
| `Amide_Coupling` | `Amide_formation` | `Amide_formation.jsonl`, `Amide_formation_db.json` |

This mapping is handled automatically by `chemtools.recommend.utils.canonical_family()`.

## Code Examples

### Loading Rule Database
```python
from chemtools import chem

# Direct loading by family name
db = chem.rules.load_database("C_N_Coupling_Pd_db.json")
```

### Programmatic File Path Construction
```python
from pathlib import Path

family = "C_N_Coupling_Cu"
dataset_dir = Path("data/reaction_dataset")
rule_db_dir = Path("data/conditionDB")

# Construct file paths
dataset_file = dataset_dir / f"{family}.jsonl"
rule_db_file = rule_db_dir / f"{family}_db.json"
drfp_file = dataset_dir / f"{family}_drfp.npz"
```

### Family Name Normalization
```python
from chemtools.recommend.utils import canonical_family

# Automatic mapping from detection names to dataset names
family = canonical_family("Ullmann_CN")  # Returns "C_N_Coupling_Cu"
family = canonical_family("Suzuki_CC")   # Returns "Suzuki"
family = canonical_family("Buchwald_CN") # Returns "C_N_Coupling_Pd"
```

## Migration from Old Names

### Old → New Mapping

| Old Rule DB Name | New Rule DB Name |
|------------------|------------------|
| `cn_coupling_cu_db.json` | `C_N_Coupling_Cu_db.json` |
| `cn_coupling_pd_db.json` | `C_N_Coupling_Pd_db.json` |
| `cn_coupling_ni.json` | `C_N_Coupling_Ni_db.json` |
| `suzuki_db.json` | `Suzuki_db.json` |
| `amide_formation_db.json` | `Amide_formation_db.json` |

### Files Updated

The following files were updated to use the new naming convention:

- ✅ `app/ui_simple.py` - RULE_DATABASES paths
- ✅ `app/ui_gradio.py` - SCDB_DB_PATH constants
- ✅ `app/README_SIMPLE_UI.md` - Documentation examples
- ✅ `chemtools/context.py` - Documentation examples
- ✅ `chemtools/recommend/utils.py` - _FAMILY_ALIASES mapping
- ✅ `chemtools/precedent/loader.py` - _dataset_family_map function

### Breaking Changes

⚠️ **None** - The system maintains backward compatibility through aliasing:

- Old lowercase names still work in code (mapped automatically)
- Rule DB files have been physically renamed
- Dataset files unchanged (already using correct names)

## Validation

Run the validation script to ensure naming is correct:

```bash
python test_unified_naming.py
```

Expected output:
```
✅ All expected rule DBs present (5 files)
✅ All expected datasets present (5 files)
🎉 All files properly aligned!
✅ All family mappings correct!
```

## Adding New Reaction Types

When adding a new reaction type, follow this naming pattern:

1. **Choose a family name** (e.g., `Heck_Coupling`)
2. **Create dataset file**: `data/reaction_dataset/Heck_Coupling.jsonl`
3. **Create rule DB file**: `data/conditionDB/Heck_Coupling_db.json`
4. **Add mapping** in `chemtools/recommend/utils.py`:
   ```python
   _FAMILY_ALIASES: Dict[str, str] = {
       "Heck_CC": "Heck_Coupling",
       "Heck_Coupling": "Heck_Coupling",
       # ... other mappings
   }
   ```
5. **Add to loader** in `chemtools/precedent/loader.py`:
   ```python
   if tl in {"heck", "heck coupling", "heck_coupling"}:
       return "Heck_Coupling"
   ```

## Directory Structure

```
data/
├── conditionDB/                   # Rule databases
│   ├── Amide_formation_db.json
│   ├── C_N_Coupling_Cu_db.json
│   ├── C_N_Coupling_Ni_db.json
│   ├── C_N_Coupling_Pd_db.json
│   └── Suzuki_db.json
│
└── reaction_dataset/              # ML precedent datasets
    ├── Amide_formation.jsonl
    ├── Amide_formation_drfp.npz
    ├── C_N_Coupling_Cu.jsonl
    ├── C_N_Coupling_Cu_drfp.npz
    ├── C_N_Coupling_Ni.jsonl
    ├── C_N_Coupling_Ni_drfp.npz
    ├── C_N_Coupling_Pd.jsonl
    ├── C_N_Coupling_Pd_drfp.npz
    ├── Suzuki.jsonl
    └── Suzuki_drfp.npz
```

## Summary

The unified naming convention ensures:
- **1:1 mapping** between rule DB and dataset files
- **Consistent capitalization** (PascalCase with underscores)
- **Easy programmatic access** (`{family}.jsonl`, `{family}_db.json`, `{family}_drfp.npz`)
- **Clear organization** with standardized patterns
- **Backward compatibility** through automatic aliasing

This makes the codebase easier to maintain and reduces the chance of naming-related bugs.
