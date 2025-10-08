# Data Processor File Organization

## Overview

The RDF data processor now saves files to **two separate locations** for better organization:

1. **Working datasets** → `data/reaction_dataset/` (used by chemtools)
2. **Archive/records** → `data-processor/original_dataset/` (for reference)

---

## File Locations

### Working Files (Used by ChemTools)

**Location**: `data/reaction_dataset/`

These files are actively used by the precedent search system:

```
data/reaction_dataset/
├── C_N_Coupling_Cu.jsonl         # Clean reaction data (7 MB)
├── C_N_Coupling_Cu_drfp.npz      # Binary fingerprints (3 MB)
├── Suzuki.jsonl                  # ~70 MB
├── Suzuki_drfp.npz               # ~12 MB
└── ...
```

**Purpose**: Direct consumption by chemtools API and precedent search

**Format**:
- `.jsonl` - Clean JSON Lines without embedded fingerprints
- `_drfp.npz` - Compressed binary DRFP fingerprints

### Archive Files (For Records)

**Location**: `data-processor/original_dataset/{family}/{year_range}/`

These files preserve the original folder structure for reference:

```
data-processor/original_dataset/
├── C_N_Coupling_Cu/
│   └── 2020-2024/
│       └── C_N_Coupling_Cu.md        # Human-readable summary
├── Suzuki/
│   ├── 2020-2022/
│   │   └── Suzuki.md
│   └── 2023-2025/
│       └── Suzuki.md
└── ...
```

**Purpose**: Documentation and record-keeping

**Format**:
- `.md` - Markdown summary with reaction details, statistics, and references

---

## Processing Flow

When you process RDF files from `C:\SciFinder\Suzuki\2023-2025\`:

1. **Extract folder structure**:
   - Parent folder: `Suzuki`
   - Current folder: `2023-2025`

2. **Save working files** to `data/reaction_dataset/`:
   ```
   data/reaction_dataset/Suzuki.jsonl
   data/reaction_dataset/Suzuki_drfp.npz
   ```

3. **Save archive file** to `data-processor/original_dataset/`:
   ```
   data-processor/original_dataset/Suzuki/2023-2025/Suzuki.md
   ```

---

## Example: Processing C_N_Coupling_Cu Dataset

### Input
```
C:\SciFinder\C_N_Coupling_Cu\2020-2024\
├── file1.rdf
├── file2.rdf
└── file3.rdf
```

### Output

**Working files** (used by chemtools):
```
data/reaction_dataset/
├── C_N_Coupling_Cu.jsonl       # 7.35 MB - Clean reaction data
└── C_N_Coupling_Cu_drfp.npz    # 2.60 MB - Binary fingerprints
```

**Archive file** (for reference):
```
data-processor/original_dataset/
└── C_N_Coupling_Cu/
    └── 2020-2024/
        └── C_N_Coupling_Cu.md  # Human-readable summary
```

---

## Benefits of This Organization

### ✅ **Clean Separation**

- **Working data**: Only what chemtools needs (JSONL + binary DRFP)
- **Archive**: Detailed documentation for human review

### ✅ **Preserves Folder Structure**

Original import structure is preserved in `original_dataset/`:
```
original_dataset/
├── Suzuki/2020-2022/
├── Suzuki/2023-2025/
└── C_N_Coupling_Cu/2020-2024/
```

### ✅ **No Duplicates**

- JSONL file saved **only once** to `data/reaction_dataset/`
- No copying or duplicate storage
- Markdown saved to archive location for records

### ✅ **Easy Cleanup**

If you need to regenerate working files:
```bash
# Delete working files
rm data/reaction_dataset/Suzuki.jsonl
rm data/reaction_dataset/Suzuki_drfp.npz

# Archive files remain intact
# Just re-process and they'll be regenerated
```

---

## File Naming Convention

### Working Files (data/reaction_dataset/)

- **Format**: `{family_name}.jsonl` and `{family_name}_drfp.npz`
- **Examples**:
  - `C_N_Coupling_Cu.jsonl`
  - `C_N_Coupling_Cu_drfp.npz`
  - `Suzuki.jsonl`
  - `Suzuki_drfp.npz`

### Archive Files (data-processor/original_dataset/)

- **Format**: `{family_name}/{year_range}/{family_name}.md`
- **Examples**:
  - `C_N_Coupling_Cu/2020-2024/C_N_Coupling_Cu.md`
  - `Suzuki/2023-2025/Suzuki.md`

---

## GUI Behavior

When you run the processor GUI:

1. **Select RDF folder**: e.g., `C:\SciFinder\Suzuki\2023-2025\`

2. **Markdown path auto-fills** to:
   ```
   data-processor\original_dataset\Suzuki\2023-2025\Suzuki.md
   ```

3. **Processing happens automatically**:
   - ✅ JSONL → `data/reaction_dataset/Suzuki.jsonl`
   - ✅ DRFP → `data/reaction_dataset/Suzuki_drfp.npz`
   - ✅ Markdown → `data-processor/original_dataset/Suzuki/2023-2025/Suzuki.md`

4. **Completion message shows all locations**:
   ```
   Successfully processed 150 RDF files with 1200 reactions.
   
   📝 Markdown (records): data-processor\original_dataset\Suzuki\2023-2025\Suzuki.md
   📊 JSONL (chemtools): data\reaction_dataset\Suzuki.jsonl
   🔬 DRFP binary: data\reaction_dataset\Suzuki_drfp.npz
   ```

---

## What Changed from Before

### Before
```
data/reaction_dataset/
├── C_N_Coupling_Cu.jsonl       # 73 MB (with embedded fingerprints)
└── C_N_Coupling_Cu.md          # Documentation

data-processor/original_dataset/
└── C_N_Coupling_Cu/
    └── 2020-2024/
        ├── C_N_Coupling_Cu.jsonl  # DUPLICATE! ❌
        └── C_N_Coupling_Cu.md
```

### After (New)
```
data/reaction_dataset/
├── C_N_Coupling_Cu.jsonl       # 7 MB (clean)
└── C_N_Coupling_Cu_drfp.npz    # 3 MB (binary)

data-processor/original_dataset/
└── C_N_Coupling_Cu/
    └── 2020-2024/
        └── C_N_Coupling_Cu.md  # Only archive ✅
```

**Benefits**:
- ✅ No duplicate JSONL files
- ✅ 86% smaller working files
- ✅ Clean separation of concerns
- ✅ Preserved folder structure for records

---

## Summary

| File Type | Location | Purpose | Format |
|-----------|----------|---------|--------|
| **JSONL** | `data/reaction_dataset/` | ChemTools API | Clean JSON without fingerprints |
| **DRFP** | `data/reaction_dataset/` | Precedent search | Compressed binary NPZ |
| **Markdown** | `data-processor/original_dataset/` | Human reference | Markdown summary |

**Result**: Clean, efficient, and well-organized! 🎯

---

**Last Updated**: October 2025
