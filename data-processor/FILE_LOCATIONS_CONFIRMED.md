# File Save Locations - Confirmation

## Current Behavior ✅

The data processor now saves files to **exactly the right locations** with **no duplicates**:

### Files Saved

| File Type | Saved To | Count | Purpose |
|-----------|----------|-------|---------|
| **JSONL** | `data/reaction_dataset/{family}.jsonl` | 1 file | ChemTools working data |
| **DRFP binary** | `data/reaction_dataset/{family}_drfp.npz` | 1 file | Binary fingerprints |
| **Markdown** | `data-processor/original_dataset/{family}/{year}/` | 1 file | Archive/records |

### Example Output

Processing `C:\SciFinder\Suzuki\2023-2025\` creates:

```
data/reaction_dataset/
├── Suzuki.jsonl           # ✅ Working data
└── Suzuki_drfp.npz        # ✅ Binary fingerprints

data-processor/original_dataset/
└── Suzuki/
    └── 2023-2025/
        └── Suzuki.md      # ✅ Archive record
```

### What's NOT Created ❌

- ❌ NO Markdown in `data/reaction_dataset/`
- ❌ NO JSONL in `data-processor/original_dataset/`
- ❌ NO duplicate files anywhere

---

## File Organization Summary

### `data/reaction_dataset/` (Working Files)

**Contains**: JSONL + Binary DRFP only  
**Used by**: ChemTools precedent search API  
**Format**: Optimized for fast loading

```
data/reaction_dataset/
├── C_N_Coupling_Cu.jsonl       # 7 MB (clean)
├── C_N_Coupling_Cu_drfp.npz    # 3 MB (binary)
├── Suzuki.jsonl                # ~70 MB
└── Suzuki_drfp.npz             # ~12 MB
```

### `data-processor/original_dataset/` (Archive)

**Contains**: Markdown only  
**Used by**: Human reference  
**Format**: Human-readable documentation

```
data-processor/original_dataset/
├── C_N_Coupling_Cu/
│   └── 2020-2024/
│       └── C_N_Coupling_Cu.md  # ✅ Only one MD file
└── Suzuki/
    ├── 2020-2022/
    │   └── Suzuki.md           # ✅ Only one MD file
    └── 2023-2025/
        └── Suzuki.md           # ✅ Only one MD file
```

---

## Processing Output Message

When processing completes, you'll see:

```
Successfully processed 150 RDF files with 1200 reactions.

📝 Markdown (records): data-processor\original_dataset\Suzuki\2023-2025\Suzuki.md
📊 JSONL (chemtools): data\reaction_dataset\Suzuki.jsonl
🔬 DRFP binary: data\reaction_dataset\Suzuki_drfp.npz
```

**Three files total** - no duplicates! ✅

---

## Confirmed: No Duplicate Markdown Files

✅ Markdown is created **only once** in `data-processor/original_dataset/`  
✅ **NOT** duplicated to `data/reaction_dataset/`  
✅ **NOT** saved to any temporary location  

**Status**: Working exactly as requested! 🎯

---

**Date**: October 2025  
**Verified**: All file paths correct, no duplicates
