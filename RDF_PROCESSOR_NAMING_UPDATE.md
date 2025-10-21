# RDF Processor Naming Update Summary

## Changes Implemented (Option 2: Replace Old with New)

### 1. RDF Processor Updates

**Modified:** `data-processor/Scifinder_rdf_processer.py`

- **Removed `_combined` suffix** from output filenames
- Output pattern: `{category}.jsonl` instead of `{category}_combined.jsonl`
- Updated UI labels to reflect simplified naming
- Files still include all subfolders (recursive scanning preserved)

**Key changes:**
```python
# OLD: output_jsonl = os.path.join(dataset_dir, final_name + "_combined.jsonl")
# NEW: output_jsonl = os.path.join(dataset_dir, final_name + ".jsonl")

# OLD: output_md = os.path.join(norm_folder, final_name + "_combined.md")
# NEW: output_md = os.path.join(norm_folder, final_name + ".md")
```

### 2. File Replacements

Renamed existing combined files to standard names:

| Old Name | New Name | Status |
|----------|----------|--------|
| `C_N_Coupling_combined.jsonl` | `C_N_Coupling.jsonl` | ✅ Replaced |
| `C_N_Coupling_combined_drfp.npz` | `C_N_Coupling_drfp.npz` | ✅ Replaced |

### 3. Current Dataset State

All dataset files now use standard naming convention:

**JSONL Files (data/reaction_dataset/):**
- `Amide_formation.jsonl` (73.64 MB)
- `C_N_Coupling.jsonl` (37.11 MB) ← **Contains combined data from all subfolders**
- `C_O_Coupling.jsonl` (9.30 MB)
- `C_S_Coupling.jsonl` (11.41 MB)
- `Suzuki.jsonl` (79.62 MB)

**DRFP Index Files:**
- `Amide_formation_drfp.npz` (2.89 MB)
- `C_N_Coupling_drfp.npz` (1.79 MB) ← **Contains DRFP fingerprints for all reactions**
- `C_O_Coupling_drfp.npz` (0.43 MB)
- `C_S_Coupling_drfp.npz` (0.55 MB)
- `Suzuki_drfp.npz` (4.04 MB)

## Compatibility Verification

### ✅ Chemtools Compatibility

1. **precedent/loader.py** - Uses `_iter_dataset_files()` to load ALL `.jsonl` files
   - ✅ Works correctly (loads by file extension, not hardcoded names)

2. **dataset_analytics.py** - Uses `get_data_dir() / f"{family}.jsonl"`
   - ✅ Works correctly (expects standard naming like `C_N_Coupling.jsonl`)

3. **DRFP similarity search** - Uses `{family}_drfp.npz` files
   - ✅ Works correctly (standard naming convention)

### ✅ No Breaking Changes

- All chemtools modules use the standard `{family}.jsonl` naming pattern
- No hardcoded references to `_combined` suffix found
- Existing code expects exactly this naming convention

## Benefits

1. **Full Chemtools Compatibility** - All analytics and search functions work without modification
2. **Simplified Naming** - Clean, standard names without suffixes
3. **Recursive Processing** - Still processes all subfolders automatically
4. **Automatic Replacement** - New processing overwrites old dataset files seamlessly
5. **Backward Compatible** - Works with both simple folders and subfolder hierarchies

## Future RDF Processing

When processing RDF files in the future:

1. **Select category folder** (e.g., `C_N_Coupling/`, `Suzuki/`)
2. **All subfolders are automatically included** (2020-2022, 2023-2025, etc.)
3. **Output files:**
   - `data/reaction_dataset/{category}.jsonl` (for chemtools)
   - `{selected_folder}/{category}.md` (for records)
4. **Old files are automatically replaced** with new combined data

## Documentation Updated

- ✅ `RDF_PROCESSOR_RECURSIVE_UPDATE.md` - Updated all naming examples
- ✅ `data-processor/Scifinder_rdf_processer.py` - Updated comments and UI labels

---

**Status:** ✅ Complete and tested  
**Date:** October 21, 2025  
**Impact:** No breaking changes - seamless transition to unified naming
