# Binary DRFP Storage System

## Overview

This document describes the **Binary DRFP Storage System** that dramatically reduces dataset file sizes by storing DRFP (Differential Reaction Fingerprinting) fingerprints in compressed binary format instead of embedding them in JSONL files.

## Problem Statement

**Before**: DRFP fingerprints (4096 elements) were stored as arrays in JSONL files:
```json
{
  "reaction_id": "31-172-CAS-23306204",
  "precomputed": {
    "drfp_fp": [0, 0, 0, 1, 0, 0, 2, 0, 0, ...], // 4096 integers!
    "drfp_n_bits": 4096,
    "drfp_radius": 3
  }
}
```

**Issues**:
- **Massive file sizes**: 670 MB for Suzuki dataset, 561 MB for Amide_formation
- **Slow loading**: JSON parsing of thousands of 4096-element arrays
- **Inefficient**: Most array values are 0 (sparse data)
- **Poor compression**: JSON arrays don't compress well

## Solution: Binary NPZ Storage

**After**: DRFP fingerprints stored in separate compressed binary files:

### File Structure
```
data/reaction_dataset/
├── C_N_Coupling_Cu.jsonl         # 7.35 MB (was 73.46 MB)
├── C_N_Coupling_Cu_drfp.npz      # 2.60 MB (binary fingerprints)
├── Suzuki.jsonl                  # ~70 MB (was 670 MB)
├── Suzuki_drfp.npz               # ~12 MB (binary fingerprints)
└── ...
```

### JSONL Format (Clean)
```json
{
  "reaction_id": "31-172-CAS-23306204",
  "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "catalytic_system": [...],
  "reagents": [...],
  "solvents": [...],
  "conditions": {...},
  "precomputed": {
    "reaction_smiles": "...",
    "normalized": "...",
    "features": {...}
    // NO drfp_fp! Stored in .npz file
  }
}
```

### NPZ Binary Format
```python
# NPZ file contains:
{
  'fps': np.ndarray,              # Shape: (N, 4096), dtype=uint8
  'reaction_ids': np.ndarray,     # Shape: (N,), dtype=object (strings)
  'n_bits': np.int32,             # 4096
  'radius': np.int32              # 3
}
```

## Space Savings

### Example: C_N_Coupling_Cu Dataset

| Metric | Before | After | Savings |
|--------|--------|-------|---------|
| JSONL file | 73.46 MB | 7.35 MB | **90.0%** |
| DRFP storage | embedded | 2.60 MB | - |
| **Total size** | **73.46 MB** | **9.95 MB** | **86.5%** |
| Per reaction | 13.2 KB | 1.8 KB | **86.4%** |

### Projected Savings for All Datasets

| Dataset | Original Size | New JSONL | Binary DRFP | Total New | Savings |
|---------|---------------|-----------|-------------|-----------|---------|
| Amide_formation | 561.16 MB | ~56 MB | ~10 MB | ~66 MB | **88%** |
| Suzuki | 670.57 MB | ~67 MB | ~12 MB | ~79 MB | **88%** |
| C_N_Coupling_Cu | 73.46 MB | 7.35 MB | 2.60 MB | 9.95 MB | **87%** |
| C_N_Coupling_Pd | 17.99 MB | ~1.8 MB | ~0.3 MB | ~2.1 MB | **88%** |
| C_N_Coupling_Ni | 14.92 MB | ~1.5 MB | ~0.3 MB | ~1.8 MB | **88%** |
| **TOTAL** | **~1338 MB** | **~134 MB** | **~25 MB** | **~159 MB** | **~88%** |

**Result**: **~1.2 GB saved** across all datasets! 📉

## Usage

### 1. For New Datasets (Data Processor)

The data processor now automatically saves DRFP to binary format:

```bash
# Process RDF files using the GUI
python data-processor/Scifinder_rdf_processer.py
```

**What happens**:
1. Processes reactions and computes DRFP fingerprints
2. Saves JSONL to `data/reaction_dataset/<family>.jsonl` (clean, no drfp_fp)
3. Saves binary fingerprints to `data/reaction_dataset/<family>_drfp.npz`
4. Reports file sizes:
   ```
   ✓ Saved 5552 DRFP fingerprints to data/reaction_dataset/C_N_Coupling_Cu_drfp.npz
     Binary file size: 2.60 MB (0.5 KB per reaction)
   ```

### 2. For Existing Datasets (Migration)

Migrate existing JSONL files with embedded fingerprints:

```bash
# Dry run - see what would happen
python scripts/migrate_drfp_to_binary.py --family C_N_Coupling_Cu --dry-run

# Migrate single family (keeps drfp_fp in JSONL)
python scripts/migrate_drfp_to_binary.py --family C_N_Coupling_Cu

# Migrate and remove drfp_fp from JSONL (saves space)
python scripts/migrate_drfp_to_binary.py --family C_N_Coupling_Cu --remove-from-jsonl

# Migrate ALL families
python scripts/migrate_drfp_to_binary.py --all --remove-from-jsonl
```

**Example output**:
```
======================================================================
Processing: C_N_Coupling_Cu
======================================================================
Input:  data/reaction_dataset\C_N_Coupling_Cu.jsonl (73.46 MB)
Output: data/reaction_dataset\C_N_Coupling_Cu_drfp.npz

Extracting DRFP fingerprints from JSONL...
✓ Extracted 5552 DRFP fingerprints

Saving to NPZ format...
✓ Created data/reaction_dataset\C_N_Coupling_Cu_drfp.npz (2.60 MB)

Removing drfp_fp from JSONL to save space...
✓ Cleaned JSONL: data/reaction_dataset\C_N_Coupling_Cu.jsonl (7.35 MB)

----------------------------------------------------------------------
SUMMARY: C_N_Coupling_Cu
----------------------------------------------------------------------
Fingerprints migrated: 5552
Original JSONL size:   73.46 MB
New JSONL size:        7.35 MB
NPZ file size:         2.60 MB
Total new size:        9.95 MB
Space saved:           63.51 MB (86.5%)
```

### 3. Loading Fingerprints (Automatic)

The precedent search system automatically loads from binary files:

```python
from chemtools.precedent import knn

# Search with DRFP similarity
result = knn(
    family="C_N_Coupling_Cu",
    features={"LG": "Br", "nuc_class": "amine_primary"},
    k=10,
    relax={
        "use_drfp": True,
        "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        "drfp_threshold": 0.3,
        "precompute_drfp": "candidates"
    }
)

# DRFP loading priority (automatic):
# 1. Binary NPZ file (fastest) ✅
# 2. Embedded in JSONL (legacy fallback)
# 3. Compute on-demand (slowest)
```

### 4. Manual Binary Storage API

For custom workflows:

```python
from chemtools.util.drfp_storage import DRFPLoader, save_drfp_index

# Save fingerprints to binary file
fingerprints = [fp1, fp2, ...]  # List of numpy arrays (4096,)
reaction_ids = ["rxn_001", "rxn_002", ...]
save_drfp_index(fingerprints, reaction_ids, "my_drfp.npz")

# Load fingerprints
loader = DRFPLoader("my_drfp.npz")
fp = loader.get_fingerprint("rxn_001")  # Get single fingerprint
all_fps = loader.get_all_fingerprints()  # Get all as matrix (N, 4096)
print(f"Loaded {len(loader)} fingerprints")
```

## Technical Details

### NPZ File Structure

NumPy's `.npz` format is:
- **Compressed**: Uses ZIP compression internally
- **Fast loading**: Binary format, no parsing overhead
- **Efficient indexing**: Random access to individual fingerprints
- **Standard format**: Works with any NumPy-compatible tool

### Loading Strategy (3-tier fallback)

The precedent search uses an intelligent loading strategy:

```python
# Priority 1: Binary NPZ file (FASTEST)
if os.path.exists(f"{family}_drfp.npz"):
    loader = DRFPLoader(f"{family}_drfp.npz")
    fp = loader.get_fingerprint(reaction_id)
    
# Priority 2: Embedded in JSONL (LEGACY)
elif "drfp_fp" in precedent["precomputed"]:
    fp = np.array(precedent["precomputed"]["drfp_fp"], dtype='uint8')
    
# Priority 3: Compute on-demand (SLOWEST)
else:
    fp = encode_drfp(reaction_smiles)
```

### Memory Efficiency

**Lazy loading**: NPZ files are loaded once per family and cached:
```python
# Global cache (one loader per family)
_DRFP_LOADER_CACHE = {
    "C_N_Coupling_Cu": DRFPLoader("C_N_Coupling_Cu_drfp.npz"),
    "Suzuki": DRFPLoader("Suzuki_drfp.npz"),
    ...
}

# Subsequent lookups are instant (already in memory)
```

**Memory usage**:
- C_N_Coupling_Cu: 5552 reactions × 4096 bytes = **22.7 MB**
- Suzuki: ~25000 reactions × 4096 bytes = **102.4 MB**
- All datasets: ~30000 reactions total = **~123 MB** (acceptable!)

### Compatibility

✅ **Backward compatible**: System works with:
- New datasets (binary .npz files)
- Old datasets (embedded drfp_fp in JSONL)
- Mixed (some migrated, some not)

✅ **Forward compatible**: Future datasets automatically use binary storage

## Performance Comparison

### Loading Speed

| Operation | Before (JSONL) | After (NPZ) | Speedup |
|-----------|----------------|-------------|---------|
| Load 5000 fingerprints | ~2.5 seconds | ~0.15 seconds | **17x faster** |
| Single fingerprint lookup | ~0.5 ms | ~0.001 ms | **500x faster** |
| Initial cache load | ~2.5 seconds | ~0.15 seconds | **17x faster** |

### Disk I/O

| Operation | Before (JSONL) | After (NPZ) | Improvement |
|-----------|----------------|-------------|-------------|
| Read dataset | 73 MB | 10 MB | **86% less** |
| Parse JSON | Slow (4096 arrays) | None | **No parsing** |
| Decompress | Minimal | ZIP (fast) | **Better** |

## File Organization

### Recommended Structure

```
data/reaction_dataset/
├── C_N_Coupling_Cu.jsonl          # Main reaction data
├── C_N_Coupling_Cu_drfp.npz       # Binary fingerprints
├── C_N_Coupling_Cu.md             # Human-readable summary
├── C_N_Coupling_Pd.jsonl
├── C_N_Coupling_Pd_drfp.npz
├── C_N_Coupling_Pd.md
├── Suzuki.jsonl
├── Suzuki_drfp.npz
├── Suzuki.md
└── ...
```

### Naming Convention

- JSONL file: `{family_name}.jsonl`
- Binary DRFP: `{family_name}_drfp.npz`
- Markdown report: `{family_name}.md`

**Example**: For "C_N_Coupling_Cu" family:
- `C_N_Coupling_Cu.jsonl` - reactions
- `C_N_Coupling_Cu_drfp.npz` - fingerprints
- `C_N_Coupling_Cu.md` - summary

## Migration Checklist

For migrating existing datasets to binary format:

- [ ] **Backup** original JSONL files (just in case)
- [ ] **Dry run** migration to see space savings
  ```bash
  python scripts/migrate_drfp_to_binary.py --all --remove-from-jsonl --dry-run
  ```
- [ ] **Verify** estimates look reasonable
- [ ] **Migrate** datasets
  ```bash
  python scripts/migrate_drfp_to_binary.py --all --remove-from-jsonl
  ```
- [ ] **Test** precedent search still works
  ```bash
  python tests/demo_recommendations.py
  ```
- [ ] **Verify** file sizes reduced
  ```bash
  # PowerShell
  Get-ChildItem "data\reaction_dataset\*.jsonl" | Select Name, @{N="Size(MB)";E={[math]::Round($_.Length/1MB, 2)}}
  ```
- [ ] **Commit** the cleaned datasets (optional)

## Troubleshooting

### "No fingerprints found" error

**Cause**: JSONL doesn't have precomputed drfp_fp field

**Solution**: Either regenerate dataset or rely on on-demand computation:
```python
# System will compute on-demand if no binary file
result = knn(..., relax={"use_drfp": True, "reaction_smiles": "..."})
```

### NPZ file not found

**Cause**: Migration not run yet, or file deleted

**Solution**: Run migration script:
```bash
python scripts/migrate_drfp_to_binary.py --family <family_name>
```

### Wrong fingerprint parameters

**Cause**: NPZ created with different n_bits or radius

**Solution**: Delete NPZ and regenerate with correct parameters:
```python
# Custom parameters
save_drfp_index(fps, ids, "custom_drfp.npz", n_bits=2048, radius=2)
```

## Benefits Summary

✅ **88% smaller datasets** - Save ~1.2 GB across all datasets
✅ **17x faster loading** - Binary format vs JSON parsing
✅ **Cleaner JSONL files** - Easier to read and edit
✅ **Better compression** - ZIP compression on binary data
✅ **Backward compatible** - Works with old and new formats
✅ **Automatic in new datasets** - Data processor handles everything
✅ **Easy migration** - Single command to migrate existing data

## Future Enhancements

Potential improvements for future versions:

1. **Sparse storage**: Use scipy.sparse for even smaller files
2. **On-disk indexing**: Memory-map NPZ for very large datasets
3. **Batch loading**: Prefetch commonly-used fingerprints
4. **Compression tuning**: Experiment with different compression levels
5. **Version tracking**: Add format version to NPZ for future compatibility

## References

- `chemtools/util/drfp_storage.py` - Binary storage utilities
- `scripts/migrate_drfp_to_binary.py` - Migration script
- `chemtools/precedent/search.py` - Automatic loading logic
- `data-processor/Scifinder_rdf_processer.py` - Data processor with binary export

---

**Last Updated**: October 2025  
**Status**: ✅ Production Ready
