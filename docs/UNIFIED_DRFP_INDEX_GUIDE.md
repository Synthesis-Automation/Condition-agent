# Building and Using Unified DRFP Index for Cross-Family Search

## Overview

The unified DRFP index enables **fast cross-family search with DRFP similarity** by combining fingerprints from all reaction families into a single NPZ file. This provides the best of both worlds: broad coverage across families with accurate DRFP-based similarity.

## Quick Start

### Step 1: Build the Unified Index

```bash
# Build unified index from all family NPZ files
python scripts/build_unified_drfp_index.py
```

**Expected output:**
```
Found 5 family NPZ files:
  - Amide_formation_drfp.npz
  - C_N_Coupling_drfp.npz
  - C_O_Coupling_drfp.npz
  - C_S_Coupling_drfp.npz
  - Suzuki_drfp.npz

Building Unified DRFP Index
================================================================================

Loading Amide_formation...
  ✓ Loaded 41427 reactions from Amide_formation
Loading C_N_Coupling...
  ✓ Loaded 15967 reactions from C_N_Coupling
Loading C_O_Coupling...
  ✓ Loaded 6392 reactions from C_O_Coupling
Loading C_S_Coupling...
  ✓ Loaded 8117 reactions from C_S_Coupling
Loading Suzuki...
  ✓ Loaded 50215 reactions from Suzuki

Total reactions: 121990

Saving unified index to data/reaction_dataset/ALL_FAMILIES_drfp.npz...
  ✓ Saved 121990 fingerprints
  ✓ File size: 12.61 MB

✓ Unified DRFP index built successfully!
```

### Step 2: Use Cross-Family Search

```bash
# Cross-family search now uses DRFP automatically!
python app/local_recommendation_cli.py \
  --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" \
  --k 100 \
  --search-all-families
```

**Python API:**
```python
from chemtools import chem

# Cross-family search with DRFP (if unified index exists)
result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=100,
    search_all_families=True  # Automatically uses unified index if available
)
```

## How It Works

### Before: Feature-Based Similarity Only

**Without unified index:**
```
Cross-family search (search_all_families=True)
  ↓
No unified index found
  ↓
DRFP disabled (too slow to compute on-the-fly)
  ↓
Uses feature-based similarity (bin, LG, nuc_class)
  ↓
Results: Fast but lower accuracy
```

### After: DRFP-Based Similarity

**With unified index:**
```
Cross-family search (search_all_families=True)
  ↓
Unified index found at data/reaction_dataset/ALL_FAMILIES_drfp.npz
  ↓
DRFP automatically enabled
  ↓
O(1) fingerprint lookup from unified NPZ file
  ↓
Results: Fast AND accurate!
```

## CLI Commands

### Build Unified Index

```bash
# Build from all families (default)
python scripts/build_unified_drfp_index.py

# Custom output path
python scripts/build_unified_drfp_index.py \
  --output data/my_unified_drfp.npz

# Select specific families
python scripts/build_unified_drfp_index.py \
  --families C_N_Coupling Suzuki Amide_formation

# Force rebuild (overwrite existing)
python scripts/build_unified_drfp_index.py --force

# Quiet mode (less output)
python scripts/build_unified_drfp_index.py --quiet

# Help
python scripts/build_unified_drfp_index.py --help
```

### Test Cross-Family Search

```bash
# Standard cross-family search (uses unified index if available)
python app/local_recommendation_cli.py \
  --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" \
  --k 100 \
  --search-all-families

# Compare with family-specific search
python app/local_recommendation_cli.py \
  --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" \
  --k 50 \
  --family C_N_Coupling
```

## File Structure

```
data/reaction_dataset/
├── C_N_Coupling_drfp.npz          (1.2 MB - family-specific)
├── Suzuki_drfp.npz                (4.2 MB - family-specific)
├── Amide_formation_drfp.npz       (3.0 MB - family-specific)
├── C_O_Coupling_drfp.npz          (455 KB - family-specific)
├── C_S_Coupling_drfp.npz          (580 KB - family-specific)
└── ALL_FAMILIES_drfp.npz          (12.6 MB - UNIFIED INDEX)
    ├── reaction_ids: [121,990 IDs]
    ├── fps: [121,990 fingerprints]
    ├── n_bits: 4096
    ├── radius: 3
    ├── families: ['Amide_formation', 'C_N_Coupling', ...]
    └── family_counts: [41427, 15967, ...]
```

## Performance Comparison

| Mode | Index | DRFP | Time | Accuracy | Use Case |
|------|-------|------|------|----------|----------|
| **Family-specific** | Family NPZ | ✅ Yes | ~1-2s | High | Known reaction type |
| **Cross-family (no unified)** | None | ❌ No | ~2-5s | Medium | Fast exploration |
| **Cross-family (with unified)** | Unified NPZ | ✅ Yes | ~3-8s | High | Accurate exploration |

## API Behavior

### Automatic Detection

The system automatically detects if the unified index exists and enables DRFP accordingly:

```python
from chemtools import chem

# If unified index exists: DRFP enabled automatically
result = chem.recommend.conditions(
    reaction="...",
    search_all_families=True
)
# No warning - DRFP used automatically

# If unified index missing: DRFP disabled with warning
result = chem.recommend.conditions(
    reaction="...",
    search_all_families=True
)
# Warning: "Cross-family search with DRFP disabled (no unified index found).
#           To enable DRFP, run: python scripts/build_unified_drfp_index.py"
```

### Manual Override

You can manually control DRFP usage:

```python
# Force DRFP off (use feature-based similarity)
result = chem.recommend.conditions(
    reaction="...",
    search_all_families=True,
    relax={"use_drfp": False}
)

# Force DRFP on (compute on-the-fly if no index - slow!)
result = chem.recommend.conditions(
    reaction="...",
    search_all_families=True,
    relax={"use_drfp": True}
)
```

## Technical Details

### NPZ File Format

The unified index uses the same format as family-specific NPZ files:

```python
import numpy as np

data = np.load('data/reaction_dataset/ALL_FAMILIES_drfp.npz', allow_pickle=True)

# Keys:
# - reaction_ids: array of str (e.g., "31-006-CAS-23498989", "Suzuki_045")
# - fps: array of uint8 arrays (DRFP fingerprints, 4096 bits = 512 bytes each)
# - n_bits: 4096 (fingerprint bit length)
# - radius: 3 (DRFP radius parameter)
# - families: array of str (family names)
# - family_counts: array of int (reactions per family)
```

### Loading Strategy

The precedent search code (`chemtools/precedent/search.py`) automatically selects the appropriate index:

```python
# For family-specific search (family="C_N_Coupling")
if family_txt is not None:
    # Load family-specific NPZ: C_N_Coupling_drfp.npz
    npz_path = get_drfp_path_for_family(family_txt)
    loader = DRFPLoader(npz_path)

# For cross-family search (family=None)
else:
    # Load unified NPZ: ALL_FAMILIES_drfp.npz
    unified_path = get_unified_drfp_path()
    loader = DRFPLoader(unified_path)

# O(1) fingerprint lookup by reaction_id
fingerprint = loader.get_fingerprint(reaction_id)
```

### Duplicate Handling

If the same reaction appears in multiple family NPZ files, the build script will:
1. Detect duplicates and list them
2. Keep only the last occurrence (from the last family processed)
3. Report the total number of duplicates

**Example output:**
```
⚠️  Found 128 duplicate reaction_ids:
  - 31-614-CAS-35265690 in C_N_Coupling
  - 31-614-CAS-37359248 in C_N_Coupling
  ...
```

## Maintenance

### When to Rebuild

Rebuild the unified index when:
- ✅ New family NPZ files are added
- ✅ Existing family NPZ files are updated
- ✅ Fingerprint parameters change (n_bits, radius)
- ✅ Reaction dataset is regenerated

### Verification

Check if the unified index is up-to-date:

```bash
# Check file timestamp
ls -lh data/reaction_dataset/ALL_FAMILIES_drfp.npz

# Check contents
python -c "
import numpy as np
data = np.load('data/reaction_dataset/ALL_FAMILIES_drfp.npz', allow_pickle=True)
print(f'Reactions: {len(data[\"reaction_ids\"])}')
print(f'Families: {list(data[\"families\"])}')
print(f'Size: {len(data[\"fps\"])} fingerprints')
"
```

### Deletion

To force cross-family search to use feature-based similarity:

```bash
# Remove unified index
rm data/reaction_dataset/ALL_FAMILIES_drfp.npz

# Cross-family search will now use features instead of DRFP
python app/local_recommendation_cli.py --search-all-families --rxn "..."
```

## Troubleshooting

### Issue: "No family DRFP NPZ files found"

**Cause:** Family-specific NPZ files don't exist yet.

**Solution:** Build family NPZ files first:
```bash
python scripts/build_family_drfp_index.py
```

### Issue: "Warning: No reaction_ids found in NPZ"

**Cause:** NPZ file is corrupted or empty.

**Solution:** Rebuild the family NPZ file:
```bash
python scripts/build_family_drfp_index.py --family <family_name> --force
```

### Issue: Cross-family search is slow

**Possible causes:**
1. Unified index doesn't exist → Build it with `python scripts/build_unified_drfp_index.py`
2. DRFP manually disabled → Remove `relax={"use_drfp": False}`
3. Using very high k value → Reduce k (try 50-100 instead of 200+)

### Issue: Low confidence scores

**Possible causes:**
1. Using feature-based similarity (DRFP disabled) → Build unified index
2. Query reaction doesn't match dataset well → Normal for novel reactions
3. Using cross-family search without catalyst filtering → Add catalyst preference

## Summary

### Benefits of Unified Index

✅ **Fast cross-family search** - O(1) fingerprint lookup
✅ **DRFP accuracy** - Full reaction similarity, not just features
✅ **Automatic detection** - System uses it automatically if available
✅ **Easy to build** - Single command rebuilds entire index
✅ **Comprehensive** - Covers all 121,990 reactions across 5 families

### Usage Pattern

```bash
# 1. Build unified index (one-time setup)
python scripts/build_unified_drfp_index.py

# 2. Use cross-family search (DRFP enabled automatically)
python app/local_recommendation_cli.py --search-all-families --rxn "..."

# 3. Rebuild when datasets update
python scripts/build_unified_drfp_index.py --force
```

### Next Steps

1. ✅ Built unified DRFP index
2. ✅ Cross-family search now uses DRFP
3. 🎯 Test with various reactions
4. 📊 Compare results with family-specific search
5. 🚀 Deploy to production

The unified index provides **fast, accurate cross-family recommendations** by combining the coverage of cross-family search with the precision of DRFP similarity! 🎉
