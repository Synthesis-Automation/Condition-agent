# Why Cross-Family Search Can't Use Precomputed DRFP Fingerprints

## Current Storage Structure (Family-Specific)

The precomputed DRFP fingerprints are stored in **separate NPZ files per family**:

```
data/reaction_dataset/
├── C_N_Coupling_drfp.npz          (1.2 MB - Cu-based C-N coupling reactions)
├── Suzuki_drfp.npz                (4.2 MB - Pd-based Suzuki reactions)
├── Amide_formation_drfp.npz       (3.0 MB - Amide formation reactions)
├── C_O_Coupling_drfp.npz          (455 KB - C-O coupling reactions)
└── C_S_Coupling_drfp.npz          (580 KB - C-S coupling reactions)
```

### How It Works (Standard Family-Specific Search)

When you search within a **specific family** (e.g., `family="C_N_Coupling"`):

```python
# Standard search - FAST ✅
result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    reaction_type="C_N_Coupling",  # Specific family
    k=50
)
```

**What happens:**
1. Load only `C_N_Coupling_drfp.npz` (~1.2 MB)
2. Load only `C_N_Coupling.jsonl` reactions (~600 reactions)
3. For each reaction, look up its fingerprint in the NPZ file (O(1) lookup by reaction_id)
4. Compute similarity: `tanimoto(query_fp, candidate_fp)` ✅ Fast!

**Performance:** ~1-2 seconds

---

## The Cross-Family Problem

When you search **across all families** (`search_all_families=True`):

```python
# Cross-family search - SLOW without fix ❌
result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    search_all_families=True,  # All families
    k=100
)
```

**The problem:**
1. Load **ALL** reaction JSONL files (~3000+ reactions total)
2. Each reaction comes from a different family (Suzuki, C-N, Amide, etc.)
3. To look up fingerprints, we need to know which NPZ file to check
4. But the code path for cross-family search doesn't track family per reaction

### Why NPZ Lookup Fails for Cross-Family

Looking at the code in `chemtools/precedent/search.py` (lines 278-299):

```python
# STRATEGY 1: Try to load from binary NPZ file
# Skip for cross-family search since NPZ files are family-specific
if _drfp_storage_available and get_drfp_path_for_family is not None and family_txt is not None:
    # ❌ Problem: family_txt is None for cross-family search!
    # We can't call get_drfp_path_for_family(None)
    # We don't know which NPZ file contains this reaction's fingerprint
    
    reaction_id = r.get("reaction_id")
    if reaction_id:
        if family_txt not in _DRFP_LOADER_CACHE:
            npz_path = get_drfp_path_for_family(family_txt)  # ❌ Requires family!
            ...
```

**The issue:**
- Each reaction has a `reaction_id` like `"C_N_Coupling_001"` or `"Suzuki_045"`
- The NPZ file path is determined by family: `get_drfp_path_for_family("C_N_Coupling")` → `"C_N_Coupling_drfp.npz"`
- For cross-family search, `family_txt=None`, so we can't determine which NPZ file to load
- We'd need to search ALL NPZ files for each reaction_id (very inefficient)

### What Happens Without Precomputed FPs

When DRFP is enabled but fingerprints aren't precomputed:

```python
# STRATEGY 3: Fall back to computing on-demand
if r_fp is None:
    r_rsmi = r.get("reaction_smiles")
    if r_rsmi:
        r_fp = rs.encode_drfp_cached(r_rsmi, n_bits=4096, radius=3)
        # ❌ This computes DRFP on-the-fly using RDKit + machine learning
        # Takes ~50-100ms PER reaction
        # For 3000 reactions: 150-300 seconds = 2.5-5 minutes! 🐌
```

**Performance:** ~30-120 seconds (too slow!)

---

## Current Solution: Disable DRFP for Cross-Family

To make cross-family search usable, we **disable DRFP by default**:

```python
# In chemtools/recommend/modules/recommender.py (lines 156-167)
if search_all_families and "use_drfp" not in relax:
    relax.setdefault("use_drfp", False)  # ✅ Disable DRFP for cross-family
    warnings.warn(
        "Cross-family search with DRFP disabled (no precomputed fingerprints available). "
        "Results will use feature-based similarity only."
    )
```

### Feature-Based Similarity (Fallback)

Instead of DRFP, we use **feature-based similarity** from `chemtools/precedent/similarity.py`:

```python
def _similarity(query_features: Dict, candidate_features: Dict) -> float:
    """
    Compute similarity based on chemical features:
    - bin: Substrate/electrophile bins (e.g., "LG:Br|NUC:aniline")
    - LG: Leaving group (Br, Cl, I, OTf, etc.)
    - nuc_class: Nucleophile class (aniline, alcohol, thiol, etc.)
    """
    # Exact bin match: score = 1.0
    # LG match: score = 0.6
    # nuc_class match: score = 0.4
    # No match: score = 0.0
```

**Pros:**
- ✅ Fast (no fingerprint computation)
- ✅ Works across all families
- ✅ Good for similar substrate classes

**Cons:**
- ❌ Lower accuracy for structurally different reactions
- ❌ Can't handle novel transformations as well
- ❌ Confidence scores may be lower

**Performance:** ~2-5 seconds ✅

---

## Future Solution: Unified Cross-Family DRFP Index

To enable fast DRFP-based cross-family search, we would need:

### Option 1: Single Unified NPZ File

Create one large NPZ file containing ALL fingerprints:

```
data/reaction_dataset/
└── ALL_FAMILIES_drfp.npz  (9.5 MB - all reactions combined)
    ├── C_N_Coupling_001 → [fingerprint]
    ├── C_N_Coupling_002 → [fingerprint]
    ├── Suzuki_001 → [fingerprint]
    ├── Suzuki_002 → [fingerprint]
    ├── Amide_formation_001 → [fingerprint]
    └── ...
```

**Implementation:**
```python
# In search.py, for cross-family search:
if family_txt is None and _drfp_storage_available:
    # Load unified index once
    if "__ALL_FAMILIES__" not in _DRFP_LOADER_CACHE:
        npz_path = "data/reaction_dataset/ALL_FAMILIES_drfp.npz"
        if os.path.exists(npz_path):
            _DRFP_LOADER_CACHE["__ALL_FAMILIES__"] = DRFPLoader(npz_path)
    
    # Lookup fingerprint (O(1) by reaction_id)
    loader = _DRFP_LOADER_CACHE.get("__ALL_FAMILIES__")
    if loader:
        r_fp = loader.get_fingerprint(reaction_id)  # ✅ Fast lookup!
```

**Performance:** ~2-5 seconds (same as feature-based, but with DRFP accuracy)

### Option 2: Load All Family NPZ Files

Load all NPZ files and track which family each reaction belongs to:

```python
# Load all family NPZ files
for family_name in ["C_N_Coupling", "Suzuki", "Amide_formation", ...]:
    npz_path = get_drfp_path_for_family(family_name)
    _DRFP_LOADER_CACHE[family_name] = DRFPLoader(npz_path)

# For each reaction, extract family from reaction_id
reaction_id = "C_N_Coupling_001"
family = reaction_id.split("_")[0] + "_" + reaction_id.split("_")[1]  # "C_N_Coupling"
loader = _DRFP_LOADER_CACHE.get(family)
r_fp = loader.get_fingerprint(reaction_id)
```

**Pros:**
- ✅ Uses existing NPZ files
- ✅ No need to rebuild unified index

**Cons:**
- ❌ More complex logic (parse reaction_id to extract family)
- ❌ Loads all NPZ files (~9.5 MB total)
- ❌ Slightly slower than single unified file

---

## Summary

### Why We Can't Use Precomputed DRFP for Cross-Family (Currently)

1. **Storage structure:** Fingerprints are stored per-family in separate NPZ files
2. **Lookup requirement:** Need to know which family/NPZ file contains each reaction
3. **Cross-family limitation:** When `family=None`, we can't determine NPZ file path
4. **Performance impact:** On-demand computation is too slow (30-120 seconds)

### Current Workaround

- **Disable DRFP by default** for cross-family search
- **Use feature-based similarity** (bin, LG, nuc_class)
- **Fast performance** (~2-5 seconds)
- **Lower accuracy** but still useful

### Future Enhancement

- **Create unified NPZ file** with all families
- **Enable DRFP for cross-family** with O(1) lookup
- **Best of both worlds:** Fast + accurate

### Current Behavior

```python
# Cross-family search (DRFP disabled by default)
result = chem.recommend.conditions(
    reaction="...",
    search_all_families=True,  # Uses feature-based similarity
    k=100
)
# Warning: "Cross-family search with DRFP disabled (no precomputed fingerprints available)"

# To enable DRFP (slow, but more accurate)
result = chem.recommend.conditions(
    reaction="...",
    search_all_families=True,
    relax={"use_drfp": True},  # Computes fingerprints on-the-fly (slow!)
    k=100
)
```

**Bottom line:** We **can't** use precomputed fingerprints because we'd need to know which NPZ file to check, and that requires knowing the family—which defeats the purpose of cross-family search. The solution is to create a unified index.
