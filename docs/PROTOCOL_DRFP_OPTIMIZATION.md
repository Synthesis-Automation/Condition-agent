# Protocol DRFP Storage Optimization

## Problem

The protocol indexer was storing DRFP fingerprints (4096 float values per protocol) directly in the JSON index file. This resulted in:

- **Large file sizes**: Each protocol required ~16-32 KB just for the DRFP fingerprint
- **Slow parsing**: Loading the entire index required parsing large JSON arrays
- **Memory inefficient**: All fingerprints loaded into memory even if not needed

## Solution

Adopted the same approach used by the reaction dataset system:

### 1. **Separate NPZ Storage**

DRFP fingerprints are now stored in a separate compressed NumPy NPZ file (`.protocol_drfp.npz`) instead of being embedded in the JSON index.

**Benefits:**
- ~90% reduction in index JSON size
- Faster JSON parsing (no large arrays)
- Fingerprints loaded on-demand only when needed

### 2. **Binary Format with uint8**

Fingerprints are stored as `uint8` arrays (0-255 values) instead of float lists.

**Benefits:**
- More compact than JSON float arrays
- Native NumPy format for fast operations
- Compressed NPZ format for additional space savings

### 3. **Lazy Loading via DRFPLoader**

The `DRFPLoader` utility (from `chemtools.util.drfp_storage`) provides efficient on-demand access:

```python
# Load index (lightweight - no DRFP data)
indexer = ProtocolIndexer.load()

# Get DRFP fingerprint only when needed
fp = indexer.get_drfp_fingerprint("protocol_001.json")
```

## File Structure

```
data/protocol_db/
  ├── .protocol_index.json        # Metadata only (~100 KB for 100 protocols)
  ├── .protocol_drfp.npz          # Binary DRFP fingerprints (~50 KB for 100 protocols)
  ├── protocol_001.json
  ├── protocol_002.json
  └── ...
```

## Size Comparison

For 100 protocols with DRFP fingerprints:

| Approach | Index Size | Total Size | Notes |
|----------|------------|------------|-------|
| **Old (embedded)** | ~2-3 MB | ~2-3 MB | All in JSON |
| **New (separate NPZ)** | ~100-200 KB | ~150-250 KB | JSON + NPZ |
| **Reduction** | **~90%** | **~90%** | 10x smaller |

## Usage

### Building Index with DRFP

```python
from chemtools.protocol.indexer import ProtocolIndexer

# Build index with DRFP fingerprints
indexer = ProtocolIndexer()
indexer.build_index(compute_drfp=True)
indexer.save()  # Saves both .protocol_index.json and .protocol_drfp.npz
```

### Loading and Using DRFP

```python
# Load index
indexer = ProtocolIndexer.load()

# Check if DRFP available
if indexer.metadata['has_drfp']:
    # Get fingerprint for similarity comparison
    fp = indexer.get_drfp_fingerprint("Suzuki_Pd_XPhos.json")
    
    if fp is not None:
        # Use for similarity calculations
        from chemtools.reaction_similarity import tanimoto
        similarity = tanimoto(query_fp, fp)
```

### CLI Usage

```bash
# Build index with DRFP (creates both JSON and NPZ files)
python -m chemtools.protocol.cli build --drfp

# Build without DRFP (smaller, faster)
python -m chemtools.protocol.cli build --no-drfp

# View stats
python -m chemtools.protocol.cli stats
```

## Implementation Details

### Data Model Changes

**ProtocolRecord (Old):**
```python
@dataclass
class ProtocolRecord:
    # ... metadata fields ...
    drfp_fingerprint: Optional[List[float]] = None  # 4096 floats in JSON
```

**ProtocolRecord (New):**
```python
@dataclass
class ProtocolRecord:
    # ... metadata fields ...
    has_drfp: bool = False  # Flag only, no data
```

### Storage Strategy

1. **During indexing**: DRFP fingerprints computed and collected in memory
2. **After processing**: All fingerprints saved to NPZ file in one batch operation
3. **JSON index**: Only stores metadata flags (`has_drfp`)
4. **On load**: NPZ file loaded lazily only when `get_drfp_fingerprint()` is called

### Backward Compatibility

The system gracefully handles:
- Missing DRFP storage utilities (falls back to no-DRFP mode)
- Missing NPZ files (returns None for fingerprint requests)
- Old indexes without DRFP support (loads normally, DRFP just unavailable)

## Performance Impact

### Index Build Time

- Slightly longer initial build (NPZ write overhead)
- But incremental updates remain fast (unchanged files skipped)

### Index Load Time

- **Much faster**: JSON parsing ~10x faster without large arrays
- DRFP data loaded only on first `get_drfp_fingerprint()` call

### Memory Usage

- **Much lower**: Index metadata fits in ~1-2 MB for 1000 protocols
- DRFP data loaded on-demand (~400 KB for 1000 protocols when accessed)

## Consistency with Dataset System

This approach mirrors the `chemtools.precedent` system:

| Feature | Precedent System | Protocol System |
|---------|------------------|-----------------|
| Index file | `.jsonl` with metadata | `.protocol_index.json` |
| DRFP storage | `{family}_drfp.npz` | `.protocol_drfp.npz` |
| Loader | `DRFPLoader` | Same `DRFPLoader` |
| Key field | `reaction_id` | `filename` |

## Migration Notes

If you have an existing protocol index with embedded DRFP data:

1. **Rebuild index**: Simply run `python -m chemtools.protocol.cli build --drfp --force`
2. **Old index**: The old `.protocol_index.json` will be replaced
3. **New files**: New `.protocol_drfp.npz` will be created

No data loss - all protocols re-processed from source JSON files.

## Future Enhancements

Possible improvements:
- [ ] Incremental DRFP updates (only recompute changed protocols)
- [ ] Multiple DRFP parameter sets (different n_bits/radius)
- [ ] Precomputed similarity matrices for small protocol collections
- [ ] Memory-mapped access for very large protocol databases
