# Protocol Database Module

## Overview

The protocol database module provides DRFP-based similarity search to find standard experimental procedures for user-supplied reactions. This system follows the same pattern as the ML-based condition recommendation but operates on curated protocol JSON files.

## Architecture

```
chemtools/protocol/
├── __init__.py           # Package exports
├── indexer.py           # Index builder with DRFP computation
├── recommend.py         # DRFP similarity-based recommendation
├── cli.py              # Command-line tools
└── matcher.py          # (Reserved for future extensions)
```

## Key Components

### 1. Protocol Indexer (`indexer.py`)

**Purpose:** Build and maintain a searchable index of protocol files

**Features:**
- Scans `data/protocol_db/*.json` files
- Extracts metadata (reaction, family, tags, source)
- Computes DRFP fingerprints for each protocol
- Builds lookup indexes (by family, tags)
- Supports incremental updates (only reindex changed files)
- Stores index as JSON for fast loading

**Data Structure:**
```python
ProtocolRecord:
    - filename: str
    - reaction_smiles: str
    - reaction_family: str
    - tags: List[str]
    - source_title, source_doi, etc.
    - drfp_fingerprint: List[float]  # Precomputed DRFP
```

### 2. Protocol Recommender (`recommend.py`)

**Purpose:** Find most similar protocols using DRFP similarity

**Features:**
- Loads precomputed DRFP fingerprints from index
- Computes DRFP for query reaction
- Finds top-k most similar protocols using cosine similarity
- Supports filtering by reaction family and tags
- Extracts detailed experimental conditions

**Similarity Metric:** Cosine similarity between DRFP vectors (same as ML recommendation)

### 3. CLI Tools (`cli.py`)

**Purpose:** Command-line interface for index management

**Commands:**
```bash
# Build index
python -m chemtools.protocol.cli build

# Show statistics
python -m chemtools.protocol.cli stats

# List families
python -m chemtools.protocol.cli list-families

# Show protocols for a family
python -m chemtools.protocol.cli show-family Suzuki
```

## Usage

### Building the Index

**First time setup:**
```bash
python -m chemtools.protocol.cli build
```

This will:
1. Scan `data/protocol_db/*.json`
2. Extract metadata from each file
3. Compute DRFP fingerprints
4. Save index to `data/protocol_db/.protocol_index.json`

**Incremental updates:**
```bash
python -m chemtools.protocol.cli build
```
- Only reindexes changed files (uses MD5 hash)
- Much faster for adding new protocols

**Force rebuild:**
```bash
python -m chemtools.protocol.cli build --force
```

### Getting Recommendations

**Python API:**
```python
from chemtools.protocol import ProtocolRecommender

# Initialize (loads index)
recommender = ProtocolRecommender()

# Get top-5 similar protocols
results = recommender.recommend(
    reaction_smiles='BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1',
    k=5
)

# Print matches
for match in results['matches']:
    print(f"{match['similarity']:.3f}: {match['source_title']}")
    print(f"  DOI: {match['source_doi']}")
```

**With filtering:**
```python
# Filter by tags
results = recommender.recommend(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=3,
    tags=['suzuki', 'palladium']
)

# Filter by reaction family
results = recommender.recommend(
    reaction_smiles='CCCI.c1ccccc1B(O)O>>CCCc1ccccc1',
    k=3,
    reaction_family='Suzuki_Cu_alkyl_halide+aryl_boron'
)
```

**With condition extraction:**
```python
# Get recommendations with extracted conditions
results = recommender.recommend_with_details(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=3
)

for match in results['matches']:
    cond = match['conditions']
    print(f"Catalyst: {cond['catalyst']}")
    print(f"Solvent: {cond['solvent']}")
    print(f"Temperature: {cond['temperature_C']} °C")
```

### Querying the Index

**Show statistics:**
```bash
python -m chemtools.protocol.cli stats
```

Output:
```
Total protocols: 16
Families: 16
Tags: 71
DRFP fingerprints: Yes

Protocols by family:
  Alkyl_Iodide_Borylation: 1
  Suzuki_Cu_alkyl_halide+aryl_boron: 1
  ...
```

**List all families:**
```bash
python -m chemtools.protocol.cli list-families
```

**Show protocols for a family:**
```bash
python -m chemtools.protocol.cli show-family "Suzuki_Cu_alkyl_halide+aryl_boron"
```

**Show protocols with a tag:**
```bash
python -m chemtools.protocol.cli show-tag palladium
```

## Protocol JSON Schema

Each protocol file follows the `Structured_Output_schema.json`:

```json
{
  "source": {
    "title": "...",
    "journal": "...",
    "year": 2025,
    "doi": "...",
    "url": "..."
  },
  "reaction": {
    "reaction_smiles": "A.B>>C",
    "family": "Suzuki_Cu_alkyl_halide+aryl_boron",
    "notes": "...",
    "tags": "Suzuki; Cu; alkyl_halide; Coupling"
  },
  "reaction_setup": [...],
  "workup_and_purification": [...],
  "original_procedure": "..."
}
```

## Index File Structure

The index file (`.protocol_index.json`) contains:

```json
{
  "metadata": {
    "version": "1.0",
    "created_at": "2025-10-10T...",
    "updated_at": "2025-10-10T...",
    "num_protocols": 16,
    "has_drfp": true
  },
  "records": {
    "Suzuki_Cu_C(sp3)-C(sp2).json": {
      "filename": "...",
      "reaction_smiles": "...",
      "reaction_family": "...",
      "tags": [...],
      "drfp_fingerprint": [0.1, 0.3, ...],  // 2048 floats
      ...
    }
  },
  "family_index": {
    "Suzuki_Cu_alkyl_halide+aryl_boron": ["Suzuki_Cu_C(sp3)-C(sp2).json"]
  },
  "tag_index": {
    "suzuki": ["Suzuki_Cu_C(sp3)-C(sp2).json", ...],
    "palladium": [...]
  }
}
```

## Performance Considerations

1. **Index Size:** ~1-2 KB per protocol (including DRFP)
   - 100 protocols: ~100-200 KB
   - 1000 protocols: ~1-2 MB

2. **Index Build Time:**
   - ~50-100 ms per protocol (DRFP computation)
   - 16 protocols: ~1 second
   - 1000 protocols: ~1 minute

3. **Recommendation Speed:**
   - Index loading: ~10-50 ms
   - DRFP computation (query): ~50 ms
   - Similarity search: ~1 ms per protocol
   - **Total for 1000 protocols: ~100-150 ms**

4. **Incremental Updates:**
   - Only recomputes DRFP for changed files
   - Uses MD5 hash for change detection
   - Typical update: <1 second

## Integration with Existing Systems

### Similar to ML Recommendation

The protocol recommendation follows the same pattern:

| Feature | ML Recommendation | Protocol Recommendation |
|---------|------------------|------------------------|
| Data source | `reaction_dataset/reactions_sample.jsonl` | `protocol_db/*.json` |
| Indexing | DRFP precomputation | DRFP precomputation |
| Similarity | Cosine similarity | Cosine similarity |
| Filtering | By family, constraints | By family, tags |
| Output | Condition recommendations | Protocol matches |

### Code Reuse

- Uses same `reaction_similarity.py` for DRFP (if available)
- Similar API design to `recommend/core.py`
- Follows ChemTools module patterns

## Adding New Protocols

**Step 1:** Create protocol JSON file

```bash
# Add new protocol to data/protocol_db/
cp my_protocol.json data/protocol_db/
```

**Step 2:** Rebuild index (incremental)

```bash
python -m chemtools.protocol.cli build
```

This will automatically detect and index the new file.

**Step 3:** Verify

```bash
python -m chemtools.protocol.cli stats
```

## Testing

Run the test suite:

```bash
python test_protocol_recommendation.py
```

Tests cover:
1. ✅ Suzuki coupling recommendation
2. ✅ Borylation recommendation
3. ✅ Filtered recommendation (by tags)
4. ✅ Recommendation with condition extraction

## Future Enhancements

### Planned Features

1. **API Integration**
   - Add `/api/v1/protocol/recommend` endpoint
   - Match the pattern of `/api/v1/recommend`

2. **Enhanced Filtering**
   - Filter by reagent availability
   - Filter by reaction conditions (temp range, time)
   - Filter by yield threshold

3. **Hybrid Recommendation**
   - Combine protocol matching with condition recommendation
   - Use protocol as starting point, then optimize with analytics

4. **Full-Text Search**
   - Index procedure text for keyword search
   - Combine with DRFP similarity

5. **Protocol Validation**
   - Check schema compliance
   - Validate SMILES and CAS numbers
   - Flag missing required fields

## Troubleshooting

### Index not found

```
❌ Index not found: .protocol_index.json
```

**Solution:** Build the index first
```bash
python -m chemtools.protocol.cli build
```

### DRFP not available

```
⚠️  DRFP package not available
```

**Solution:** Install DRFP
```bash
pip install drfp
```

### No matches found

**Possible causes:**
1. Query reaction too different from protocols
2. Filters too restrictive (family/tags)
3. Index out of date

**Solution:**
- Try without filters
- Broaden tags
- Rebuild index

## Summary

The protocol database module provides:

✅ **Fast indexing** with incremental updates  
✅ **DRFP-based similarity** (same as ML recommendation)  
✅ **Flexible filtering** by family and tags  
✅ **CLI tools** for index management  
✅ **Clean API** matching existing ChemTools patterns  
✅ **Scalable** to thousands of protocols  

Ready for integration with the web API and recommendation workflow!
