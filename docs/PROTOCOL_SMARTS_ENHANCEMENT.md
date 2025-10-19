# Protocol SMARTS-Based Recommendation Enhancement

## Summary

Successfully updated the protocol-based recommendation system to use `reaction_SMARTS` patterns from protocol JSON files for intelligent structural matching of input reactions.

## Changes Made

### 1. **Indexer Module** (`chemtools/protocol/indexer.py`)
- Added `reaction_smarts` field to `ProtocolRecord` dataclass
- Updated `_process_protocol_file()` to extract `reaction_SMARTS` from protocol JSON files
- Maintains backward compatibility with protocols that don't have SMARTS patterns

### 2. **Recommender Module** (`chemtools/protocol/recommend.py`)
- Added `match_reaction_smarts()` function for structural matching:
  - Uses RDKit's reaction SMARTS matching when available
  - Checks if input reaction matches any of the protocol's SMARTS patterns
  - Graceful fallback when RDKit is unavailable (permissive mode)
  
- Added `_filter_by_smarts()` method to filter candidates by structural match

- Updated `recommend()` method:
  - New parameter: `use_smarts_filter` (default: `True`)
  - Applies SMARTS-based pre-filtering before DRFP similarity scoring
  - Tracks filtering statistics in metadata
  
- Updated all related methods to support the new parameter:
  - `recommend_with_details()`
  - `recommend_protocol()` standalone function

- Enhanced output to include `reaction_smarts` in protocol information

### 3. **Matcher Module** (`chemtools/protocol/matcher.py`)
- Updated for consistency (marked as legacy)
- Added `reaction_smarts` field to `ProtocolMetadata` and `ProtocolMatch`
- Updated `_extract_metadata()` to extract SMARTS patterns

### 4. **Testing**
- Created `test_protocol_smarts.py` test script
- Verified SMARTS matching works correctly:
  - ✅ Matches structurally similar reactions
  - ✅ Rejects non-matching reactions
  - ✅ Integrates with full recommendation pipeline

## How It Works

### Workflow:
1. **Index Building**: Extract `reaction_SMARTS` patterns from protocol JSON files
2. **Query Processing**: User provides reaction SMILES
3. **SMARTS Filtering** (optional):
   - For each protocol, check if query matches any of its SMARTS patterns
   - Uses RDKit's substructure matching on both reactants and products
   - Filters out protocols that don't structurally match
4. **DRFP Similarity**: Compute similarity scores on remaining candidates
5. **Ranking**: Return top-k most similar protocols

### Benefits:
- **More Precise**: Filters out structurally incompatible protocols before similarity scoring
- **Better Recommendations**: Focuses on protocols that match the reaction type
- **Configurable**: Can enable/disable SMARTS filtering as needed
- **Backward Compatible**: Protocols without SMARTS patterns are still included (permissive)

## Example Usage

```python
from chemtools.protocol.recommend import ProtocolRecommender

recommender = ProtocolRecommender()

# Enable SMARTS filtering (default)
results = recommender.recommend(
    reaction_smiles='O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>O=C(C)c1ccc(C2C(CCC2)=O)cc1',
    k=5,
    use_smarts_filter=True  # Enable structural matching
)

# Disable SMARTS filtering (use only DRFP similarity)
results = recommender.recommend(
    reaction_smiles='O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>O=C(C)c1ccc(C2C(CCC2)=O)cc1',
    k=5,
    use_smarts_filter=False  # Disable, use only fingerprint similarity
)
```

## Test Results

### Example: Alpha-Arylation Query
**Input**: `O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>O=C(C)c1ccc(C2C(CCC2)=O)cc1`

**With SMARTS filtering enabled:**
- 1 protocol matched (perfectly matches the structural pattern)
- Confidence: 1.0000
- Protocol: α-Arylation of Cyclopentanones by Palladium/Enamine Cooperative Catalysis

**Without SMARTS filtering:**
- 3 protocols returned (based on DRFP similarity only)
- Includes structurally unrelated protocols with low similarity scores

## Protocol JSON Format

Protocols should include `reaction_SMARTS` field with structural patterns:

```json
{
  "reaction": {
    "reaction_smiles": "O=C1CCCC1.Brc2ccc(C(C)=O)cc2>>O=C(C)c3ccc(C4C(CCC4)=O)cc3",
    "reaction_SMARTS": ["O=C1CCCC1.Br[a]>>[a]C2C(CCC2)=O"],
    "family": "Pd/Enamine_α-Arylation_C(sp3)-C(sp2)",
    ...
  }
}
```

## Rebuilding the Index

After adding or updating `reaction_SMARTS` in protocol JSON files, rebuild the index:

```bash
python -m chemtools.protocol.cli build
```

This will:
- Scan all protocol JSON files
- Extract reaction_SMARTS patterns
- Compute DRFP fingerprints
- Build searchable index with SMARTS patterns included

## Technical Details

### SMARTS Matching Algorithm
1. Parse input reaction SMILES using RDKit
2. Parse each SMARTS pattern as a reaction template
3. Check if all pattern reactants are substructures of input reactants
4. Check if all pattern products are substructures of input products
5. Return true if both reactants and products match

### Fallback Behavior
- If RDKit is unavailable: Returns `True` (permissive, allows all protocols)
- If protocol has no SMARTS patterns: Included in results (permissive)
- If SMARTS parsing fails: Logs warning and continues (robust)

## Future Enhancements

Possible improvements:
1. Support multiple SMARTS patterns per protocol (already implemented!)
2. Add SMARTS pattern confidence/weight
3. Generate SMARTS patterns automatically from reaction SMILES
4. Support fuzzy SMARTS matching with tolerance
5. Combine SMARTS match score with DRFP similarity score

## Files Modified

1. `chemtools/protocol/indexer.py` - Index building with SMARTS extraction
2. `chemtools/protocol/recommend.py` - Recommendation with SMARTS filtering
3. `chemtools/protocol/matcher.py` - Legacy matcher updated for consistency
4. `test_protocol_smarts.py` - Test script for verification

## Compatibility

- ✅ Backward compatible with existing code
- ✅ Works with or without RDKit
- ✅ Protocols without SMARTS patterns still work
- ✅ Can enable/disable SMARTS filtering as needed
- ✅ No breaking changes to API
