# UnifiedRecommender: DRFP-Based Reaction Condition Recommendations

## Overview

The **UnifiedRecommender** provides intelligent reaction condition recommendations by searching a unified index of both **experimental protocols** from literature and **general reaction rules** using DRFP (Differentiable Reaction Fingerprints) similarity matching.

### Key Features

- **Unified Search**: Single interface for both protocols and rules
- **DRFP Similarity**: Reaction-level fingerprinting for accurate structural matching  
- **Flexible Filtering**: Filter by source type, similarity threshold, and result count
- **Rich Metadata**: Includes reaction family, tags, versions, and source attribution
- **Fast Performance**: Pre-computed fingerprints enable sub-second searches

### What's Included

**Protocols (18 sources)**:
- Full experimental procedures from Organic Syntheses and literature
- Complete reagent lists, conditions, workup steps, and notes
- Specific to individual reactions with detailed procedures

**Rules (9 sources)**:
- General reaction guidelines and best practices
- Covers major cross-coupling families (Suzuki, Buchwald-Hartwig, etc.)
- Multiple reference reactions with condition recommendations

---

## Installation & Setup

### Requirements

```bash
pip install drfp numpy
```

### Pre-built Index

The system uses a pre-built unified index located at:
```
build/unified_index_complete/
├── index.json          # Metadata for all sources
├── fingerprints.npz    # Pre-computed DRFP vectors
└── stats.json          # Build statistics
```

To rebuild or update the index:
```bash
python data-processor/Build_unified_drfp_index.py \
    --protocol-dir data/protocol_db_v2 \
    --rule-dir data/rule_db_v2 \
    --output build/unified_index_complete
```

---

## Usage

### Python API

#### Basic Usage

```python
from chemtools.recommend import UnifiedRecommender

# Initialize (uses default index)
recommender = UnifiedRecommender()

# Get recommendations
results = recommender.recommend(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    top_k=5,
    min_similarity=0.3
)

# Display results
for r in results:
    print(f"[{r.rank}] {r.name}")
    print(f"    Similarity: {r.similarity:.3f}")
    print(f"    Family: {r.family}")
    print(f"    Type: {r.source_type}")
    print()
```

Output:
```
[1] Buchwald–Hartwig C–N Coupling
    Similarity: 1.000
    Family: Buchwald–Hartwig_C–N
    Type: rule

[2] Ullmann–Goldberg (Cu) C–N Coupling
    Similarity: 0.375
    Family: Ullmann_C–N
    Type: rule
...
```

#### Filter by Source Type

```python
# Protocols only (full procedures)
protocols = recommender.recommend(
    reaction_smiles="C=CCN(C(=O)OC(C)(C)C)CC=C>>CC(C)(C)OC(=O)N1CCC=C1",
    top_k=3,
    source_types=["protocol"]
)

# Rules only (general guidelines)
rules = recommender.recommend(
    reaction_smiles="CC(=O)c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>CC(=O)c1ccc(-c2ccccc2)cc1",
    top_k=5,
    source_types=["rule"],
    min_similarity=0.5
)
```

#### Load Full Source Data

```python
# Get top recommendation
results = recommender.recommend(reaction_smiles, top_k=1)
best_match = results[0]

# Load complete protocol/rule data
full_data = recommender.get_source_details(best_match.id)

# Access detailed information
if full_data:
    print(json.dumps(full_data, indent=2))
```

#### Get Index Statistics

```python
stats = recommender.get_statistics()
print(f"Protocols: {stats['protocols']['count']}")
print(f"Rules: {stats['rules']['count']}")
print(f"Total Fingerprints: {stats['drfp']['computed']}")
```

### Command-Line Interface

#### Basic Query

```bash
python scripts/unified_recommend_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
```

#### Advanced Options

```bash
# Top 3 results with minimum similarity 0.5
python scripts/unified_recommend_cli.py "reaction.smi" -k 3 --min-sim 0.5

# Protocol-only search
python scripts/unified_recommend_cli.py "reaction.smi" --type protocol

# Show index statistics
python scripts/unified_recommend_cli.py "reaction.smi" --stats

# JSON output for scripting
python scripts/unified_recommend_cli.py "reaction.smi" --json

# Show detailed metadata
python scripts/unified_recommend_cli.py "reaction.smi" --details
```

### LangChain Agent Integration

The UnifiedRecommender is available as a LangChain tool in `chem_assistant`:

```python
from chem_assistant.chemtools_wrapper import unified_recommender_tool

# Use in LangChain agent
result = unified_recommender_tool.invoke({
    "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "top_k": 5,
    "min_similarity": 0.3,
    "source_type": "protocol"
})

print(result)
```

---

## API Reference

### `UnifiedRecommender`

Main class for DRFP-based recommendation.

#### `__init__(index_dir=None)`

Initialize recommender with unified index.

**Parameters:**
- `index_dir` (str | Path | None): Path to index directory. If None, uses default `build/unified_index_complete`.

**Raises:**
- `FileNotFoundError`: If index files missing
- `ImportError`: If DRFP library not installed

#### `recommend(reaction_smiles, top_k=5, min_similarity=0.0, source_types=None)`

Get ranked recommendations for a reaction.

**Parameters:**
- `reaction_smiles` (str): Reaction SMILES (reactants>>products)
- `top_k` (int): Number of results to return (default: 5)
- `min_similarity` (float): Minimum similarity 0.0-1.0 (default: 0.0)
- `source_types` (List[str] | None): Filter by ["protocol"], ["rule"], or None (both)

**Returns:**
- `List[RecommendationResult]`: Ranked results with similarity scores

**Example:**
```python
results = recommender.recommend(
    "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    top_k=3,
    min_similarity=0.5,
    source_types=["rule"]
)
```

#### `get_source_details(source_id)`

Load full source data from file.

**Parameters:**
- `source_id` (str): ID of protocol or rule

**Returns:**
- `Dict | None`: Complete source data or None if not found

**Example:**
```python
details = recommender.get_source_details("suzuki_v2")
print(details["metadata"]["name"])
print(details["default_rule"]["conditions"])
```

#### `get_statistics()`

Get index statistics.

**Returns:**
- `Dict`: Statistics including protocol/rule counts, fingerprints, families

**Example:**
```python
stats = recommender.get_statistics()
print(f"Version: {stats['build_info']['version']}")
print(f"Protocols: {stats['protocols']['count']}")
```

### `RecommendationResult`

Dataclass representing a single recommendation.

**Attributes:**
- `id` (str): Unique source identifier
- `name` (str): Human-readable name
- `source_type` (str): "protocol" or "rule"
- `family` (str): Reaction family
- `similarity` (float): DRFP similarity score (0.0-1.0)
- `rank` (int): Position in results (1-indexed)
- `tags` (List[str]): Descriptive tags
- `version` (str): Source version
- `source_file` (str): Path to source file
- `full_data` (Dict | None): Complete source data (loaded on demand)

**Methods:**
- `to_dict()`: Convert to dictionary representation

---

## Understanding DRFP Similarity

### What is DRFP?

DRFP (Differentiable Reaction Fingerprints) encodes entire chemical reactions as binary vectors by comparing molecular fingerprints of reactants and products. This captures structural changes rather than individual molecules.

### Similarity Interpretation

- **1.000**: Perfect match (identical reaction)
- **0.80-0.99**: Very similar (same reaction class, minor structural variations)
- **0.50-0.79**: Similar (related mechanism, different substituents)
- **0.30-0.49**: Somewhat similar (same general transformation type)
- **< 0.30**: Dissimilar (different reaction types)

### Recommended Thresholds

- **Exact matches**: `min_similarity=0.95`
- **Similar reactions**: `min_similarity=0.5`
- **Exploratory search**: `min_similarity=0.0` (default)

---

## Use Cases

### 1. Find Literature Precedents for New Reaction

```python
# Looking for experimental protocols for a Suzuki coupling
results = recommender.recommend(
    "CC(=O)c1ccc(I)cc1.c1ccc(B(O)O)cc1>>CC(=O)c1ccc(-c2ccccc2)cc1",
    top_k=5,
    source_types=["protocol"],
    min_similarity=0.7
)

# Get full procedure for top match
if results:
    details = recommender.get_source_details(results[0].id)
    print(details["original_procedure"])
```

### 2. Compare Protocol vs Rule Recommendations

```python
# Get both types
all_results = recommender.recommend(reaction, top_k=10)

# Separate by type
protocols = [r for r in all_results if r.source_type == "protocol"]
rules = [r for r in all_results if r.source_type == "rule"]

print(f"Protocols found: {len(protocols)}")
print(f"Rules found: {len(rules)}")
```

### 3. Explore Reaction Families

```python
# Find all sources for Buchwald-Hartwig reactions
results = recommender.recommend(
    "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    top_k=20,
    min_similarity=0.0
)

# Group by family
families = {}
for r in results:
    families.setdefault(r.family, []).append(r)

for family, matches in families.items():
    print(f"{family}: {len(matches)} matches")
```

### 4. Validate Planned Conditions

```python
# Check if planned conditions align with precedents
planned_reaction = "custom_reaction>>product"
results = recommender.recommend(planned_reaction, top_k=5)

if results and results[0].similarity > 0.8:
    print(f"High similarity to: {results[0].name}")
    print(f"Consider reviewing their conditions")
else:
    print("No close precedents found - proceed with caution")
```

---

## Troubleshooting

### No Results Returned

**Problem**: `recommender.recommend()` returns empty list

**Solutions**:
1. Lower `min_similarity` threshold (try 0.0)
2. Remove `source_types` filter (search both protocols and rules)
3. Verify reaction SMILES is valid (check for >> separator)
4. Increase `top_k` to see more distant matches

### Low Similarity Scores

**Problem**: All results have similarity < 0.3

**Solutions**:
1. Your reaction may be novel/underrepresented in database
2. Try searching for component reactions (e.g., just the C-C coupling step)
3. Use rule-based recommendations instead: `source_types=["rule"]`
4. Check if reaction family is represented: `recommender.get_statistics()`

### Import Error

**Problem**: `ImportError: drfp library is required`

**Solution**:
```bash
pip install drfp
```

### FileNotFoundError

**Problem**: `FileNotFoundError: Index file not found`

**Solutions**:
1. Ensure unified index exists at `build/unified_index_complete/`
2. Rebuild index: `python data-processor/Build_unified_drfp_index.py`
3. Specify custom index: `UnifiedRecommender(index_dir="path/to/index")`

---

## Technical Details

### Index Structure

```
build/unified_index_complete/
├── index.json
│   ├── version: "2.0"
│   ├── protocols: [{id, name, family, file_path, ...}, ...]
│   └── rules: [{id, name, family, file_path, ...}, ...]
├── fingerprints.npz
│   ├── protocol_fps: (18, 2048) uint8 array
│   ├── rule_fps: (50, 2048) uint8 array
│   ├── protocol_ids: array of protocol IDs
│   ├── rule_ids: array of rule IDs
│   └── rule_fp_indices: mapping of rule IDs to fingerprint ranges
└── stats.json
    └── build_info, protocols, rules, drfp statistics
```

### DRFP Computation

1. **Parse reaction**: Separate reactants and products from SMILES
2. **Generate molecular fingerprints**: Morgan/ECFP for each molecule
3. **Compute difference**: Subtract reactant FPs from product FPs
4. **Hash to fixed length**: 2048-bit binary vector

### Similarity Calculation

Tanimoto similarity for binary vectors:

```
similarity = (A ∩ B) / (A ∪ B)
           = count(A AND B) / (count(A) + count(B) - count(A AND B))
```

Where A and B are DRFP bit vectors.

### Performance

- **Index loading**: ~0.1s (one-time per session)
- **Single query**: ~0.01-0.05s
- **Memory footprint**: ~10MB (index + fingerprints)

---

## Future Enhancements

### Planned Features

1. **Dynamic indexing**: Add new protocols/rules without rebuild
2. **Multi-step reactions**: Handle reaction sequences
3. **Condition extraction**: Auto-extract conditions from protocols
4. **Yield prediction**: ML model for outcome prediction
5. **Interactive web UI**: Browser-based search interface

### Contributing

To add new protocols or rules:

1. Add JSON file to `data/protocol_db_v2/` or `data/rule_db_v2/`
2. Follow v2.0 schema (see `chemtools/schema/`)
3. Rebuild index: `python data-processor/Build_unified_drfp_index.py`
4. Validate: `pytest tests/test_unified_recommender.py`

---

## Related Tools

- **`protocol_recommendation_tool`**: Protocol-only search (legacy)
- **`rule_based_conditions_tool`**: Deterministic rule engine
- **`recommend_conditions_tool`**: ML-based recommendations
- **`search_precedents_tool`**: Precedent database search

The UnifiedRecommender provides a modern, similarity-based alternative that combines the best of protocols and rules in a single query.

---

## References

- **DRFP Paper**: Probst et al., "Reaction classification and yield prediction using the differential reaction fingerprint", Nat. Mach. Intell. (2022)
- **Unified Index Design**: Condition-agent project documentation
- **Protocol Database**: data/protocol_db_v2/ (Organic Syntheses and literature)
- **Rule Database**: data/rule_db_v2/ (HTE-validated guidelines)

---

## Support

For issues, questions, or contributions:
- File an issue on GitHub
- Check existing tests: `tests/test_unified_recommender.py`
- Review demo script: `demo_unified_recommender.py`
