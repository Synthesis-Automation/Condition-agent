# Cross-Family ML Recommendation Search

## Overview

The ML-based recommendation system now supports **cross-family search**, which allows you to search across ALL reaction datasets without filtering by reaction type. This is useful when:

- The reaction type is uncertain or ambiguous
- You want to explore conditions from multiple reaction families
- You're working with novel transformations that don't fit standard categories
- You want to discover similar reactions regardless of classification

## How It Works

### Standard Search (Family-Specific)
```python
from chemtools import chem

# Searches only within detected/specified family (e.g., "C_N_Coupling")
result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50,
    search_all_families=False  # Default
)
```

- ✅ Detects reaction family (auto or user-specified)
- ✅ Searches only within that family's dataset
- ✅ Uses family-specific reranking rules
- ✅ Faster (smaller dataset)
- ✅ Higher precision for well-defined reactions

### Cross-Family Search (All Datasets)
```python
from chemtools import chem

# Searches across ALL reaction families
result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=100,  # Use more precedents for broader search
    search_all_families=True  # NEW PARAMETER
)
```

- ✅ Skips reaction family detection
- ✅ Searches ALL datasets (Suzuki, C-N Coupling, Amide, etc.)
- ✅ Uses DRFP similarity only (no family-specific rules)
- ✅ Broader coverage (finds diverse precedents)
- ✅ Better for uncertain or novel reactions

## Key Differences

| Feature | Standard Search | Cross-Family Search |
|---------|----------------|---------------------|
| **Family Detection** | Yes (auto or user) | Skipped |
| **Dataset Scope** | Single family | All families |
| **Reranking** | Rule or analytics | DRFP similarity only |
| **Performance** | Faster | Slower (more data) |
| **Use Case** | Known reaction type | Unknown/uncertain type |
| **Precision** | High (focused) | Lower (broader) |
| **Discovery** | Limited to family | Cross-family insights |

## API Reference

### ChemTools API

```python
from chemtools import chem

result = chem.recommend.conditions(
    reaction: str,                    # Reaction SMILES
    reaction_type: Optional[str] = None,  # Ignored if search_all_families=True
    k: int = 50,                      # Number of precedents
    limit: int = 5,                   # Max recommendations
    search_all_families: bool = False,  # NEW: Enable cross-family search
    rerank_strategy: str = 'rule',    # 'rule', 'analytics', or 'none'
    filter_unknown_reagents: bool = False,
)
```

### CLI Interface

```bash
# Standard search (auto-detect family)
python app/local_recommendation_cli.py --rxn "Brc1ccccc1.c1ccccc1>>c1ccccc1Nc1ccccc1" --k 50

# Cross-family search (search all datasets)
python app/local_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"  --k 100 --search-all-families
```

### FastAPI Endpoint

```python
import requests

# Standard search
response = requests.post("http://localhost:8000/api/v1/recommend/conditions", json={
    "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "k": 50,
    "search_all_families": False
})

# Cross-family search
response = requests.post("http://localhost:8000/api/v1/recommend/conditions", json={
    "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "k": 100,
    "search_all_families": True
})
```

## Output Format

The output format is identical for both search modes, with additional metadata:

```json
{
  "meta": {
    "status": "success",
    "model": "drfp_similarity",
    "result_count": 3
  },
  "detection": {
    "detected_reaction_type": "C_N_Coupling",  // or "All" for cross-family
    "source": "auto",  // or "cross_family_search"
    "search_mode": "family_specific"  // or "cross_family"
  },
  "recommendations": [
    {
      "rank": 1,
      "confidence": 0.87,
      "chemicals": [...],
      "conditions": {...}
    }
  ],
  "precedents_used": {
    "total_count": 50,
    "top_precedents": [...]
  }
}
```

## Performance Considerations

### Data Loading
- **Standard**: Loads only requested family dataset (~1-5 MB)
- **Cross-family**: Loads ALL family datasets (~10-20 MB)
- **Recommendation**: Use `k=100-200` for cross-family (vs `k=25-50` for standard)

### DRFP Similarity

**IMPORTANT: Cross-family search currently has DRFP disabled by default due to performance.**

- **Standard**: Can use binary NPZ files (fastest, if available) + DRFP enabled by default
- **Cross-family**: DRFP disabled by default (no precomputed fingerprints across families)
  - Uses feature-based similarity only (bin, LG, nuc_class)
  - DRFP would require on-demand computation for all reactions (very slow)
  - Can be enabled with `relax={"use_drfp": True}` but will be significantly slower
- **Future**: Precomputed cross-family DRFP index for faster similarity search

**Performance comparison:**
- Standard (DRFP enabled): ~1-2 seconds
- Cross-family (DRFP disabled): ~2-5 seconds
- Cross-family (DRFP enabled): ~30-120 seconds (not recommended without precomputed FPs)

### Caching
- Results are cached by family + features + relax parameters
- Cross-family search uses special cache key `__ALL__`
- Cache hit rate may be lower for cross-family search

## Best Practices

### When to Use Cross-Family Search

✅ **Use cross-family search when:**
- Reaction type is uncertain or ambiguous
- Working with novel transformations
- Exploring diverse reaction conditions
- Need broad coverage over precision
- Have sufficient computational resources

❌ **Don't use cross-family search when:**
- Reaction type is well-known and clear
- Need fast response times (API/UI)
- Working with standard transformations
- Precision is more important than coverage

### Recommended Parameters

```python
# Cross-family search: DRFP disabled by default for performance
result = chem.recommend.conditions(
    reaction="...",
    k=100,  # More precedents for better coverage
    limit=5,  # Standard limit
    search_all_families=True,
    rerank_strategy='none',  # Feature-based similarity only (DRFP disabled)
    filter_unknown_reagents=True,  # Optional: filter unknown reagents
)

# To enable DRFP (slower, but more accurate for dissimilar structures)
result = chem.recommend.conditions(
    reaction="...",
    k=100,
    search_all_families=True,
    relax={"use_drfp": True},  # Enable DRFP (slow without precomputed FPs)
)
```

## Implementation Details

### Modified Modules

1. **`chemtools/precedent/search.py`**
   - Updated `knn()` to accept `family=None` or `family="All"`
   - Modified `_candidate_pool()` to skip family filtering when `family_txt=None`
   - Added special cache key `__ALL__` for cross-family searches

2. **`chemtools/recommend/modules/recommender.py`**
   - Added `search_all_families` parameter to `recommend_from_reaction()`
   - Skips family detection when `search_all_families=True`
   - Disables rule-based reranking for cross-family search

3. **`chemtools/recommend/modules/structured.py`**
   - Updated `recommend_conditions_structured()` to pass `search_all_families`

4. **`chemtools/context.py`**
   - Updated `RecommendNamespace.conditions()` to support new parameter

5. **`app/local_recommendation_cli.py`**
   - Added `--search-all-families` CLI flag

### Technical Notes

- Cross-family search automatically disables selective loading
- Rule-based and analytics reranking are skipped (family-specific)
- Only DRFP similarity is used for ranking
- Binary NPZ fingerprint files are not used (family-specific optimization)
- Detection metadata includes `search_mode: "cross_family"` marker

## Examples

### Example 1: Standard vs Cross-Family

```python
from chemtools import chem

reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

# Standard: Auto-detect as C-N coupling, search only C-N dataset
standard = chem.recommend.conditions(
    reaction=reaction,
    k=50,
    search_all_families=False
)
print(f"Family: {standard['detection']['detected_reaction_type']}")
# Output: Family: C_N_Coupling

# Cross-family: Search all datasets
cross = chem.recommend.conditions(
    reaction=reaction,
    k=100,
    search_all_families=True
)
print(f"Family: {cross['detection']['detected_reaction_type']}")
# Output: Family: All
```

### Example 2: Uncertain Reaction Type

```python
from chemtools import chem

# Novel/ambiguous transformation
reaction = "BrC1=CC=CC=C1.N2C=NC=C2>>C1=CC=CC=C1N2C=NC=C2"

# Use cross-family search for exploration
result = chem.recommend.conditions(
    reaction=reaction,
    k=150,  # Cast a wide net
    limit=10,  # Get more options
    search_all_families=True
)

# Examine which families contributed to recommendations
precedents = result.get('precedents_used', {}).get('top_precedents', [])
families = {p.get('rxn_type', 'Unknown') for p in precedents}
print(f"Found precedents from families: {families}")
```

### Example 3: CLI Usage

```bash
# Test cross-family search with example reaction
python app/local_recommendation_cli.py \
  --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" \
  --k 100 \
  --limit 5 \
  --search-all-families \
  --strategy ml \
  --save-dir results/cross_family_test

# Compare with standard search
python app/local_recommendation_cli.py \
  --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" \
  --k 50 \
  --limit 5 \
  --family C_N_Coupling \
  --strategy ml \
  --save-dir results/standard_test
```

## Testing

Run the test script to verify the implementation:

```bash
python test_cross_family_search.py
```

This will compare standard vs cross-family search results and show:
- Detected families
- Number of recommendations
- Search modes
- Top recommendations from each method

## Future Enhancements

Potential improvements for cross-family search:

1. **Precomputed cross-family DRFP index** - Store all fingerprints in a single NPZ file
2. **Hybrid reranking** - Use reaction-type-agnostic rules for cross-family search
3. **Family clustering** - Group similar families for faster search
4. **Adaptive k-value** - Automatically adjust k based on dataset size
5. **Cross-family analytics** - Track which families contribute most to recommendations

## Summary

The cross-family search feature provides a powerful way to explore reaction conditions without being constrained by reaction family classifications. It's particularly valuable for:

- ✅ Uncertain or novel reactions
- ✅ Broad exploration and discovery
- ✅ Finding diverse precedents
- ✅ Cross-family insights

Use `search_all_families=True` to enable this feature in any ML recommendation API call.
