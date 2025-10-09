# ML-Based Precedent Filtering by Reagent Database

## Overview

The ML-based recommendation system now includes **reagent database filtering** as a default feature. This ensures that recommended precedents only include reagents that are available in the curated reagent database, making recommendations more practical and actionable.

## Feature Description

### What It Does

When enabled (default), the precedent search filters out reactions where any reagent (catalyst, ligand, base, solvent, additive, etc.) is **not found** in the reagent database. This ensures:

1. **Database Completeness**: All recommended chemicals have full metadata (CAS numbers, SMILES, aliases)
2. **Practical Recommendations**: Only suggest reagents that are well-documented and available
3. **Quality Control**: Filter out obscure or proprietary reagents without database entries

### Default Behavior

✅ **Filtering is ENABLED by default**

This is controlled by the `filter_by_reagent_database` parameter in the `relax` dict passed to `precedent.knn()`.

## Usage

### Default (Filtering Enabled)

```python
from chemtools import precedent

# Get precedents - filtering enabled by default
pack = precedent.knn(
    family="C_N_Coupling_Pd",
    features={},
    k=25,
    relax={
        "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    }
)

# Only precedents with all reagents in database are returned
print(f"Support: {pack['support']} precedents")  # e.g., 14 precedents
```

### Disable Filtering (Get All Precedents)

```python
# Disable filtering to get all precedents regardless of database availability
pack = precedent.knn(
    family="C_N_Coupling_Pd",
    features={},
    k=25,
    relax={
        "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        "filter_by_reagent_database": False,  # Disable filtering
    }
)

print(f"Support: {pack['support']} precedents")  # e.g., 1307 precedents
```

### Explicit Enable (Same as Default)

```python
pack = precedent.knn(
    family="C_N_Coupling_Pd",
    features={},
    k=25,
    relax={
        "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        "filter_by_reagent_database": True,  # Explicitly enable (same as default)
    }
)
```

## How It Works

### 1. Reagent Extraction

For each precedent, the system extracts:

- **Catalytic System**: Metal precursors and ligands from `catalytic_system` array
- **Reagents**: Bases, additives, acids, oxidants, reductants from `reagents` array  
- **Solvents**: Solvents from `solvents` array

### 2. Database Lookup

Each reagent name is checked against the appropriate database:

| Role | Database File |
|------|---------------|
| Metal precursor / Catalyst | `metal_precursor.json` |
| Ligand | `ligand.json` |
| Base | `base.json` |
| Solvent | `solvent.json` |
| Additive | `additive.json` |
| Acid | `acid.json` |
| Oxidant | `oxidant.json` |
| Reductant | `reductant.json` |

### 3. Filtering Decision

- **Keep**: If ALL reagents are found in their respective databases
- **Discard**: If ANY reagent is missing from the database

### 4. Re-ranking

After filtering, precedents are re-ranked by similarity, and the top results are returned.

## Example Results

### Test Case: Buchwald-Hartwig Amination

**Reaction**: `Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1`

**Without Filtering**:
```
Support: 1307 precedents
Top precedent: Pd (yield 99%)
  - Missing reagent: '2599846-83-8' (not in database)
```

**With Filtering (Default)**:
```
Support: 14 precedents (filtered out 1293 = 98.9%)
Top precedent: Pd/Tri-tert-butylphosphinetetrafluoroborate (yield 100%)
  - All reagents found:
    ✅ Tris(dibenzylideneacetone)dipalladium(0)
    ✅ Tri-tert-butylphosphine tetrafluoroborate
    ✅ Sodium tert-butoxide
    ✅ Toluene
```

## Checking Reagent Availability

You can manually check if precedent reagents are in the database:

```python
from chemtools.reagent_lookup import check_precedent_reagents_in_database

# Check a precedent
precedent = pack["precedents"][0]
result = check_precedent_reagents_in_database(precedent)

print(f"All found: {result['all_found']}")
print(f"Found: {result['found_count']}/{result['total_count']}")
print(f"Missing: {result['missing']}")

# Output:
# All found: True
# Found: 4/4
# Missing: []
```

### Result Structure

```python
{
    "all_found": bool,          # True if all reagents found
    "missing": [                # List of missing reagents
        {"name": str, "type": str},
        ...
    ],
    "found_count": int,         # Number found
    "total_count": int,         # Total number checked
}
```

## Integration with Recommendation Pipeline

The filtering is automatically applied in:

1. **`chemtools.precedent.knn()`**: Core k-NN search function
2. **`chemtools.recommend.recommend_from_reaction()`**: Main recommendation API
3. **API endpoints**: `/api/recommend` and `/api/ml-recommend`

### API Usage

The FastAPI endpoints inherit the default filtering behavior:

```bash
# POST /api/ml-recommend
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "k": 25
}

# Returns precedents with all reagents in database (default)
```

To disable filtering via API:

```bash
# POST /api/ml-recommend
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "k": 25,
  "relax": {
    "filter_by_reagent_database": false
  }
}
```

## Performance Impact

- **Minimal overhead**: Database lookups are cached via `@lru_cache`
- **Typical filtering time**: <10ms for 1000+ precedents
- **First call**: Slightly slower due to loading reagent databases (~50-100ms)
- **Subsequent calls**: Very fast due to caching

## Configuration Options

All configuration is done via the `relax` parameter:

```python
relax = {
    # Enable/disable filtering (default: True)
    "filter_by_reagent_database": True,
    
    # Other precedent search options
    "use_drfp": True,
    "selective_loading": True,
    "debug_timing": False,
}
```

## Benefits

### 1. **Higher Quality Recommendations**
- Only suggest well-documented reagents
- Avoid obscure chemicals without metadata

### 2. **Better User Experience**
- All recommended reagents have CAS numbers
- Full chemical names and abbreviations available
- SMILES strings for structure visualization

### 3. **Practical Guidance**
- Recommendations are actionable
- Reagents are likely commercially available
- Reduced risk of suggesting proprietary compounds

### 4. **Database Consistency**
- Ensures reagent enrichment works properly
- No null values for important metadata
- Consistent data quality across recommendations

## Testing

Run the test script to verify filtering:

```bash
python test_ml_filtering.py
```

Expected output:
```
Without filtering: 1307 precedents
With filtering:    14 precedents
Filtered out:      1293 precedents (98.9%)

✅ Filtering is working!
```

## Migration Notes

### For Existing Code

**No changes required** - filtering is enabled by default. Your existing code will automatically benefit from this feature.

### To Preserve Old Behavior

If you need the old behavior (all precedents regardless of database availability):

```python
# Add this to your relax dict
relax["filter_by_reagent_database"] = False
```

## Future Enhancements

Potential improvements:

- **Partial filtering**: Keep precedents with ≥80% reagents in database
- **Substitution suggestions**: Suggest similar reagents from database
- **Coverage metrics**: Report database coverage per reaction family
- **Custom reagent lists**: Filter by user-provided reagent inventory

## Related Documentation

- `docs/EQUIVALENTS_CALCULATION.md` - Equivalents calculation from rule files
- `docs/RULE_FORMAT_UPDATE.md` - Standardized rule format
- `chemtools/reagent_lookup.py` - Reagent database utilities
- `chemtools/precedent/search.py` - Precedent search implementation
