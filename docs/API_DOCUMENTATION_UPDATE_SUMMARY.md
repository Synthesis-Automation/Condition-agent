# API Documentation Update Summary

## Overview
Updated `API_DOCUMENTATION.md` to reflect new recommendation endpoint features and improvements implemented in 2024.

---

## What Was Updated

### 1. ✨ New "What's New (2024)" Section

Added prominent section at top of document highlighting:

**Specific Catalyst Preservation:**
- Shows exact catalyst complexes (e.g., "Pd(dppf)Cl2·DCM adduct" with CAS number)
- No more generic "Palladium" or "Copper" in recommendations
- Preserves specific catalyst names from top precedents

**New Reranking Strategies:**
- `rerank_strategy='none'` - Pure similarity (fastest)
- `rerank_strategy='rule'` - Reaction-specific rules (recommended)
- `rerank_strategy='analytics'` - Dataset analytics (novel reactions)

**Unknown Reagent Filtering:**
- `filter_unknown_reagents=True` to exclude recommendations with "[Unknown ...]" components
- Ensures all returned reagents are identified in the database

---

### 2. 📝 Enhanced Use Case 1 (Primary Recommendation Endpoint)

**What Changed:**
- Added comprehensive examples showing new parameters
- Included before/after catalyst specificity examples
- Added reranking strategy comparison table
- Showed filter_unknown_reagents usage

**Example Code Added:**

```python
# Advanced example with new parameters
response = requests.post(
    'http://localhost:8000/api/v1/recommend',
    json={
        'reaction': 'c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1c1ccccc1',  # Suzuki
        'reaction_family': 'suzuki',
        'k': 100,
        'rerank_strategy': 'rule',  # NEW
        'filter_unknown_reagents': True  # NEW
    }
)

# Now shows specific catalysts:
# "Dichloro(1,1'-bis(diphenylphosphino)ferrocene)palladium(II) dichloromethane adduct"
# CAS: 95464-05-4
```

**Reranking Strategy Guide Table:**

| Strategy | When to Use | Effect |
|----------|-------------|--------|
| `'none'` | Fast lookups, well-known reactions | Pure DRFP similarity ranking |
| `'rule'` | Standard reactions (Suzuki, Buchwald, etc.) | Rerank using reaction-specific rules ⭐ **Recommended** |
| `'analytics'` | Novel/complex reactions | Rerank using dataset analytics |

---

### 3. ⚠️ Deprecated Fusion Endpoint

**What Changed:**
- Marked `/api/v1/recommend/fusion` as **DEPRECATED**
- Added clear warning banner
- Provided migration guide to new unified endpoint

**Migration Example:**

```python
# OLD (deprecated)
response = requests.post(
    'http://localhost:8000/api/v1/recommend/fusion',
    json={'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1', 'k': 50}
)

# NEW (recommended)
response = requests.post(
    'http://localhost:8000/api/v1/recommend',
    json={
        'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
        'k': 50,
        'rerank_strategy': 'analytics'  # Replaces fusion endpoint
    }
)
```

---

### 4. 📋 Updated Endpoint List & Parameter Table

**Endpoint List Changes:**
- Marked `/api/v1/recommend` as ⭐ **Primary** endpoint
- Updated description to mention "supports reranking, filtering"
- Marked `/api/v1/recommend/fusion` as ⚠️ **DEPRECATED**

**New Parameter Reference Table:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reaction` | string | *required* | SMILES reaction string |
| `reaction_family` | string | `null` | Reaction family (e.g., "suzuki", "buchwald") |
| `k` | integer | `50` | Number of precedents to retrieve |
| `rerank_strategy` | string | `'rule'` | Reranking method: `'none'`, `'rule'`, `'analytics'` |
| `filter_unknown_reagents` | boolean | `false` | Filter out recommendations with unknown reagents |
| `max_variants` | integer | `5` | Maximum condition variants to return |

---

## Key Improvements

### For Users

1. **Clear Feature Discovery**: "What's New" section immediately shows latest capabilities
2. **Practical Examples**: Real code showing how to use new parameters
3. **Migration Path**: Easy upgrade from deprecated fusion endpoint
4. **Decision Support**: Table showing when to use each reranking strategy

### For API Quality

1. **Catalyst Specificity**: Documents the fix that preserves exact catalyst complexes
2. **Backward Compatibility**: All examples work with both old and new code
3. **Best Practices**: Recommends `rerank_strategy='rule'` for standard reactions
4. **Deprecation Warnings**: Clear communication about fusion endpoint deprecation

---

## Testing

To verify the documentation examples work:

```powershell
# Start the server
uvicorn app.main:app --reload --port 8000

# Test basic recommendation (from Use Case 1)
python -c "import requests; r = requests.post('http://localhost:8000/api/v1/recommend', json={'reaction': 'CCBr.CCCO>>CCOCC', 'k': 50}); print(f'Found {len(r.json()[\"formatted\"][\"recommended_conditions\"])} recommendations')"

# Test with new parameters (from Advanced Example)
python -c "import requests; r = requests.post('http://localhost:8000/api/v1/recommend', json={'reaction': 'c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1c1ccccc1', 'reaction_family': 'suzuki', 'k': 100, 'rerank_strategy': 'rule', 'filter_unknown_reagents': True}); print(r.json()['formatted']['recommended_conditions'][0]['catalyst'])"
```

---

## Files Modified

- `docs/API_DOCUMENTATION.md` - Main API documentation (updated)

---

## Related Documentation

- `CATALYST_SPECIFICITY_FIX.md` - Technical details of catalyst preservation
- `WEB_CLI_UPDATE.md` - CLI feature parity documentation
- `TASKS_1_2_4_COMPLETE.md` - Previous session completion summary

---

**Summary**: Documentation now clearly shows how to use the improved recommendation endpoint with specific catalyst preservation, reranking strategies, and unknown reagent filtering. Users can quickly understand and adopt the new features.
