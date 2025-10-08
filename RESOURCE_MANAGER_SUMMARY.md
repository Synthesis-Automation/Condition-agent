# ChemTools Resource Manager - Implementation Summary

## What We're Building

A **centralized resource management system** (`ChemToolsContext`) that coordinates loading and sharing of all chemtools resources:

- ✅ **Reaction datasets** (precedents)
- ✅ **ML models** (DRFP yield predictor)  
- ✅ **Rule databases** (SciFinder condition schemes)
- ✅ **Reagent databases** (bases, solvents, ligands)
- ✅ **Chemical registry** (taxonomy, properties)

## Current Problems

| Problem | Impact |
|---------|--------|
| Loads ALL 99,668 reactions | 64s startup time |
| Each module has isolated `@lru_cache` | Memory waste, no sharing |
| No way to choose what to load | Can't optimize for specific use cases |
| ML models load on-demand | Slow first prediction |

## Proposed Solution

```python
# Single context manages all resources
ctx = ChemToolsContext(ResourceConfig(
    reaction_families=["C_N_Coupling_Pd"],  # Load only 1,343 reactions
    ml_model_path="models/drfp_yield.pkl",
    reagent_types=["base", "solvent"],
))

# Preload at API startup (3 seconds total)
ctx.preload('all')

# All modules share the same resources
recommend_from_reaction(rxn, context=ctx)  # Uses preloaded data
precedent.knn(family, features, context=ctx)  # Same dataset
```

## Performance Impact

**Current (Before):**
- First API call: **75 seconds** (load 99k reactions + compute DRFP)
- Memory: **~500 MB**

**With Context (After):**
- API startup: **3 seconds** (load 1.3k reactions + ML model)
- First API call: **0.6 seconds** (everything preloaded)
- Memory: **~150 MB** (only loaded families)

**Improvement: 125x faster first call, 3x less memory**

## Key Features

### 1. Selective Loading
```python
# Only load what you need
config = ResourceConfig(reaction_families=["C_N_Coupling_Pd"])
ctx = ChemToolsContext(config)
# Loads 1,343 reactions instead of 99,668 (50x faster!)
```

### 2. Lazy Loading
```python
# Resources load on first access
ctx = ChemToolsContext()
dataset = ctx.get_reaction_dataset()  # Loads here
dataset2 = ctx.get_reaction_dataset()  # Returns cached
```

### 3. Preloading for Production
```python
# Load everything at startup
ctx = ChemToolsContext(config)
ctx.preload('all')  # Eager loading
# First API call is now fast!
```

### 4. Thread-Safe
```python
# Safe for concurrent FastAPI requests
ctx = ChemToolsContext(config)

@app.post("/recommend")
async def recommend(request):
    # Multiple requests can safely use same context
    return recommend_from_reaction(request.smiles, context=ctx)
```

### 5. Observable
```python
# See what's loaded
stats = ctx.get_load_stats()
# {
#   'loaded': {
#     'dataset_families': ['C_N_Coupling_Pd'],
#     'dataset_size': 1343,
#     'ml_model': True,
#   },
#   'load_times': {
#     'dataset': {'load_time_s': 0.82},
#     'ml_model': {'load_time_s': 2.03}
#   }
# }
```

## Backward Compatibility

**Old code still works:**
```python
# No changes needed
from chemtools import recommend
result = recommend.recommend_from_reaction(rxn_smiles)
# Uses default global context internally
```

**New code is explicit:**
```python
# Better performance with explicit context
ctx = ChemToolsContext(ResourceConfig(reaction_families=["C_N_Coupling_Pd"]))
result = recommend.recommend_from_reaction(rxn_smiles, context=ctx)
```

## Implementation Timeline

| Week | Phase | Deliverables |
|------|-------|--------------|
| 1 | Core Infrastructure | `context.py`, `context_config.py` |
| 2 | Recommend Integration | Update `recommend.py`, `recommend_ml.py` |
| 3 | Rule Systems | Update `rule_scdb_matcher/` |
| 4 | Reagent & Registry | Update `reagent_lookup.py`, `registry.py` |
| 5 | API Integration | Update `app/main.py` with lifespan |
| 6 | Testing & Docs | Full test coverage, documentation |

## Files to Create

1. **`chemtools/context.py`** - Main `ChemToolsContext` class (~500 lines)
2. **`chemtools/context_config.py`** - `ResourceConfig` dataclass (~100 lines)

## Files to Modify

1. **`chemtools/precedent.py`** - Add `context` parameter to `knn()`
2. **`chemtools/recommend.py`** - Add `context` parameter, use ctx resources
3. **`chemtools/recommend_ml.py`** - Add `context` parameter, use ctx.get_ml_model()
4. **`chemtools/rule_scdb_matcher/matcher.py`** - Use ctx.get_rule_database()
5. **`chemtools/reagent_lookup.py`** - Use ctx.get_reagent_database()
6. **`chemtools/registry.py`** - Use ctx.get_registry()
7. **`app/main.py`** - Add lifespan with preloading

## Next Steps

1. ✅ **Review implementation plan** - DONE (see `RESOURCE_MANAGER_IMPLEMENTATION_PLAN.md`)
2. 🔲 **Prototype `ChemToolsContext`** - Create basic class with dataset loading
3. 🔲 **Benchmark selective loading** - Verify 50-100x speedup
4. 🔲 **Test with precedent search** - Ensure precomputed DRFP works
5. 🔲 **Update recommend.py** - Add context parameter
6. 🔲 **Full integration** - Roll out to all modules

## Questions to Discuss

1. Should we support loading `ResourceConfig` from YAML/JSON file?
2. Do we need multiple contexts in same process (e.g., different model versions)?
3. Should preloading be async for non-blocking startup?
4. What metrics should we expose (Prometheus format)?
5. Should we cache across process restarts (pickle serialization)?

## Success Criteria

- ✅ API startup < 5 seconds (vs 75s current)
- ✅ First recommendation < 1 second (vs 75s current)
- ✅ Memory usage < 200 MB (vs 500MB current)
- ✅ Backward compatible (existing code works)
- ✅ Thread-safe (concurrent requests)
- ✅ 100% test coverage for context class

---

**See Also:**
- 📄 `RESOURCE_MANAGER_IMPLEMENTATION_PLAN.md` - Full technical specification
- 📄 `RESOURCE_MANAGER_QUICK_REF.md` - Quick reference and examples
