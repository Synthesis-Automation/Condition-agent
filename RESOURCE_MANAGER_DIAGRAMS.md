# ChemTools Resource Manager - Architecture Diagrams

## Current Architecture (Without Resource Manager)

```
┌─────────────────────────────────────────────────────────────────┐
│                     FastAPI Application                         │
│                         (app/main.py)                           │
└────────────┬────────────────────────┬───────────────────────────┘
             │                        │
             ▼                        ▼
    ┌────────────────┐      ┌────────────────┐
    │  recommend.py  │      │ precedent.py   │
    │                │      │                │
    │ @lru_cache(1)  │      │ @lru_cache(1)  │
    │ _load()        │      │ _load()        │
    └────────┬───────┘      └────────┬───────┘
             │                       │
             │ Load ALL datasets     │ Load ALL datasets
             │ (99,668 reactions)    │ (99,668 reactions)
             │ 64 seconds            │ 64 seconds
             ▼                       ▼
    ┌────────────────────────────────────────┐
    │      data/reaction_dataset/            │
    │  ├─ C_N_Coupling_Pd.jsonl (1,343)      │
    │  ├─ C_N_Coupling_Cu.jsonl (10,234)     │
    │  ├─ Suzuki.jsonl (45,123)              │
    │  ├─ Amide_formation.jsonl (25,456)     │
    │  └─ C_N_Coupling_Ni.jsonl (17,512)     │
    └────────────────────────────────────────┘

    ┌────────────────┐      ┌────────────────┐
    │ recommend_ml.py│      │reagent_lookup  │
    │                │      │                │
    │ Loads ML model │      │ @lru_cache(32) │
    │ on-demand      │      │ load_reagent() │
    │ 2 seconds      │      │                │
    └────────┬───────┘      └────────┬───────┘
             │                       │
             ▼                       ▼
    ┌──────────────┐       ┌──────────────┐
    │ ML Models    │       │ Reagent DBs  │
    │ (DRFP)       │       │ (JSON files) │
    └──────────────┘       └──────────────┘

PROBLEMS:
❌ Each module loads dataset independently (duplicate work)
❌ Can't choose which families to load (all or nothing)
❌ First API call takes 75+ seconds (cold start)
❌ Memory waste from duplicate caches
❌ No coordination between modules
```

## Proposed Architecture (With Resource Manager)

```
┌─────────────────────────────────────────────────────────────────┐
│                     FastAPI Application                         │
│                         (app/main.py)                           │
│                                                                 │
│  Startup Event:                                                 │
│  ┌────────────────────────────────────────────┐                │
│  │ config = ResourceConfig(                   │                │
│  │   reaction_families=["C_N_Coupling_Pd"],   │                │
│  │   ml_model_path="models/drfp_yield.pkl",   │                │
│  │   reagent_types=["base", "solvent"],       │                │
│  │ )                                           │                │
│  │                                             │                │
│  │ ctx = ChemToolsContext(config)             │                │
│  │ ctx.preload('all')  # 3 seconds total      │                │
│  └────────────────────────────────────────────┘                │
│                                                                 │
│  API Endpoints:                                                 │
│  ┌────────────────────────────────────────────┐                │
│  │ @app.post("/recommend")                    │                │
│  │ async def recommend(request):              │                │
│  │     return recommend_from_reaction(        │                │
│  │         request.smiles,                    │                │
│  │         context=ctx  # Use shared context  │                │
│  │     )                                       │                │
│  └────────────────────────────────────────────┘                │
└─────────────────────────────┬───────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│              ChemToolsContext (Resource Manager)                │
│                    (chemtools/context.py)                       │
│                                                                 │
│  Configuration:                                                 │
│  ┌────────────────────────────────────────────┐                │
│  │ reaction_families: ["C_N_Coupling_Pd"]     │                │
│  │ ml_model_path: "models/drfp_yield.pkl"     │                │
│  │ load_drfp_precomputed: True                │                │
│  │ reagent_types: ["base", "solvent"]         │                │
│  └────────────────────────────────────────────┘                │
│                                                                 │
│  Cached Resources (Thread-Safe with RLock):                    │
│  ┌────────────────────────────────────────────┐                │
│  │ _reaction_dataset: List[Dict] (1,343)      │                │
│  │ _ml_model: DRFPYieldPredictor              │                │
│  │ _rule_databases: Dict[str, RuleDB]         │                │
│  │ _reagent_databases: Dict[str, List[Dict]]  │                │
│  │ _registry: RegistryIndex                   │                │
│  └────────────────────────────────────────────┘                │
│                                                                 │
│  Public API:                                                    │
│  ┌────────────────────────────────────────────┐                │
│  │ get_reaction_dataset(families) → List[Dict]│                │
│  │ get_ml_model() → DRFPYieldPredictor        │                │
│  │ get_rule_database(rxn_type) → RuleDB       │                │
│  │ get_reagent_database(type) → List[Dict]    │                │
│  │ get_registry() → RegistryIndex             │                │
│  │ preload(*resources) → None                 │                │
│  │ get_load_stats() → Dict                    │                │
│  │ clear() → None                             │                │
│  └────────────────────────────────────────────┘                │
└────┬────────────┬────────────┬────────────┬──────────┬─────────┘
     │            │            │            │          │
     ▼            ▼            ▼            ▼          ▼
┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐
│recommend │ │precedent │ │recommend │ │  rule_   │ │ reagent_ │
│   .py    │ │   .py    │ │  _ml.py  │ │  scdb_   │ │ lookup.py│
│          │ │          │ │          │ │  matcher │ │          │
│ Uses ctx │ │ Uses ctx │ │ Uses ctx │ │ Uses ctx │ │ Uses ctx │
│ No cache │ │ No cache │ │ No cache │ │ No cache │ │ No cache │
└────┬─────┘ └────┬─────┘ └────┬─────┘ └────┬─────┘ └────┬─────┘
     │            │            │            │            │
     └────────────┴────────────┴────────────┴────────────┘
                              │
                    All modules share same
                    resources from context

BENEFITS:
✅ Load only what you need (1,343 vs 99,668 reactions)
✅ Load once, share everywhere (no duplicate work)
✅ Fast API startup (3s preload vs 75s first call)
✅ Thread-safe (concurrent requests)
✅ Observable (load stats, metrics)
```

## Data Flow - Before vs After

### BEFORE: First API Request (Cold Start)

```
Time: 0s
┌────────────────────┐
│ API Request        │
│ /recommend         │
└─────────┬──────────┘
          │
          ▼
Time: 0s-64s
┌────────────────────┐
│ precedent.py       │
│ _load()            │
│ Loading ALL        │
│ datasets...        │
│ (99,668 reactions) │
└─────────┬──────────┘
          │
          ▼
Time: 64s-73s
┌────────────────────┐
│ precedent.knn()    │
│ Computing DRFP     │
│ on-demand for      │
│ 1,343 reactions... │
│ (9 seconds)        │
└─────────┬──────────┘
          │
          ▼
Time: 73s-75s
┌────────────────────┐
│ recommend.py       │
│ Processing results │
│ (2 seconds)        │
└─────────┬──────────┘
          │
          ▼
Time: 75s
┌────────────────────┐
│ Response           │
│ (75 seconds total) │
└────────────────────┘
```

### AFTER: API Startup + First Request

```
Time: 0s (API Startup)
┌────────────────────┐
│ API Startup        │
│ lifespan event     │
└─────────┬──────────┘
          │
          ▼
Time: 0s-0.8s
┌────────────────────┐
│ ChemToolsContext   │
│ Loading C_N_Pd     │
│ dataset only       │
│ (1,343 reactions)  │
│ with precomputed   │
│ DRFP fingerprints  │
└─────────┬──────────┘
          │
          ▼
Time: 0.8s-2.8s
┌────────────────────┐
│ Loading ML model   │
│ DRFP predictor     │
│ (2 seconds)        │
└─────────┬──────────┘
          │
          ▼
Time: 2.8s-3.0s
┌────────────────────┐
│ Loading reagent    │
│ databases          │
│ (0.2 seconds)      │
└─────────┬──────────┘
          │
          ▼
Time: 3.0s
┌────────────────────┐
│ API Ready!         │
│ (3 seconds total)  │
└────────────────────┘

═══════════════════════════════════════

Time: 3.0s
┌────────────────────┐
│ First API Request  │
│ /recommend         │
└─────────┬──────────┘
          │
          ▼
Time: 3.0s-3.1s
┌────────────────────┐
│ precedent.knn()    │
│ Dataset already    │
│ in memory!         │
│ DRFP precomputed!  │
│ (0.1 seconds)      │
└─────────┬──────────┘
          │
          ▼
Time: 3.1s-3.6s
┌────────────────────┐
│ ML prediction      │
│ Model already      │
│ loaded!            │
│ (0.5 seconds)      │
└─────────┬──────────┘
          │
          ▼
Time: 3.6s
┌────────────────────┐
│ Response           │
│ (0.6 seconds)      │
└────────────────────┘

Total time from startup to response: 3.6s
(vs 75s without context)
```

## Memory Comparison

### BEFORE: Multiple Caches

```
┌────────────────────────────────────┐
│ Python Process Memory: ~500 MB     │
├────────────────────────────────────┤
│                                    │
│ precedent.py @lru_cache:           │
│ └─ All datasets: 250 MB            │
│                                    │
│ recommend.py @lru_cache:           │
│ └─ All datasets: 250 MB (duplicate)│
│                                    │
│ ML model cache: 50 MB              │
│                                    │
│ DRFP fingerprint cache: 100 MB     │
│                                    │
│ Reagent DB cache: 10 MB            │
│                                    │
└────────────────────────────────────┘
```

### AFTER: Shared Resources

```
┌────────────────────────────────────┐
│ Python Process Memory: ~150 MB     │
├────────────────────────────────────┤
│                                    │
│ ChemToolsContext:                  │
│ ├─ C_N_Pd dataset: 80 MB (only 1!) │
│ ├─ ML model: 50 MB                 │
│ ├─ DRFP cache: 10 MB (smaller)     │
│ └─ Reagent DBs: 10 MB              │
│                                    │
│ Module code: ~20 MB                │
│                                    │
└────────────────────────────────────┘

Memory savings: 350 MB (70% reduction)
```

## Thread Safety Design

```
┌─────────────────────────────────────────────────┐
│          ChemToolsContext                       │
│                                                 │
│  ┌─────────────────────────────────────────┐   │
│  │ self._lock = RLock()                    │   │
│  │ (Reentrant Lock)                        │   │
│  └─────────────────────────────────────────┘   │
│                                                 │
│  def get_reaction_dataset(self, families):     │
│      with self._lock:  ← Thread-safe           │
│          if self._reaction_dataset is None:    │
│              # Load from disk                  │
│              self._reaction_dataset = load()   │
│          return self._reaction_dataset         │
└─────────────────────────────────────────────────┘
                          │
          ┌───────────────┼───────────────┐
          │               │               │
          ▼               ▼               ▼
    ┌──────────┐    ┌──────────┐    ┌──────────┐
    │ Thread 1 │    │ Thread 2 │    │ Thread 3 │
    │ Request  │    │ Request  │    │ Request  │
    └────┬─────┘    └────┬─────┘    └────┬─────┘
         │               │               │
         └───────────────┴───────────────┘
                         │
         All threads see same dataset
         No race conditions
         Only loaded once
```

## Configuration Options

```
ResourceConfig
├─ reaction_families: Optional[List[str]]
│  ├─ None = Load all families (99,668 reactions)
│  ├─ [] = Load no families (dataset disabled)
│  └─ ["C_N_Coupling_Pd"] = Load specific (1,343 reactions)
│
├─ load_drfp_precomputed: bool = True
│  └─ Use precomputed DRFP from dataset (50x faster)
│
├─ ml_model_path: Optional[str] = None
│  ├─ None = Disable ML predictions
│  └─ "models/drfp.pkl" = Enable with model
│
├─ drfp_fingerprint_path: Optional[str] = None
│  └─ Path to precomputed .npz fingerprints
│
├─ rule_db_paths: Optional[List[str]] = None
│  └─ Custom paths to rule database JSON files
│
├─ reagent_types: Optional[List[str]] = None
│  ├─ None = Load on-demand
│  └─ ["base", "solvent"] = Preload specific
│
├─ load_registry: bool = True
│  └─ Chemical taxonomy (names, CAS, properties)
│
└─ enable_caching: bool = True
   └─ DRFP encoding cache, KNN cache, etc.
```

## Common Usage Patterns

### Pattern 1: API Production (Preload Everything)

```
┌───────────────────────────────────────┐
│ FastAPI Startup                       │
│                                       │
│ config = ResourceConfig(              │
│   reaction_families=["C_N_Pd"],       │
│   ml_model_path="models/drfp.pkl",    │
│   reagent_types=["base", "solvent"],  │
│ )                                     │
│                                       │
│ ctx = ChemToolsContext(config)        │
│ ctx.preload('all')  ← Load now        │
└─────────────┬─────────────────────────┘
              │
              ▼
┌───────────────────────────────────────┐
│ Every API Request                     │
│                                       │
│ result = recommend_from_reaction(     │
│     rxn_smiles,                       │
│     context=ctx  ← Use preloaded      │
│ )                                     │
│                                       │
│ Fast! (0.6s)                          │
└───────────────────────────────────────┘
```

### Pattern 2: Jupyter Notebook (Lazy Loading)

```
┌───────────────────────────────────────┐
│ Cell 1: Import                        │
│                                       │
│ from chemtools.context import         │
│     ChemToolsContext                  │
│                                       │
│ ctx = ChemToolsContext()              │
│ # Nothing loaded yet (instant)        │
└───────────────────────────────────────┘

┌───────────────────────────────────────┐
│ Cell 2: First Use                     │
│                                       │
│ result = recommend_from_reaction(     │
│     "Brc1ccccc1.Nc1ccccc1>>...",      │
│     context=ctx                       │
│ )                                     │
│                                       │
│ # Dataset loads here (0.8s)           │
│ # Then cached for future cells        │
└───────────────────────────────────────┘

┌───────────────────────────────────────┐
│ Cell 3: Subsequent Use                │
│                                       │
│ result2 = recommend_from_reaction(    │
│     "another reaction",               │
│     context=ctx                       │
│ )                                     │
│                                       │
│ # Uses cached dataset (fast!)         │
└───────────────────────────────────────┘
```

### Pattern 3: CLI Script (Minimal Loading)

```
┌───────────────────────────────────────┐
│ script.py                             │
│                                       │
│ config = ResourceConfig(              │
│   reaction_families=["C_N_Pd"],       │
│   ml_model_path=None,  ← No ML        │
│   load_drfp_precomputed=True,         │
│ )                                     │
│                                       │
│ ctx = ChemToolsContext(config)        │
│                                       │
│ for rxn in batch_reactions:           │
│     result = recommend(rxn, ctx)      │
│     save(result)                      │
│                                       │
│ ctx.clear()  ← Free memory            │
└───────────────────────────────────────┘
```

---

**Document**: Architecture Diagrams  
**Version**: 1.0  
**Date**: 2025-10-07

**See Also:**
- `RESOURCE_MANAGER_IMPLEMENTATION_PLAN.md` - Full technical spec
- `RESOURCE_MANAGER_QUICK_REF.md` - Quick reference
- `RESOURCE_MANAGER_SUMMARY.md` - Executive summary
