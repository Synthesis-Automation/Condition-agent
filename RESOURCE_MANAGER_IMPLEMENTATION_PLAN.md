# ChemTools Resource Manager Implementation Plan

## Executive Summary
Create a centralized resource management system (`ChemToolsContext`) to coordinate loading, caching, and sharing of datasets, ML models, and rule databases across all chemtools modules.

**Current Problems:**
- Each module loads resources independently with isolated `@lru_cache`
- `precedent.py` loads ALL 99,668 reactions (64s) even when only 1,343 needed
- No shared caching between `recommend.py`, `recommend_ml.py`, and `precedent.py`
- ML models loaded on-demand causing slow first predictions
- Rule databases loaded separately by `recommend.py` and `rule_scdb_matcher/`
- Memory duplication across modules

**Target Benefits:**
- **50-100x faster startup**: Load only needed families (1.3k vs 99k reactions)
- **Unified caching**: Share datasets/models across all modules
- **Lazy loading**: Load resources only when needed
- **Memory efficient**: Single instance of each resource
- **Explicit control**: API users can preload or defer loading
- **Thread-safe**: Support concurrent API requests

---

## Architecture Design

### 1. Core Class: `ChemToolsContext`

```python
# chemtools/context.py

from typing import Dict, List, Any, Optional, Set
from dataclasses import dataclass, field
from threading import RLock
from pathlib import Path

@dataclass
class ResourceConfig:
    """Configuration for what resources to load."""
    # Datasets
    reaction_families: Optional[List[str]] = None  # None = all, [] = none
    load_drfp_precomputed: bool = True
    
    # ML Models
    ml_model_path: Optional[str] = None
    drfp_fingerprint_path: Optional[str] = None
    
    # Rule databases
    rule_db_paths: Optional[List[str]] = None
    
    # Reagent databases
    reagent_types: Optional[List[str]] = None  # ["base", "solvent", "ligand"]
    
    # Registry (chemical taxonomy)
    load_registry: bool = True
    
    # Performance
    enable_caching: bool = True
    cache_sizes: Dict[str, int] = field(default_factory=lambda: {
        "drfp_encode": 200000,
        "knn_search": 512,
        "reagent_lookup": 32,
    })


class ChemToolsContext:
    """
    Centralized resource manager for chemtools.
    
    Manages loading, caching, and lifecycle of:
    - Reaction datasets (precedents)
    - ML models (DRFP yield predictor)
    - Rule databases (SciFinder condition schemes, selectors)
    - Reagent databases (bases, solvents, ligands, catalysts)
    - Chemical registry (taxonomy, properties)
    
    Thread-safe for concurrent API requests.
    """
    
    def __init__(self, config: Optional[ResourceConfig] = None):
        self.config = config or ResourceConfig()
        self._lock = RLock()
        
        # Cached resources (lazy-loaded)
        self._reaction_dataset: Optional[List[Dict[str, Any]]] = None
        self._ml_model: Optional[Any] = None  # DRFPYieldPredictor
        self._rule_databases: Dict[str, Any] = {}
        self._reagent_databases: Dict[str, List[Dict]] = {}
        self._registry: Optional[Any] = None
        
        # Tracking
        self._loaded_families: Set[str] = set()
        self._load_stats: Dict[str, Any] = {}
    
    # === Dataset Management ===
    
    def get_reaction_dataset(
        self, 
        families: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Get reaction dataset, loading on first access.
        
        Args:
            families: Specific families to load. If None, uses config.
        
        Returns:
            List of reaction records with precomputed features/DRFP.
        """
        with self._lock:
            target_families = families or self.config.reaction_families
            
            # Check if already loaded
            if self._reaction_dataset is not None:
                if target_families is None or self._loaded_families >= set(target_families):
                    return self._reaction_dataset
            
            # Load from disk
            import time
            start = time.perf_counter()
            
            from .precedent import _load_selective
            self._reaction_dataset = _load_selective(target_families)
            
            if target_families:
                self._loaded_families.update(target_families)
            
            elapsed = time.perf_counter() - start
            self._load_stats['dataset'] = {
                'families': target_families or 'all',
                'n_reactions': len(self._reaction_dataset),
                'load_time_s': elapsed,
            }
            
            return self._reaction_dataset
    
    def get_ml_model(self) -> Optional[Any]:
        """Get ML yield prediction model, loading on first access."""
        with self._lock:
            if self._ml_model is not None:
                return self._ml_model
            
            if not self.config.ml_model_path:
                return None
            
            import time
            start = time.perf_counter()
            
            from .ml.drfp_predictor import DRFPYieldPredictor
            self._ml_model = DRFPYieldPredictor.load(self.config.ml_model_path)
            
            elapsed = time.perf_counter() - start
            self._load_stats['ml_model'] = {
                'path': self.config.ml_model_path,
                'load_time_s': elapsed,
            }
            
            return self._ml_model
    
    def get_rule_database(self, reaction_type: str) -> Optional[Any]:
        """Get rule database for specific reaction type."""
        with self._lock:
            if reaction_type in self._rule_databases:
                return self._rule_databases[reaction_type]
            
            # Load rule DB
            from .rule_scdb_matcher.loader import load_db
            from .rule_scdb_matcher.matcher import find_matching_db
            
            # Find matching DB file
            db_path = find_matching_db(reaction_type, self.config.rule_db_paths)
            if not db_path:
                return None
            
            db = load_db(db_path)
            self._rule_databases[reaction_type] = db
            
            return db
    
    def get_reagent_database(self, reagent_type: str) -> List[Dict[str, Any]]:
        """Get reagent database (base, solvent, ligand, etc.)."""
        with self._lock:
            if reagent_type in self._reagent_databases:
                return self._reagent_databases[reagent_type]
            
            from .reagent_lookup import load_reagent_database
            db = load_reagent_database(reagent_type)
            self._reagent_databases[reagent_type] = db
            
            return db
    
    def get_registry(self) -> Any:
        """Get chemical registry (taxonomy)."""
        with self._lock:
            if self._registry is not None:
                return self._registry
            
            if not self.config.load_registry:
                return None
            
            from .registry import _load_registry
            self._registry = _load_registry()
            
            return self._registry
    
    # === Preloading ===
    
    def preload(self, *resource_types: str):
        """
        Eagerly load specified resources.
        
        Args:
            resource_types: One or more of:
                'dataset', 'ml_model', 'rule_databases', 
                'reagent_databases', 'registry', 'all'
        """
        if 'all' in resource_types:
            resource_types = ('dataset', 'ml_model', 'rule_databases', 
                             'reagent_databases', 'registry')
        
        for rt in resource_types:
            if rt == 'dataset':
                self.get_reaction_dataset()
            elif rt == 'ml_model':
                self.get_ml_model()
            elif rt == 'registry':
                self.get_registry()
            elif rt == 'reagent_databases' and self.config.reagent_types:
                for rtype in self.config.reagent_types:
                    self.get_reagent_database(rtype)
    
    # === Statistics ===
    
    def get_load_stats(self) -> Dict[str, Any]:
        """Get statistics about loaded resources."""
        return {
            'config': {
                'reaction_families': self.config.reaction_families,
                'ml_model_enabled': self.config.ml_model_path is not None,
            },
            'loaded': {
                'dataset_families': list(self._loaded_families),
                'dataset_size': len(self._reaction_dataset) if self._reaction_dataset else 0,
                'ml_model': self._ml_model is not None,
                'rule_databases': list(self._rule_databases.keys()),
                'reagent_databases': list(self._reagent_databases.keys()),
                'registry': self._registry is not None,
            },
            'load_times': self._load_stats,
        }
    
    # === Cleanup ===
    
    def clear(self):
        """Clear all cached resources to free memory."""
        with self._lock:
            self._reaction_dataset = None
            self._ml_model = None
            self._rule_databases.clear()
            self._reagent_databases.clear()
            self._registry = None
            self._loaded_families.clear()
            self._load_stats.clear()


# Global singleton (optional - can also pass context explicitly)
_default_context: Optional[ChemToolsContext] = None


def get_default_context() -> ChemToolsContext:
    """Get or create the default global context."""
    global _default_context
    if _default_context is None:
        _default_context = ChemToolsContext()
    return _default_context


def set_default_context(context: ChemToolsContext):
    """Set the default global context."""
    global _default_context
    _default_context = context


def reset_default_context():
    """Clear the default global context."""
    global _default_context
    if _default_context:
        _default_context.clear()
    _default_context = None
```

---

## 2. Module Integration Plan

### Phase 1: Core Infrastructure (Week 1)

**Files to create:**
- `chemtools/context.py` - ResourceManager core class
- `chemtools/context_config.py` - Configuration dataclasses

**Files to modify:**
- `chemtools/__init__.py` - Export context classes
- `chemtools/precedent.py` - Add context parameter to functions

**Changes:**

```python
# chemtools/precedent.py

def knn(
    family: str, 
    features: Dict[str, Any], 
    k: int = 50, 
    relax: Dict[str, Any] | None = None,
    context: Optional[ChemToolsContext] = None,  # NEW
) -> Dict[str, Any]:
    """KNN precedent search with optional context."""
    ctx = context or get_default_context()
    
    # Use context to get dataset
    rows = ctx.get_reaction_dataset(families=[_family_text(family)])
    
    # Rest of logic unchanged...
```

### Phase 2: Recommend Module (Week 2)

**Files to modify:**
- `chemtools/recommend.py`
- `chemtools/recommend_ml.py`

**Key changes:**
```python
# chemtools/recommend.py

def recommend_from_reaction(
    reaction_smiles: str,
    k: int = 10,
    context: Optional[ChemToolsContext] = None,  # NEW
    **kwargs
) -> Dict[str, Any]:
    ctx = context or get_default_context()
    
    # Use context for:
    # - Precedent search (via precedent.knn with context)
    # - Rule databases (ctx.get_rule_database)
    # - Reagent lookup (ctx.get_reagent_database)
```

```python
# chemtools/recommend_ml.py

def hybrid_recommend(
    reaction_smiles: str,
    context: Optional[ChemToolsContext] = None,  # NEW
    **kwargs,
) -> Dict[str, Any]:
    ctx = context or get_default_context()
    
    # Use context for:
    # - ML model (ctx.get_ml_model())
    # - Precedent search (via recommend_from_reaction with context)
```

### Phase 3: Rule-Based Systems (Week 3)

**Files to modify:**
- `chemtools/rule_scdb_matcher/matcher.py`
- `chemtools/rule_scdb_matcher/loader.py`

**Add helper:**
```python
# chemtools/rule_scdb_matcher/matcher.py

def find_matching_db(
    reaction_type: str,
    search_paths: Optional[List[str]] = None,
) -> Optional[Path]:
    """Find rule database file for reaction type."""
    # Logic to locate DB files
```

### Phase 4: Reagent & Registry (Week 4)

**Files to modify:**
- `chemtools/reagent_lookup.py` - Use context for caching
- `chemtools/registry.py` - Use context for registry

**Remove module-level `@lru_cache`, use context instead**

### Phase 5: API Integration (Week 5)

**Files to modify:**
- `app/main.py` - FastAPI application

**Add startup event:**
```python
# app/main.py

from chemtools.context import ChemToolsContext, ResourceConfig
from contextlib import asynccontextmanager

# Global context for API
api_context: Optional[ChemToolsContext] = None

@asynccontextmanager
async def lifespan(app: FastAPI):
    global api_context
    
    # Startup: Initialize context
    config = ResourceConfig(
        reaction_families=["C_N_Coupling_Pd", "Suzuki"],  # Only load what we need
        load_drfp_precomputed=True,
        ml_model_path="models/drfp_yield_v1.pkl",
        reagent_types=["base", "solvent", "ligand"],
    )
    
    api_context = ChemToolsContext(config)
    
    # Preload critical resources
    api_context.preload('dataset', 'ml_model', 'reagent_databases')
    
    print(f"Loaded resources: {api_context.get_load_stats()}")
    
    yield
    
    # Shutdown: Cleanup
    api_context.clear()

app = FastAPI(lifespan=lifespan)

@app.post("/recommend")
async def recommend_endpoint(request: RecommendRequest):
    # Use global context
    result = recommend_from_reaction(
        request.reaction_smiles,
        context=api_context,
    )
    return result
```

---

## 3. Migration Strategy

### Backward Compatibility

**Keep existing APIs working without context:**

```python
# Old code still works
from chemtools import precedent
results = precedent.knn("C_N_Coupling_Pd", features, k=10)

# New code with context
from chemtools.context import ChemToolsContext, ResourceConfig
ctx = ChemToolsContext(ResourceConfig(reaction_families=["C_N_Coupling_Pd"]))
results = precedent.knn("C_N_Coupling_Pd", features, k=10, context=ctx)
```

**Default behavior:**
- If `context=None`, use `get_default_context()` (auto-creates singleton)
- Default context loads resources lazily on first use
- Existing `@lru_cache` decorators stay until context is widely adopted

### Deprecation Path

1. **Phase 1-4**: Add context parameters, keep old behavior as default
2. **Phase 5**: Update API to use explicit context
3. **Phase 6** (future): Deprecate module-level caches, require context

---

## 4. Testing Plan

### Unit Tests

```python
# tests/chemtools/test_context.py

def test_context_loads_dataset_on_demand():
    config = ResourceConfig(reaction_families=["C_N_Coupling_Pd"])
    ctx = ChemToolsContext(config)
    
    assert ctx._reaction_dataset is None  # Not loaded yet
    
    dataset = ctx.get_reaction_dataset()
    assert len(dataset) > 0
    assert ctx._reaction_dataset is not None  # Now cached


def test_context_selective_loading_performance():
    # Load only one family
    config1 = ResourceConfig(reaction_families=["C_N_Coupling_Pd"])
    ctx1 = ChemToolsContext(config1)
    
    start = time.perf_counter()
    ds1 = ctx1.get_reaction_dataset()
    time1 = time.perf_counter() - start
    
    # Load all families
    config2 = ResourceConfig(reaction_families=None)
    ctx2 = ChemToolsContext(config2)
    
    start = time.perf_counter()
    ds2 = ctx2.get_reaction_dataset()
    time2 = time.perf_counter() - start
    
    assert len(ds1) < len(ds2)
    assert time1 < time2 / 10  # At least 10x faster


def test_context_thread_safety():
    ctx = ChemToolsContext()
    
    def worker():
        dataset = ctx.get_reaction_dataset()
        return len(dataset)
    
    import concurrent.futures
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(worker) for _ in range(100)]
        results = [f.result() for f in futures]
    
    assert all(r == results[0] for r in results)  # All same length
```

### Integration Tests

```python
# tests/integration/test_context_recommend.py

def test_recommend_with_context():
    config = ResourceConfig(
        reaction_families=["C_N_Coupling_Pd"],
        ml_model_path="models/drfp_yield_v1.pkl",
    )
    ctx = ChemToolsContext(config)
    ctx.preload('dataset', 'ml_model')
    
    result = recommend_from_reaction(
        "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        context=ctx,
    )
    
    assert 'recommended_conditions' in result
    assert ctx.get_load_stats()['loaded']['dataset_size'] > 0
```

---

## 5. Performance Benchmarks

### Before (Current System)

```
First API call (cold start):
- Load all datasets: 64s (99,668 reactions)
- Load ML model: 2s
- First recommendation: 75s total

Subsequent calls:
- Cached precedent search: 0.3s
- ML prediction: 0.5s
```

### After (With Context)

```
API startup with preloading:
- Load C_N_Coupling_Pd only: 0.8s (1,343 reactions)
- Load ML model: 2s
- Startup total: 3s

First API call:
- Precedent search: 0.1s (already loaded)
- ML prediction: 0.5s (already loaded)
- Total: 0.6s

Subsequent calls:
- Same as first: 0.6s
```

**Improvement: 125x faster first call, 3s startup vs 75s first-call**

---

## 6. Implementation Timeline

| Week | Phase | Deliverables |
|------|-------|--------------|
| 1 | Core Infrastructure | `context.py`, `context_config.py` |
| 2 | Recommend Integration | Update `recommend.py`, `recommend_ml.py` |
| 3 | Rule Systems | Update `rule_scdb_matcher/` |
| 4 | Reagent & Registry | Update `reagent_lookup.py`, `registry.py` |
| 5 | API Integration | Update `app/main.py` with preloading |
| 6 | Testing & Docs | Full test coverage, documentation |

---

## 7. Benefits Summary

### Performance
- ✅ **50-100x faster startup**: Load only needed families
- ✅ **Sub-second API responses**: Resources preloaded
- ✅ **Memory efficient**: Single instance per resource

### Architecture
- ✅ **Centralized control**: One place to manage all resources
- ✅ **Thread-safe**: Support concurrent requests
- ✅ **Explicit lifecycle**: Clear startup/shutdown

### Developer Experience
- ✅ **Easy configuration**: Simple config object
- ✅ **Backward compatible**: Existing code still works
- ✅ **Testable**: Easy to mock context in tests

### Production Readiness
- ✅ **Observability**: Load stats and metrics
- ✅ **Resource cleanup**: Explicit `clear()` method
- ✅ **Flexible deployment**: API vs CLI vs notebook

---

## 8. Example Usage

### FastAPI Production

```python
# app/main.py

from chemtools.context import ChemToolsContext, ResourceConfig

# Create context at startup
config = ResourceConfig(
    reaction_families=["C_N_Coupling_Pd", "Suzuki"],
    ml_model_path="models/drfp_yield_v1.pkl",
    reagent_types=["base", "solvent", "ligand"],
)

api_context = ChemToolsContext(config)
api_context.preload('all')  # Eager load everything

@app.post("/recommend")
async def recommend(request: RecommendRequest):
    return recommend_from_reaction(
        request.reaction_smiles,
        context=api_context,
    )
```

### Jupyter Notebook

```python
# notebook.ipynb

from chemtools.context import ChemToolsContext, ResourceConfig

# Lazy loading for interactive work
ctx = ChemToolsContext()

# Resources load on first use
results = recommend_from_reaction(
    "Brc1ccccc1.Nc1ccccc1>>...",
    context=ctx,
)

# Check what's loaded
print(ctx.get_load_stats())
```

### CLI Script

```python
# scripts/batch_recommend.py

from chemtools.context import ChemToolsContext, ResourceConfig

def main():
    # Load only what we need
    config = ResourceConfig(
        reaction_families=["C_N_Coupling_Pd"],
        ml_model_path=None,  # Disable ML
    )
    
    ctx = ChemToolsContext(config)
    
    for rxn in reactions:
        result = recommend_from_reaction(rxn, context=ctx)
        print(result)
    
    # Cleanup
    ctx.clear()
```

---

## 9. Open Questions

1. **Caching strategy**: Should context use same `@lru_cache` sizes as current modules?
2. **Config file**: Support loading `ResourceConfig` from YAML/JSON?
3. **Multiple contexts**: Support multiple contexts in same process (e.g., different model versions)?
4. **Async loading**: Should `preload()` support async for non-blocking startup?
5. **Metrics**: Should context expose Prometheus metrics for monitoring?

---

## 10. Next Steps

1. **Review this plan** with team
2. **Prototype `context.py`** with basic dataset loading
3. **Benchmark** selective loading vs full loading
4. **Update `precedent.py`** to accept context parameter
5. **Create test harness** for context functionality
6. **Iterate** based on feedback

---

**Document Version**: 1.0  
**Author**: AI Assistant  
**Date**: 2025-10-07  
**Status**: Proposed
