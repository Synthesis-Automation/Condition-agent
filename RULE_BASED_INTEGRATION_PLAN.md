# Rule-Based System Integration Plan

**Problem**: Rule-based recommendation system (`rule_scdb_matcher/`) is currently NOT integrated into ChemTools v2.0, creating inconsistency with other subsystems.

**Solution**: Integrate it properly with the same pattern as ML recommendations, precedent search, and reagent lookup.

---

## Current Issues

### 1. Import Path Mismatch ❌

**Directory name**: `chemtools/rule_scdb_matcher/`  
**App imports**: `from chemtools.scdb_matcher import ...`

This mismatch means the imports are **broken** in:
- `app/main.py` (line 41)
- `app/ui_simple.py` (line 51)
- `app/ui_gradio.py` (line 59)

### 2. No ChemTools Integration ❌

- NOT exported from `chemtools/__init__.py`
- NO `RuleNamespace` class in `context.py`
- NOT accessible via `chem.*` namespace
- Manual database loading required

### 3. Inconsistent Architecture ❌

| System | Integration Status |
|--------|-------------------|
| SMILES normalization | ✅ `chem.smiles.*` |
| Router | ✅ `chem.router.*` |
| Precedent search | ✅ `chem.precedent.*` |
| ML recommendations | ✅ `chem.recommend.*` |
| Reagent lookup | ✅ `chem.reagent.*` |
| **Rule-based** | ❌ **NOT integrated** |

---

## Integration Solution

### Phase 1: Fix Import Path

**Option A**: Rename directory to match imports (breaking change)
```powershell
# Rename rule_scdb_matcher → scdb_matcher
mv chemtools/rule_scdb_matcher chemtools/scdb_matcher
```

**Option B**: Create alias in `__init__.py` (non-breaking)
```python
# chemtools/__init__.py
from . import rule_scdb_matcher as scdb_matcher
```

**Recommendation**: Option B (non-breaking, backward compatible)

### Phase 2: Create RuleNamespace Class

**Add to `chemtools/context.py`**:

```python
class RuleNamespace:
    """
    Rule-based scheme matching for reaction conditions.
    
    Provides access to SMARTS-based deterministic matching against
    curated reaction schemes.
    
    Example:
        >>> from chemtools import chem
        >>> db = chem.rules.load_database("cn_coupling_pd_db.json")
        >>> result = chem.rules.match(db, reaction)
        >>> print(result.conditions)
    """
    
    def __init__(self):
        """Initialize rule namespace with lazy database loading."""
        self._databases = {}  # Cache for loaded databases
    
    def load_database(
        self, 
        path: str,
        cache: bool = True
    ) -> Any:
        """
        Load a scheme database from JSON file.
        
        Args:
            path: Path to database JSON file (absolute or relative to data/conditionDB/)
            cache: Whether to cache the loaded database (default: True)
            
        Returns:
            SchemeConditionDB or SelectorRuleDB instance
            
        Example:
            >>> db = chem.rules.load_database("cn_coupling_pd_db.json")
            >>> db = chem.rules.load_database("data/conditionDB/suzuki_db.json")
        """
        from .rule_scdb_matcher import load_db
        from pathlib import Path
        
        # Normalize path
        p = Path(path)
        if not p.is_absolute():
            # Try relative to data/conditionDB/
            candidate = Path("data/conditionDB") / path
            if candidate.exists():
                p = candidate.resolve()
            else:
                p = Path(path).resolve()
        else:
            p = p.resolve()
        
        path_str = str(p)
        
        # Check cache
        if cache and path_str in self._databases:
            return self._databases[path_str]
        
        # Load database
        db = load_db(path_str)
        
        # Cache if requested
        if cache:
            self._databases[path_str] = db
        
        return db
    
    def match(
        self,
        db: Any,
        reaction: str = None,
        *,
        features: Dict[str, Any] = None
    ) -> Any:
        """
        Match a reaction against a rule database.
        
        Args:
            db: Loaded database from load_database()
            reaction: Reaction SMILES (for SchemeConditionDB)
            features: Feature dict (for SelectorRuleDB)
            
        Returns:
            MatchResult with conditions, trace, and metadata
            
        Example:
            >>> db = chem.rules.load_database("cn_coupling_pd_db.json")
            >>> result = chem.rules.match(db, "c1ccccc1Br.Nc1ccccc1>>c1ccccc1Nc1ccccc1")
            >>> print(result.conditions)
        """
        from .rule_scdb_matcher import match
        return match(db, reaction, features=features)
    
    def list_databases(self) -> List[str]:
        """
        List available databases in data/conditionDB/.
        
        Returns:
            List of database filenames
        """
        from pathlib import Path
        db_dir = Path("data/conditionDB")
        if db_dir.exists():
            return sorted([f.name for f in db_dir.glob("*.json")])
        return []
    
    def clear_cache(self):
        """Clear cached databases to free memory."""
        self._databases.clear()
```

### Phase 3: Add to ChemTools Class

**Update `ChemTools.__init__()` in `context.py`**:

```python
class ChemTools:
    def __init__(self, config: ResourceConfig = None):
        # ... existing code ...
        
        # Add rule-based namespace
        self.rules = RuleNamespace()
```

### Phase 4: Export from Package

**Update `chemtools/__init__.py`**:

```python
# Add alias for backward compatibility
from . import rule_scdb_matcher as scdb_matcher

# ChemTools now includes chem.rules.*
from .context import ChemTools, chem, ResourceConfig

__all__ = [
    "chem",
    "ChemTools", 
    "ResourceConfig",
    "scdb_matcher",  # Legacy import path
]
```

---

## Usage After Integration

### Before (Broken/Inconsistent)

```python
# Import path mismatch - BROKEN
from chemtools.scdb_matcher import load_db, match  # Module doesn't exist!

# Manual database loading
db = load_db("data/conditionDB/cn_coupling_pd_db.json")
result = match(db, reaction)
```

### After (Clean & Consistent)

```python
from chemtools import chem

# Integrated API - matches ML/precedent/reagent patterns
db = chem.rules.load_database("cn_coupling_pd_db.json")
result = chem.rules.match(db, reaction)

# Or list available databases
print(chem.rules.list_databases())
# ['amide_formation_db.json', 'cn_coupling_cu_db.json', 'cn_coupling_ni.json', ...]
```

### Hybrid Recommendation Workflow

```python
from chemtools import chem

reaction = "c1ccccc1Br.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

# 1. Try rule-based first (fast, deterministic)
db = chem.rules.load_database("cn_coupling_pd_db.json")
rule_result = chem.rules.match(db, reaction)

if rule_result.match_type == "scheme" and rule_result.priority == 0:
    # High-confidence scheme match
    print(f"Rule-based match: {rule_result.entry_name}")
    conditions = rule_result.conditions
else:
    # Fall back to ML-based recommendations
    print("Using ML-based recommendations")
    ml_result = chem.recommend.hybrid_recommend(reaction)
    conditions = ml_result.conditions

# 3. Enrich with reagent details
enriched = chem.reagent.enrich_conditions(conditions)
```

---

## Migration for Existing Code

### App Files to Update

1. **`app/main.py`** (line 41):
```python
# OLD (broken)
from chemtools.scdb_matcher import load_db as scdb_load_db, match as scdb_match

# NEW (working)
from chemtools import chem
db = chem.rules.load_database(path)
result = chem.rules.match(db, reaction)
```

2. **`app/ui_simple.py`** (line 51):
```python
# OLD
from chemtools.scdb_matcher import load_db as scdb_load_db, match as scdb_match
scdb_db = scdb_load_db("data/conditionDB/buchwald_scheme_db.json")

# NEW
from chemtools import chem
scdb_db = chem.rules.load_database("cn_coupling_pd_db.json")
```

3. **`app/ui_gradio.py`** (line 59):
Same pattern as `ui_simple.py`

### Backward Compatibility

For legacy code that uses direct imports:

```python
# Legacy import still works (via alias)
from chemtools.scdb_matcher import load_db, match

# Or direct import
from chemtools.rule_scdb_matcher import load_db, match
```

---

## Benefits of Integration

### 1. **Architectural Consistency**
All recommendation engines accessible via unified API:
- ✅ `chem.recommend.*` - ML-based
- ✅ `chem.precedent.*` - Precedent search
- ✅ `chem.rules.*` - Rule-based (NEW!)

### 2. **Resource Management**
- Database caching handled by ChemTools
- Memory-efficient (load once, reuse)
- Thread-safe access
- Clear cache on demand

### 3. **Developer Experience**
- Single import: `from chemtools import chem`
- Discoverable API (all features in one place)
- Consistent patterns across subsystems
- Better IDE autocomplete

### 4. **Hybrid Workflows**
Easy to combine multiple recommendation strategies:
```python
# Try rule-based → fall back to ML → enrich with reagent data
rule_result = chem.rules.match(db, reaction)
ml_result = chem.recommend.hybrid_recommend(reaction)
enriched = chem.reagent.enrich_conditions(best_conditions)
```

### 5. **Path Normalization**
Automatic handling of:
- Relative vs. absolute paths
- `data/conditionDB/` prefix
- Path resolution and validation

---

## Implementation Checklist

- [ ] **Phase 1**: Fix import path mismatch
  - [ ] Add alias in `chemtools/__init__.py`: `scdb_matcher = rule_scdb_matcher`
  - [ ] Test existing imports still work

- [ ] **Phase 2**: Create `RuleNamespace` class
  - [ ] Add class to `chemtools/context.py`
  - [ ] Implement `load_database()` with caching
  - [ ] Implement `match()` with proper forwarding
  - [ ] Add `list_databases()` helper
  - [ ] Add `clear_cache()` for memory management

- [ ] **Phase 3**: Integrate into ChemTools
  - [ ] Add `self.rules = RuleNamespace()` to `ChemTools.__init__()`
  - [ ] Test `chem.rules.*` API works
  - [ ] Verify caching behavior

- [ ] **Phase 4**: Update application code
  - [ ] Fix `app/main.py` to use `chem.rules.*`
  - [ ] Fix `app/ui_simple.py` to use `chem.rules.*`
  - [ ] Fix `app/ui_gradio.py` to use `chem.rules.*`

- [ ] **Phase 5**: Documentation & Testing
  - [ ] Update `CHEMTOOLS_QUICKSTART.md` with rule-based examples
  - [ ] Create `test_rules_integration.py`
  - [ ] Add hybrid workflow examples
  - [ ] Update API documentation

---

## Summary

**Question**: Why should rule-based be standalone?

**Answer**: **It shouldn't!** The current standalone status is an **oversight**, not a design decision. 

**Evidence**:
- Import path mismatch (broken imports in app/*.py)
- Inconsistent with ML/precedent/reagent integration patterns
- Missing from ChemTools v2.0 API
- No resource management or caching

**Solution**: Integrate it properly following the same pattern as other subsystems, providing a unified `chem.rules.*` API that's consistent with `chem.recommend.*`, `chem.precedent.*`, and `chem.reagent.*`.

**Priority**: HIGH - This fixes broken imports and completes the ChemTools v2.0 architecture.
