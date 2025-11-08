# SMARTS Cache Quick Reference

**TL;DR:** Always use `compile_smarts()` instead of `Chem.MolFromSmarts()` - it provides automatic caching across the entire codebase.

## Centralized SMARTS Source

- All canonical SMARTS (functional groups, coupling handles, nucleophiles, etc.) now live in `chemtools/featurizers/calculable_features.json`. Update that file once and every consumer (router, detection, featurizers, substrate classifier) picks up the change automatically.
- Use `chemtools.util.functional_groups.detect_all()` / `detect_any()` to reuse cached detections sourced from the JSON spec instead of defining ad-hoc SMARTS dictionaries.
- Need the raw SMARTS for mapping atoms? Import `FUNCTIONAL_GROUP_SMARTS` from `chemtools.util.functional_groups` for a read-only `{name: tuple(smarts)}` view that is always in sync with the spec.

## Quick Start

### Basic Usage

```python
from chemtools.util.smarts_cache import compile_smarts

# Compile a pattern (automatically cached)
pattern = compile_smarts("[CX4][Cl,Br,I]", validate=False)

if pattern and mol.HasSubstructMatch(pattern):
    print("Alkyl halide found!")
```

### Module-Level Pattern Definitions

**✅ Recommended Pattern:**

```python
# Define patterns as strings at module level
_PATTERNS = {
    "aryl_halide": "[$(c[Cl,Br,I]),$(c-[Cl,Br,I])]",
    "alcohol": "[OX2H]",
    "amine": "[NX3;H1,H2]",
}

def detect_features(mol):
    # Compile lazily when needed
    pattern = compile_smarts(_PATTERNS["aryl_halide"], validate=False)
    if pattern:
        return mol.HasSubstructMatch(pattern)
    return False
```

**❌ Anti-Pattern (Don't Do This):**

```python
from rdkit import Chem

# DON'T: Eager compilation at module load
_PATTERN = Chem.MolFromSmarts("[CX4][Cl,Br,I]")  # ❌ Not cached globally!

def detect_features(mol):
    # DON'T: Inline compilation without caching
    pattern = Chem.MolFromSmarts("[OX2H]")  # ❌ Recompiled every call!
    return mol.HasSubstructMatch(pattern)
```

## API Reference

### `compile_smarts(smarts, validate=False)`

Compile a SMARTS pattern with automatic LRU caching.

**Parameters:**
- `smarts` (str): SMARTS pattern string
- `validate` (bool): If True, raise ValueError on invalid pattern

**Returns:**
- Compiled RDKit Mol object, or None if invalid/RDKit unavailable

**Examples:**

```python
# Production: Skip validation for speed
pattern = compile_smarts("[CX4][Br,I]", validate=False)

# Development: Validate new patterns
try:
    pattern = compile_smarts(user_input, validate=True)
except ValueError as e:
    print(f"Invalid SMARTS: {e}")
```

### `compile_smarts_batch(patterns, validate=False, skip_invalid=True)`

Compile multiple patterns efficiently.

**Parameters:**
- `patterns` (dict): Dict mapping names to SMARTS strings
- `validate` (bool): Raise on first invalid pattern
- `skip_invalid` (bool): Return None for invalid patterns instead of raising

**Returns:**
- Dict mapping names to compiled patterns (or None for invalid)

**Examples:**

```python
from chemtools.util.smarts_cache import compile_smarts_batch

patterns = {
    "halide": "[Cl,Br,I]",
    "alcohol": "[OX2H]",
    "amine": "[NX3;H1,H2]",
}

# Compile all at once (still benefits from cache)
compiled = compile_smarts_batch(patterns, skip_invalid=True)

for name, pattern in compiled.items():
    if pattern and mol.HasSubstructMatch(pattern):
        print(f"Found {name}")
```

### `get_cache_info()`

Get cache statistics (hits, misses, size, hit rate).

**Returns:**
- Dict with keys: `hits`, `misses`, `size`, `maxsize`, `hit_rate`, `total_calls`

**Example:**

```python
from chemtools.util.smarts_cache import get_cache_info

info = get_cache_info()
print(f"Cache hit rate: {info['hit_rate']:.1%}")
print(f"Patterns cached: {info['size']}/{info['maxsize']}")
```

### `clear_smarts_cache()`

Clear all cached patterns (useful for testing).

**Example:**

```python
from chemtools.util.smarts_cache import clear_smarts_cache

# Clear cache before test
clear_smarts_cache()
```

## Common Patterns

### Pattern 1: Simple Feature Detection

```python
from chemtools.util.smarts_cache import compile_smarts

def has_alkyl_halide(mol):
    """Check if molecule contains alkyl halide."""
    pattern = compile_smarts("[CX4][Cl,Br,I]", validate=False)
    return pattern and mol.HasSubstructMatch(pattern)
```

### Pattern 2: Module-Level Pattern Library

```python
from chemtools.util.smarts_cache import compile_smarts

# Define patterns as constants
_FUNCTIONAL_GROUPS = {
    "alcohol": "[OX2H]",
    "aldehyde": "[CX3H](=O)",
    "ketone": "[CX3](=O)[C]",
    "ester": "[CX3](=O)[OX2][C,H]",
}

def detect_functional_groups(mol):
    """Detect all functional groups in molecule."""
    results = {}
    for name, smarts in _FUNCTIONAL_GROUPS.items():
        pattern = compile_smarts(smarts, validate=False)
        if pattern:
            results[name] = mol.HasSubstructMatch(pattern)
        else:
            results[name] = False
    return results
```

### Pattern 3: Dynamic Pattern Generation

```python
from chemtools.util.smarts_cache import compile_smarts

def has_halogen(mol, halogen="Br"):
    """Check for specific halogen attached to carbon."""
    smarts = f"[C][{halogen}]"
    pattern = compile_smarts(smarts, validate=False)
    return pattern and mol.HasSubstructMatch(pattern)
```

### Pattern 4: Counting Matches

```python
from chemtools.util.smarts_cache import compile_smarts

def count_alcohols(mol):
    """Count number of alcohol groups."""
    pattern = compile_smarts("[OX2H]", validate=False)
    if not pattern:
        return 0
    matches = mol.GetSubstructMatches(pattern)
    return len(matches)
```

## Performance Tips

1. **Use `validate=False` in production**: Skip validation for 2-3× faster compilation
2. **Prefer module-level pattern defs**: Easier to maintain and validate in bulk
3. **Use batch compilation for startup**: If loading many patterns at once
4. **Monitor cache statistics**: Use `get_cache_info()` to track hit rates
5. **Clear cache in tests**: Use `clear_smarts_cache()` for isolated test runs

## Cache Configuration

- **Default max size**: 1024 patterns (configurable via `@lru_cache(maxsize=...)`)
- **Eviction policy**: LRU (Least Recently Used)
- **Thread safety**: Yes (`functools.lru_cache` is thread-safe)
- **Scope**: Global (shared across all modules)

## Migration from Direct RDKit Usage

### Before (No Caching):

```python
from rdkit import Chem

def detect_features(mol):
    # Compiled every time function is called
    pattern = Chem.MolFromSmarts("[CX4][Br,I]")
    if pattern:
        return mol.HasSubstructMatch(pattern)
    return False
```

### After (Centralized Cache):

```python
from chemtools.util.smarts_cache import compile_smarts

# Pattern string as module constant
_ALKYL_HALIDE_SMARTS = "[CX4][Br,I]"

def detect_features(mol):
    # Compiled once, cached globally
    pattern = compile_smarts(_ALKYL_HALIDE_SMARTS, validate=False)
    if pattern:
        return mol.HasSubstructMatch(pattern)
    return False
```

**Benefits:**
- First call: ~1-5ms compilation time
- Subsequent calls: ~0.001ms (cached lookup)
- **1000-5000× speedup** for repeated patterns

## Troubleshooting

### Pattern Returns None

**Symptom:** `compile_smarts()` returns `None`

**Causes:**
1. Invalid SMARTS syntax
2. RDKit not available

**Solution:**
```python
# Use validate=True to get error message
try:
    pattern = compile_smarts(smarts, validate=True)
except ValueError as e:
    print(f"Invalid pattern: {e}")
```

### Low Cache Hit Rate

**Symptom:** `get_cache_info()['hit_rate'] < 0.1`

**Causes:**
1. Too many unique patterns
2. Patterns with variables/f-strings

**Solutions:**
1. Consolidate similar patterns
2. Use pattern libraries instead of dynamic generation
3. Consider increasing cache size (edit `smarts_cache.py`)

### Module Import Slow

**Symptom:** Module takes long to import

**Cause:** Eager pattern compilation at module level

**Solution:** Use lazy compilation (see "Module-Level Pattern Definitions" above)

## See Also

- `SMARTS_COMPILATION_CONSOLIDATION.md` - Original design document
- `SMARTS_PHASE1_COMPLETE.md` - Phase 1 implementation details
- `SMARTS_PHASE2_COMPLETE.md` - Phase 2 implementation details  
- `scripts/validate_smarts_patterns.py` - Pattern validation tool
- `tests/test_smarts_cache.py` - Comprehensive test suite

---

**Last Updated:** November 7, 2025  
**Module:** `chemtools/util/smarts_cache.py`  
**Global Cache Size:** 1024 patterns
