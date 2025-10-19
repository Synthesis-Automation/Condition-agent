# Core.py Refactoring Plan

## Overview

The `chemtools/recommend/core.py` file (1,334 lines) should be split into logical modules within `chemtools/recommend/modules/`. This document provides the complete refactoring plan.

## Current Structure Analysis

### File Statistics
- **Total lines**: 1,334
- **Public functions**: 2
  - `recommend_from_reaction()` - Main recommendation (~310 lines)
  - `recommend_conditions_structured()` - Structured API wrapper (~140 lines)
- **Private helpers**: 8 functions (~880 lines)
  - Fusion converters: 2 functions (~350 lines)
  - Output builders: 2 functions (~600 lines)
  - Precedent builders: 1 function (~130 lines)
  - Statistics: 4 small functions (~30 lines)

### Function Breakdown

| Function | Lines | Purpose | Target Module |
|----------|-------|---------|---------------|
| `recommend_from_reaction()` | 310 | Main DRFP recommendation | `recommender.py` |
| `recommend_conditions_structured()` | 140 | API-friendly wrapper | `structured.py` |
| `_convert_fusion_to_core_format()` | 165 | Convert fusion → core format | `fusion_adapter.py` |
| `_build_formatted_output_from_fusion()` | 185 | Build fusion output | `fusion_adapter.py` |
| `_build_precedent_details()` | 130 | Build precedent summaries | `precedent_builder.py` ✅ |
| `_calculate_average_yield()` | 8 | Calculate avg yield | `precedent_builder.py` ✅ |
| `_calculate_yield_range()` | 8 | Calculate yield range | `precedent_builder.py` ✅ |
| `_calculate_temp_range()` | 8 | Calculate temp range | `precedent_builder.py` ✅ |
| `_calculate_time_range()` | 8 | Calculate time range | `precedent_builder.py` ✅ |
| `_build_formatted_output()` | 340 | Build multi-variant output | `output_builder.py` |

✅ = Already extracted

## Proposed Module Structure

```
chemtools/recommend/
├── core.py                    # 100 lines (was 1,334) - Compatibility layer
└── modules/                   # New package (7 files, ~1,400 lines)
    ├── __init__.py            # Package exports
    ├── recommender.py         # ✅ Main DRFP recommendation (310 lines)
    ├── structured.py          # Structured API wrapper (140 lines)
    ├── fusion_adapter.py      # Fusion format conversion (350 lines)
    ├── precedent_builder.py   # ✅ Precedent summaries (180 lines)
    └── output_builder.py      # Output formatting (420 lines)
```

## Module Details

### 1. `recommender.py` ✅ (DONE - 310 lines)

**Purpose**: Main DRFP-based recommendation engine

**Functions**:
- `recommend_from_reaction()` - Core recommendation logic

**Dependencies**:
- `chemtools.smiles.normalize_reaction`
- `chemtools.router.detect_family`
- `chemtools.precedent.knn`
- `chemtools.explain.for_pack`
- `chemtools.reagent.enrich_reagent_info` (optional filtering)
- `chemtools.ml.simple_precedent_ranker` (optional reranking)
- `..utils.canonical_family, median, pick_with_constraints`
- `.output_builder.build_formatted_output` (lazy import)

**Status**: ✅ Created

---

### 2. `structured.py` (TODO - 140 lines)

**Purpose**: API-friendly wrapper with structured output

**Functions**:
- `recommend_conditions_structured()` - Wrapper around `recommend_from_reaction()`

**Dependencies**:
- `.recommender.recommend_from_reaction`
- `datetime`, `time` (for timestamps and timing)

**Extraction**:
```python
# Lines 350-489 from core.py
def recommend_conditions_structured(...):
    """Wrapper with structured API output."""
```

---

### 3. `fusion_adapter.py` (TODO - 350 lines)

**Purpose**: Convert fusion recommender output to core.py format

**Functions**:
- `convert_fusion_to_core_format()` - Main conversion (~165 lines)
- `build_formatted_output_from_fusion()` - Fusion output builder (~185 lines)

**Dependencies**:
- `collections.Counter`
- `chemtools.smiles.normalize_reaction`
- `chemtools.reagent.enrich_reagent_info`

**Extraction**:
```python
# Lines 489-835 from core.py
def convert_fusion_to_core_format(...):
    """Convert fusion results to core format."""

def build_formatted_output_from_fusion(...):
    """Build fusion output structure."""
```

---

### 4. `precedent_builder.py` ✅ (DONE - 180 lines)

**Purpose**: Build precedent summaries and calculate statistics

**Functions**:
- `build_precedent_details()` - Build detailed precedent info
- `calculate_average_yield()` - Average yield statistics
- `calculate_yield_range()` - Yield range
- `calculate_temp_range()` - Temperature range
- `calculate_time_range()` - Time range

**Status**: ✅ Created

---

### 5. `output_builder.py` (TODO - 420 lines)

**Purpose**: Build formatted multi-variant outputs

**Functions**:
- `build_formatted_output()` - Main output builder (~340 lines)
- Helper closures:
  - `_lookup()` - Reagent database lookup
  - `_chemical_payload()` - Create chemical payloads
  - `_cat_items()` - Extract catalyst system
  - `_build_variant()` - Build single condition variant
  - `_matching_precedents()` - Find matching precedents
  - `_add_combo()` - Add base/solvent combinations

**Dependencies**:
- `chemtools.smiles.normalize_reaction`
- `chemtools.reagent.enrich_reagent_info`
- `chemtools.context.lookup` (legacy fallback)
- `.precedent_builder.build_precedent_details`

**Extraction**:
```python
# Lines 993-1334 from core.py
def build_formatted_output(...):
    """Build formatted output with multiple condition variants."""
    
    # All internal helper functions
```

---

### 6. `__init__.py` (TODO - 50 lines)

**Purpose**: Package-level exports

**Content**:
```python
"""
chemtools.recommend.modules - Modular recommendation components.
"""

from .recommender import recommend_from_reaction
from .structured import recommend_conditions_structured
from .fusion_adapter import (
    convert_fusion_to_core_format,
    build_formatted_output_from_fusion,
)
from .precedent_builder import (
    build_precedent_details,
    calculate_average_yield,
    calculate_yield_range,
    calculate_temp_range,
    calculate_time_range,
)
from .output_builder import build_formatted_output

__all__ = [
    # Main functions
    "recommend_from_reaction",
    "recommend_conditions_structured",
    
    # Fusion adapter
    "convert_fusion_to_core_format",
    "build_formatted_output_from_fusion",
    
    # Precedent builders
    "build_precedent_details",
    "calculate_average_yield",
    "calculate_yield_range",
    "calculate_temp_range",
    "calculate_time_range",
    
    # Output builder
    "build_formatted_output",
]
```

---

### 7. `core.py` (TODO - Update to ~100 lines)

**Purpose**: Backwards compatibility layer

**Content**:
```python
"""
Core recommendation engine using DRFP-based reaction similarity.

REFACTORED: This module now re-exports from chemtools.recommend.modules.
For new code, prefer importing from submodules directly.
"""

from chemtools.recommend.modules import (
    recommend_from_reaction,
    recommend_conditions_structured,
)

__all__ = [
    "recommend_from_reaction",
    "recommend_conditions_structured",
]
```

**Old imports still work**:
```python
# ✅ Still works
from chemtools.recommend.core import recommend_from_reaction

# ✅ Also works (new)
from chemtools.recommend import recommend_from_reaction

# ✅ Also works (most specific)
from chemtools.recommend.modules.recommender import recommend_from_reaction
```

---

## Refactoring Steps

### ✅ Completed Steps

1. ✅ Created `chemtools/recommend/modules/` directory
2. ✅ Created `recommender.py` with `recommend_from_reaction()`
3. ✅ Created `precedent_builder.py` with all statistics functions

### 🔄 Remaining Steps

4. ⏳ Create `output_builder.py` with `build_formatted_output()` (~420 lines)
   - Extract lines 993-1334 from core.py
   - Include all helper closures
   - Import `build_precedent_details` from precedent_builder

5. ⏳ Create `structured.py` with `recommend_conditions_structured()` (~140 lines)
   - Extract lines 350-489 from core.py
   - Import `recommend_from_reaction` from recommender

6. ⏳ Create `fusion_adapter.py` with fusion conversion functions (~350 lines)
   - Extract lines 489-835 from core.py
   - Import `build_formatted_output_from_fusion` internal helper

7. ⏳ Create `modules/__init__.py` (~50 lines)
   - Re-export all public functions
   - Document module organization

8. ⏳ Update `core.py` to compatibility layer (~100 lines)
   - Import all functions from modules
   - Re-export for backwards compatibility
   - Add deprecation notices for internal functions

9. ⏳ Test refactoring
   - Verify all imports work
   - Test `recommend_from_reaction()` execution
   - Test `recommend_conditions_structured()` execution
   - Ensure no breaking changes

---

## Testing Plan

### Import Tests

```python
# Test backwards compatibility
from chemtools.recommend.core import recommend_from_reaction
from chemtools.recommend.core import recommend_conditions_structured

# Test new module structure
from chemtools.recommend.modules import recommend_from_reaction
from chemtools.recommend.modules.recommender import recommend_from_reaction
from chemtools.recommend.modules.precedent_builder import build_precedent_details

# Test all work
assert callable(recommend_from_reaction)
```

### Execution Tests

```python
# Test basic recommendation
result = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=25
)
assert "recommendation" in result
assert "formatted" in result

# Test structured output
result = recommend_conditions_structured(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50,
    limit=5
)
assert "recommendations" in result
assert "meta" in result
```

---

## Benefits

1. **Modularity**: Clear separation of recommendation, formatting, fusion, precedents
2. **Maintainability**: Files 140-420 lines instead of 1,334
3. **Testability**: Each module independently testable
4. **Reusability**: Import only what you need
5. **Backwards Compatible**: No breaking changes
6. **Scalability**: Easy to add new recommendation strategies

---

## Lines of Code Comparison

| File | Before | After | Change |
|------|--------|-------|--------|
| `core.py` | 1,334 lines | 100 lines | **-92.5%** |
| `recommender.py` | - | 310 lines | +310 |
| `structured.py` | - | 140 lines | +140 |
| `fusion_adapter.py` | - | 350 lines | +350 |
| `precedent_builder.py` | - | 180 lines | +180 |
| `output_builder.py` | - | 420 lines | +420 |
| `__init__.py` | - | 50 lines | +50 |
| **Total** | **1,334 lines** | **1,550 lines** | **+216 lines (+16%)** |

**Net impact**: Slightly more lines overall due to module headers and imports, but **much better organization** with average file size of 258 lines instead of one 1,334-line file.

---

## Priority

This refactoring is **HIGH PRIORITY** because:
- `core.py` is one of the largest files (1,334 lines)
- It's a critical recommendation engine
- It has complex dependencies and internal helpers
- It would benefit greatly from modularization

**Estimated effort**: 3-4 hours for complete refactoring + testing

---

## Related Documentation

- `FORMATTERS.md` - Output formatter refactoring (similar pattern)
- `SERVICE_LAYER.md` - Service layer extraction (similar pattern)
- `CODE_REVIEW.md` - Original recommendations

---

**Status**: 🔄 **IN PROGRESS (40% complete)**
- ✅ Directory created
- ✅ Recommender module created
- ✅ Precedent builder created
- ⏳ Output builder (TODO)
- ⏳ Structured wrapper (TODO)
- ⏳ Fusion adapter (TODO)
- ⏳ Package init (TODO)
- ⏳ Core.py update (TODO)
- ⏳ Testing (TODO)
