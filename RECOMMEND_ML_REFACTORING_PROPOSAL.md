# Recommendation & ML Module Refactoring Proposal

**Date**: October 7, 2025  
**Status**: PROPOSAL - Awaiting approval  
**Scope**: `recommend.py`, `recommend_ml.py`, `/ml` folder reorganization

---

## 📊 Current State Analysis

### Module Sizes & Complexity

| Module | Lines | Responsibilities | Issues |
|--------|-------|-----------------|--------|
| `recommend.py` | 1,454 | • Family detection<br>• Substrate analysis<br>• Amide-specific logic (538 lines)<br>• C-N coupling features<br>• Precedent search<br>• Plate design<br>• Output formatting | ⚠️ **Too large**<br>⚠️ Multiple responsibilities<br>⚠️ Amide logic embedded |
| `recommend_ml.py` | 226 | • ML/k-NN hybrid<br>• Yield prediction integration<br>• Re-ranking logic | ✅ Reasonably sized<br>⚠️ Duplicates some logic |
| `ml/drfp_predictor.py` | 370 | • DRFP encoding<br>• LightGBM training<br>• Yield prediction | ✅ Well-focused |
| `ml/evaluation.py` | 250 | • Metrics computation<br>• Plotting<br>• Model comparison | ✅ Well-focused |

### Current Architecture Issues

#### 1. **recommend.py is a Monolith** (1,454 lines)
```python
# All in one file:
- Amide-specific analysis (538 lines) ← Should be separate
- C-N coupling features           ← Should be separate  
- Suzuki/Sonogashira features      ← Should be separate
- General recommendation logic     ← Core
- Plate design                     ← Could be separate
- Role featurization               ← Utility
```

#### 2. **Duplicate Logic Between recommend.py and recommend_ml.py**
- Both parse reaction SMILES
- Both call precedent search
- Both format outputs
- recommend_ml.py wraps recommend.py but adds complexity

#### 3. **ML Module Lacks Clear Integration Point**
- `ml/drfp_predictor.py` is standalone
- `recommend_ml.py` bridges the gap awkwardly
- No clear separation between:
  - k-NN recommendations (rule-based)
  - ML recommendations (learned)
  - Hybrid recommendations

#### 4. **Family-Specific Logic is Embedded**
```python
# In recommend.py - hard to maintain/extend:
def _amide_rule_feature_builder():  # 538 lines!
    # Acid classification
    # Amine classification
    # Functional groups
    # Water tolerance
    # ...
```

---

## 🎯 Proposed Refactoring

### Strategy: Split by Responsibility

Following the successful `precedent/` refactoring pattern, we'll create focused modules:

```
chemtools/
├── recommend/                    # ✨ NEW: Package structure
│   ├── __init__.py               # Public API exports
│   ├── core.py                   # Main recommendation engine
│   ├── families/                 # Family-specific logic
│   │   ├── __init__.py
│   │   ├── amide.py              # Amide formation (538 lines)
│   │   ├── cn_coupling.py        # C-N coupling (Pd/Cu/Ni)
│   │   ├── cc_coupling.py        # Suzuki, Sonogashira
│   │   └── base.py               # Base family interface
│   ├── substrate_analysis.py     # Substrate classification
│   ├── plate_design.py           # Plate layout generation
│   └── utils.py                  # Helper functions
│
├── ml/                           # ✨ ENHANCED: Clear ML package
│   ├── __init__.py               # Unified ML API
│   ├── predictors/               # ✨ NEW: Predictor modules
│   │   ├── __init__.py
│   │   ├── drfp_predictor.py     # DRFP baseline (moved)
│   │   └── neural_predictor.py   # Future: Neural models
│   ├── hybrid.py                 # ✨ NEW: ML+k-NN integration
│   └── evaluation.py             # (unchanged)
│
└── recommend_ml.py               # ⚠️ DEPRECATED: Remove after migration
```

---

## 📦 Detailed Module Breakdown

### 1. `recommend/` Package

#### **`recommend/__init__.py`** (~30 lines)
**Purpose**: Public API exports for backward compatibility

```python
"""
Reaction condition recommendation system.

Provides k-NN precedent-based and ML-augmented recommendations.
"""

# Main recommendation functions (from core.py)
from .core import (
    recommend_from_reaction,
    recommend_conditions_structured,
)

# Plate design (from plate_design.py)
from .plate_design import design_plate_from_reaction

# Family-specific builders (from families/)
from .families import get_family_analyzer

__all__ = [
    'recommend_from_reaction',
    'recommend_conditions_structured', 
    'design_plate_from_reaction',
    'get_family_analyzer',
]
```

#### **`recommend/core.py`** (~300-400 lines)
**Purpose**: Main recommendation engine (family-agnostic)

**Responsibilities**:
- Reaction normalization and parsing
- Family detection integration
- Precedent search coordination
- Output formatting coordination
- Constraint application

**Key Functions**:
```python
def recommend_from_reaction(
    reaction: str,
    k: int = 25,
    constraints: Dict = None,
    **kwargs
) -> Dict[str, Any]:
    """Main recommendation entry point."""
    
def _build_recommendation(
    family: str,
    precedent_pack: Dict,
    features: Dict,
    constraints: Dict,
) -> Dict[str, Any]:
    """Build final recommendation from precedents."""

def _format_output(
    recommendation: Dict,
    precedent_pack: Dict,
    **kwargs
) -> Dict[str, Any]:
    """Format final output structure."""
```

**No family-specific logic** - delegates to `families/` modules.

---

#### **`recommend/families/base.py`** (~100 lines)
**Purpose**: Base interface for family-specific analyzers

```python
from abc import ABC, abstractmethod
from typing import Dict, Any, List

class FamilyAnalyzer(ABC):
    """Base class for reaction family-specific analysis."""
    
    @abstractmethod
    def analyze_substrates(
        self,
        reactants: List[str],
        features: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Analyze substrates for this family."""
        pass
    
    @abstractmethod
    def build_features(
        self,
        role_pack: Dict[str, Any],
        features: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Build family-specific features for precedent search."""
        pass
    
    def get_default_conditions(self) -> Dict[str, Any]:
        """Get default conditions for this family."""
        return {}
```

---

#### **`recommend/families/amide.py`** (~550 lines)
**Purpose**: Amide formation specific logic (extracted from recommend.py)

**Extracted from**: Lines 297-836 of current `recommend.py`

**Key Functions** (all prefixed with `_` → now in dedicated module):
```python
# Substrate analysis (currently 297-507)
def analyze_carboxylic_acid(smiles: str) -> Dict[str, Any]:
    """Analyze carboxylic acid substrate."""
    # Aromatic vs aliphatic
    # α-amino, β-hydroxy detection
    # Protected amino acids
    # Steric hindrance
    
def analyze_amine(smiles: str, features: Dict) -> Dict[str, Any]:
    """Analyze amine substrate."""
    # Primary/secondary/aniline
    # Benzylic, α-amino detection
    # Nucleophilicity assessment
    # Steric hindrance

# Classification (currently 507-596)
def classify_acid(acid_info: Dict) -> Dict[str, Any]:
    """Classify acid type for rule matching."""
    
def classify_amine(amine_info: Dict, features: Dict) -> Dict[str, Any]:
    """Classify amine type for rule matching."""

# Feature building (currently 729-794)
def build_amide_features(
    role_pack: Dict,
    features: Dict,
) -> Dict[str, Any]:
    """Build amide-specific features for precedent search."""
    # Combines all above logic
    # Water tolerance inference
    # Green chemistry flags

# Analyzer class
class AmideAnalyzer(FamilyAnalyzer):
    """Amide formation analyzer."""
    
    def analyze_substrates(self, reactants, features):
        acid = self._find_acid(reactants)
        amine = self._find_amine(reactants)
        return {
            'acid': analyze_carboxylic_acid(acid),
            'amine': analyze_amine(amine, features),
        }
    
    def build_features(self, role_pack, features):
        return build_amide_features(role_pack, features)
```

---

#### **`recommend/families/cn_coupling.py`** (~200 lines)
**Purpose**: C-N coupling (Pd/Cu/Ni) specific logic

**Extracted from**: Default feature builder + C-N specific code

```python
class CNCouplingAnalyzer(FamilyAnalyzer):
    """C-N coupling analyzer (Pd, Cu, Ni variants)."""
    
    def __init__(self, metal: str = 'Pd'):
        self.metal = metal  # 'Pd', 'Cu', or 'Ni'
    
    def analyze_substrates(self, reactants, features):
        elec = self._find_electrophile(reactants)
        nuc = self._find_nucleophile(reactants)
        return {
            'electrophile': self._analyze_aryl_halide(elec),
            'nucleophile': self._analyze_amine(nuc),
            'metal': self.metal,
        }
    
    def build_features(self, role_pack, features):
        return {
            'LG': self._get_leaving_group(features),
            'nuc_class': self._get_nucleophile_class(features),
            'ortho': features.get('ortho_subst', 0),
            'para_EWG': features.get('para_EWG', False),
            'metal': self.metal,
        }
```

---

#### **`recommend/families/cc_coupling.py`** (~150 lines)
**Purpose**: C-C coupling (Suzuki, Sonogashira) specific logic

```python
class SuzukiAnalyzer(FamilyAnalyzer):
    """Suzuki coupling analyzer."""
    
class SonogashiraAnalyzer(FamilyAnalyzer):
    """Sonogashira coupling analyzer."""
```

---

#### **`recommend/families/__init__.py`** (~50 lines)
**Purpose**: Family registry and factory

```python
from .base import FamilyAnalyzer
from .amide import AmideAnalyzer
from .cn_coupling import CNCouplingAnalyzer
from .cc_coupling import SuzukiAnalyzer, SonogashiraAnalyzer

FAMILY_REGISTRY = {
    'Amide_Formation': AmideAnalyzer,
    'C_N_Coupling_Pd': lambda: CNCouplingAnalyzer(metal='Pd'),
    'C_N_Coupling_Cu': lambda: CNCouplingAnalyzer(metal='Cu'),
    'C_N_Coupling_Ni': lambda: CNCouplingAnalyzer(metal='Ni'),
    'Suzuki_Coupling': SuzukiAnalyzer,
    'Sonogashira': SonogashiraAnalyzer,
}

def get_family_analyzer(family: str) -> FamilyAnalyzer | None:
    """Get analyzer for reaction family."""
    factory = FAMILY_REGISTRY.get(family)
    if factory:
        return factory() if callable(factory) else factory
    return None
```

---

#### **`recommend/substrate_analysis.py`** (~150 lines)
**Purpose**: Generic substrate analysis utilities

**Extracted from**: Lines 462-650 of current `recommend.py`

```python
def has_free_alcohol(smiles: str) -> bool:
    """Check for free alcohol groups."""
    
def has_phenol(smiles: str) -> bool:
    """Check for phenol groups."""

def detect_functional_groups(smiles: str) -> Dict[str, bool]:
    """Detect common functional groups."""
    
def assess_steric_hindrance(smiles: str) -> str:
    """Assess steric hindrance (low/medium/high)."""
```

---

#### **`recommend/plate_design.py`** (~200 lines)
**Purpose**: Experimental plate layout generation

**Extracted from**: Lines 1339-1454 of current `recommend.py`

```python
def design_plate_from_reaction(
    reaction: str,
    plate_size: int = 24,
    **kwargs
) -> Dict[str, Any]:
    """Generate experimental plate layout."""
    
def _well_ids(n: int) -> List[str]:
    """Generate well IDs (A1, A2, ...)."""
    
def _assign_conditions_to_wells(
    conditions: List[Dict],
    plate_size: int,
) -> List[Dict]:
    """Assign conditions to plate wells."""
```

---

#### **`recommend/utils.py`** (~100 lines)
**Purpose**: Helper functions

**Extracted from**: Small utility functions scattered in recommend.py

```python
def canonical_family(family: str | None) -> str:
    """Normalize family name."""
    
def pick_electrophile_nucleophile(reactants: List[str]) -> Tuple[str, str]:
    """Pick electrophile and nucleophile from reactants."""
    
def median(vals: List[float]) -> float | None:
    """Compute median value."""
    
def pick_with_constraints(
    candidates: List[str],
    constraints: Dict,
) -> Tuple[str | None, Dict]:
    """Apply constraints to candidate list."""
```

---

### 2. `ml/` Package Enhancement

#### **`ml/__init__.py`** (updated, ~40 lines)
**Purpose**: Unified ML API with clear exports

```python
"""
Machine learning models for reaction condition prediction.

Provides:
- Yield predictors (DRFP-based, neural)
- Hybrid ML + k-NN recommenders
- Model evaluation utilities
"""

# Predictors
from .predictors import DRFPYieldPredictor

# Hybrid recommendation
from .hybrid import (
    hybrid_recommend,
    recommend_with_yield_prediction,
)

# Evaluation
from .evaluation import (
    evaluate_yield_predictor,
    compute_metrics,
    compare_models,
)

__all__ = [
    # Predictors
    'DRFPYieldPredictor',
    
    # Hybrid
    'hybrid_recommend',
    'recommend_with_yield_prediction',
    
    # Evaluation
    'evaluate_yield_predictor',
    'compute_metrics',
    'compare_models',
]
```

---

#### **`ml/predictors/__init__.py`** (~20 lines)
**Purpose**: Predictor exports

```python
from .drfp_predictor import DRFPYieldPredictor

__all__ = ['DRFPYieldPredictor']
```

---

#### **`ml/predictors/drfp_predictor.py`** (moved, unchanged)
**Purpose**: DRFP-based yield predictor (already well-structured)

---

#### **`ml/hybrid.py`** (~250 lines)
**Purpose**: ML + k-NN hybrid recommendation (moved from recommend_ml.py)

**Changes from current `recommend_ml.py`**:
- Import from `chemtools.recommend` instead of `chemtools.recommend` (same)
- Cleaner integration with new recommend package structure
- Better error handling
- Unified output format

```python
"""
Hybrid ML + k-NN recommender.

Integrates yield prediction with precedent-based recommendations.
"""

from chemtools.recommend import recommend_from_reaction
from .predictors import DRFPYieldPredictor

def hybrid_recommend(
    reaction_smiles: str,
    model_path: str | None = None,
    ml_threshold: int = 50,
    k: int = 10,
    **kwargs,
) -> Dict[str, Any]:
    """Hybrid ML + k-NN recommendation."""
    # (same logic as current recommend_ml.py)
    
def recommend_with_yield_prediction(
    reaction_smiles: str,
    model: Any | None = None,
    model_path: str | None = None,
    top_n: int = 5,
    **knn_kwargs,
) -> List[Dict[str, Any]]:
    """Simplified API for top-N with yield predictions."""
    # (same logic as current recommend_ml.py)
```

---

## 🔄 Migration Plan

### Phase 1: Create New Structure (No Breaking Changes)

**Step 1**: Create `recommend/` package
- Create directory structure
- Move amide logic to `families/amide.py`
- Move C-N coupling logic to `families/cn_coupling.py`
- Extract plate design to `plate_design.py`
- Create core recommendation engine in `core.py`
- Create public API in `__init__.py`

**Step 2**: Enhance `ml/` package
- Create `ml/predictors/` subdirectory
- Move `drfp_predictor.py` to `ml/predictors/`
- Move `recommend_ml.py` logic to `ml/hybrid.py`
- Update `ml/__init__.py` with unified exports

**Step 3**: Update imports in `recommend/` package
- Ensure `recommend/__init__.py` exports same API as current `recommend.py`
- **100% backward compatible**: `from chemtools.recommend import recommend_from_reaction` still works

### Phase 2: Update Consumers

**Step 4**: Update `context.py`
```python
class RecommendNamespace:
    def hybrid_recommend(self, ...):
        # OLD: from . import recommend_ml
        # NEW: from .ml import hybrid
        from .ml.hybrid import hybrid_recommend
        return hybrid_recommend(...)
```

**Step 5**: Update tests
- Update import paths if needed
- Add new tests for family-specific analyzers

### Phase 3: Cleanup

**Step 6**: Deprecate old files
- Mark `recommend_ml.py` as deprecated (keep shim for 1 version)
- Add deprecation warnings

**Step 7**: Remove old files (next major version)
- Delete `recommend.py` (replaced by `recommend/`)
- Delete `recommend_ml.py` (replaced by `ml/hybrid.py`)

---

## ✅ Benefits

### 1. **Maintainability**
- ✅ Each module has single responsibility
- ✅ Family-specific logic isolated
- ✅ Easy to add new families (just add new analyzer)
- ✅ Amide logic (538 lines) no longer embedded

### 2. **Testability**
- ✅ Can test amide analyzer independently
- ✅ Can test C-N coupling independently
- ✅ Can mock family analyzers in core tests
- ✅ Smaller, focused test files

### 3. **Extensibility**
- ✅ Adding new family = new file in `families/`
- ✅ No need to modify core recommendation logic
- ✅ Plugin-style architecture

### 4. **Clarity**
- ✅ Clear separation: k-NN vs ML vs Hybrid
- ✅ `recommend/` = rule-based precedent search
- ✅ `ml/` = machine learning models
- ✅ `ml/hybrid.py` = integration layer

### 5. **Performance**
- ✅ Can lazy-load family analyzers
- ✅ No performance regression (same logic, better organized)

---

## 📊 Comparison: Before vs After

### Before (Current State)

```
chemtools/
├── recommend.py              (1,454 lines) ❌ Monolith
│   ├── Amide logic           (538 lines) - embedded
│   ├── C-N coupling logic    (~150 lines) - embedded
│   ├── Suzuki logic          (~100 lines) - embedded
│   ├── Core logic            (~300 lines) - mixed
│   ├── Plate design          (~200 lines) - mixed
│   └── Utilities             (~166 lines) - scattered
├── recommend_ml.py           (226 lines) ⚠️ Duplicate logic
└── ml/
    ├── drfp_predictor.py     (370 lines) ✅ Good
    └── evaluation.py         (250 lines) ✅ Good
```

**Issues**:
- ❌ 1,454-line monolith
- ❌ Family logic embedded (hard to extend)
- ❌ ML integration unclear
- ❌ Duplicate logic between recommend.py and recommend_ml.py

### After (Proposed State)

```
chemtools/
├── recommend/                         ✅ Package
│   ├── __init__.py               (30 lines)   - Public API
│   ├── core.py                   (350 lines)  - Main engine
│   ├── families/                              - Family-specific
│   │   ├── __init__.py           (50 lines)   - Registry
│   │   ├── base.py               (100 lines)  - Interface
│   │   ├── amide.py              (550 lines)  - Amide logic
│   │   ├── cn_coupling.py        (200 lines)  - C-N logic
│   │   └── cc_coupling.py        (150 lines)  - C-C logic
│   ├── substrate_analysis.py     (150 lines)  - Utilities
│   ├── plate_design.py           (200 lines)  - Plate layouts
│   └── utils.py                  (100 lines)  - Helpers
│
└── ml/                                        ✅ Enhanced
    ├── __init__.py               (40 lines)   - Unified API
    ├── predictors/                            - Predictor modules
    │   ├── __init__.py           (20 lines)
    │   └── drfp_predictor.py     (370 lines)  - Moved here
    ├── hybrid.py                 (250 lines)  - ML+k-NN integration
    └── evaluation.py             (250 lines)  - (unchanged)
```

**Improvements**:
- ✅ No file > 550 lines (vs 1,454-line monolith)
- ✅ Clear family separation (easy to extend)
- ✅ ML integration point clear (`ml/hybrid.py`)
- ✅ No duplicate logic
- ✅ 100% backward compatible

---

## 📋 File-by-File Changes Summary

| Current File | Action | New Location | Lines | Status |
|--------------|--------|--------------|-------|--------|
| `recommend.py` | ✂️ **SPLIT** | `recommend/` package | 1,454 → ~1,730 total (8 files) | Cleaner |
| ├─ Main logic | → | `recommend/core.py` | 350 | ✅ Focused |
| ├─ Amide logic | → | `recommend/families/amide.py` | 550 | ✅ Isolated |
| ├─ C-N logic | → | `recommend/families/cn_coupling.py` | 200 | ✅ Isolated |
| ├─ C-C logic | → | `recommend/families/cc_coupling.py` | 150 | ✅ Isolated |
| ├─ Substrate analysis | → | `recommend/substrate_analysis.py` | 150 | ✅ Reusable |
| ├─ Plate design | → | `recommend/plate_design.py` | 200 | ✅ Separate |
| └─ Utilities | → | `recommend/utils.py` | 100 | ✅ Clean |
| `recommend_ml.py` | 🔄 **MOVE** | `ml/hybrid.py` | 226 | ✅ Better location |
| `ml/drfp_predictor.py` | 🔄 **MOVE** | `ml/predictors/drfp_predictor.py` | 370 | ✅ Better organized |
| `ml/evaluation.py` | ✅ **KEEP** | (unchanged) | 250 | ✅ Already good |

**Net Effect**:
- Total lines: ~2,000 (slightly more due to interfaces/exports)
- Average module size: ~200 lines (vs 1,454-line monolith)
- Maintainability: ⭐⭐⭐⭐⭐ (vs ⭐⭐ current)

---

## 🎯 Recommendations

### ✅ **RECOMMENDED: Proceed with Refactoring**

**Reasoning**:
1. ✅ Proven pattern (successful `precedent/` refactoring)
2. ✅ Clear benefits (maintainability, testability, extensibility)
3. ✅ 100% backward compatible migration path
4. ✅ No performance regression
5. ✅ Addresses real pain points (1,454-line monolith)

### 📅 **Suggested Timeline**

**Week 1**: Create `recommend/` package structure
- Day 1-2: Create base structure and `core.py`
- Day 3-4: Extract `families/amide.py` (largest module)
- Day 5: Extract other family modules + tests

**Week 2**: Enhance `ml/` package
- Day 1: Create `ml/predictors/` and move files
- Day 2: Create `ml/hybrid.py` from `recommend_ml.py`
- Day 3: Update `ml/__init__.py` and tests

**Week 3**: Integration & testing
- Day 1-2: Update `context.py` integration
- Day 3: Update all tests
- Day 4: Integration testing
- Day 5: Documentation updates

### 🚀 **Next Steps**

1. **Review this proposal** - Get team feedback
2. **Create feature branch** - `feature/recommend-ml-refactoring`
3. **Implement Phase 1** - Create new structure (backward compatible)
4. **Test thoroughly** - Ensure no regressions
5. **Merge with deprecation warnings** - Keep old files temporarily
6. **Remove old files** - In next major version

---

## 📚 Additional Recommendations

### Future Enhancements (After Refactoring)

1. **Add unit tests for each family analyzer**
   ```python
   tests/
   ├── test_recommend_core.py
   ├── test_amide_analyzer.py
   ├── test_cn_coupling_analyzer.py
   └── test_hybrid_ml.py
   ```

2. **Add family analyzer configuration**
   ```python
   # Allow users to customize family analyzers
   from chemtools.recommend.families import AmideAnalyzer
   
   custom_analyzer = AmideAnalyzer(
       water_tolerance='strict',
       green_chemistry=True,
   )
   ```

3. **Add plugin system for custom families**
   ```python
   # Register custom family analyzer
   from chemtools.recommend.families import register_family
   
   @register_family('Custom_Reaction')
   class CustomAnalyzer(FamilyAnalyzer):
       ...
   ```

4. **Add ML model registry**
   ```python
   # ml/model_registry.py
   MODELS = {
       'drfp_v1': 'models/drfp_yield_v1.pkl',
       'neural_v1': 'models/neural_yield_v1.pt',
   }
   ```

---

## ❓ Questions for Discussion

1. **Naming**: Happy with `recommend/` and `ml/hybrid.py`? Alternatives?
2. **Backward compatibility**: Keep `recommend_ml.py` shim for how long? (suggest: 1 major version)
3. **Family interface**: Is `FamilyAnalyzer` base class sufficient? Need more methods?
4. **Testing strategy**: Add tests during refactoring or after?
5. **Documentation**: Update README.md during or after refactoring?

---

## 📝 Summary

**Current Issues**:
- ❌ `recommend.py` is 1,454-line monolith
- ❌ Family-specific logic embedded (hard to extend)
- ❌ ML integration unclear
- ❌ Duplicate logic between files

**Proposed Solution**:
- ✅ Split `recommend.py` into focused `recommend/` package
- ✅ Isolate family logic in `recommend/families/`
- ✅ Enhance `ml/` with clear integration point
- ✅ 100% backward compatible migration

**Benefits**:
- ✅ Average module size: ~200 lines (vs 1,454)
- ✅ Easy to extend (add new families)
- ✅ Better testing (isolated modules)
- ✅ Clearer architecture (k-NN vs ML vs Hybrid)

**Recommendation**: ✅ **PROCEED with refactoring** (proven pattern, clear benefits)
