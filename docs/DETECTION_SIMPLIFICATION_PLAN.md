# Detection System Simplification Plan

## 🎯 TL;DR - Where Does rxn-insight Go?

**Answer**: The entire `reaction_type_detector.py` file (~304 lines) gets consolidated into a **single method** inside the new `_DetectionEngine` class:

```python
# chemtools/detection.py (NEW unified file)

class _DetectionEngine:
    def _ml_detection(self) -> dict:
        """
        rxn-insight ML detection - ALL ML logic goes here.

        Consolidates from reaction_type_detector.py:
        - _call_insight() → calls rxn-insight library
        - _extract_fields() → parses ML response
        - _map_to_family() → maps to unified taxonomy (via resolve_to_taxonomy())
        - _refine_cn_family() → catalyst-aware refinements

        ALL outputs validated against chemtools/taxonomy/data/reaction_types.json
        """
        # ~150 lines of consolidated ML logic
        # Optional layer - gracefully degrades if unavailable
        # Uses TaxonomyRegistry for automatic alias resolution
```

**File Migration:**

```
BEFORE:
chemtools/router.py (548 lines) ────────┐
chemtools/reaction_type_detector.py (304)│──► chemtools/detection.py (400 lines)
                                        │    ├─ Uses TaxonomyRegistry
                                        │    ├─ _DetectionEngine._ml_detection()
                                        │    └─ resolve_to_taxonomy() for mapping
                                        └──► DELETE these files

AFTER:
chemtools/detection.py (400 lines total)
├─ One class, one public function: detect_reaction()
└─ Uses chemtools/taxonomy/ as single source of truth
   └─ 80+ canonical reaction types automatically available
```

**Key Insight: Taxonomy Alignment**

The NEW system uses your existing `chemtools/taxonomy/data/reaction_types.json` (80+ reaction types) as the **single source of truth**. All detection methods—SMARTS, ML, catalyst overrides—automatically map to canonical taxonomy IDs:

```python
# Taxonomy registry handles ALL aliases
resolve_to_taxonomy("Suzuki") → "suzuki_miyaura"
resolve_to_taxonomy("Buchwald-Hartwig") → "buchwald_hartwig_c_n"
resolve_to_taxonomy("C-N Coupling", catalysts={"Pd"}) → "buchwald_hartwig_c_n"

# No manual mapping tables needed!
# See TAXONOMY_ALIGNMENT_STRATEGY.md for details
```

---

## 🎯 Problem Statement

The current reaction detection system is **overly complex** with:

- **5 different detection methods** across multiple modules
- **3 separate entry points** that users must choose between
- **Confusing naming** (`detect_family` vs `detect_family_from_reaction` vs `detect_reaction_type`)
- **Duplicate logic** for catalyst detection, family mapping, and confidence scoring
- **Unclear API surface** - users don't know which method to use

### Current Complexity

```
User wants to detect reaction type
    ↓
Which function should I use?
    ├─ detect_family(reactants)?
    ├─ detect_family_from_reaction(reaction)?
    ├─ detect_reaction_type(reaction)?
    ├─ classify_reactants_with_context(reaction)?
    └─ What's the difference? 🤷
```

### Code Duplication

**Catalyst detection** - duplicated in:

- `router.py::_detect_agent_metals()`
- `reaction_type_detector.py::detect_reaction_type()` (indirectly via imports)
- `router.py::detect_family_from_reaction()` (calls \_detect_agent_metals)

**Family mapping** - duplicated in:

- `reaction_type_detector.py::_map_to_family()`
- `router.py::resolve_reaction_family()`
- `analysis/reactions.py::resolve_reaction_family()`

**Confidence calculation** - scattered across:

- `router.py::detect_family()` (hardcoded values)
- `reaction_type_detector.py::_refine_cn_family()` (override logic)
- `router.py::detect_family_from_reaction()` (merging logic)

---

## ✨ Proposed Solution: Single Unified API

### New Simple API

```python
from chemtools import detect_reaction

# ONE function for all use cases
result = detect_reaction(
    "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    use_ml=True  # optional: use rxn-insight ML (default: True if available)
)

# Returns consistent schema:
{
    "family": "buchwald_cn",          # Canonical family label
    "confidence": 0.95,                # 0.0-1.0 confidence score
    "method": "rule+catalyst",         # How it was detected
    "details": {                       # Full analysis details
        "reactants": [...],
        "catalysts": ["Pd"],
        "functional_groups": {...},
        "ml_prediction": {...}         # If ML was used
    }
}
```

### Benefits

✅ **One function to learn** - `detect_reaction()`  
✅ **Smart defaults** - automatically uses best available method  
✅ **Clear output** - consistent schema always  
✅ **Backwards compatible** - old functions still work, deprecated gradually  
✅ **Less code** - consolidate duplicated logic

---

## 📋 Implementation Plan

### Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    User Code (Simple)                           │
│  from chemtools import detect_reaction                          │
│  result = detect_reaction(rxn_smiles, use_ml=True)              │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│               chemtools/detection.py (NEW)                      │
│  ┌───────────────────────────────────────────────────────────┐ │
│  │ _DetectionEngine                                          │ │
│  │                                                           │ │
│  │  1. _normalize() ────────────► Reaction SMILES parsing   │ │
│  │  2. _detect_catalysts() ────► Pd, Cu, Ni, Co extraction  │ │
│  │  3. _detect_functional_groups() ─► SMARTS patterns       │ │
│  │  4. _rule_based_detection() ─► Deterministic rules       │ │
│  │  5. _ml_detection() ──────────► rxn-insight (optional)   │ │  ← ML HERE
│  │  6. _apply_catalyst_overrides() ─► Metal-based tweaks    │ │
│  │  7. detect() ─────────────────► Orchestrate all above    │ │
│  └───────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│         Consolidated Detection Methods (Internal)               │
│                                                                 │
│  Rule-Based (router.py) ────────► _rule_based_detection()      │
│  ML-Based (reaction_type_detector.py) ──► _ml_detection()      │  ← ENTIRE FILE
│  Catalyst Logic (router.py) ────► _apply_catalyst_overrides()  │     GOES HERE
└─────────────────────────────────────────────────────────────────┘
```

### Phase 1: Create Unified Core (Week 1)

**New file: `chemtools/detection.py`**

```python
"""
Unified reaction detection API.

Public Functions:
    detect_reaction(reaction_smiles, use_ml=True) -> dict

Internal:
    _DetectionEngine class (consolidates all detection logic)
"""

class _DetectionEngine:
    """Internal engine - consolidates all detection logic."""

    def __init__(self, reaction_smiles: str):
        self.reaction = reaction_smiles
        self.normalized = None
        self.reactants = []
        self.agents = []
        self.catalysts = set()
        self._normalize()

    def _normalize(self):
        """Normalize reaction and extract components."""
        # Consolidate normalization logic

    def _detect_catalysts(self):
        """Detect catalyst metals - single implementation."""
        # Move from router._detect_agent_metals()

    def _detect_functional_groups(self):
        """Detect functional groups in reactants."""
        # Move from router._rule_hits()

    def _rule_based_detection(self) -> dict:
        """Rule-based SMARTS detection."""
        # Consolidate router.detect_family() logic

    def _ml_detection(self) -> dict:
        """
        ML-based detection via rxn-insight (optional).

        This wraps the rxn-insight library to provide:
        - Broad reaction class (e.g., "C-C Coupling")
        - Specific reaction name (e.g., "Suzuki coupling with boronic acids")
        - Confidence score from ML model
        - Automatic mapping to ChemTools family taxonomy

        Returns:
            {
                "available": bool,      # Is rxn-insight installed?
                "family": str | None,   # Mapped family label
                "rxn_class": str | None,
                "rxn_name": str | None,
                "confidence": float | None
            }
        """
        # Consolidate ALL rxn-insight logic from reaction_type_detector.py here:
        # - _call_insight() - calls rxn-insight API
        # - _extract_fields() - parses ML response
        # - _map_to_family() - maps to ChemTools taxonomy
        # - _refine_cn_family() - applies catalyst-aware refinements

    def _apply_catalyst_overrides(self, family: str) -> str:
        """
        Apply catalyst-based family overrides.

        For C-N coupling reactions:
            Pd catalyst → buchwald_cn
            Cu catalyst → ullmann_cn

        For C-C couplings:
            Preserves specific family (suzuki, sonogashira, etc.)
        """
        # Consolidate override logic from router.py

    def detect(self, use_ml: bool = True) -> dict:
        """
        Main detection orchestrator - combines all detection methods.

        Detection Pipeline:
        ┌──────────────────────────────────────┐
        │ 1. Rule-Based Detection (ALWAYS)     │ ← SMARTS patterns, fast
        │    - Functional groups               │
        │    - Deterministic rules             │
        │    - Confidence: 0.7-0.9             │
        └──────────────────────────────────────┘
                        ↓
        ┌──────────────────────────────────────┐
        │ 2. Catalyst Detection (ALWAYS)       │ ← Metal extraction
        │    - Extract Pd, Cu, Ni, Co          │
        │    - From agents in reaction         │
        └──────────────────────────────────────┘
                        ↓
        ┌──────────────────────────────────────┐
        │ 3. ML Detection (OPTIONAL)           │ ← rxn-insight
        │    - Only if use_ml=True             │   (goes here!)
        │    - Only if rxn-insight installed   │
        │    - Broad class + specific name     │
        │    - Confidence: ML model score      │
        └──────────────────────────────────────┘
                        ↓
        ┌──────────────────────────────────────┐
        │ 4. Catalyst Overrides (ALWAYS)       │ ← Smart refinement
        │    - Apply metal-specific families   │
        │    - Pd → buchwald_cn                │
        │    - Cu → ullmann_cn                 │
        │    - Boost confidence to 0.95        │
        └──────────────────────────────────────┘
                        ↓
        ┌──────────────────────────────────────┐
        │ 5. Result Merging & Ranking          │ ← Choose best
        │    - Catalyst override > ML > Rule   │
        │    - Return highest confidence       │
        │    - Include all metadata            │
        └──────────────────────────────────────┘

        Args:
            use_ml: Use ML detection if available (default: True)

        Returns:
            Unified detection result with all metadata
        """
        # Implementation consolidates detect_family_from_reaction() logic


def detect_reaction(reaction_smiles: str, use_ml: bool = True) -> dict:
    """
    Detect reaction family from reaction SMILES.

    This is the MAIN entry point for all reaction detection.

    Args:
        reaction_smiles: Full reaction SMILES (reactants>>products)
        use_ml: Use ML-based detection if available (default: True)

    Returns:
        {
            "family": str,           # Canonical family label
            "confidence": float,     # 0.0-1.0
            "method": str,           # "rule", "ml", "rule+catalyst", "ml+catalyst"
            "details": {
                "reactants": [...],
                "catalysts": [...],
                "functional_groups": {...},
                "ml_prediction": {...}  # if ML used
            }
        }
    """
    engine = _DetectionEngine(reaction_smiles)
    return engine.detect(use_ml=use_ml)
```

**Consolidation Map:**

| Current Location                          | New Location                                   | Purpose                        |
| ----------------------------------------- | ---------------------------------------------- | ------------------------------ |
| `router._detect_agent_metals()`           | `_DetectionEngine._detect_catalysts()`         | Metal extraction from agents   |
| `router._rule_hits()`                     | `_DetectionEngine._detect_functional_groups()` | SMARTS pattern matching        |
| `router.detect_family()`                  | `_DetectionEngine._rule_based_detection()`     | Deterministic rule engine      |
| `reaction_type_detector.py` (entire file) | `_DetectionEngine._ml_detection()`             | **rxn-insight ML integration** |
| `router.detect_family_from_reaction()`    | `_DetectionEngine.detect()`                    | Main orchestration             |

### 🤖 Where rxn-insight ML Goes

**Current situation (confusing):**

```
chemtools/
├── router.py                    ← Rule-based detection
├── reaction_type_detector.py   ← ML detection (separate file)
│                                  ❌ User must know about this file
│                                  ❌ Separate import path
│                                  ❌ Different API
│                                  ❌ Manual mapping to taxonomy
└── taxonomy/                    ← Unified reaction taxonomy
    └── data/
        └── reaction_types.json  ← ~80 canonical reaction types
                                  ⚠️ Not automatically used by detection!
```

**New unified approach (clear):**

```
chemtools/
├── detection.py                 ← EVERYTHING in one file
│   ├── _DetectionEngine
│   │   ├── _rule_based_detection()     ← SMARTS rules (always runs)
│   │   ├── _ml_detection()             ← rxn-insight goes HERE
│   │   │                                  ✅ Optional layer
│   │   │                                  ✅ Graceful fallback
│   │   │                                  ✅ Transparent to user
│   │   │                                  ✅ Auto-maps to taxonomy
│   │   └── detect()                    ← Combines both intelligently
│   └── detect_reaction()               ← Public API (hides complexity)
└── taxonomy/                    ← **SINGLE SOURCE OF TRUTH**
    ├── registry.py              ← TaxonomyRegistry with 80+ types
    └── data/
        └── reaction_types.json  ← suzuki_miyaura, buchwald_hartwig_c_n, etc.
                                  ✅ All detection methods use this!
```

**Key Design Principles:**

1. **Unified Taxonomy is the Single Source of Truth**

   - All detection methods MUST return IDs from `chemtools/taxonomy/data/reaction_types.json`
   - 80+ canonical reaction types (suzuki_miyaura, buchwald_hartwig_c_n, heck, etc.)
   - Aliases automatically resolved via `TaxonomyRegistry.resolve_alias()`
   - No hardcoded mappings scattered across files

2. **rxn-insight is an ENHANCEMENT, not a replacement**

   - Rule-based detection ALWAYS runs (fast, deterministic baseline)
   - ML detection adds refinement when available
   - System works perfectly without rxn-insight installed
   - **ML predictions mapped to taxonomy IDs using `canonical_family_label()`**

3. **Transparent integration**

   ```python
   # User doesn't need to know about rxn-insight
   result = detect_reaction(rxn)  # Uses ML if available, rules otherwise

   # Power users can control it
   result = detect_reaction(rxn, use_ml=True)   # Prefer ML
   result = detect_reaction(rxn, use_ml=False)  # Rules only

   # ALWAYS returns canonical taxonomy IDs
   result["family"]  # e.g., "buchwald_hartwig_c_n" (not "Buchwald_CN")
   ```

4. **All ML logic in `_ml_detection()` method**

   - Consolidates entire `reaction_type_detector.py` (~304 lines)
   - Includes: `_call_insight()`, `_extract_fields()`, `_map_to_taxonomy()`, `_refine_cn_family()`
   - **Uses `chemtools.featurizers.analysis.reactions.canonical_family_label()` for mapping**
   - Single responsibility: call rxn-insight and map to unified taxonomy

5. **Detection pipeline priority**

   ```
   Priority 1 (Highest): Catalyst Override
   ├─ If Pd + C-N → "buchwald_hartwig_c_n" (conf: 0.95)
   └─ If Cu + C-N → "ullmann_cn" (conf: 0.90)

   Priority 2: ML Detection (if available)
   └─ rxn-insight → mapped via taxonomy registry

   Priority 3 (Baseline): Rule-Based
   └─ SMARTS → mapped via taxonomy registry (conf: 0.7-0.9)

   ALL outputs validated against reaction_types.json
   ```

**Example: ML detection method implementation**

```python
class _DetectionEngine:
    def _ml_detection(self) -> dict:
        """
        Optional ML enhancement via rxn-insight.

        IMPORTANT: rxn-insight is an ML model and returns UNPREDICTABLE names:
        - "Suzuki coupling" vs "Suzuki-Miyaura" vs "Cross-coupling"
        - "Buchwald-Hartwig" vs "Pd-catalyzed amination" vs "C-N coupling"
        - Names vary between runs and versions

        Solution: Robust mapping with context + conservative fallback
        """

        # Check if rxn-insight is available
        try:
            import rxn_insight
        except ImportError:
            return {
                "available": False,
                "family": None,
                "confidence": None
            }

        # Call rxn-insight (consolidated from reaction_type_detector.py)
        try:
            raw = self._call_rxn_insight()  # API wrapper
            rxn_class, rxn_name, ml_conf = self._extract_ml_fields(raw)

            # CRITICAL: rxn_name is unpredictable!
            # Could be: "Suzuki coupling", "Cross-coupling reaction",
            #           "Pd-catalyzed C-C bond formation", etc.

            # Use robust mapping with functional group + catalyst context
            family = resolve_to_taxonomy(
                rxn_name or rxn_class,
                catalysts=self.catalysts,
                is_cn_coupling="c-n" in (rxn_class or "").lower(),
                functional_groups=self.functional_groups  # Context helps!
            )

            # Adjust confidence based on mapping certainty
            # (registry exact match vs keyword fallback vs failed)
            confidence = self._calculate_ml_confidence(
                rxn_name or rxn_class,
                family,
                ml_conf or 0.5
            )

            return {
                "available": True,
                "family": family,  # May be None if unmappable
                "rxn_class": rxn_class,
                "rxn_name": rxn_name,  # Keep original for debugging
                "confidence": confidence,
                "raw": raw,
                "mapping_method": self._get_mapping_method(rxn_name, family)
            }
        except Exception as e:
            # Graceful fallback on ML failure
            logger.warning(f"ML detection failed: {e}")
            return {
                "available": True,
                "family": None,
                "confidence": None,
                "error": str(e)
            }

    def _get_mapping_method(self, ml_name: str, mapped_family: Optional[str]) -> str:
        """Track HOW we mapped unpredictable ML name."""
        if not mapped_family:
            return "failed"
        if registry.resolve_alias(ml_name) == mapped_family:
            return "exact_alias"  # High confidence
        if mapped_family in ml_name.lower():
            return "keyword_match"  # Medium confidence
        return "context_inference"  # Lower confidence (used catalysts/FGs)
```

**Benefits of this approach:**

✅ **No separate file to maintain** - `reaction_type_detector.py` becomes obsolete  
✅ **No import confusion** - one entry point: `detect_reaction()`  
✅ **Graceful degradation** - works without ML, enhanced with ML  
✅ **Clear control** - `use_ml` parameter is obvious  
✅ **Easy testing** - can test with/without rxn-insight in same test suite  
✅ **Handles ML variability** - robust mapping strategy for unpredictable names  
✅ **Context-aware** - uses catalysts + functional groups for disambiguation  
✅ **Conservative fallback** - returns None if uncertain (better than wrong!)  
✅ **Logging & tracking** - logs unmapped predictions for taxonomy improvement

### Phase 2: Deprecate Old APIs (Week 2)

**Update `router.py`:**

```python
import warnings
from .detection import detect_reaction

def detect_family_from_reaction(reaction_smiles: str, *, use_rxn_insight: bool = True):
    """
    DEPRECATED: Use chemtools.detect_reaction() instead.

    This function will be removed in v2.0.
    """
    warnings.warn(
        "detect_family_from_reaction() is deprecated. Use detect_reaction() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    result = detect_reaction(reaction_smiles, use_ml=use_rxn_insight)
    # Convert to old schema for backwards compatibility
    return {
        "family": result["family"],
        "confidence": result["confidence"],
        "hits": result["details"]["functional_groups"],
        "auto": result["details"].get("ml_prediction")
    }

def detect_family(reactants: List[str]):
    """
    DEPRECATED: Use chemtools.detect_reaction() instead.

    This function will be removed in v2.0.
    """
    warnings.warn(
        "detect_family() is deprecated. Use detect_reaction() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    # Convert reactants list to pseudo-reaction
    reaction = ".".join(reactants) + ">>"
    result = detect_reaction(reaction, use_ml=False)
    return {
        "family": result["family"],
        "confidence": result["confidence"],
        "hits": result["details"]["functional_groups"]
    }
```

**Update `reaction_type_detector.py`:**

```python
import warnings
from .detection import detect_reaction

def detect_reaction_type(reaction_smiles: str):
    """
    DEPRECATED: Use chemtools.detect_reaction() instead.

    This function will be removed in v2.0.
    """
    warnings.warn(
        "detect_reaction_type() is deprecated. Use detect_reaction() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    result = detect_reaction(reaction_smiles, use_ml=True)
    # Convert to old schema
    ml_pred = result["details"].get("ml_prediction", {})
    return {
        "available": ml_pred.get("available", False),
        "success": result["confidence"] > 0.5,
        "rxn_class": ml_pred.get("rxn_class"),
        "rxn_name": ml_pred.get("rxn_name"),
        "mapped_family": result["family"],
        "confidence": result["confidence"],
        "raw": ml_pred.get("raw"),
        "catalysts": result["details"]["catalysts"]
    }
```

### Phase 3: Update Consumers (Week 3)

**Update API endpoints (`app/main.py`):**

```python
from chemtools import detect_reaction

@app.post("/api/v1/detection")
def detect_family_endpoint(req: DetectionRequest):
    """Unified detection endpoint."""
    return detect_reaction(req.reaction_smiles, use_ml=True)

# Keep old endpoints with deprecation warnings for backwards compatibility
@app.post("/api/v1/router/detect-family")
@deprecated("Use /api/v1/detection instead")
def old_detect_family(req: DetectionRequest):
    return router.detect_family_from_reaction(req.reaction_smiles)
```

**Update LangChain tools (`lang_chain/chemtools_wrapper.py`):**

```python
from chemtools import detect_reaction

@tool
def detect_reaction_family_tool(reaction_smiles: str) -> str:
    """
    Detect reaction family from reaction SMILES.

    Args:
        reaction_smiles: Full reaction SMILES string

    Returns:
        JSON with family, confidence, and detection details
    """
    result = detect_reaction(reaction_smiles, use_ml=True)
    return json.dumps(result, indent=2)
```

**Update CLIs:**

```python
from chemtools import detect_reaction

# In local_recommendation_cli.py:
detection = detect_reaction(reaction, use_ml=False)  # Rule-based only
family = detection["family"]
confidence = detection["confidence"]

# In cli_AI_recommend.py:
detection = detect_reaction(reaction_smiles, use_ml=True)  # ML-enabled
```

### Phase 4: Update Documentation (Week 3)

- Update `REACTION_DETECTION_METHODS.md` to focus on `detect_reaction()`
- Update `README.md` examples
- Create migration guide for users

### Phase 5: Testing & Validation (Week 4)

- Write comprehensive tests for `detect_reaction()`
- Ensure backwards compatibility tests pass
- Performance benchmarking (should be faster with consolidated code)
- Integration testing with all consumers

---

## 📊 Before & After Comparison

### Before (Complex)

```python
# User confusion - which function?
from chemtools.router import detect_family, detect_family_from_reaction
from chemtools.reaction_type_detector import detect_reaction_type

# 3 different functions, different signatures, different outputs
result1 = detect_family(["Brc1ccccc1", "Nc1ccccc1"])
result2 = detect_family_from_reaction("Brc1ccccc1.Nc1ccccc1>>...")
result3 = detect_reaction_type("Brc1ccccc1.Nc1ccccc1>>...")

# Different output schemas - hard to use
family1 = result1["family"]
family2 = result2["family"]
family3 = result3["mapped_family"]  # Why different key? 🤔
```

### After (Simple)

```python
# One import, one function
from chemtools import detect_reaction

# Simple API - works for all cases
result = detect_reaction("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")

# Consistent schema
family = result["family"]
confidence = result["confidence"]
method = result["method"]  # "rule+catalyst"
```

---

## 🎓 Migration Guide for Users

### Quick Migration

**Old code:**

```python
from chemtools.router import detect_family_from_reaction
result = detect_family_from_reaction(reaction, use_rxn_insight=True)
family = result["family"]
```

**New code:**

```python
from chemtools import detect_reaction
result = detect_reaction(reaction, use_ml=True)
family = result["family"]
```

**Changes:**

- Import from `chemtools` (top-level) instead of `chemtools.router`
- Rename `use_rxn_insight` → `use_ml` (clearer name)
- Same output schema for `family` and `confidence` keys

### Detailed Migration Table

| Old Function                                             | New Function                                                        | Notes                                    |
| -------------------------------------------------------- | ------------------------------------------------------------------- | ---------------------------------------- |
| `detect_family(reactants)`                               | `detect_reaction(".".join(reactants) + ">>")`                       | Convert reactant list to pseudo-reaction |
| `detect_family_from_reaction(rxn, use_rxn_insight=True)` | `detect_reaction(rxn, use_ml=True)`                                 | Direct replacement                       |
| `detect_reaction_type(rxn)`                              | `detect_reaction(rxn, use_ml=True)`                                 | Same function, different output schema   |
| `classify_reactants_with_context(rxn)`                   | `detect_reaction(rxn)` then access `result["details"]["reactants"]` | Details now in nested structure          |

### Deprecation Timeline

- **v1.8 (Current)**: Old functions work with deprecation warnings
- **v1.9**: Old functions still work but log warnings prominently
- **v2.0**: Old functions removed, only `detect_reaction()` available

---

## 📈 Expected Benefits

### Code Reduction

- **Before**: ~850 lines across 3 files
- **After**: ~400 lines in 1 file
- **Savings**: 53% less code to maintain

### Performance

- **No redundant normalization** (currently done multiple times)
- **Cached functional group detection**
- **Single pass through reaction components**
- **Estimated speedup**: 20-30% faster

### User Experience

- **1 function to learn** vs 5
- **1 output schema** vs 4 different schemas
- **Clear documentation** - one page, not scattered
- **Easier debugging** - single code path

### Maintainability

- **Single source of truth** for detection logic
- **Easier to add new families** (one place to update)
- **Simpler testing** (one function to test thoroughly)
- **Less duplication** (DRY principle)

---

## ⚠️ Risks & Mitigation

### Risk: Breaking Changes

**Mitigation:**

- Keep old functions working with deprecation warnings
- Provide clear migration guide
- Version bump to v2.0 for removals
- 6-month deprecation period

### Risk: Performance Regression

**Mitigation:**

- Benchmark before/after
- Profile code paths
- Optimize critical sections
- Maintain performance tests

### Risk: Missing Edge Cases

**Mitigation:**

- Comprehensive test suite (100+ test cases)
- Test against existing production data
- Gradual rollout (opt-in beta period)
- Keep old code available for fallback

---

## 🚀 Implementation Priority

### Must Have (Phase 1-2)

✅ Create `chemtools/detection.py` with `detect_reaction()`  
✅ Consolidate all detection logic into `_DetectionEngine`  
✅ Add deprecation warnings to old functions  
✅ Ensure backwards compatibility

### Should Have (Phase 3)

✅ Update all internal consumers (API, CLI, LangChain)  
✅ Update documentation  
✅ Write comprehensive tests

### Nice to Have (Phase 4)

✅ Performance optimizations (caching, vectorization)  
✅ Enhanced ML integration (multiple models)  
✅ Interactive detection visualizer

---

## 📝 Next Steps

1. **Review this plan** - Get stakeholder approval
2. **Create feature branch** - `feature/unified-detection-api`
3. **Implement Phase 1** - Core `detect_reaction()` function
4. **Write tests** - Ensure 100% coverage of new code
5. **Gradual rollout** - Beta testing with select users
6. **Full deployment** - Update all consumers
7. **Monitor metrics** - Performance, error rates, user feedback

---

## 🙋 Questions for Review

1. **API Design**: Is `detect_reaction(reaction_smiles, use_ml=True)` the right signature?
2. **Output Schema**: Should we add more fields to the response (e.g., `alternatives`, `warnings`)?
3. **Deprecation Timeline**: Is 6 months enough notice before removing old functions?
4. **Performance**: Should we add caching layer for repeated detections?
5. **ML Integration**: Should we support multiple ML models (not just rxn-insight)?

---

**Status**: 📋 **Proposal - Pending Approval**  
**Est. Implementation Time**: 4 weeks  
**Breaking Changes**: Yes (in v2.0, with 6-month deprecation)  
**Backwards Compatibility**: Yes (until v2.0)

---

## 📚 Related Documentation

- **`TAXONOMY_ALIGNMENT_STRATEGY.md`** - How to align detection with unified taxonomy
- **`RXN_INSIGHT_NAMING_GUIDE.md`** - ⚠️ Critical: Handling rxn-insight's unpredictable naming
- **`REACTION_DETECTION_METHODS.md`** - Current detection methods overview

---

## 📌 Quick Reference: Where Everything Goes

### File Consolidation Map

| What                    | From (Current)                            | To (New)                                       | Lines      |
| ----------------------- | ----------------------------------------- | ---------------------------------------------- | ---------- |
| **Public API**          | `router.detect_family_from_reaction()`    | `detection.detect_reaction()`                  | Main entry |
| **Rule Detection**      | `router.detect_family()`                  | `_DetectionEngine._rule_based_detection()`     | ~150       |
| **SMARTS Patterns**     | `router._rule_hits()`                     | `_DetectionEngine._detect_functional_groups()` | ~80        |
| **Catalyst Extraction** | `router._detect_agent_metals()`           | `_DetectionEngine._detect_catalysts()`         | ~60        |
| **🤖 ML Detection**     | `reaction_type_detector.py` (entire file) | `_DetectionEngine._ml_detection()`             | **~150**   |
| **Catalyst Overrides**  | `router.detect_family_from_reaction()`    | `_DetectionEngine._apply_catalyst_overrides()` | ~40        |
| **Normalization**       | `router.detect_family_from_reaction()`    | `_DetectionEngine._normalize()`                | ~30        |

**Total consolidation**: ~850 lines → ~400 lines (53% reduction)

### Method Hierarchy (New Structure)

```
detect_reaction(rxn, use_ml=True)          ← PUBLIC API
    ↓
_DetectionEngine(rxn)                      ← Internal class
    ├─ __init__()                          ← Initialize
    ├─ _normalize()                        ← Parse reaction SMILES
    ├─ _detect_catalysts()                 ← Extract Pd/Cu/Ni/Co
    ├─ _detect_functional_groups()         ← SMARTS matching
    ├─ _rule_based_detection()             ← Deterministic rules
    ├─ _ml_detection()                     ← 🤖 rxn-insight (OPTIONAL)
    ├─ _apply_catalyst_overrides()         ← Metal-based tweaks
    └─ detect(use_ml)                      ← Orchestrate all above
```

### rxn-insight Integration Details

**All ML logic consolidated into `_ml_detection()` method:**

```python
def _ml_detection(self) -> dict:
    """
    Optional ML detection using rxn-insight library.

    Consolidates from reaction_type_detector.py:
    ───────────────────────────────────────────────
    Line 74-97:   _call_insight() → Try different API entrypoints
    Line 122-169: _extract_fields() → Parse ML response
    Line 172-203: _map_to_family() → Map to taxonomy
    Line 206-235: _refine_cn_family() → Catalyst refinements
    Line 238-304: detect_reaction_type() → Main orchestrator

    Returns:
    ───────
    {
        "available": bool,        # Is rxn-insight installed?
        "family": str | None,     # Mapped to ChemTools taxonomy
        "rxn_class": str | None,  # Broad ML class
        "rxn_name": str | None,   # Specific ML name
        "confidence": float,      # ML model score
        "raw": dict              # Raw rxn-insight output
    }
    """
    # Implementation here (~150 lines)
```

**Key design decisions:**

1. ✅ **Optional dependency** - system works without rxn-insight
2. ✅ **Graceful fallback** - catches import/runtime errors
3. ✅ **Transparent** - user doesn't need to know about ML
4. ✅ **Controlled** - `use_ml` parameter gives explicit control
5. ✅ **Testable** - can mock ML responses in tests

---

Would you like to proceed with this simplification? Any concerns or suggestions?
