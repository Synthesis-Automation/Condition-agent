# Taxonomy Alignment Strategy for Detection System

## 🎯 Problem: Detection Not Using Unified Taxonomy

### Current State (Broken)

```
Detection Methods                   Unified Taxonomy
─────────────────                   ────────────────
router.py                           taxonomy/
├─ Returns: "suzuki_miyaura"    ┐   ├─ registry.py
├─ Returns: "Suzuki_CC"         │   │   └─ TaxonomyRegistry
├─ Returns: "Suzuki"            │   │       ├─ 80+ reaction types
                                ├──X│       ├─ Alias resolution
reaction_type_detector.py       │   │       └─ Validation
├─ Returns: "Suzuki coupling"   │   │
├─ Manual mapping table         │   └─ data/
├─ Hardcoded overrides          │       └─ reaction_types.json
                                │           ├─ suzuki_miyaura
❌ No centralized validation    ┘           ├─ buchwald_hartwig_c_n
❌ Inconsistent IDs returned                ├─ ullmann_cn
❌ Taxonomy not used as source of truth     └─ ... 77 more types
```

**Problems:**

1. **3 different sources of truth**
   - `router.py` has hardcoded family names
   - `reaction_type_detector.py` has separate mapping table
   - `taxonomy/reaction_types.json` exists but isn't used

2. **No validation**
   - Detection can return IDs not in taxonomy
   - Aliases handled inconsistently
   - No guarantee output is valid

3. **Duplicate mapping logic**
   - `reaction_type_detector._map_to_family()` - manual ML→taxonomy mapping
   - `router.detect_family()` - hardcoded SMARTS→family
   - `analysis/reactions.canonical_family_label()` - taxonomy lookup
   - **All should use the same taxonomy!**

---

## ✨ Solution: Taxonomy-First Detection

### New Architecture

```
┌──────────────────────────────────────────────────────────┐
│  chemtools/taxonomy/ (SINGLE SOURCE OF TRUTH)            │
│  ┌────────────────────────────────────────────────────┐  │
│  │ TaxonomyRegistry                                   │  │
│  │ ├─ 80+ canonical reaction IDs                      │  │
│  │ ├─ Alias resolution (Suzuki → suzuki_miyaura)      │  │
│  │ ├─ Category mapping                                │  │
│  │ └─ Validation (is ID valid?)                       │  │
│  └────────────────────────────────────────────────────┘  │
└──────────────────────────────────────────────────────────┘
                          ↑
                          │ ALL detection uses this
                          │
┌──────────────────────────────────────────────────────────┐
│  chemtools/detection.py (NEW unified detector)           │
│  ┌────────────────────────────────────────────────────┐  │
│  │ _DetectionEngine                                   │  │
│  │                                                    │  │
│  │ _rule_based_detection():                          │  │
│  │   ├─ SMARTS match: "aryl_halide + boron"          │  │
│  │   ├─ Lookup via registry.resolve_alias("Suzuki")  │  │
│  │   └─ Returns: "suzuki_miyaura" ✅                  │  │
│  │                                                    │  │
│  │ _ml_detection():                                  │  │
│  │   ├─ rxn-insight: "Suzuki coupling"              │  │
│  │   ├─ Lookup via registry.resolve_alias(...)      │  │
│  │   └─ Returns: "suzuki_miyaura" ✅                  │  │
│  │                                                    │  │
│  │ _apply_catalyst_overrides():                      │  │
│  │   ├─ Pd + C-N detected                            │  │
│  │   ├─ Lookup via registry.get_reaction_type(...)  │  │
│  │   └─ Returns: "buchwald_hartwig_c_n" ✅           │  │
│  │                                                    │  │
│  │ ALL OUTPUTS VALIDATED AGAINST TAXONOMY            │  │
│  └────────────────────────────────────────────────────┘  │
└──────────────────────────────────────────────────────────┘
```

---

## 🔧 Implementation Steps

### Step 1: Create Taxonomy Mapping Helper

**New file: `chemtools/detection_mapper.py`**

```python
"""
Centralized taxonomy mapping for all detection methods.
All detection logic must use these functions to ensure taxonomy alignment.
"""

from typing import Optional, Set
from .taxonomy import load_registry
from .taxonomy.registry import TaxonomyRegistry
from .analysis.reactions import canonical_family_label

_registry: Optional[TaxonomyRegistry] = None

def get_taxonomy_registry() -> TaxonomyRegistry:
    """Get cached taxonomy registry."""
    global _registry
    if _registry is None:
        _registry = load_registry()
    return _registry


def resolve_to_taxonomy(
    raw_prediction: str,
    *,
    catalysts: Optional[Set[str]] = None,
    is_cn_coupling: bool = False
) -> Optional[str]:
    """
    Resolve ANY detection result to canonical taxonomy ID.
    
    This is the ONLY function that should be used to convert
    detection results to taxonomy IDs.
    
    Args:
        raw_prediction: Raw string from SMARTS, ML, or heuristic
        catalysts: Set of detected catalyst metals (Pd, Cu, Ni, Co)
        is_cn_coupling: Whether C-N coupling signature detected
        
    Returns:
        Canonical taxonomy ID from reaction_types.json or None
        
    Examples:
        >>> resolve_to_taxonomy("Suzuki")
        "suzuki_miyaura"
        
        >>> resolve_to_taxonomy("C-N Coupling", catalysts={"Pd"})
        "buchwald_hartwig_c_n"
        
        >>> resolve_to_taxonomy("Suzuki coupling with boronic acids")
        "suzuki_miyaura"
    """
    if not raw_prediction:
        return None
    
    # 1. Try direct lookup in taxonomy
    registry = get_taxonomy_registry()
    if registry.get_reaction_type(raw_prediction):
        return raw_prediction
    
    # 2. Try alias resolution
    canonical = canonical_family_label(raw_prediction)
    if canonical:
        # 3. Apply catalyst overrides for C-N couplings
        if catalysts and (is_cn_coupling or canonical in CN_FAMILIES):
            canonical = apply_catalyst_override_with_taxonomy(
                canonical, 
                catalysts,
                is_cn_coupling=is_cn_coupling
            )
        return canonical
    
    # 4. No match found
    return None


def apply_catalyst_override_with_taxonomy(
    family: str,
    catalysts: Set[str],
    *,
    is_cn_coupling: bool
) -> str:
    """
    Apply catalyst-based overrides using taxonomy validation.
    
    Ensures returned ID exists in taxonomy.
    """
    from .analysis.reactions import apply_catalyst_override, CN_FAMILIES_CANONICAL
    
    override = apply_catalyst_override(family, catalysts, is_cn_coupling=is_cn_coupling)
    
    # Validate override exists in taxonomy
    registry = get_taxonomy_registry()
    if registry.get_reaction_type(override):
        return override
    
    # Fallback to original if override invalid
    if registry.get_reaction_type(family):
        return family
    
    return override  # Return anyway, will be caught by validation


def validate_taxonomy_id(reaction_id: str) -> bool:
    """Check if reaction_id is valid in taxonomy."""
    registry = get_taxonomy_registry()
    return registry.get_reaction_type(reaction_id) is not None


def get_all_taxonomy_ids() -> Set[str]:
    """Get all valid reaction type IDs from taxonomy."""
    registry = get_taxonomy_registry()
    return set(registry.reaction_types.keys())


# Constants from taxonomy
CN_FAMILIES = {
    "cn_coupling",
    "buchwald_hartwig_c_n",
    "ullmann_cn",
    "chan_lam",
    "snar_cn"
}
```

### Step 2: Update Detection to Use Taxonomy

**In `chemtools/detection.py`:**

```python
from .detection_mapper import resolve_to_taxonomy, validate_taxonomy_id

class _DetectionEngine:
    def _rule_based_detection(self) -> dict:
        """Rule-based SMARTS detection - returns taxonomy IDs."""
        
        # Existing SMARTS matching logic...
        hits = self._detect_functional_groups()
        
        # Determine raw family based on SMARTS
        if hits.get("aryl_halide") and hits.get("boron"):
            raw_family = "Suzuki"  # Raw prediction
        elif hits.get("aryl_halide") and hits.get("nucleophile_n"):
            raw_family = "C-N Coupling"
        # ... more rules
        
        # RESOLVE TO TAXONOMY (centralized)
        family = resolve_to_taxonomy(
            raw_family,
            catalysts=self.catalysts,
            is_cn_coupling=self._is_cn_signature()
        )
        
        # Validate
        if family and not validate_taxonomy_id(family):
            raise ValueError(f"Detection returned invalid taxonomy ID: {family}")
        
        return {
            "family": family or "Unknown",
            "confidence": self._calculate_confidence(hits),
            "hits": hits,
            "method": "rule"
        }
    
    def _ml_detection(self) -> dict:
        """ML-based detection via rxn-insight - returns taxonomy IDs."""
        
        # Call rxn-insight
        raw_result = self._call_rxn_insight()
        rxn_class, rxn_name, conf = self._extract_ml_fields(raw_result)
        
        # Raw ML prediction might be:
        # - "Suzuki coupling with boronic acids"
        # - "C-N Coupling"
        # - "Heteroatom Alkylation"
        
        # RESOLVE TO TAXONOMY (centralized)
        family = resolve_to_taxonomy(
            rxn_name or rxn_class,  # Try specific name first
            catalysts=self.catalysts,
            is_cn_coupling="C-N" in (rxn_class or "")
        )
        
        # Validate
        if family and not validate_taxonomy_id(family):
            raise ValueError(f"ML detection returned invalid taxonomy ID: {family}")
        
        return {
            "family": family,
            "rxn_class": rxn_class,
            "rxn_name": rxn_name,
            "confidence": conf,
            "method": "ml"
        }
```

### Step 3: Remove Duplicate Mapping Logic

**Delete these hardcoded mappings:**

1. ❌ `reaction_type_detector._map_to_family()` - replaced by `resolve_to_taxonomy()`
2. ❌ `router.detect_family()` direct family strings - use taxonomy lookups
3. ❌ Hardcoded family constants scattered across files

**Keep only:**
✅ `chemtools/taxonomy/data/reaction_types.json` - canonical source
✅ `chemtools/analysis/reactions.canonical_family_label()` - uses taxonomy
✅ `chemtools/detection_mapper.resolve_to_taxonomy()` - centralized resolver

---

## 📊 Mapping Examples

### Rule-Based Detection

```python
# BEFORE (hardcoded)
if hits["aryl_halide"] and hits["boron"]:
    family = "suzuki_miyaura"  # Hardcoded string
    
# AFTER (taxonomy-aligned)
if hits["aryl_halide"] and hits["boron"]:
    raw = "Suzuki"  # Can be any alias
    family = resolve_to_taxonomy(raw, catalysts=catalysts)
    # Returns: "suzuki_miyaura" (validated against taxonomy)
```

### ML Detection

```python
# BEFORE (manual mapping)
if "suzuki" in rxn_name.lower():
    family = "suzuki_miyaura"  # Manual if/else chain
elif "buchwald" in rxn_name.lower():
    family = "buchwald_hartwig_c_n"
# ... 20 more manual mappings
    
# AFTER (taxonomy-aligned)
family = resolve_to_taxonomy(rxn_name, catalysts=catalysts)
# Taxonomy registry handles ALL aliases automatically
# "Suzuki coupling" → "suzuki_miyaura"
# "Buchwald-Hartwig" → "buchwald_hartwig_c_n"
# No manual mapping needed!
```

### Catalyst Overrides

```python
# BEFORE (scattered logic)
if "Pd" in catalysts and family in CN_FAMILIES:
    family = "buchwald_hartwig_c_n"  # Hardcoded
    
# AFTER (taxonomy-aligned)
family = resolve_to_taxonomy(
    "C-N Coupling",
    catalysts={"Pd"},
    is_cn_coupling=True
)
# Returns: "buchwald_hartwig_c_n" (validated)
```

---

## ✅ Validation Strategy

### Runtime Validation

```python
def detect(self, use_ml: bool = True) -> dict:
    """Main detection - validates all outputs."""
    
    # Run detection
    rule_result = self._rule_based_detection()
    ml_result = self._ml_detection() if use_ml else None
    
    # Validate ALL results against taxonomy
    for result in [rule_result, ml_result]:
        if result and result.get("family"):
            family = result["family"]
            if not validate_taxonomy_id(family):
                raise ValueError(
                    f"Detection produced invalid taxonomy ID: {family}. "
                    f"Valid IDs: {get_all_taxonomy_ids()}"
                )
    
    # Merge results...
    return final_result
```

### Test Suite

```python
# tests/test_taxonomy_alignment.py

def test_all_detections_use_taxonomy():
    """Ensure all detection methods return valid taxonomy IDs."""
    
    registry = load_registry()
    valid_ids = set(registry.reaction_types.keys())
    
    test_reactions = [
        "Brc1ccccc1.c1ccc(B(O)O)cc1>>...",  # Suzuki
        "Brc1ccccc1.Nc1ccccc1>>...",         # C-N coupling
        # ... 50 more test cases
    ]
    
    for rxn in test_reactions:
        result = detect_reaction(rxn)
        family = result["family"]
        
        assert family in valid_ids, (
            f"Detection returned '{family}' which is not in taxonomy. "
            f"Valid IDs: {sorted(valid_ids)}"
        )


def test_alias_resolution():
    """Test that aliases resolve to canonical IDs."""
    
    test_cases = [
        ("Suzuki", "suzuki_miyaura"),
        ("Buchwald-Hartwig", "buchwald_hartwig_c_n"),
        ("Ullmann CN", "ullmann_cn"),
        ("C-N Coupling", "cn_coupling"),
        ("Heck", "heck"),
    ]
    
    for alias, expected_id in test_cases:
        result = resolve_to_taxonomy(alias)
        assert result == expected_id


def test_ml_predictions_map_correctly():
    """Test rxn-insight predictions map to taxonomy."""
    
    ml_predictions = [
        ("Suzuki coupling with boronic acids", "suzuki_miyaura"),
        ("Buchwald-Hartwig amination", "buchwald_hartwig_c_n"),
        ("Heteroatom Alkylation", "cn_coupling"),  # With catalyst override
    ]
    
    for ml_name, expected_id in ml_predictions:
        result = resolve_to_taxonomy(ml_name, catalysts={"Pd"})
        assert result == expected_id
```

---

## 🎯 Benefits

### Before (Broken Alignment)

```python
# Different sources of truth
router.detect_family() → "suzuki_miyaura"       # Hardcoded
detect_reaction_type() → "Suzuki coupling"      # ML output
taxonomy.json          → "suzuki_miyaura"       # Canonical

# No validation
family = "suzuki_cc"  # Typo! But no error raised
recommend_conditions(family)  # Fails silently
```

### After (Perfect Alignment)

```python
# Single source of truth
detect_reaction() → "suzuki_miyaura"  # From taxonomy
taxonomy.json     → "suzuki_miyaura"  # Canonical

# Automatic validation
family = "suzuki_cc"  # Typo!
resolve_to_taxonomy(family)  # Returns None (not in taxonomy)
detect_reaction() # Raises: ValueError("Invalid taxonomy ID")

# Aliases work automatically
resolve_to_taxonomy("Suzuki") → "suzuki_miyaura" ✅
resolve_to_taxonomy("Suzuki coupling") → "suzuki_miyaura" ✅
resolve_to_taxonomy("Buchwald-Hartwig") → "buchwald_hartwig_c_n" ✅
```

---

## 📝 Migration Checklist

### Phase 1: Infrastructure
- [ ] Create `chemtools/detection_mapper.py` with `resolve_to_taxonomy()`
- [ ] Add validation functions using `TaxonomyRegistry`
- [ ] Write comprehensive test suite for taxonomy alignment

### Phase 2: Detection Updates
- [ ] Update `_rule_based_detection()` to use `resolve_to_taxonomy()`
- [ ] Update `_ml_detection()` to use `resolve_to_taxonomy()`
- [ ] Update catalyst override logic to use taxonomy validation
- [ ] Remove hardcoded family constants from `router.py`

### Phase 3: Cleanup
- [ ] Delete `reaction_type_detector._map_to_family()` (replaced)
- [ ] Delete manual mapping tables
- [ ] Consolidate all alias resolution to taxonomy registry

### Phase 4: Validation
- [ ] Run full test suite with taxonomy validation enabled
- [ ] Test against production data to ensure no regressions
- [ ] Verify all 80+ reaction types can be detected

---

## 🚀 Result

**Before:**
- 3 sources of truth (router, detector, taxonomy)
- Manual mapping tables scattered across files
- No validation of detection outputs
- Inconsistent IDs returned

**After:**
- **1 source of truth**: `chemtools/taxonomy/data/reaction_types.json`
- **Automatic alias resolution** via `TaxonomyRegistry`
- **Runtime validation** ensures only valid IDs returned
- **100% alignment** between detection and taxonomy

---

## ⚠️ Important: rxn-insight Naming Variability

### The Problem

**rxn-insight is an ML model and returns unpredictable, varying names**:

```python
# Same reaction, different names returned:
"Brc1ccccc1.c1ccc(B(O)O)cc1>>..."

# Run 1: "Suzuki coupling"
# Run 2: "Suzuki-Miyaura coupling with boronic acids"
# Run 3: "Cross-coupling reaction"
# Run 4: "C-C bond formation"
# Run 5: "Pd-catalyzed coupling"

❌ Inconsistent
❌ Not guaranteed to match our taxonomy
❌ Can change between versions
```

### Why This Matters

1. **Cannot rely on exact string matching**
   ```python
   # ❌ WRONG - brittle
   if ml_result == "Suzuki coupling":
       family = "suzuki_miyaura"
   ```

2. **Need fuzzy/semantic mapping**
   ```python
   # ✅ CORRECT - flexible
   family = resolve_to_taxonomy(ml_result, catalysts=catalysts)
   # Handles: "Suzuki", "Suzuki coupling", "Suzuki-Miyaura", 
   #          "Cross-coupling", etc.
   ```

3. **May return names not in our taxonomy**
   - rxn-insight might say: "Heterocyclic arylation"
   - Our taxonomy has: "buchwald_hartwig_c_n"
   - Need intelligent mapping, not exact match

### Solution: Robust Mapping Strategy

```python
def resolve_to_taxonomy(
    raw_prediction: str,
    *,
    catalysts: Optional[Set[str]] = None,
    is_cn_coupling: bool = False,
    functional_groups: Optional[Dict[str, bool]] = None
) -> Optional[str]:
    """
    Resolve unpredictable ML predictions to canonical taxonomy IDs.
    
    Handles rxn-insight's variable naming by:
    1. Alias resolution (registry handles common variations)
    2. Keyword/pattern matching (fallback for new variations)
    3. Context-aware mapping (uses catalysts + functional groups)
    4. Conservative fallback (returns None if uncertain)
    """
    
    # Strategy 1: Try registry alias resolution (handles known variations)
    canonical = registry.resolve_alias(raw_prediction)
    if canonical:
        return canonical
    
    # Strategy 2: Keyword-based fuzzy matching
    # rxn-insight might say "Pd-catalyzed amination" instead of "Buchwald-Hartwig"
    lower = raw_prediction.lower()
    
    # C-C Couplings
    if any(kw in lower for kw in ["suzuki", "boronic", "boron coupling"]):
        if functional_groups and functional_groups.get("boron"):
            return "suzuki_miyaura"
    
    if any(kw in lower for kw in ["sonogashira", "alkyne coupling"]):
        return "sonogashira"
    
    if any(kw in lower for kw in ["heck", "vinyl coupling", "olefin arylation"]):
        return "heck"
    
    # C-N Couplings (context-aware!)
    if any(kw in lower for kw in ["amination", "c-n coupling", "amine arylation", 
                                   "heterocyclic arylation", "buchwald", "hartwig"]):
        # Use catalyst context to disambiguate
        if catalysts:
            if "Pd" in catalysts:
                return "buchwald_hartwig_c_n"
            elif "Cu" in catalysts:
                return "ullmann_cn"
        # Generic if no catalyst info
        return "cn_coupling"
    
    if any(kw in lower for kw in ["chan-lam", "chan lam", "oxidative coupling"]):
        return "chan_lam"
    
    # Strategy 3: Category-level matching
    # If specific name unclear, map to general category
    if "cross-coupling" in lower or "c-c bond formation" in lower:
        # Use functional groups to disambiguate
        if functional_groups:
            if functional_groups.get("boron"):
                return "suzuki_miyaura"
            elif functional_groups.get("terminal_alkyne"):
                return "sonogashira"
        return None  # Too ambiguous
    
    # Strategy 4: Conservative fallback
    # If none match, return None (better than wrong classification)
    return None


def _calculate_ml_confidence(
    ml_prediction: str,
    mapped_family: Optional[str],
    original_confidence: float
) -> float:
    """
    Adjust confidence based on mapping certainty.
    
    rxn-insight confidence doesn't account for our taxonomy mapping.
    """
    if not mapped_family:
        return 0.0  # No mapping = no confidence
    
    # Exact alias match = high confidence
    if registry.resolve_alias(ml_prediction) == mapped_family:
        return original_confidence * 0.95
    
    # Keyword-based match = moderate confidence
    if mapped_family in ml_prediction.lower():
        return original_confidence * 0.85
    
    # Context-based inference = lower confidence
    return original_confidence * 0.70
```

### Updated ML Detection Method

```python
class _DetectionEngine:
    def _ml_detection(self) -> dict:
        """
        ML-based detection handling rxn-insight's unpredictable naming.
        """
        try:
            import rxn_insight
        except ImportError:
            return {"available": False, "family": None}
        
        # Get ML prediction (may be unpredictable!)
        raw_result = self._call_rxn_insight()
        rxn_class, rxn_name, ml_conf = self._extract_ml_fields(raw_result)
        
        # Examples of what rxn-insight might return:
        # - "Suzuki coupling with boronic acids"
        # - "Pd-catalyzed cross-coupling"
        # - "Heterocyclic C-N bond formation"
        # - "Buchwald-Hartwig amination"
        # All different ways to describe same/similar reactions!
        
        # Robust mapping with context
        family = resolve_to_taxonomy(
            rxn_name or rxn_class,
            catalysts=self.catalysts,
            is_cn_coupling="c-n" in (rxn_class or "").lower(),
            functional_groups=self.functional_groups
        )
        
        # Adjust confidence based on mapping certainty
        confidence = _calculate_ml_confidence(
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
            "method": "ml",
            "mapping_method": "exact" if registry.resolve_alias(rxn_name) 
                            else "keyword" if family else "failed"
        }
```

### Fallback Strategy

```python
def detect(self, use_ml: bool = True) -> dict:
    """
    Main detection with robust ML fallback.
    """
    # ALWAYS run rule-based (reliable baseline)
    rule_result = self._rule_based_detection()
    
    # Try ML if requested
    ml_result = None
    if use_ml:
        ml_result = self._ml_detection()
    
    # Decision logic:
    # 1. If ML mapped successfully with high confidence → use ML
    # 2. If ML failed to map → use rule-based
    # 3. If both available → prefer higher confidence
    
    if ml_result and ml_result.get("family"):
        ml_conf = ml_result["confidence"]
        rule_conf = rule_result["confidence"]
        
        # ML only wins if significantly more confident
        if ml_conf > rule_conf + 0.1:
            return {
                "family": ml_result["family"],
                "confidence": ml_conf,
                "method": "ml",
                "fallback": rule_result  # Keep for debugging
            }
    
    # Default to rule-based (more reliable)
    return {
        "family": rule_result["family"],
        "confidence": rule_result["confidence"],
        "method": "rule",
        "ml_attempted": ml_result is not None,
        "ml_mapping_failed": ml_result and not ml_result.get("family")
    }
```

### Logging & Debugging

```python
import logging

logger = logging.getLogger(__name__)

def resolve_to_taxonomy(raw_prediction: str, **kwargs) -> Optional[str]:
    """Resolve with logging for unknown ML outputs."""
    
    canonical = _try_resolve(raw_prediction, **kwargs)
    
    if not canonical:
        logger.warning(
            f"ML prediction '{raw_prediction}' could not be mapped to taxonomy. "
            f"Consider adding alias to taxonomy/data/aliases.json"
        )
        # Log for future taxonomy expansion
        _log_unmapped_prediction(raw_prediction, **kwargs)
    
    return canonical


def _log_unmapped_prediction(prediction: str, **kwargs):
    """Track unmapped predictions to improve taxonomy coverage."""
    unmapped_log = Path("data/unmapped_ml_predictions.jsonl")
    with open(unmapped_log, "a") as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "ml_prediction": prediction,
            "catalysts": list(kwargs.get("catalysts", [])),
            "functional_groups": kwargs.get("functional_groups", {})
        }, f)
        f.write("\n")
```

---

## 🙋 Questions?

1. **What if rxn-insight returns a reaction not in our taxonomy?**
   - `resolve_to_taxonomy()` returns `None` (conservative)
   - Falls back to rule-based detection (reliable)
   - Logs warning with prediction for future taxonomy expansion
   - **Better to admit uncertainty than misclassify!**

2. **How to handle new ML naming variations?**
   - Check `data/unmapped_ml_predictions.jsonl` for patterns
   - Add new aliases to `taxonomy/data/aliases.json`
   - Update keyword matching in `resolve_to_taxonomy()`
   - No code changes to core detection logic needed

3. **What if ML and rule-based disagree?**
   - ML only wins if confidence > rule-based + 0.1
   - Prefer rule-based (deterministic) over uncertain ML
   - Log disagreements for analysis
   - Return both results for debugging

4. **How to add new reaction types?**
   - Add to `taxonomy/data/reaction_types.json`
   - Add common variations to `aliases.json`
   - Update keyword patterns if needed
   - Test with real rxn-insight outputs

5. **What about legacy code using old family names?**
   - Taxonomy aliases handle migration
   - `"Suzuki_CC"` → resolves to `"suzuki_miyaura"`
   - Deprecation warnings guide updates

---

**Status**: 📋 **Strategy Document**  
**Dependencies**: DETECTION_SIMPLIFICATION_PLAN.md  
**Implementation**: Integrated into unified detection system
