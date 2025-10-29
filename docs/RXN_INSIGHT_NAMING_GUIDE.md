# rxn-insight ML Naming Variability - Quick Reference

## ⚠️ Critical Issue

**rxn-insight is an ML model** → returns **UNPREDICTABLE, VARYING names**

### Same Reaction, Different Names

```python
reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"

# rxn-insight might return any of these:
- "Suzuki coupling"
- "Suzuki-Miyaura coupling"
- "Suzuki coupling with boronic acids"
- "Cross-coupling reaction"
- "Pd-catalyzed cross-coupling"
- "C-C bond formation"
- "Aryl-aryl coupling"

❌ Cannot rely on exact string matching!
```

---

## ❌ WRONG Approach (Brittle)

```python
# Don't do this - will break!
def map_ml_to_family(ml_name: str) -> str:
    if ml_name == "Suzuki coupling":  # Only matches ONE variation
        return "suzuki_miyaura"
    elif ml_name == "Buchwald-Hartwig":  # Will miss "Pd-catalyzed amination"
        return "buchwald_hartwig_c_n"
    # Breaks when rxn-insight changes naming!
```

---

## ✅ CORRECT Approach (Robust)

### Strategy 1: Alias Resolution (Taxonomy Registry)

```python
# Let taxonomy handle known variations
family = registry.resolve_alias(ml_name)

# Handles:
# "Suzuki" → "suzuki_miyaura"
# "Suzuki coupling" → "suzuki_miyaura"
# "Suzuki-Miyaura" → "suzuki_miyaura"
# "Buchwald-Hartwig" → "buchwald_hartwig_c_n"
# "BH coupling" → "buchwald_hartwig_c_n"
```

### Strategy 2: Keyword Matching (Fallback)

```python
# For new/unknown variations
lower = ml_name.lower()

if any(kw in lower for kw in ["suzuki", "boronic", "boron coupling"]):
    return "suzuki_miyaura"

if any(kw in lower for kw in ["amination", "buchwald", "hartwig", "c-n coupling"]):
    # Use catalyst context!
    if "Pd" in catalysts:
        return "buchwald_hartwig_c_n"
    elif "Cu" in catalysts:
        return "ullmann_cn"
    return "cn_coupling"
```

### Strategy 3: Context-Aware Mapping

```python
def resolve_to_taxonomy(
    ml_name: str,
    catalysts: Set[str],
    functional_groups: Dict[str, bool]
) -> Optional[str]:
    """
    Use ALL available context to disambiguate.
    """
    
    # Vague ML output: "Cross-coupling reaction"
    # Use functional groups to narrow down:
    
    if "cross-coupling" in ml_name.lower():
        if functional_groups.get("boron"):
            return "suzuki_miyaura"
        elif functional_groups.get("terminal_alkyne"):
            return "sonogashira"
        elif functional_groups.get("nucleophile_n"):
            # Use catalyst context
            if "Pd" in catalysts:
                return "buchwald_hartwig_c_n"
        # Too ambiguous - return None
        return None
```

### Strategy 4: Conservative Fallback

```python
# If uncertain, admit it!
family = resolve_to_taxonomy(ml_name, catalysts, functional_groups)

if not family:
    # Fall back to rule-based detection (more reliable)
    return rule_based_result
    
# Log for future improvement
logger.warning(f"Could not map ML prediction: {ml_name}")
```

---

## 🎯 Complete Implementation

```python
def resolve_to_taxonomy(
    ml_prediction: str,
    *,
    catalysts: Optional[Set[str]] = None,
    functional_groups: Optional[Dict[str, bool]] = None
) -> Optional[str]:
    """
    Robust mapping handling rxn-insight's unpredictable naming.
    
    Returns:
        Canonical taxonomy ID or None if uncertain
    """
    
    # 1. Try exact alias (handles common variations)
    canonical = registry.resolve_alias(ml_prediction)
    if canonical:
        return canonical  # High confidence
    
    # 2. Try lowercase alias
    canonical = registry.resolve_alias(ml_prediction.lower())
    if canonical:
        return canonical
    
    # 3. Keyword-based fuzzy matching
    lower = ml_prediction.lower()
    
    # Suzuki family
    if any(kw in lower for kw in ["suzuki", "boronic", "boron coupling"]):
        if functional_groups and functional_groups.get("boron"):
            return "suzuki_miyaura"
    
    # Sonogashira family
    if any(kw in lower for kw in ["sonogashira", "alkyne coupling", "terminal alkyne"]):
        return "sonogashira"
    
    # Heck family
    if any(kw in lower for kw in ["heck", "mizoroki", "vinyl arylation"]):
        return "heck"
    
    # C-N coupling family (context-dependent!)
    if any(kw in lower for kw in [
        "amination", "c-n", "buchwald", "hartwig", "ullmann",
        "heteroatom alkylation", "amine arylation", "chan-lam"
    ]):
        # Disambiguate using catalyst
        if catalysts:
            if "Pd" in catalysts:
                return "buchwald_hartwig_c_n"
            elif "Cu" in catalysts:
                # Check for Chan-Lam (oxidative)
                if "oxidative" in lower or "chan" in lower:
                    return "chan_lam"
                return "ullmann_cn"
        # Generic fallback
        return "cn_coupling"
    
    # 4. Category-level (requires context)
    if "cross-coupling" in lower or "c-c bond" in lower:
        # Use functional groups
        if functional_groups:
            if functional_groups.get("boron"):
                return "suzuki_miyaura"
            elif functional_groups.get("terminal_alkyne"):
                return "sonogashira"
            elif functional_groups.get("grignard"):
                return "kumada"
        # Too vague
        return None
    
    # 5. Conservative: admit uncertainty
    return None
```

---

## 📊 Confidence Adjustment

```python
def calculate_ml_confidence(
    ml_name: str,
    mapped_family: Optional[str],
    original_ml_confidence: float
) -> float:
    """
    Adjust confidence based on HOW we mapped the name.
    """
    
    if not mapped_family:
        return 0.0  # Failed to map
    
    # Exact alias match = trust ML confidence
    if registry.resolve_alias(ml_name) == mapped_family:
        return original_ml_confidence * 0.95
    
    # Keyword match = moderate confidence
    if mapped_family in ml_name.lower():
        return original_ml_confidence * 0.85
    
    # Context-based inference = lower confidence
    return original_ml_confidence * 0.70
```

---

## 🔍 Logging & Improvement

```python
def resolve_to_taxonomy(ml_prediction: str, **kwargs) -> Optional[str]:
    """Resolve with tracking for unknown variations."""
    
    canonical = _try_resolve(ml_prediction, **kwargs)
    
    if not canonical:
        # Log for future taxonomy improvement
        logger.warning(
            f"Unmapped ML prediction: '{ml_prediction}'. "
            f"Consider adding to taxonomy/data/aliases.json"
        )
        
        # Track for analysis
        with open("data/unmapped_ml_predictions.jsonl", "a") as f:
            json.dump({
                "timestamp": datetime.now().isoformat(),
                "ml_prediction": ml_prediction,
                "catalysts": list(kwargs.get("catalysts", [])),
                "functional_groups": kwargs.get("functional_groups", {})
            }, f)
            f.write("\n")
    
    return canonical
```

### Periodic Review

```bash
# Review unmapped predictions
cat data/unmapped_ml_predictions.jsonl | jq -r '.ml_prediction' | sort | uniq -c | sort -rn

# Example output:
# 15 "Pd-catalyzed amination"       ← Add alias!
# 12 "Aryl-aryl cross-coupling"     ← Add alias!
#  8 "Heterocyclic C-N coupling"    ← Add alias!
#  3 "Oxidative N-arylation"        ← Add alias!

# Add to taxonomy/data/aliases.json:
{
  "text": "Pd-catalyzed amination",
  "entity_type": "reaction_type",
  "entity_id": "buchwald_hartwig_c_n"
}
```

---

## 🎯 Key Principles

1. **Never use exact string matching with ML outputs**
   - ML naming is unpredictable
   - Changes between versions
   - Use fuzzy/semantic matching

2. **Always use context (catalysts + functional groups)**
   - ML might say "amination" (vague)
   - Pd catalyst → Buchwald-Hartwig
   - Cu catalyst → Ullmann

3. **Be conservative**
   - If uncertain, return None
   - Fall back to rule-based (reliable)
   - **Better to admit uncertainty than misclassify**

4. **Track and improve**
   - Log unmapped predictions
   - Review periodically
   - Add new aliases to taxonomy
   - No code changes needed!

5. **Validate everything**
   - Check mapped ID exists in taxonomy
   - Raise error if invalid
   - Prevent downstream failures

---

## 🧪 Test Coverage

```python
# Test known variations
@pytest.mark.parametrize("ml_name,expected", [
    ("Suzuki coupling", "suzuki_miyaura"),
    ("Suzuki-Miyaura coupling with boronic acids", "suzuki_miyaura"),
    ("Cross-coupling reaction", None),  # Too vague without context
    ("Buchwald-Hartwig amination", "buchwald_hartwig_c_n"),
    ("Pd-catalyzed C-N coupling", "buchwald_hartwig_c_n"),  # With Pd catalyst
])
def test_ml_name_variations(ml_name, expected):
    result = resolve_to_taxonomy(ml_name, catalysts={"Pd"})
    assert result == expected


# Test context-dependent disambiguation
def test_context_aware_mapping():
    # Same vague name, different context
    
    # With boron → Suzuki
    assert resolve_to_taxonomy(
        "Cross-coupling",
        functional_groups={"boron": True}
    ) == "suzuki_miyaura"
    
    # With alkyne → Sonogashira
    assert resolve_to_taxonomy(
        "Cross-coupling",
        functional_groups={"terminal_alkyne": True}
    ) == "sonogashira"
    
    # With amine + Pd → Buchwald-Hartwig
    assert resolve_to_taxonomy(
        "Amination",
        catalysts={"Pd"},
        functional_groups={"nucleophile_n": True}
    ) == "buchwald_hartwig_c_n"
```

---

**See also:**
- `TAXONOMY_ALIGNMENT_STRATEGY.md` - Full implementation guide
- `DETECTION_SIMPLIFICATION_PLAN.md` - Unified detection architecture
- `chemtools/taxonomy/data/aliases.json` - Canonical alias definitions
