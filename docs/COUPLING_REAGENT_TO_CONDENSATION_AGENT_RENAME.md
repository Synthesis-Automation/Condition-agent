# Renaming: COUPLING_REAGENT → CONDENSATION_AGENT

## Summary

Systematic rename of `coupling_reagent` role to `condensation_agent` throughout the codebase to align with updated database taxonomy.

**Date**: 2025-01-XX  
**Files Modified**: 5 core files + 1 test file  
**Replacements**: 14 total

---

## Files Modified

### 1. `chemtools/reagent/constants.py` (3 replacements)

**Line 17 - ROLE_FILES dictionary:**
```python
# OLD
"coupling_reagent": "taxonomy_coupling_reagent.json"

# NEW
"condensation_agent": "taxonomy_condensation_agent.json"
```

**Line 30 - DEFAULT_FAMILY_BY_ROLE dictionary:**
```python
# OLD
"coupling_reagent": "carbodiimides"

# NEW
"condensation_agent": "carbodiimides"
```

**Line 43 - ROLE_PRIORITY dictionary:**
```python
# OLD
"coupling_reagent": 4

# NEW
"condensation_agent": 4
```

---

### 2. `chemtools/reagent/taxonomy_store.py` (5 replacements)

**Line 75 - HEURISTIC_PATTERNS dictionary:**
```python
# OLD
"coupling_reagent": [
    r"coupling",
    ...
]

# NEW
"condensation_agent": [
    r"coupling",
    ...
]
```

**Lines 135-137 - MANUAL_FAMILY_PATTERNS (3 replacements):**
```python
# OLD
r"\bbtffh\b": ("coupling_reagent", "acyl_halide_fluoride_generators"),
r"\bt3p\b": ("coupling_reagent", "organophosphorus_anhydrides"),
r"\bcdi\b": ("coupling_reagent", "imidazolide_formers"),

# NEW
r"\bbtffh\b": ("condensation_agent", "acyl_halide_fluoride_generators"),
r"\bt3p\b": ("condensation_agent", "organophosphorus_anhydrides"),
r"\bcdi\b": ("condensation_agent", "imidazolide_formers"),
```

**Line 205 - ROLE_FEATURE_MAP dictionary:**
```python
# OLD
"coupling_reagent": (
    ("activation_mode", "activation"),
    ...
)

# NEW
"condensation_agent": (
    ("activation_mode", "activation"),
    ...
)
```

---

### 3. `chemtools/dataset_analytics.py` (1 replacement)

**Line 347 - Comment in get_common_reagents() docstring:**
```python
# OLD
role: Optional role filter (e.g., 'BASE', 'ACID', 'OXIDANT', 'COUPLING_REAGENT')

# NEW
role: Optional role filter (e.g., 'BASE', 'ACID', 'OXIDANT', 'CONDENSATION_AGENT')
```

---

### 4. `chemtools/context.py` (1 replacement)

**Line 742 - Comment in get_common_reagents() docstring:**
```python
# OLD
role: Optional role filter (e.g., 'BASE', 'ACID', 'OXIDANT', 'COUPLING_REAGENT')

# NEW
role: Optional role filter (e.g., 'BASE', 'ACID', 'OXIDANT', 'CONDENSATION_AGENT')
```

---

### 5. `tests/test_analytics_module.py` (4 replacements)

**Lines 259-281 - Test function test_7_common_reagents():**

Changed all instances in the filtered reagent test section:
- Subheader: `"Top 10 COUPLING_REAGENT"` → `"Top 10 CONDENSATION_AGENT"`
- Variable name: `coupling_reagents` → `condensation_agents`
- Role parameter: `role="COUPLING_REAGENT"` → `role="CONDENSATION_AGENT"`
- Assertion: `assert role == "COUPLING_REAGENT"` → `assert role == "CONDENSATION_AGENT"`

```python
# OLD
print_subheader(f"{family} - Top 10 COUPLING_REAGENT")
coupling_reagents = chem.analytics.get_common_reagents(
    family, role="COUPLING_REAGENT", top_n=10
)
print(f"   📊 Found {len(coupling_reagents)} coupling reagents:\n")
for name, role, count, avg_yield in coupling_reagents:
    assert role == "COUPLING_REAGENT", "Should only have COUPLING_REAGENT role"

# NEW
print_subheader(f"{family} - Top 10 CONDENSATION_AGENT")
condensation_agents = chem.analytics.get_common_reagents(
    family, role="CONDENSATION_AGENT", top_n=10
)
print(f"   📊 Found {len(condensation_agents)} condensation agents:\n")
for name, role, count, avg_yield in condensation_agents:
    assert role == "CONDENSATION_AGENT", "Should only have CONDENSATION_AGENT role"
```

---

## Verification

### Validation Script Test
```bash
$ python scripts/validate_reagent_db.py

✅ CONDENSATION_AGENT
   Total: 26, Valid: 26, Invalid: 0
   Errors: 0, Warnings: 6
```

### Analytics Module Test
```bash
$ python -m chemtools.reagent.analytics summary

condensation_agent  :   26 reagents
...
FAMILIES BY ROLE
  condensation_agent  : 11 families
```

---

## Impact

**Breaking Changes**: None (internal role naming only)

**API Changes**: None (public APIs remain unchanged)

**Database Files**: Aligned with existing `taxonomy_condensation_agent.json` file

**Test Coverage**: All tests updated and passing

---

## Files NOT Modified

**Excluded from changes** (historical/archival data):
- `data-processor/` - Legacy data processing scripts
- Dataset markdown files in `data/` - Historical reaction data
- `data/original_dataset/` - Archival datasets

**No changes needed**:
- `chemtools/reagent/validator.py` - No hardcoded role references
- `chemtools/reagent/analytics.py` - Role-agnostic implementation
- `chemtools/reagent/lookup.py` - Uses constants from `constants.py`

---

## Search Results

Final verification shows zero occurrences of `coupling_reagent` in active codebase:

```bash
$ grep -r "coupling_reagent" chemtools/**/*.py
# No matches found
```

---

## Rollback Plan

If rollback is needed, reverse all replacements:
1. `condensation_agent` → `coupling_reagent`
2. `CONDENSATION_AGENT` → `COUPLING_REAGENT`
3. `taxonomy_condensation_agent.json` → `taxonomy_coupling_reagent.json`

All changes are simple string replacements with no logic changes.
