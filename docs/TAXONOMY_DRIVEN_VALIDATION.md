# Taxonomy-Driven Validation Refactor

## Overview

Refactored `chemtools/featurizers/formatters/detection_validation.py` to be fully **taxonomy-driven** instead of hardcoding motif patterns. This eliminates duplication and technical debt.

## What Changed

### Before (Hardcoded Patterns)
```python
def _has_aryl_halide(motifs: Set[str]) -> bool:
    """Check for aryl/heteroaryl halide electrophiles."""
    aryl_halides = ["Ar-Br", "Ar-Cl", "Ar-I", "Ar-F"]  # Hardcoded!
    heteroaryl_halides = ["HeteroAr-Br", "HeteroAr-Cl", "HeteroAr-I", "HeteroAr-F"]
    alkenyl_halides = ["Alkenyl-Br", "Alkenyl-Cl", "Alkenyl-I", "Alkenyl-F"]
    return any(m in motifs for m in aryl_halides + heteroaryl_halides + alkenyl_halides)
```

Had **6 hardcoded validation rules**:
1. Suzuki: `if _has_organoboron() and _has_aryl_halide():`
2. C-N Coupling: `if _has_aryl_halide() and _has_amine():`
3. Stille: `if _has_organotin() and _has_aryl_halide():`
4. Negishi: `if _has_organozinc() and _has_aryl_halide():`
5. Heck: `if _has_aryl_halide() and "Alkenyl" in reacted_set:`
6. Sonogashira: `if _has_aryl_halide() and any(m in reacted_set for m in ["Alkynyl", ...]):`

**Problems:**
- ❌ Duplication: Same patterns in taxonomy files AND validation code (two sources of truth)
- ❌ Maintenance burden: Bug fixes require updating multiple places (e.g., HeteroAr-I bug)
- ❌ Incomplete: Only covered 6 reactions when taxonomy has 9+ coupling reactions
- ❌ Manual work: Each new reaction required adding hardcoded rules

### After (Taxonomy-Driven)
```python
def validate_detection_with_reacted_motifs(...) -> Dict[str, Any]:
    # Load taxonomy
    definitions, _ = _get_catalog()
    
    # Check each reaction definition in the catalog for pattern matches
    for reaction_id, defn in definitions.items():
        if not defn.reactants:
            continue
            
        # Check if all reactant slots match the reacted motifs
        all_slots_match = True
        matched_slots = []
        
        for slot_name, slot_req in defn.reactants.items():
            if not slot_req.allowed:
                continue
                
            if _motifs_match_slot(reacted_set, slot_req.allowed):
                matched_slots.append(slot_name)
            else:
                all_slots_match = False
                break
        
        # If all required reactant slots matched, check product formation
        if all_slots_match and len(matched_slots) >= 2:
            # Check if expected products are formed
            product_match = True
            if defn.products:
                product_match = False
                for slot_name, slot_req in defn.products.items():
                    if slot_req.allowed and _motifs_match_slot(formed_set, slot_req.allowed):
                        product_match = True
                        break
            
            # If pattern matches and it's different from initial detection, correct it
            if product_match and reaction_id != initial_detection:
                return {
                    "reaction_type": reaction_id,
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": f"Taxonomy pattern: {' + '.join(matched_slots)} → {defn.name}",
                }
```

**Benefits:**
- ✅ Single source of truth: Patterns only in taxonomy files
- ✅ Automatic coverage: Supports ALL reactions in taxonomy (not just 6)
- ✅ Zero maintenance: New reactions work automatically (no code changes needed)
- ✅ Bug fixes propagate: Fix taxonomy = fix validation instantly
- ✅ Future-proof: Schema evolution doesn't require validation code changes

## Taxonomy Updates

Also updated `chemtools/taxonomy/data/compound_logic.json` to include **HeteroAr variants** (proper fix vs previous bandaid):

```json
"sp2_electrophiles": {
  "members": [
    "Ar-I", "Ar-Br", "Ar-Cl", "Ar-F",
    "HeteroAr-I", "HeteroAr-Br", "HeteroAr-Cl", "HeteroAr-F",  // Added
    "Alkenyl-I", "Alkenyl-Br", "Alkenyl-Cl", "Alkenyl-F"
  ]
},
"organoboron": {
  "members": [
    "Ar-B(OH)2", "Ar-Bpin", "Ar-BF3K", "Ar-B(OR)2",
    "HeteroAr-B(OH)2", "HeteroAr-Bpin", "HeteroAr-BF3K", "HeteroAr-B(OR)2",  // Added
    "Alkenyl-B(OH)2", "Alkenyl-Bpin", "Alkenyl-BF3K", "Alkenyl-B(OR)2"
  ]
}
```

This is the **root cause fix** for the HeteroAr-I detection bug.

## Legacy Functions

Old hardcoded helper functions are now **deprecated but kept for backward compatibility**:
- `_has_organoboron()`
- `_has_aryl_halide()`
- `_has_amine()`
- `_has_organotin()`
- `_has_organozinc()`
- `_has_organometallic_nucleophile()`

All marked with `[DEPRECATED]` docstrings.

## Testing

Created comprehensive test suite (`test_taxonomy_validation.py`):

```
Testing taxonomy-driven validation...
============================================================
✓ Suzuki with HeteroAr-I: PASSED
✓ Suzuki with HeteroAr-B(OH)2: PASSED
✓ C-N coupling: PASSED
✓ C-N coupling with HeteroAr-Br: PASSED
✓ No correction needed: PASSED
✓ Organometallic exclusion: PASSED
============================================================
✅ All taxonomy-driven validation tests PASSED!
```

User's original Suzuki reaction now works perfectly:
```
SMILES: Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1
Reaction Type:
  - reaction_type: Suzuki_miyaura
  - confidence: 0.9500
```

## Architecture Impact

### Before: Two Sources of Truth
```
taxonomy/data/reaction_types.v4.0.json  ← Definition
       └─> "@sp2_electrophiles"
       
taxonomy/data/compound_logic.json      ← Expansion
       └─> "Ar-I", "Ar-Br", ... (missing HeteroAr)
       
featurizers/formatters/detection_validation.py  ← Validation (DUPLICATED!)
       └─> ["Ar-Br", "Ar-Cl", "Ar-I", "HeteroAr-Br", ...]  ← Had to patch manually
```

**Result:** Bug in compound_logic.json → had to bandaid in validation code

### After: Single Source of Truth
```
taxonomy/data/reaction_types.v4.0.json  ← Definition
       └─> "@sp2_electrophiles"
       
taxonomy/data/compound_logic.json      ← Single source
       └─> "Ar-I", "Ar-Br", "HeteroAr-I", "HeteroAr-Br", ...  ← Fixed once
       
featurizers/formatters/detection_validation.py  ← Reads from taxonomy
       └─> load_reaction_catalog() → uses taxonomy patterns
```

**Result:** Bug fix in compound_logic.json → automatically fixed everywhere

## Files Modified

1. **chemtools/featurizers/formatters/detection_validation.py**
   - Added taxonomy imports
   - Replaced 6 hardcoded validation rules with single taxonomy-driven loop
   - Added `_get_catalog()`, `_motifs_match_slot()`, `_has_organometallic_nucleophile_taxonomy()`
   - Deprecated old helper functions (kept for compatibility)

2. **chemtools/taxonomy/data/compound_logic.json**
   - Added HeteroAr-I, HeteroAr-Br, HeteroAr-Cl, HeteroAr-F to `sp2_electrophiles`
   - Added HeteroAr-B(OH)2, HeteroAr-Bpin, HeteroAr-BF3K, HeteroAr-B(OR)2 to `organoboron`

3. **test_taxonomy_validation.py** (new)
   - Comprehensive test suite for taxonomy-driven validation
   - 6 test cases covering Suzuki, C-N coupling, HeteroAr variants, exclusion rules

## Migration Notes

### For Developers
- ✅ **No breaking changes**: Public API remains identical
- ✅ **Drop-in replacement**: Old code keeps working
- ✅ **Performance**: Cached catalog loading (LRU cache)
- ⚠️ **New reactions**: Just add to taxonomy - validation works automatically

### For Maintainers
- 🎯 **Single edit point**: Fix patterns in compound_logic.json only
- 🎯 **No code changes**: Validation logic is now schema-driven
- 🎯 **Better coverage**: ALL reactions (not just 6) are validated
- 🎯 **Test once**: Add reaction to taxonomy → automatically covered

## Future Work

1. **Remove deprecated functions** (in next major version):
   - Mark `_has_organoboron()` etc. for removal
   - Add deprecation warnings

2. **Extend validation** (now trivial with taxonomy):
   - Add product chirality checks from taxonomy metadata
   - Add stereochemistry validation rules
   - Add regioselectivity patterns

3. **Performance optimization** (if needed):
   - Pre-compile all SMARTS patterns at startup
   - Add reaction-specific validation caches

## Conclusion

This refactor achieves the goal of making validation **taxonomy-driven** with:
- ✅ Zero duplication (single source of truth)
- ✅ Automatic coverage (all reactions)
- ✅ Easy maintenance (edit taxonomy, not code)
- ✅ Proper fix for HeteroAr detection (in taxonomy)
- ✅ Full backward compatibility (no breaking changes)

The system is now **schema-driven** and **future-proof** - adding new reactions or motifs no longer requires code changes in validation logic.
