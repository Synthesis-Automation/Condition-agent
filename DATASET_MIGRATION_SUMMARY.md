# Dataset Structure Migration Summary

## Overview

Successfully migrated the reaction dataset structure from metal-specific C-N coupling datasets to a unified approach:

**Old Structure:**
- `C_N_Coupling_Pd` (Buchwald-Hartwig, Palladium-catalyzed)
- `C_N_Coupling_Cu` (Ullmann, Copper-catalyzed)  
- `C_N_Coupling_Ni` (Nickel-catalyzed)

**New Structure:**
- `C_N_Coupling` (Unified - all metals)
- `C_O_Coupling` (NEW - C-O coupling reactions)
- `C_S_Coupling` (NEW - C-S coupling reactions)

## Motivation

1. **Simplified Dataset Management**: One unified C-N coupling dataset instead of three metal-specific ones
2. **Metal Preference via Constraints**: Metal selection now handled through the recommendation engine's constraint system (more flexible)
3. **Expanded Coverage**: Added support for C-O and C-S coupling reactions
4. **Cleaner Architecture**: Separation of concerns - dataset structure vs. recommendation filtering

## Files Modified

### 1. Core Chemistry Tools

#### `chemtools/recommend/utils.py`
**Changes:**
- Updated `_FAMILY_ALIASES` dictionary to map all C-N variants to unified `C_N_Coupling`
- Added legacy aliases for backward compatibility:
  - `C_N_Coupling_Cu` → `C_N_Coupling`
  - `C_N_Coupling_Pd` → `C_N_Coupling`
  - `C_N_Coupling_Ni` → `C_N_Coupling`
  - `Ullmann_CN` → `C_N_Coupling`
  - `Buchwald_CN` → `C_N_Coupling`
- Added new coupling types:
  - `Ullmann_CO` → `C_O_Coupling`
  - `Ullmann_CS` → `C_S_Coupling`

**Impact:** All code using `canonical_family()` automatically gets correct mapping.

#### `chemtools/router.py`
**Changes:**
1. Added S-nucleophile detection to SMARTS patterns:
   ```python
   "nucleophile_s": Chem.MolFromSmarts("[SX2H]")
   ```

2. Updated `detect_family()` function:
   - C-N Coupling: Returns unified `C_N_Coupling` (was `Ullmann_CN`)
   - C-O Coupling: Returns `C_O_Coupling` (was `Ullmann_O` or ignored)
   - C-S Coupling: Returns `C_S_Coupling` (NEW)
   - Suzuki still returns `Suzuki_CC`
   - Amide still returns `Amide_Coupling`

3. Updated `_apply_catalyst_override()`:
   - Deprecated metal-specific routing
   - Now returns unified `C_N_Coupling` regardless of metal
   - Added deprecation comment explaining metal preference is now via constraints

4. Updated `detect_family_from_reaction()`:
   - Added intelligent preference logic for C-O and C-S coupling
   - Prefers rule-based detection over rxn_insight for heteroatom couplings
   - Maintains backward compatibility with rxn_insight for other reactions

**Impact:** All reaction detection automatically uses new unified structure.

### 2. CLI Application

#### `app/cli_recommend.py`
**Changes:**
1. Removed metal-specific routing logic in `determine_final_reaction_type()`
   - Deleted code that routed C_N_Coupling + Cu → C_N_Coupling_Cu
   - Deleted code that routed C_N_Coupling + Pd → C_N_Coupling_Pd
   - Deleted code that routed C_N_Coupling + Ni → C_N_Coupling_Ni

2. Updated JSON schema description:
   - Old: "e.g., Suzuki, Buchwald, Ullmann, C_N_Coupling"
   - New: "e.g., Suzuki, C_N_Coupling, C_O_Coupling, C_S_Coupling"

3. Metal preference now handled via constraints:
   - `metal_preference` field in constraints
   - `required_reagents` array includes metal catalyst specifications
   - Recommendation engine filters/ranks based on these constraints

**Impact:** CLI routes to unified dataset, metal preference via constraints.

### 3. Data Processing

#### `data-processor/Scifinder_rdf_processer.py`
**Changes:**
1. Updated SciFinder reaction type mapping:
   ```python
   # Old
   "buchwald": "C_N_Coupling_Pd",
   "ullmann": "C_N_Coupling_Cu",
   
   # New
   "buchwald": "C_N_Coupling",
   "ullmann": "C_N_Coupling",
   "c-o coupling": "C_O_Coupling",
   "c-s coupling": "C_S_Coupling",
   ```

2. Updated filename comment:
   - Old: `(e.g., "Suzuki", "C_N_Coupling_Cu")`
   - New: `(e.g., "Suzuki", "C_N_Coupling", "C_O_Coupling")`

**Impact:** All newly processed reactions go to correct unified datasets.

### 4. Documentation & Examples

#### `chemtools/__init__.py`
- Updated example from `C_N_Coupling_Pd` to `C_N_Coupling`

#### `chemtools/util/drfp_storage.py`
- Updated docstring examples from `C_N_Coupling_Cu` to `C_N_Coupling`
- Updated `get_drfp_path_for_family()` docstring

#### `chemtools/precedent/search.py`
- Updated example from `C_N_Coupling_Pd` to `C_N_Coupling`

#### `chemtools/dataset_analytics.py`
- Updated docstring example from `C_N_Coupling_Pd` to `C_N_Coupling`

#### `chemtools/schema/output_core_format.json`
- Updated example JSON from `C_N_Coupling_Cu` to `C_N_Coupling`

## Dataset Files

### Current Structure (Verified)
```
data/reaction_dataset/
├── Amide_formation.jsonl
├── Amide_formation_drfp.npz
├── C_N_Coupling.jsonl          ← UNIFIED (was 3 separate files)
├── C_N_Coupling_drfp.npz
├── C_O_Coupling.jsonl          ← NEW
├── C_O_Coupling_drfp.npz
├── C_S_Coupling.jsonl          ← NEW
├── C_S_Coupling_drfp.npz
├── Suzuki.jsonl
└── Suzuki_drfp.npz
```

### Migration Notes
- Old `C_N_Coupling_Cu.jsonl`, `C_N_Coupling_Pd.jsonl`, `C_N_Coupling_Ni.jsonl` were combined into `C_N_Coupling.jsonl`
- DRFP fingerprints were recombined into single `C_N_Coupling_drfp.npz`
- Metal-specific information preserved in reaction metadata
- New C-O and C-S coupling datasets added

## Backward Compatibility

### Aliasing System
All old family names automatically map to new structure via `canonical_family()`:

```python
canonical_family("C_N_Coupling_Cu")   # → "C_N_Coupling"
canonical_family("C_N_Coupling_Pd")   # → "C_N_Coupling"
canonical_family("C_N_Coupling_Ni")   # → "C_N_Coupling"
canonical_family("Ullmann_CN")        # → "C_N_Coupling"
canonical_family("Buchwald_CN")       # → "C_N_Coupling"
canonical_family("Ullmann_CO")        # → "C_O_Coupling"
canonical_family("Ullmann_CS")        # → "C_S_Coupling"
```

### API Compatibility
- Existing API calls with old family names work correctly
- Family names automatically normalized via `canonical_family()`
- No breaking changes to API contracts

### Metal Preference Handling
**Old Way (DEPRECATED):**
```python
# System automatically routed to metal-specific datasets
reaction_type = "C_N_Coupling_Cu"  # Forced copper-based results
```

**New Way (RECOMMENDED):**
```python
# Unified dataset with metal preference via constraints
{
  "reaction_type": "C_N_Coupling",
  "constraints": {
    "metal_preference": "Cu",
    "required_reagents": ["copper catalyst"]
  }
}
```

## Testing

### Test File: `tests/test_dataset_update.py`

**Test Coverage:**
1. ✅ C-N Coupling detection → Returns `C_N_Coupling`
2. ✅ C-O Coupling detection → Returns `C_O_Coupling`
3. ✅ C-S Coupling detection → Returns `C_S_Coupling`
4. ✅ Suzuki detection → Returns `Suzuki_CC` → Canonical `Suzuki`
5. ✅ All legacy aliases map correctly to new structure

**Test Results:**
```
C-N Coupling:
  Detected: C_N_Coupling (confidence: 0.90)
  Canonical: C_N_Coupling

C-O Coupling:
  Detected: C_O_Coupling (confidence: 0.85)
  Canonical: C_O_Coupling

C-S Coupling:
  Detected: C_S_Coupling (confidence: 0.85)
  Canonical: C_S_Coupling

Suzuki:
  Detected: Suzuki_CC (confidence: 0.90)
  Canonical: Suzuki
```

## Benefits

### 1. Simplified Data Management
- **Before**: 3 separate C-N coupling datasets to maintain
- **After**: 1 unified C-N coupling dataset
- **Impact**: Easier updates, consistent versioning, simpler deployment

### 2. More Flexible Metal Preference
- **Before**: Hard-coded routing based on detected/specified family
- **After**: Dynamic filtering via constraint system
- **Impact**: Can prefer Cu but still see Pd options, better recommendations

### 3. Expanded Reaction Coverage
- **Before**: Only C-N coupling (3 variants) + Suzuki + Amide
- **After**: C-N + C-O + C-S + Suzuki + Amide
- **Impact**: Broader reaction types supported

### 4. Cleaner Architecture
- **Before**: Dataset selection mixed with recommendation logic
- **After**: Clear separation - dataset structure vs. recommendation filtering
- **Impact**: More maintainable, easier to extend

### 5. Better Heteroatom Coupling Support
- **Before**: C-O coupling poorly supported (called "Ullmann_O", low priority)
- **After**: First-class support for C-O and C-S coupling
- **Impact**: Better recommendations for these important reaction types

## Migration Checklist

- [x] Update family alias mappings
- [x] Update SMARTS patterns for S-nucleophile
- [x] Update detection logic for C-O and C-S
- [x] Remove metal-specific routing from CLI
- [x] Update SciFinder processor mappings
- [x] Update all docstring examples
- [x] Update schema examples
- [x] Create and run comprehensive tests
- [x] Verify backward compatibility
- [x] Document all changes

## Future Work

### Recommendations
1. **Update Documentation**: README files should reflect new structure
2. **API Migration Guide**: Create guide for external users migrating from old family names
3. **Performance Testing**: Verify unified dataset doesn't impact search performance
4. **Constraint Optimization**: Tune metal preference filtering in recommendation engine

### Considerations
1. **Dataset Size**: Monitor C_N_Coupling.jsonl size as it combines 3 datasets
2. **DRFP Index**: Ensure fingerprint loading performance is acceptable
3. **Metal Distribution**: Analyze metal distribution in unified dataset
4. **Recommendation Quality**: A/B test recommendation quality vs. old system

## Rollback Plan

If issues arise, rollback is straightforward:

1. **Revert Code Changes**: Git revert all commits
2. **Restore Datasets**: Keep old metal-specific datasets as backup
3. **Update Aliases**: Temporarily map new names to old structure
4. **Gradual Migration**: Support both old and new simultaneously

## Summary

✅ **Successfully migrated from metal-specific to unified C-N coupling datasets**
✅ **Added C-O and C-S coupling support**  
✅ **Maintained 100% backward compatibility via aliasing**
✅ **Improved architecture with constraint-based metal preference**
✅ **All tests passing**

The codebase is now ready for the new unified dataset structure!
