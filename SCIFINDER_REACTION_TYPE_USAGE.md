# Using SciFinder Reaction Type in Preprocessing

## Enhancement Summary

**Key Insight**: SciFinder RDF files already contain reaction type metadata. We should use this as the primary source for family detection instead of re-detecting from SMILES patterns.

## Changes Made

### 1. Priority Order for Family Detection

**New approach** (in `generate_jsonl_export()`):

```python
# 1. Try SciFinder reaction type first (exact match)
scifinder_type = row.get("ReactionType", "").strip()
if scifinder_type.lower() in scifinder_map:
    detected_family = scifinder_map[scifinder_type.lower()]
    family_source = "scifinder"
    family_confidence = 1.0  # High confidence

# 2. Try partial match (e.g., "Buchwald-Hartwig" contains "Buchwald")
elif scifinder_type:
    for key, family in scifinder_map.items():
        if key in scifinder_lower or scifinder_lower in key:
            detected_family = family
            family_source = "scifinder_partial"
            family_confidence = 0.8
            break

# 3. Fallback: SMARTS detection from reactant patterns
if detected_family == "Unknown":
    fam_result = router.detect_family(reactants_normalized)
    detected_family = fam_result.get("family", "Unknown")
    family_confidence = fam_result.get("confidence", 0.0)
    family_source = "smarts_detection"
```

### 2. SciFinder Reaction Type Mapping

**Mapping table** for common SciFinder names:

```python
scifinder_map = {
    # C-N Coupling variations
    "buchwald": "C_N_Coupling_Pd",
    "buchwald-hartwig": "C_N_Coupling_Pd",
    "buchwald c-n": "C_N_Coupling_Pd",
    "buchwald_cn": "C_N_Coupling_Pd",
    "ullmann": "C_N_Coupling_Cu",
    "ullmann c-n": "C_N_Coupling_Cu",
    "ullmann_cn": "C_N_Coupling_Cu",
    "goldberg": "C_N_Coupling_Cu",
    
    # Suzuki variations
    "suzuki": "Suzuki_CC",
    "suzuki-miyaura": "Suzuki_CC",
    "suzuki coupling": "Suzuki_CC",
    
    # Others
    "sonogashira": "Sonogashira_CC",
    "amide": "Amide_Coupling",
    "amide coupling": "Amide_Coupling",
    "amide formation": "Amide_Coupling",
}
```

### 3. Enhanced Precomputed Data

**New fields in dataset:**

```json
{
  "precomputed": {
    "reaction_smiles": "c1ccccc1Br.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "normalized": "c1ccccc1Br.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "reactants_normalized": ["c1ccccc1Br", "Nc1ccccc1"],
    "detected_family": "C_N_Coupling_Pd",
    "family_confidence": 1.0,
    "family_source": "scifinder",           // ← NEW: Track source
    "scifinder_type": "Buchwald-Hartwig",    // ← NEW: Keep original name
    "features": {...}
  }
}
```

### 4. Statistics Reporting

**Enhanced output** showing family detection sources:

```
============================================================
PREPROCESSING STATISTICS
============================================================
Total reactions:           1000
Successfully preprocessed: 985 (98.5%)
Failed:                    10 (1.0%)
Skipped:                   5 (0.5%)

Family Detection Sources:
  SciFinder exact match:   920 (93.4%)  ← Most reliable!
  SciFinder partial match: 45 (4.6%)
  SMARTS detection:        15 (1.5%)
  Unknown:                 5 (0.5%)
============================================================

✓ Dataset saved with precomputed normalization, family detection, and features!
✓ Family names prioritized from SciFinder metadata (more reliable)
✓ This will significantly speed up precedent searches at runtime.
```

## Benefits

### 1. Higher Accuracy
- ✅ **SciFinder metadata is curated** - expert chemists assigned reaction types
- ✅ **SMARTS heuristics are fallible** - may misclassify edge cases
- ✅ **Consistency** - same naming as SciFinder literature

### 2. Better Coverage
- ✅ **Handles complex cases** - multi-step reactions, unusual conditions
- ✅ **Recognizes variants** - "Buchwald-Hartwig", "Buchwald C-N", etc.
- ✅ **Fallback still works** - SMARTS detection for unlabeled data

### 3. Traceability
- ✅ **`family_source` field** - know where family came from
- ✅ **`scifinder_type` preserved** - original metadata available
- ✅ **Confidence scores** - different levels for exact/partial/detected

### 4. Performance
- ✅ **Skip SMARTS when not needed** - faster preprocessing
- ✅ **One-time cost** - family stored in dataset
- ✅ **No runtime overhead** - already decided during creation

## Example Mappings

### Exact Matches (confidence = 1.0)
| SciFinder Type | Mapped Family | Source |
|----------------|---------------|--------|
| `Buchwald` | `C_N_Coupling_Pd` | scifinder |
| `Ullmann` | `C_N_Coupling_Cu` | scifinder |
| `Suzuki-Miyaura` | `Suzuki_CC` | scifinder |
| `Sonogashira` | `Sonogashira_CC` | scifinder |

### Partial Matches (confidence = 0.8)
| SciFinder Type | Matched Key | Mapped Family | Source |
|----------------|-------------|---------------|--------|
| `Buchwald-Hartwig amination` | `buchwald` | `C_N_Coupling_Pd` | scifinder_partial |
| `Modified Ullmann` | `ullmann` | `C_N_Coupling_Cu` | scifinder_partial |
| `Suzuki cross-coupling` | `suzuki` | `Suzuki_CC` | scifinder_partial |

### SMARTS Fallback (confidence = variable)
| SciFinder Type | Detected Family | Confidence | Source |
|----------------|-----------------|------------|--------|
| `(empty)` | `C_N_Coupling_Pd` | 0.90 | smarts_detection |
| `Unclassified` | `Suzuki_CC` | 0.85 | smarts_detection |
| `Custom` | `Unknown` | 0.30 | smarts_detection |

## Edge Cases Handled

### Case 1: SciFinder type with different spelling
```
SciFinder: "Buchwald Hartwig"  (space, no hyphen)
Partial match: "buchwald" → C_N_Coupling_Pd ✓
```

### Case 2: Descriptive SciFinder type
```
SciFinder: "Pd-catalyzed Buchwald aryl amination"
Partial match: "buchwald" → C_N_Coupling_Pd ✓
```

### Case 3: No SciFinder type metadata
```
SciFinder: "" (empty)
Fallback: SMARTS detection → C_N_Coupling_Pd (confidence 0.90) ✓
```

### Case 4: Unknown SciFinder type
```
SciFinder: "Novel coupling method"
Partial match: (none)
Fallback: SMARTS detection → Suzuki_CC (confidence 0.85) ✓
```

## Backward Compatibility

### Old datasets (no precomputed field)
```python
# Runtime detection in precedent.py still works
features = precomputed.get("features")
if not features:
    # Fallback: compute on-the-fly
    features = feat_molecular.featurize(elec, nuc)
```

### New datasets (with precomputed + family_source)
```python
# Can check source for quality assessment
if precomputed.get("family_source") == "scifinder":
    # High confidence - use directly
elif precomputed.get("family_source") == "smarts_detection":
    # May want to review or add user override
```

## Adding New SciFinder Mappings

To support new reaction types from SciFinder:

1. **Identify SciFinder name**
   ```
   Look at RDF files: <ReactionType>Heck</ReactionType>
   ```

2. **Add to mapping**
   ```python
   scifinder_map = {
       ...
       "heck": "Heck_CC",
       "heck-mizoroki": "Heck_CC",
       "heck coupling": "Heck_CC",
   }
   ```

3. **Regenerate datasets**
   ```bash
   python data-processor/Scifinder_rdf_processer.py
   ```

4. **Verify statistics**
   ```
   Family Detection Sources:
     SciFinder exact match:   950 (96.5%)  ← Should increase!
   ```

## Migration Path

### For New RDF Imports
✅ Automatic - mapping applies immediately

### For Existing Datasets
1. Keep original RDF files
2. Re-run processor with updated mapping
3. Old datasets continue working (fallback)

### For Custom Reaction Types
1. Add to `scifinder_map` dictionary
2. Or rely on SMARTS detection fallback
3. Can manually edit JSONL if needed

---

**Date**: 2025-10-07  
**Impact**: Higher accuracy family detection from SciFinder metadata  
**Breaking Changes**: None (additive enhancement)
