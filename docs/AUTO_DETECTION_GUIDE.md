# Auto-Detection Implementation Guide

## Overview

The Condition Agent now supports **automatic reaction type detection** when users select "Auto-detect" as the reaction type. This feature uses a two-tier detection system:

1. **Primary**: `rxn-insight` ML-based classifier (if available)
2. **Fallback**: Rule-based functional group detection via `chemtools.router`

## Implementation Summary

### Files Modified

1. **`app/ui_simple.py`** (Lines 35-48, 323-476, 531-571, 793-825)
   - Added `detect_reaction_type` import from `chemtools.reaction_type_detector`
   - Created `detect_and_map_reaction_type()` helper function
   - Updated `get_ml_recommendations()` to use auto-detection
   - Updated `get_rule_recommendations()` to use auto-detection

### Key Features

✅ **Multi-method detection**: rxn-insight → router fallback  
✅ **Intelligent mapping**: Maps detected families to both ML and rule-based naming conventions  
✅ **Catalyst awareness**: Distinguishes Pd/Cu/Ni-catalyzed C-N couplings  
✅ **Confidence scoring**: Returns detection confidence for transparency  
✅ **Graceful degradation**: Works without rxn-insight installed  

## How It Works

### Detection Flow

```
User Input: SMILES + "Auto-detect"
         ↓
detect_and_map_reaction_type()
         ↓
Try rxn-insight detector
  ├─ Success → Map to ML/Rule families
  └─ Fail → Try router.detect_family()
         ├─ Success → Map to ML/Rule families
         └─ Fail → Return error message
```

### Family Name Mapping

The system handles three different naming conventions:

| Internal Family | ML API Name | Rule Database Name |
|----------------|-------------|-------------------|
| `Buchwald_CN` | `C_N_Coupling_Pd` | `C-N Coupling (Pd)` |
| `Ullmann_CN` | `C_N_Coupling_Cu` | `C-N Coupling (Cu)` |
| `C_N_Coupling_Ni` | `C_N_Coupling_Ni` | `C-N Coupling (Ni)` |
| `Suzuki_CC` | `Suzuki_CC` | `Suzuki Coupling` |
| `Amide_Coupling` | `Amide_Coupling` | `Amide Formation` |

### Catalyst Detection

The system examines:
- **SMILES agents** (middle part of `reactants>agents>products`)
- **Functional groups** in reactants
- **rxn-insight classification** output

Example:
```python
# Cu-catalyzed (Ullmann)
"Brc1ccccc1.Nc1ccccc1>Cu>c1ccc(Nc2ccccc2)cc1"
→ Detected: Ullmann_CN (Cu-catalyzed)

# Pd-catalyzed (Buchwald-Hartwig)
"Brc1ccccc1.Nc1ccccc1>Pd>c1ccc(Nc2ccccc2)cc1"
→ Detected: Buchwald_CN (Pd-catalyzed)

# Without catalyst (relies on rxn-insight)
"Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
→ Detected: Heteroatom Alkylation and Arylation
→ Mapped: Ullmann_CN (default for C-N without catalyst info)
```

## Usage Examples

### Example 1: ML Recommendations with Auto-Detection

```python
from app.ui_simple import get_ml_recommendations

reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
json_output, table = get_ml_recommendations(
    reaction_smiles=reaction,
    reaction_type="Auto-detect",  # ← Auto-detect enabled
    top_k=3
)

# Output includes detection info:
# **Auto-detected (rxn-insight):** Suzuki_CC
#   Class: C-C Coupling | Name: Suzuki coupling with boronic acids
#   Confidence: 80.00%
```

### Example 2: Rule-Based Recommendations with Auto-Detection

```python
from app.ui_simple import get_rule_recommendations

reaction = "Brc1ccccc1.Nc1ccccc1>Pd>c1ccc(Nc2ccccc2)cc1"
json_output, table = get_rule_recommendations(
    reaction_smiles=reaction,
    reaction_type="Auto-detect"  # ← Auto-detect enabled
)

# Output includes detection info:
# **Auto-detected (rxn-insight):** Buchwald_CN (Pd-catalyzed)
#   Class: Heteroatom Alkylation and Arylation
#   Catalysts: Pd
```

### Example 3: Direct Detection Testing

```python
from app.ui_simple import detect_and_map_reaction_type

result = detect_and_map_reaction_type(
    reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    requested_type="Auto-detect"
)

print(result)
# {
#     'detected_family': 'Suzuki_CC',
#     'ml_family': 'Suzuki_CC',
#     'rule_db_name': 'Suzuki Coupling',
#     'confidence': 0.8,
#     'method': 'rxn_insight',
#     'success': True,
#     'message': '**Auto-detected (rxn-insight):** Suzuki_CC...'
# }
```

## Detection Methods

### Method 1: rxn-insight (Primary)

**Advantages:**
- ML-based, more accurate
- Handles complex reaction patterns
- Provides confidence scores
- Extracts catalyst information

**Requirements:**
- `rxn_insight` package installed
- RDKit available

**Detection Process:**
1. Call `chemtools.reaction_type_detector.detect_reaction_type()`
2. Extract `rxn_class` (e.g., "Heteroatom Alkylation and Arylation")
3. Extract `rxn_name` (e.g., "Suzuki coupling with boronic acids")
4. Map to chemtools family via `_map_to_family()`
5. Refine based on catalysts via `_refine_cn_family()`

### Method 2: Router-Based (Fallback)

**Advantages:**
- Always available (no external dependencies)
- Fast functional group matching
- Explicit rules

**Limitations:**
- Less accurate than ML
- May miss complex patterns
- Lower confidence scores

**Detection Process:**
1. Call `chemtools.router.detect_family()`
2. Check functional group hits (aryl halide, nucleophile, etc.)
3. Apply confidence thresholds
4. Map to family names

## Error Handling

The system provides helpful error messages when detection fails:

### No ML Model Available
```
**Auto-detected (rxn-insight):** Custom_Reaction

**No ML model available** for detected family: Custom_Reaction

Available ML families:
  - C-N Coupling (Cu) (C_N_Coupling_Cu)
  - C-N Coupling (Pd) (C_N_Coupling_Pd)
  - Suzuki Coupling (Suzuki_CC)
  ...
```

### No Rule Database Available
```
**Auto-detected (rxn-insight):** Custom_Reaction

**No rule database available** for detected family: Custom_Reaction

Available rule databases:
  - C-N Coupling (Cu)
  - C-N Coupling (Pd)
  - Suzuki Coupling
  ...
```

### Detection Failed
```
**Auto-detection failed:** [Error details]
```

## Configuration

### Adding New Reaction Types

To add support for a new reaction type:

1. **Add to ML_FAMILY_MAP** (in `app/ui_simple.py`):
```python
ML_FAMILY_MAP = {
    "Auto-detect": None,
    "New Reaction": "New_Reaction_Family",  # ← Add here
    ...
}
```

2. **Add to RULE_DATABASES** (if you have a rule DB):
```python
RULE_DATABASES = {
    "New Reaction": str(SCDB_DIR / "new_reaction_db.json"),  # ← Add here
    ...
}
```

3. **Add mapping in `detect_and_map_reaction_type()`**:
```python
# In the rxn-insight success branch:
elif mapped_family == "New_Reaction_Family":
    ml_family = "New_Reaction_Family"
    rule_db_name = "New Reaction"
    detected_family = "New_Reaction_Family"
```

4. **Add to router fallback mappings**:
```python
family_to_ml = {
    "New_Reaction_Family": "New_Reaction_Family",  # ← Add here
    ...
}

family_to_db = {
    "New_Reaction_Family": "New Reaction",  # ← Add here
    ...
}
```

## Testing

### Test Script

Run the comprehensive test:
```powershell
$env:PYTHONIOENCODING='utf-8'
python scripts/test_auto_detection.py
```

### Manual Testing

1. **Start the UI**:
```powershell
python app/ui_simple.py
```

2. **Enter a reaction SMILES**:
```
Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
```

3. **Select "Auto-detect" from dropdown**

4. **Click "Get ML Recommendations" or "Get Rule Recommendations"**

5. **Verify detection message appears**:
```
**Auto-detected (rxn-insight):** Suzuki_CC
  Class: C-C Coupling | Name: Suzuki coupling with boronic acids
  Confidence: 80.00%
```

## Performance Notes

### Detection Speed

- **rxn-insight**: ~0.5-2 seconds (first call), ~0.1-0.5s (cached)
- **Router**: <0.1 seconds

The detection overhead is negligible compared to recommendation generation (70-110s for ML).

### Caching

- rxn-insight caches models after first load
- Router has no caching (fast enough not to need it)

## Troubleshooting

### rxn-insight Not Available

**Symptom**: Always falls back to router detection

**Solution**: Install rxn-insight:
```bash
pip install rxn-insight
```

### Low Confidence Detections

**Symptom**: Confidence < 50%

**Possible Causes:**
1. Reaction type not in training data
2. Unusual functional groups
3. Incomplete SMILES

**Solution**: Manually select reaction type instead of auto-detect

### Wrong Catalyst Detection

**Symptom**: Pd-catalyzed detected as Cu or vice versa

**Cause**: Catalyst not specified in SMILES agents section

**Solution**: Use proper SMILES format:
```
# Correct (with catalyst)
Brc1ccccc1.Nc1ccccc1>Pd>product

# Without catalyst (may default to wrong type)
Brc1ccccc1.Nc1ccccc1>>product
```

## API Reference

### `detect_and_map_reaction_type()`

```python
def detect_and_map_reaction_type(
    reaction_smiles: str,
    requested_type: str
) -> Dict[str, Any]:
    """
    Detect reaction type and map to both ML and rule-based naming conventions.
    
    Args:
        reaction_smiles: Reaction SMILES string
        requested_type: User-requested type or "Auto-detect"
        
    Returns:
        Dict with:
            - detected_family: str - Internal chemtools family
            - ml_family: str | None - ML API family name
            - rule_db_name: str | None - Rule database name
            - confidence: float - Detection confidence (0-1)
            - method: str - Detection method used
            - success: bool - Whether detection succeeded
            - message: str - Human-readable message
            - raw_detection: dict - Raw detection result
    """
```

## Future Improvements

### Planned Enhancements

1. **Confidence Thresholds**: Warn users when confidence < 60%
2. **Multi-Family Support**: Handle reactions matching multiple families
3. **User Feedback Loop**: Learn from user corrections
4. **Batch Detection**: Optimize for multiple reactions
5. **Extended Catalyst Support**: Detect Ni, Cu(I) vs Cu(II), ligands

### Known Limitations

- Cannot distinguish between similar C-N coupling protocols without catalyst info
- May misclassify unusual or hybrid reaction types
- Confidence scores from router are heuristic-based
- No support for reaction types not in training data

## Summary

The auto-detection system provides:

✅ **Seamless Integration**: Works with both ML and rule-based systems  
✅ **High Accuracy**: ML-based primary detection  
✅ **Robustness**: Fallback to rule-based detection  
✅ **User-Friendly**: Clear detection messages and confidence scores  
✅ **Flexible**: Easy to add new reaction types  

The system significantly improves user experience by eliminating the need to manually select reaction types while maintaining high accuracy through the two-tier detection approach.
