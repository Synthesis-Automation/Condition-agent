# Analysis Module Integration Summary

## Overview

Successfully integrated the **new taxonomy/analysis module** (`chemtools.analysis.reaction_context`) into the recommendation pipeline (`chemtools.recommend.modules.recommender`).

## Integration Details

### Changes Made

1. **`chemtools/recommend/modules/recommender.py`**:
   - Added import for `classify_reactants_with_context` from analysis module
   - Modified reaction type detection to use Three-Tier Priority System:
     - **Priority 1**: New analysis module (Two-Pass Approach with SMARTS matching)
     - **Priority 2**: rxn-insight ML model
     - **Priority 3**: Rule-based detection (fallback)
   - Added `analysis_module_used` flag to detection metadata
   - Added `reactant_classification` to detection metadata with full Two-Pass results

2. **`chemtools/recommend/modules/structured.py`**:
   - Enhanced detection section to include reactant classification info
   - Added `analysis_module_used` flag to API output

3. **Configuration**:
   - Analysis module enabled by default
   - Can be disabled via `relax={"use_analysis_module": False}`

### Detection Priority

```
User Override → Analysis Module → rxn-insight ML → Rule-based → Unknown
     ↓                ↓                ↓               ↓
  (1.0 conf)     (ML conf)        (ML conf)        (low conf)
```

### Family Mapping

The analysis module returns specific reaction taxonomy IDs (e.g., `ullmann_cn`, `suzuki_miyaura`), which are automatically mapped to precedent database families using `canonical_family()`:

- `ullmann_cn` → `C_N_Coupling`
- `buchwald_hartwig_c_n` → `C_N_Coupling`
- `suzuki_miyaura` → `Suzuki`
- `chan_lam` → `C_N_Coupling`
- `ullmann_ether` → `C_O_Coupling`
- etc.

## API Response Format

The recommendation API now includes enhanced detection metadata:

```json
{
  "detection": {
    "detected_reaction_type": "C_N_Coupling",
    "confidence": 0.85,
    "source": "auto",
    "analysis_module_used": true,
    "reactant_classification": {
      "reaction_type": "ullmann_cn",
      "confidence": 0.85,
      "detection_method": "ml_detected",
      "num_reactants": 2,
      "has_multi_functional_substrates": false,
      "reactants": [
        {
          "position": 0,
          "category": "ArX*",
          "member_type": "ArBr",
          "name": "aryl bromide",
          "role": "electrophile",
          "is_expected": true,
          "has_alternatives": false
        },
        {
          "position": 1,
          "category": "Aniline-type",
          "member_type": "ArNH2",
          "name": "aniline (primary)",
          "role": null,
          "is_expected": false,
          "has_alternatives": false
        }
      ]
    }
  }
}
```

## Benefits

### 1. Improved Reaction Type Detection
- **Context-Aware**: Uses Two-Pass Approach to classify reactants based on detected reaction type
- **Higher Accuracy**: Combines ML detection (rxn-insight) with SMARTS pattern matching
- **Multi-Functional Substrates**: Correctly identifies primary functional group based on reaction context

### 2. Enhanced Metadata
- **Reactant Roles**: Identifies electrophile, nucleophile, coupling_partner, substrate
- **Functional Groups**: Detailed classification (ArBr, ArB(OH)2, ArNH2, etc.)
- **Alternative Matches**: Tracks other potential functional groups in multi-functional substrates
- **Confidence Scores**: ML-based confidence for both reaction type and reactant classification

### 3. Backward Compatibility
- Existing code continues to work unchanged
- Analysis module can be disabled if needed
- Falls back to rxn-insight/rule-based detection if analysis module fails

## Usage Examples

### Basic Usage (Analysis Module Enabled by Default)

```python
from chemtools import chem

result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    k=20
)

# Check if analysis module was used
if result['detection']['analysis_module_used']:
    print("Analysis module detected:", 
          result['detection']['reactant_classification']['reaction_type'])
```

### Disable Analysis Module

```python
result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    k=20,
    relax={"use_analysis_module": False}  # Use old rxn-insight/rule-based only
)
```

### Access Reactant Classification

```python
reactants = result['detection']['reactant_classification']['reactants']
for r in reactants:
    print(f"Reactant {r['position']}: {r['member_type']} ({r['category']})")
    if r['role']:
        print(f"  Role: {r['role']}")
    print(f"  Expected: {r['is_expected']}")
```

## Testing

Run integration tests:

```bash
python test_recommendation_integration.py
```

Test output shows:
- ✅ Analysis module successfully integrated
- ✅ Reaction type detection working (ullmann_cn → C_N_Coupling)
- ✅ Reactant classification working (ArBr as electrophile, ArNH2 as aniline)
- ✅ Backward compatibility maintained

## Files Modified

1. `chemtools/recommend/modules/recommender.py` - Core integration
2. `chemtools/recommend/modules/structured.py` - API output enhancement
3. `test_recommendation_integration.py` - Integration tests (NEW)

## Related Modules

- `chemtools/analysis/reaction_context.py` - Two-Pass Approach implementation
- `chemtools/analysis/_registry.py` - Taxonomy registry access
- `chemtools/analysis/reactants.py` - Reactant type detection
- `chemtools/taxonomy/data/` - Taxonomy data files

## Next Steps (Optional)

- [ ] Add confidence calibration for reactant classification
- [ ] Use reactant roles for smarter precedent filtering
- [ ] Add reactant classification to formatted output display
- [ ] Create visualization of detected functional groups
- [ ] Add batch processing support

## Conclusion

The new analysis module is now fully integrated into the recommendation pipeline, providing enhanced reaction type detection and detailed reactant classification while maintaining full backward compatibility.
