# Generalized Coupling Confirmation System

## Summary of Changes

The reaction detection system has been **generalized** to support all major coupling reactions, removing Suzuki-specific parameters and making the API more consistent.

## What Changed

### Before (Suzuki-specific):
```python
# Old API - Suzuki-specific parameter
featurize_reaction(
    smiles,
    options={"confirm_suzuki_products": True}  # Only for Suzuki!
)
```

### After (General):
```python
# New API - works for 9+ coupling reaction types
featurize_reaction(
    smiles,
    options={"confirm_coupling_products": True}  # Default is now True
)
```

## Supported Coupling Reactions

The general coupling confirmation system now validates:

1. **Suzuki-Miyaura** (Ar-B(OH)₂ + Ar-X → Ar-Ar)
2. **Negishi** (Ar-ZnX + Ar-X → Ar-Ar)
3. **Stille** (Ar-SnR₃ + Ar-X → Ar-Ar)
4. **Kumada** (Ar-MgX + Ar-X → Ar-Ar)
5. **Hiyama** (Ar-SiR₃ + Ar-X → Ar-Ar)
6. **Sonogashira** (Ar-X + R-C≡C-H → Ar-C≡C-R)
7. **Buchwald-Hartwig (C-N)** (Ar-X + R₂NH → Ar-NR₂)
8. **C-O Coupling** (Ar-X + R-OH → Ar-OR)
9. **C-S Coupling** (Ar-X + R-SH → Ar-SR)

## How It Works

### 1. Motif-Based Detection
Detects reactant and product motifs (e.g., `Ar-B(OH)2`, `Ar-I`, `Ar-Ar`)

### 2. Reacted Motifs Validation
Validates using consumption patterns:
- **Reacted motifs**: Consumed in reaction (`Ar-B(OH)2`, `Ar-I`)
- **Formed motifs**: Created in product (`Ar-Ar`)
- **Pattern matching**: Organoboron + aryl halide → biaryl = Suzuki

### 3. Product Confirmation (Optional but Default)
Uses bond formation analysis to confirm coupling products by checking:
- Nucleophile attachment points
- Electrophile attachment points
- Product structure validation

## API Changes

### `featurize_reaction()`

**Old Parameters (DEPRECATED):**
```python
confirm_suzuki_products: bool = False  # Suzuki-specific ❌
```

**New Parameters:**
```python
confirm_coupling_products: bool = True  # General, supports 9+ reactions ✓
```

### `detect_reaction_types()`

**Old Signature:**
```python
def detect_reaction_types(
    reaction_smiles: str,
    confirm_suzuki_products: bool = False,      # ❌ Removed
    confirm_coupling_products: bool = False,
)
```

**New Signature:**
```python
def detect_reaction_types(
    reaction_smiles: str,
    confirm_coupling_products: bool = True,     # ✓ Now default True
)
```

## Backward Compatibility

The old `confirm_suzuki_products` parameter still works for backward compatibility:

```python
# Old code still works - maps to general parameter
featurize_reaction(
    smiles,
    options={"confirm_suzuki_products": True}  # ✓ Still works
)

# But new code is preferred
featurize_reaction(
    smiles,
    options={"confirm_coupling_products": True}  # ✓ Recommended
)
```

## Examples

### Example 1: Suzuki Coupling (auto-detected)
```python
from chemtools.featurizers.unified import featurize_reaction

result = featurize_reaction(
    "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1"
)

print(result["reaction_type"])  # "Suzuki_miyaura"
print(result["confidence"])     # 0.95
```

### Example 2: Buchwald-Hartwig C-N Coupling
```python
result = featurize_reaction(
    "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
)

print(result["reaction_type"])  # "C_N_Coupling"
print(result["confidence"])     # 0.95
```

### Example 3: Disable Coupling Confirmation
```python
# For faster detection without product validation
result = featurize_reaction(
    smiles,
    options={"confirm_coupling_products": False}
)
```

## Benefits

1. **Consistency**: Single parameter for all coupling reactions
2. **Extensibility**: Easy to add new coupling types
3. **Better Defaults**: Coupling confirmation now enabled by default
4. **Cleaner API**: Removed reaction-specific parameters
5. **Better Detection**: Multi-pass validation catches more edge cases

## Migration Guide

### If you were using `confirm_suzuki_products`:

**Before:**
```python
result = featurize_reaction(
    smiles,
    options={
        "confirm_suzuki_products": True,
        "detailed": True
    }
)
```

**After:**
```python
result = featurize_reaction(
    smiles,
    options={
        "confirm_coupling_products": True,  # Now default True
        "detailed": True
    }
)

# Or just rely on defaults:
result = featurize_reaction(smiles, options={"detailed": True})
```

### If you have custom detection logic:

The system now has two validation layers:

1. **Slot-based detection** (checks reactant motifs)
2. **Reacted motifs validation** (checks consumption patterns)
3. **Product confirmation** (validates bond formation) - Optional, default True

You can disable product confirmation for speed:
```python
result = featurize_reaction(
    smiles,
    options={"confirm_coupling_products": False}
)
```

## Technical Details

### Validation Priority

1. **Bond change analysis** (if enabled via `use_bond_changes=True`)
2. **Slot-based motif matching** (default)
3. **Reacted motifs pattern validation** (always applied)
4. **Product confirmation** (applied if `confirm_coupling_products=True`)

### Confidence Scores

- **Slot-based match**: confidence = matched_slots / required_slots
- **Validated by pattern**: confidence = 0.95 (high confidence correction)
- **Product confirmed**: adds metadata but doesn't change confidence

## Files Changed

1. `chemtools/featurizers/formatters/reaction.py`
   - Removed `confirm_suzuki_products` parameter
   - Added backward compatibility mapping
   - Changed default to `confirm_coupling_products=True`

2. `chemtools/featurizers/detection/core.py`
   - Removed `confirm_suzuki_products` from function signatures
   - Simplified coupling confirmation logic
   - Updated documentation

3. `chemtools/featurizers/formatters/detection_validation.py`
   - General validation patterns (already implemented)
   - Supports 9+ coupling reaction types

## Testing

Run the test suite to verify:
```bash
python test_general_coupling_detection.py
```

Expected output:
```
✓ PASS: Suzuki_miyaura
✓ PASS: C_N_Coupling
✓ PASS: Other coupling types
```
