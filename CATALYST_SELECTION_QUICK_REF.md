# Catalyst Selection - Quick Reference

## ✅ Feature Added: Interactive Catalyst Selection for C-N Coupling

### What Changed?

Users can now select their preferred catalyst when running C-N coupling recommendations:
- **Option 1**: No preference (any catalyst) - default
- **Option 2**: Palladium (Pd) - Buchwald-Hartwig conditions
- **Option 3**: Copper (Cu) - Ullmann-type conditions  
- **Option 4**: Nickel (Ni) - Ni-catalyzed conditions

### How to Use

#### Interactive Mode
```bash
python scripts/local_recommendation_cli.py

# Prompts will appear:
Reaction Type Options:
  1) Auto-detect (server decides) (default)
  2) Suzuki Coupling
  3) C–N Coupling (unified)
  4) Amide Formation
Select reaction type [1]: 3

Catalyst Preference:
  1) No preference (any catalyst) (default)
  2) Palladium (Pd)
  3) Copper (Cu)
  4) Nickel (Ni)
Select catalyst preference [1]: 3
```

#### Command-Line Mode
```bash
# Specify catalyst directly
python scripts/local_recommendation_cli.py \
  --family C_N_Coupling \
  --catalyst Cu \
  --strategy ml

# Options: --catalyst Pd | Cu | Ni
# Omit --catalyst for "no preference" (all catalysts)
```

#### Programmatic Mode
```python
from scripts.local_recommendation_cli import local_ml_recommendation

result = local_ml_recommendation(
    reaction="Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1",
    reaction_type="C_N_Coupling",
    k_value=50,
    limit=5,
    catalyst_preference="Cu"  # or "Pd", "Ni", None
)
```

### When Does It Prompt?

✅ **Prompts for catalyst** when:
- Reaction type is C-N Coupling (unified)
- `--catalyst` argument NOT provided
- Running in interactive mode

❌ **Does NOT prompt** when:
- Using other reaction types (Suzuki, Amide)
- `--catalyst` already specified
- `--family` explicitly set and `--catalyst` omitted (defaults to None)

### How It Works

1. User selects C-N Coupling → System prompts for catalyst
2. Catalyst preference passed to ML recommendation
3. Uses `relax={"catalyst_class": "Cu"}` parameter
4. Precedent search filters by catalyst class
5. Returns only reactions using specified catalyst

### Files Modified

- ✅ `scripts/recommendation_cli_utils.py` - Added `CATALYST_CHOICES` and `choose_catalyst()`
- ✅ `scripts/local_recommendation_cli.py` - Added `--catalyst` arg, catalyst prompt, and filtering

### Example Output

```
Enter reaction SMILES: Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1
Select reaction type [1]: 3
Selected reaction type: C–N Coupling (unified)
Select catalyst preference [1]: 3
Catalyst preference: Copper (Cu)

Running local pipelines...

ML detected type: C_N_Coupling (confidence 0.85)
Top ML recommendation:
  Catalyst: CuI
  Ligand: L-Proline
  Base: K2CO3
  Solvent: DMSO
  Temperature: 90°C
  Precedent count: 15
```

### Verification

```python
# Test catalyst choices are available
from scripts.recommendation_cli_utils import CATALYST_CHOICES
print(len(CATALYST_CHOICES))  # Should be 4

# Test ML function accepts catalyst_preference
from scripts.local_recommendation_cli import local_ml_recommendation
import inspect
sig = inspect.signature(local_ml_recommendation)
print('catalyst_preference' in sig.parameters)  # Should be True
```

---

**Status**: ✅ IMPLEMENTED  
**Date**: 2025-10-16  
**Feature**: Catalyst selection for C-N coupling  
**Usage**: Interactive prompt or `--catalyst Pd|Cu|Ni` argument
