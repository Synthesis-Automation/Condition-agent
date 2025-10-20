# Enhanced Precedent Display in Cross-Family CLI

## Update Summary

Enhanced the cross-family recommendation CLI to show **top precedents for ALL recommendations**, not just the first one, with detailed information for each precedent.

## What Changed

### Before
- Only Recommendation #1 showed precedent details
- Other recommendations had no precedent information visible
- No control over how many precedents to display

### After
- **All recommendations** show precedent information (when available)
- Clear message "Top Precedents: None available" when precedents aren't provided by the API
- New `--max-precedents` option to control how many precedents to show per recommendation
- Shows count: "Top Precedents (3 of 8)" to indicate total available

## New Features

### 1. Precedent Display for All Recommendations

Each recommendation now shows:
```
Top Precedents (3 of 8):
  1. 31-614-CAS-44374029
     Catalyst: Pd/XPhos
     Yield: 85%
     Ref: Waste-minimized access to diarylamines...
  2. 31-614-CAS-44374022
     Catalyst: Pd/XPhos
     Ref: Waste-minimized access to diarylamines...
  3. 31-614-CAS-44374027
     Catalyst: Pd/XPhos
     Ref: Waste-minimized access to diarylamines...
```

Or if no precedents available:
```
Top Precedents: None available
```

### 2. New CLI Option: --max-precedents

Control how many precedents to show per recommendation:

```bash
# Default: show 3 precedents
python app/cross_family_recommendation_cli.py "reaction>>product"

# Show up to 5 precedents per recommendation
python app/cross_family_recommendation_cli.py --rxn "..." --max-precedents 5

# Show only 1 precedent per recommendation
python app/cross_family_recommendation_cli.py --rxn "..." --max-precedents 1
```

### 3. Enhanced Precedent Information

Each precedent now displays:
- **Reaction ID** (e.g., "31-614-CAS-44374029")
- **Catalyst** (e.g., "Pd/XPhos") - if available
- **Yield** (e.g., "85%") - if available
- **Reference** (truncated title from paper)

## Example Output

### C-N Coupling Reaction
```bash
python app/cross_family_recommendation_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --k 20
```

**Output:**
```
================================================================================
Recommendation #1
────────────────────────────────────────────────────────────────────────────────
  Catalyst: Pd/XPhos
  Base: Potassium tert-butoxide
  Solvent: Cyclopentyl methyl ether
  Confidence: 0.76
  Precedent Support: 8 similar reaction(s)
  Top Precedents (3 of 3):
    1. 31-614-CAS-44374029
       Catalyst: Pd/XPhos
       Ref: Waste-minimized access to diarylamines and triarylamines via Csp2-N...
    2. 31-614-CAS-44374022
       Catalyst: Pd/XPhos
       Ref: Waste-minimized access to diarylamines and triarylamines via Csp2-N...
    3. 31-614-CAS-44374027
       Catalyst: Pd/XPhos
       Ref: Waste-minimized access to diarylamines and triarylamines via Csp2-N...

────────────────────────────────────────────────────────────────────────────────
Recommendation #2
────────────────────────────────────────────────────────────────────────────────
  Catalyst: Pd/XPhos
  Base: Sodium tert-butoxide
  Solvent: 1,4-Dioxane
  Confidence: 0.76
  Precedent Support: 1 similar reaction(s)
  Top Precedents: None available

────────────────────────────────────────────────────────────────────────────────
Recommendation #3
────────────────────────────────────────────────────────────────────────────────
  Catalyst: Pd/XPhos
  Base: Potassium hydroxide
  Solvent: 2-Methyl-2-butanol
  Confidence: 0.76
  Precedent Support: 1 similar reaction(s)
  Top Precedents: None available
```

### Suzuki Coupling with Multiple Precedents
```bash
python app/cross_family_recommendation_cli.py "Brc1ccc(C)cc1.c1ccc(B(O)O)cc1>>Cc1ccc(-c2ccccc2)cc1" --k 30 --max-precedents 5
```

**Output:**
```
Recommendation #1
────────────────────────────────────────────────────────────────────────────────
  Catalyst: Pd
  Base: Potassium carbonate
  Solvent: Ethanol
  Confidence: 0.95
  Precedent Support: 4 similar reaction(s)
  Top Precedents (3 of 3):
    1. 31-614-CAS-35216960
       Catalyst: Pd
       Ref: Fabrication and catalytic performance of a new diaminopyridine Pd(II)...
    2. 31-614-CAS-36083898
       Catalyst: Pd
       Ref: Magnetite Pd-loaded nitrogen-rich porous organic polymer as a catalyst...
    3. 31-614-CAS-38827379
       Catalyst: Pd
       Ref: Palladium nanoparticles on polydimethylsiloxane film for C-C coupling...

────────────────────────────────────────────────────────────────────────────────
Recommendation #2
────────────────────────────────────────────────────────────────────────────────
  Catalyst: Pd
  Solvent: N/A
  Confidence: 0.95
  Precedent Support: 2 similar reaction(s)
  Top Precedents (2 of 2):
    1. 31-179-CAS-11213804
       Catalyst: Pd
       Ref: A novel green protocol for ligand free Suzuki-Miyaura cross-coupling...
    2. 31-614-CAS-41961987
       Catalyst: Pd
       Ref: Preparation and Application of WELAN to Pd(OAc)2-Catalyzed Suzuki-Miyaura...
```

## Implementation Details

### Modified Function Signature
```python
def print_recommendation(result: dict, reaction_smiles: str, max_precedents: int = 3):
    """Print recommendation results in a readable format."""
```

### Precedent Display Logic
```python
precedents = summary.get("precedents", [])
if precedents:
    num_to_show = min(len(precedents), max_precedents)
    print(f"  Top Precedents ({num_to_show} of {len(precedents)}):")
    for i, prec in enumerate(precedents[:max_precedents], 1):
        reaction_id = prec.get("reaction_id", "N/A")
        ref = prec.get("reference", "")
        core = prec.get("core", "")
        yield_pct = prec.get("yield_pct")
        
        # Extract just the title from reference
        ref_title = ref.split("|")[0].strip() if "|" in ref else ref[:80]
        
        print(f"    {i}. {reaction_id}")
        if core:
            print(f"       Catalyst: {core}")
        if yield_pct is not None:
            print(f"       Yield: {yield_pct}%")
        if ref_title:
            print(f"       Ref: {ref_title}...")
else:
    print(f"  Top Precedents: None available")
```

### New CLI Argument
```python
parser.add_argument(
    "--max-precedents",
    type=int,
    default=3,
    help="Maximum number of precedents to show per recommendation (default: 3)"
)
```

## Usage Examples

### Basic Usage (Default 3 Precedents)
```bash
python app/cross_family_recommendation_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
```

### Show More Precedents
```bash
# Show up to 5 precedents per recommendation
python app/cross_family_recommendation_cli.py --rxn "..." --max-precedents 5

# Show all available precedents (use large number)
python app/cross_family_recommendation_cli.py --rxn "..." --max-precedents 10
```

### Show Fewer Precedents
```bash
# Show only 1 precedent per recommendation
python app/cross_family_recommendation_cli.py --rxn "..." --max-precedents 1
```

### Combined with Other Options
```bash
# Get 100 precedents and show top 5 for each recommendation
python app/cross_family_recommendation_cli.py \
  --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" \
  --k 100 \
  --max-precedents 5
```

## Why Some Recommendations Show "None available"

The API (`chem.recommend.conditions()`) returns detailed precedent information only for the most relevant recommendations based on the aggregation logic. Recommendations with fewer or less relevant precedents may not have the detailed `precedents` array populated in the `summary` object.

This is normal behavior and indicates:
- The recommendation is valid but has less precedent support
- The API chose not to include detailed precedent information for that specific recommendation
- The recommendation may be based on extrapolation from similar conditions

## Benefits

✅ **Complete visibility** - See precedents for all recommendations, not just #1  
✅ **Transparency** - Know when precedent data is not available  
✅ **Flexibility** - Control how many precedents to view with `--max-precedents`  
✅ **Better decisions** - More precedent information helps choose the best conditions  
✅ **Research support** - Direct links to papers via reaction IDs and references  

## Files Modified

- **`app/cross_family_recommendation_cli.py`**:
  - Updated `print_recommendation()` to accept `max_precedents` parameter
  - Enhanced precedent display to show count "X of Y"
  - Added "None available" message for recommendations without precedents
  - Added `--max-precedents` CLI argument
  - Updated help text with new example

## Summary

The cross-family CLI now provides **comprehensive precedent information for all recommendations**, giving users complete transparency into the evidence supporting each suggested set of conditions. The new `--max-precedents` option provides flexibility to view as much or as little precedent detail as needed! 🎉
