# Catalyst Selection Feature - Implementation Summary

## Feature Overview

Added interactive catalyst selection for C-N coupling reactions, allowing users to filter recommendations by preferred metal catalyst (Pd, Cu, Ni, or no preference).

## Changes Implemented

### 1. New Catalyst Choices (`recommendation_cli_utils.py`)

Added catalyst selection constants and prompt function:

```python
# Catalyst options for C-N coupling reactions
CATALYST_CHOICES: Tuple[Tuple[str, Optional[str]], ...] = (
    ("No preference (any catalyst)", None),
    ("Palladium (Pd)", "Pd"),
    ("Copper (Cu)", "Cu"),
    ("Nickel (Ni)", "Ni"),
)

def choose_catalyst() -> Tuple[str, Optional[str]]:
    """Allow the user to choose a preferred catalyst for C-N coupling."""
    print("\nCatalyst Preference:")
    for idx, (label, _) in enumerate(CATALYST_CHOICES, start=1):
        default_marker = " (default)" if idx == 1 else ""
        print(f"  {idx}) {label}{default_marker}")

    while True:
        choice = input("Select catalyst preference [1]: ").strip()
        if not choice:
            return CATALYST_CHOICES[0]
        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(CATALYST_CHOICES):
                return CATALYST_CHOICES[idx - 1]
        print(f"Please enter a number between 1 and {len(CATALYST_CHOICES)}.")
```

### 2. Updated ML Recommendation Function (`local_recommendation_cli.py`)

Added `catalyst_preference` parameter:

```python
def local_ml_recommendation(
    reaction: str,
    reaction_type: Optional[str],
    k_value: int,
    limit: int,
    rerank_strategy: str = 'rule',
    filter_unknown_reagents: bool = False,
    catalyst_preference: Optional[str] = None,  # NEW
) -> Dict[str, Any]:
    """Replicate the /api/v1/recommend/conditions endpoint locally.
    
    Args:
        catalyst_preference: Preferred catalyst class (e.g., 'Pd', 'Cu', 'Ni')
    """
    # Build relax constraints for catalyst filtering
    relax = {}
    if catalyst_preference:
        relax["catalyst_class"] = catalyst_preference
    
    try:
        raw_data = chem.recommend.conditions(
            reaction=reaction,
            reaction_type=reaction_type,
            k=k_value,
            limit=limit,
            relax=relax,  # Passes catalyst filter
            constraints={},
            rerank=rerank_strategy,
            filter_unknown_reagents=filter_unknown_reagents
        )
```

### 3. Added Command-Line Argument

```python
parser.add_argument(
    "--catalyst",
    type=str,
    default=None,
    choices=[None, "Pd", "Cu", "Ni"],
    help="Catalyst preference for C-N coupling (Pd, Cu, Ni). If not provided and C-N coupling "
         "is selected, will prompt interactively. Uses relax parameter for filtering."
)
```

### 4. Interactive Prompt Logic

```python
# Prompt for catalyst preference if C-N coupling is selected
catalyst_preference = None
if reaction_type == "C_N_Coupling" and not args.catalyst:
    catalyst_label, catalyst_preference = choose_catalyst()
    if catalyst_preference:
        print(f"Catalyst preference: {catalyst_label}")
    else:
        print("Catalyst preference: No preference (any catalyst)")
elif args.catalyst:
    catalyst_preference = args.catalyst
    print(f"Catalyst preference: {catalyst_preference}")
```

### 5. Updated Function Call

```python
if run_ml:
    ml_result = local_ml_recommendation(
        reaction, 
        reaction_type, 
        k_value, 
        limit_value,
        rerank_strategy=args.rerank,
        filter_unknown_reagents=args.filter_unknown,
        catalyst_preference=catalyst_preference  # NEW
    )
```

## Usage Examples

### Interactive Mode

```bash
python scripts/local_recommendation_cli.py

# User interaction:
Enter reaction SMILES: Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1

Reaction Type Options:
  1) Auto-detect (server decides) (default)
  2) Suzuki Coupling
  3) C–N Coupling (unified)
  4) Amide Formation
Select reaction type [1]: 3
Selected reaction type: C–N Coupling (unified)

Catalyst Preference:
  1) No preference (any catalyst) (default)
  2) Palladium (Pd)
  3) Copper (Cu)
  4) Nickel (Ni)
Select catalyst preference [1]: 3
Catalyst preference: Copper (Cu)
```

### Command-Line Mode

```bash
# Specify catalyst via command-line
python scripts/local_recommendation_cli.py \
  --family C_N_Coupling \
  --catalyst Cu \
  --strategy ml \
  --k 50 \
  --limit 5

# No catalyst preference (all catalysts)
python scripts/local_recommendation_cli.py \
  --family C_N_Coupling \
  --strategy ml
```

### Programmatic Usage

```python
from scripts.local_recommendation_cli import local_ml_recommendation

rsmi = "Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1"

# Get Cu-catalyzed recommendations
result_cu = local_ml_recommendation(
    reaction=rsmi,
    reaction_type="C_N_Coupling",
    k_value=50,
    limit=5,
    rerank_strategy='rule',
    catalyst_preference="Cu"
)

# Get Pd-catalyzed recommendations
result_pd = local_ml_recommendation(
    reaction=rsmi,
    reaction_type="C_N_Coupling",
    k_value=50,
    limit=5,
    rerank_strategy='rule',
    catalyst_preference="Pd"
)

# Get all recommendations (no catalyst filter)
result_all = local_ml_recommendation(
    reaction=rsmi,
    reaction_type="C_N_Coupling",
    k_value=50,
    limit=5,
    rerank_strategy='rule',
    catalyst_preference=None
)
```

## Technical Details

### Catalyst Filtering Mechanism

Uses the existing `relax` parameter in the precedent search:

```python
relax = {"catalyst_class": "Pd"}  # or "Cu", "Ni"
```

This filters the precedent pool to only include reactions using the specified catalyst class before DRFP similarity ranking.

### Catalyst Classification

Catalyst classes are determined by `chemtools/precedent/catalyst.py`:
- **Pd**: Palladium-based catalysts (Buchwald-Hartwig)
- **Cu**: Copper-based catalysts (Ullmann-type)
- **Ni**: Nickel-based catalysts
- **other**: Non-catalyzed or other metal catalysts

### Workflow

1. **User selects C-N coupling** → System prompts for catalyst preference
2. **Catalyst preference stored** → Passed as `catalyst_preference` parameter
3. **relax constraint built** → `{"catalyst_class": "Cu"}` if Cu selected
4. **Precedent search filters** → Only matches specified catalyst class
5. **DRFP ranking** → Top-k most similar reactions returned
6. **Reranking applied** → Rule-based or analytics boost
7. **Recommendations returned** → Filtered by catalyst preference

## Related Documentation

- `CATALYST_FILTERING_GUIDE.md` - Detailed guide on catalyst filtering
- `CATALYST_CLASSIFICATION_FIX.md` - Correct chemistry classification
- `CATALYST_DISTRIBUTION_COMPLETE.md` - Catalyst statistics
- `chemtools/precedent/search.py` - DRFP-based search with relax parameter

## Files Modified

1. ✅ `scripts/recommendation_cli_utils.py`
   - Added `CATALYST_CHOICES` constant
   - Added `choose_catalyst()` function

2. ✅ `scripts/local_recommendation_cli.py`
   - Updated imports to include `choose_catalyst`
   - Added `--catalyst` CLI argument
   - Updated `local_ml_recommendation()` signature
   - Added catalyst prompt logic
   - Updated ML recommendation call

## Benefits

### 1. **User Control**
Users can explicitly request:
- Pd-catalyzed conditions (Buchwald-Hartwig)
- Cu-catalyzed conditions (Ullmann-type)
- Ni-catalyzed conditions
- Any catalyst (no preference)

### 2. **Improved Relevance**
Recommendations match user's catalyst preference:
- Lab has Cu expertise → Get Cu conditions
- Pd-free required → Get Cu/Ni conditions
- Screening study → Get all catalysts

### 3. **Consistent Interface**
Same interaction pattern as reaction type selection:
- Interactive prompts with numbered choices
- Command-line argument for automation
- Programmatic API for scripting

### 4. **Backward Compatible**
- Default: No preference (all catalysts)
- Only prompts for C-N coupling reactions
- Other reaction types unaffected

## Testing

### Manual Test

```bash
cd c:\Git-softwares\Condition-agent
python scripts/local_recommendation_cli.py

# Test inputs:
Reaction: Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1
Type: 3 (C-N Coupling)
Catalyst: 3 (Cu)
```

**Expected**:
- Prompt for catalyst preference appears
- Cu selection is stored
- ML recommendations filtered to Cu-catalyzed only

### Automated Test

```python
from scripts.recommendation_cli_utils import CATALYST_CHOICES

# Verify 4 catalyst options available
assert len(CATALYST_CHOICES) == 4
assert CATALYST_CHOICES[0][1] is None  # No preference
assert CATALYST_CHOICES[1][1] == "Pd"
assert CATALYST_CHOICES[2][1] == "Cu"
assert CATALYST_CHOICES[3][1] == "Ni"
```

## Future Enhancements

### 1. Extend to Other Reaction Types
- Suzuki: Pd vs Ni catalysts
- C-O Coupling: Cu vs Pd catalysts

### 2. Multi-Catalyst Selection
```python
catalyst_preference=["Pd", "Cu"]  # Accept either Pd or Cu
```

### 3. Catalyst Property Filters
```python
catalyst_preference={
    "metal": "Pd",
    "ligand_type": "phosphine"
}
```

### 4. Rule-Based Catalyst Selection
Integrate with rule-based matching for catalyst filtering.

---

**Status**: ✅ **IMPLEMENTED**
**Date**: 2025-10-16
**Feature**: Interactive catalyst selection for C-N coupling reactions
**Method**: Uses existing `relax={"catalyst_class": "Pd/Cu/Ni"}` parameter
