# HTE Agent Integration Update

## Summary

Updated the HTE recommendation tool to support reaction SMILES input format, enabling the agent to answer queries like:

```
"use HTE system, recommend conditions for reactions: Brc1ccco1.Nc1ccccc1>>c1ccccc1Nc1ccco1, 
its a C-N coupling reaction, and the catalyst is copper"
```

## Changes Made

### 1. Updated `HTERecommendInput` Schema

Added `reaction_smiles` parameter:
- Users can now provide either individual reactants OR complete reaction SMILES
- Reactants are automatically extracted from reaction SMILES format

```python
class HTERecommendInput(BaseModel):
    reactant_a_smiles: Optional[str]  # Individual reactant (OR use reaction_smiles)
    reactant_b_smiles: Optional[str]  # Individual reactant (optional)
    reaction_smiles: Optional[str]    # NEW: Complete reaction SMILES
    reaction_type_filter: Optional[str]  # e.g., "C_N_Coupling"
    catalyst_filter: Optional[str]       # e.g., "Cu", "copper"
```

### 2. Updated `hte_recommend_tool` Function

**New Capabilities:**
- Parses reaction SMILES using `_split_reaction_smiles`
- Extracts reactants from `reactants>>products` format
- Handles dot notation for multiple reactants (`Brc1ccco1.Nc1ccccc1`)
- Returns z-scores (primary ranking metric)

**Auto-extraction Logic:**
```python
if reaction_smiles:
    parts = _split_reaction_smiles(reaction_smiles.strip())
    reactants_str = parts[0]
    reactant_list = reactants_str.split(".")
    reactant_a_smiles = reactant_list[0]
    reactant_b_smiles = reactant_list[1] if len(reactant_list) >= 2 else None
```

### 3. Enhanced Error Messages

When no data is found, provides:
- Detected reactant types (e.g., "furan-diene + ArNH2")
- Applied filters (reaction type, catalyst)
- Actionable suggestions:
  - Remove filters to broaden search
  - Use `hte_analytics_tool` to explore available pairs
  - Check if reactant types are common in database

### 4. Added Z-Score to Output

Recommendations now include `avg_z_score` (primary ranking metric):

```python
{
    "rank": 1,
    "avg_z_score": 2.61,  # NEW: Primary ranking metric
    "catalyst": "CuI",
    "ligand": "PPBO",
    "base": "NaOtBu",
    "solvent": "tAmOH",
    "avg_yield": 49.8,
    "num_experiments": 1
}
```

## Example Usage

### Agent Query Format

User: "use HTE system, recommend conditions for reactions: Brc1ccco1.Nc1ccccc1>>c1ccccc1Nc1ccco1, its a C-N coupling reaction, and the catalyst is copper"

### Tool Invocation

```python
hte_recommend_tool(
    reaction_smiles="Brc1ccco1.Nc1ccccc1>>c1ccccc1Nc1ccco1",
    reaction_type_filter="C_N_Coupling",
    catalyst_filter="Cu",
    top_k=10,
    min_experiments=1
)
```

### Expected Response

If data exists:
```json
{
    "success": true,
    "data": {
        "reactant_a_type": "ArBr",
        "reactant_b_type": "ArNH2",
        "reactant_a_smiles": "Brc1ccco1",
        "reactant_b_smiles": "Nc1ccccc1",
        "predicted_reaction_type": "C_N_Coupling",
        "matching_experiments": 112,
        "recommendations": [
            {
                "rank": 1,
                "avg_z_score": 2.61,
                "catalyst": "CuI",
                "ligand": "PPBO",
                "base": "NaOtBu",
                "solvent": "tAmOH",
                ...
            }
        ]
    }
}
```

If no data:
```json
{
    "success": false,
    "error": "No HTE data found for reactant combination: furan-diene + ArNH2 with reaction type 'C_N_Coupling' and catalyst 'Cu'",
    "data": {
        "reactant_a_type": "furan-diene",
        "reactant_b_type": "ArNH2",
        "matching_experiments": 0,
        "suggestion": "The reactants were successfully classified... Try: (1) removing filters...",
        "filters_applied": {
            "reaction_type": "C_N_Coupling",
            "catalyst": "Cu"
        }
    }
}
```

## Testing

### Test Case 1: Reaction SMILES with Data

```python
# ArBr + ArNH2 + Cu has 112 experiments
hte_recommend_tool(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    reaction_type_filter="C_N_Coupling",
    catalyst_filter="Cu",
    min_experiments=1
)
# ✅ Returns 112 conditions ranked by z-score
```

### Test Case 2: Reaction SMILES with No Data

```python
# furan-diene + ArNH2 + Cu has 0 experiments
hte_recommend_tool(
    reaction_smiles="Brc1ccco1.Nc1ccccc1>>c1ccccc1Nc1ccco1",
    reaction_type_filter="C_N_Coupling",
    catalyst_filter="Cu"
)
# ✅ Returns helpful error with suggestions
```

### Test Case 3: Individual Reactants (Original Format)

```python
hte_recommend_tool(
    reactant_a_smiles="Brc1ccccc1",
    reactant_b_smiles="Nc1ccccc1",
    catalyst_filter="Cu"
)
# ✅ Works as before, backward compatible
```

## Benefits

1. **Natural Language Support**: Agent can handle reaction SMILES from user queries
2. **Automatic Parsing**: No manual reactant extraction needed
3. **Better Error Messages**: Users understand why no data was found
4. **Actionable Guidance**: Suggestions help users refine their queries
5. **Z-Score Ranking**: Results prioritize statistically successful conditions
6. **Backward Compatible**: Original reactant input format still works

## Available C-N Coupling Data with Copper

Based on HTE database:
- ArI + Carbamate: 736 experiments
- ArBr + RNH2-a-branch: 672 experiments
- ArBr + arom-NH: 671 experiments
- ArI + arom-NH: 240 experiments
- ArBr + ArNH2: **112 experiments** ← Closest to user's query
- ArBr + Lactam: 160 experiments
- And 9 more combinations...

## Files Modified

- `chem_assistant/chemtools_wrapper.py`:
  - `HTERecommendInput` schema
  - `hte_recommend_tool` function
  - Error handling and messages

## Next Steps

The agent can now:
1. ✅ Parse reaction SMILES from queries
2. ✅ Extract and classify reactants automatically
3. ✅ Apply reaction type and catalyst filters
4. ✅ Return z-score ranked recommendations
5. ✅ Provide helpful suggestions when no data is found

To improve further:
- Add fuzzy matching for similar reactant types
- Suggest closest available combinations
- Pre-fetch common combinations for faster response
