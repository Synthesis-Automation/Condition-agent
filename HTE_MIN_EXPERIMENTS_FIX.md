# HTE Agent Fix: min_experiments Parameter

## Problem

Agent was returning "no data found" for queries like:
```
use HTE system, recommend conditions for reactions: Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1, its a C-N coupling reaction, and the catalyst is copper
```

Even though:
- HTERecommender found 112 matching experiments
- Reactants correctly classified as ArBr + ArNH2
- Database contains relevant data

## Root Cause

The `hte_recommend_tool` had `min_experiments=2` as default parameter:
```python
def hte_recommend_tool(
    ...
    min_experiments: int = 2,  # ❌ Too strict!
    ...
)
```

HTERecommender groups experiments by condition (catalyst, ligand, base, solvent) and **filters out any group with fewer than min_experiments**. In this case:
- 112 total experiments found
- BUT many individual conditions had only 1 experiment each
- With `min_experiments=2`, all conditions with 1 experiment were filtered out
- Result: `num_recommendations=0` despite `total_matching_experiments=112`

## Solution

Changed default to `min_experiments=1`:
```python
def hte_recommend_tool(
    ...
    min_experiments: int = 1,  # ✅ More permissive
    ...
)
```

## Results

**Before Fix:**
```
total_matching_experiments: 112
num_recommendations: 0
Agent: "no data found"
```

**After Fix:**
```
total_matching_experiments: 112
num_recommendations: 5
Agent: Returns top 5 recommendations with CuI/PPBO/NaOtBu/tAmOH as top hit
```

## Testing

Run test script:
```bash
python test_agent_hte.py
```

Expected output: Agent returns 5 recommendations with z-scores, yields, and confidence scores.

## Impact

- Users can now get recommendations even for less-studied substrate combinations
- Still maintains data quality (z-score based ranking)
- Users can override with higher `min_experiments` if they want more robust statistics
- More permissive default aligns with exploratory use case

## Files Changed

- `chem_assistant/chemtools_wrapper.py`: Changed default `min_experiments` from 2 to 1
- `chem_assistant/gui/main_window.py`: Enabled verbose mode for debugging (optional, can revert)

## Related Issues

This fix completes the HTE agent integration along with previous fixes:
1. ✅ SMARTS pattern fix for benzylic halides (Bn-Br classification)
2. ✅ Added HTE tools to system prompt
3. ✅ **min_experiments parameter fix (this document)**
