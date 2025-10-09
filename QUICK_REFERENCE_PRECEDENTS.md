# Quick Reference: ML Recommendation Output with Precedent Details

## ✅ Your Request is Complete!

The ML-based recommendation output JSON now includes **all precedent information** you requested.

## What You Asked For ✅

| Item | Status | Location in JSON |
|------|--------|------------------|
| CAS reaction number | ✅ Complete | `precedents_used.top_precedents[].reaction_id` |
| Reaction SMILES | ✅ Complete | `precedents_used.top_precedents[].reaction_smiles` |
| Reagents (catalysts, bases, solvents) | ✅ Complete | `precedents_used.top_precedents[].catalysts/reagents/solvents` |
| Simple literature reference | ✅ Complete | `precedents_used.top_precedents[].reference` |

## Example Precedent Entry

```json
{
  "rank": 1,
  "reaction_id": "31-179-CAS-9031327",
  "reaction_smiles": "OB(O)c1ccoc1.Cc1ccc(...)>>...",
  "core": "Pd/SPhos",
  "yield": 71.0,
  "catalysts": [
    {
      "name": "Palladium(II) acetate",
      "cas": "3375-31-3",
      "role": "catalyst/ligand"
    },
    {
      "name": "SPhos",
      "cas": "657408-07-6",
      "role": "catalyst/ligand"
    }
  ],
  "reagents": [
    {
      "name": "Tripotassium phosphate",
      "cas": "7778-53-2",
      "role": "base"
    }
  ],
  "solvents": [
    {
      "name": "Toluene",
      "cas": "108-88-3"
    }
  ],
  "conditions": {
    "temperature_C": null,
    "time_h": null,
    "yield_pct": 71.0
  },
  "reference": "Full citation here..."
}
```

## How to Use Right Now

### Option 1: Run your existing CLI script
```bash
python scripts/local_recommendation_cli.py
```

The output JSON files will now include `precedents_used` section.

### Option 2: Quick test
```bash
python test_cli_precedents_used.py
```

### Option 3: Generate sample JSON
```bash
python generate_ml_with_precedents.py
```

This creates: `results/ml_recommendation_with_precedents_cli_test.json`

## What Changed

✅ **`chemtools/recommend/core.py`**: Added precedent details extraction  
✅ **`scripts/local_recommendation_cli.py`**: Preserves precedents_used in output  
✅ **All ML recommendations**: Automatically include precedent details  

## Test Results

```
✅ precedents_used section found in output!
   Total precedents: 4
   Top precedents count: 4

   First precedent:
   - Reaction ID: 31-179-CAS-9031327
   - Core: Pd/SPhos
   - Yield: 71.0%
   - Catalysts: ['Palladium(II) acetate', 'SPhos']
   - Reagents: ['Tripotassium phosphate']
   - Solvents: ['Toluene']

   Statistics:
   - Average yield: 83.8%
   - Yield range: [71.0, 98.0]

✅ Test PASSED
```

## Next Time You Run

When you run `local_recommendation_cli.py` or any ML recommendation, the output will automatically include:

- **Top 10 precedents** with full details
- **CAS reaction numbers** for each precedent
- **Complete reagent information** (catalysts, bases, solvents with CAS numbers)
- **Literature references** for each precedent
- **Statistics** (average yield, ranges)

---

**Everything is ready to use!** 🎉

Just run your normal recommendation commands and the JSON output will include all the precedent information you requested.
