# How to Test the Simple Condition Recommender

## Quick Start - 3 Ways to Test

### 1. **Quick Test (Recommended for First Time)**
Runs 5 example reactions showing different scenarios:

```powershell
python data-processor/other_data/quick_test.py
```

**What you'll see:**
- Common reactions with exact matches (Buchwald-Hartwig, Suzuki)
- Less common combinations that fall back to category matching
- Full output with catalysts, ligands, bases, solvents
- Z-scores and success rates for each recommendation

**Run time:** ~10 seconds

---

### 2. **Interactive Testing**
Test your own specific reactions:

```powershell
python data-processor/other_data/test_recommender.py
```

Then select **Option 3: Interactive testing**

**You can test:**
- Any reaction type (Buchwald-Hartwig, Suzuki, Negishi, etc.)
- Any electrophile (ArBr, ArCl, Alkyl-Br, etc.)
- Any nucleophile (RNH2, ArB(OH)2, ROH, etc.)

**Example session:**
```
Reaction Type: Buchwald-Hartwig
Electrophile Type: ArBr  
Nucleophile Type: RNH2

→ Returns top 5 recommendations with full details
```

---

### 3. **Comprehensive Testing**
Run full test suite (all automated tests):

```powershell
python data-processor/other_data/test_recommender.py
```

Then select **Option 6: Run all tests**

**Tests include:**
1. Known combinations validation
2. Coverage by reaction type
3. Edge case handling
4. Recommendation quality analysis

**Run time:** ~1-2 minutes

---

## Test Files Overview

| File | Purpose | When to Use |
|------|---------|-------------|
| `quick_test.py` | 5 quick examples | First time testing, quick validation |
| `test_recommender.py` | Full test framework (7 options) | Comprehensive testing, validation |
| `simple_condition_recommender.py` | The actual recommender | Import into your code |

---

## What to Look For in Test Results

### ✅ **Good Signs**

**Match Levels:**
- `exact` for common reactions (Buchwald-Hartwig, Suzuki, etc.)
- `category` for less common substrates
- `reaction_type` only for truly novel combinations

**Precedent Counts:**
- \>100 precedents = Excellent coverage
- 20-100 precedents = Good coverage
- 5-20 precedents = Moderate coverage
- <5 precedents = Poor coverage, be cautious

**Z-Scores:**
- Top recommendation z-score > 2.0 = Excellent
- Top recommendation z-score 1.5-2.0 = Good
- Top recommendation z-score 1.0-1.5 = Moderate
- Top recommendation z-score < 1.0 = Poor (shouldn't happen with default threshold)

**Confidence Scores:**
- \>0.4 = High confidence (clear winner)
- 0.3-0.4 = Good confidence
- 0.2-0.3 = Moderate (multiple viable options)
- <0.2 = Low (many competing conditions)

### ⚠️ **Warning Signs**

- `reaction_type` match for supposedly common reactions → Data coverage issue
- Very few precedents (<5) → Recommendations less reliable
- All recommendations have similar confidence scores → No clear winner
- Z-score < 1.5 for all recommendations → No clearly successful conditions

### 🚩 **Red Flags**

- No recommendations returned → Reaction type not in database
- Recommended conditions don't make chemical sense → Check data quality
- Very high z-score variance (std > 2.0) → Inconsistent results

---

## Example Test Session

```powershell
PS> python data-processor/other_data/quick_test.py
```

**Output you'll see:**

```
Loading recommender...
Total reactions: 66,308
High performers (z > 1.0): 10,724 (16.2%)

*** TEST 1: Common Buchwald-Hartwig ***
✓ Found 112 exact matches (high performers)

Reaction: Buchwald-Hartwig
Substrate: ArBr + RNH2
Match Level: EXACT
Total Precedents: 112

Top 3 Recommendations:
#1 (Confidence: 0.504)
  Catalyst: dtbpfPdCl2
  Ligand:   dtbpf
  Base:     KHMDS
  Solvent:  tAmOH
  Evidence:
    - Success cases: 1 / 112 (0.9%)
    - Avg z-score: 8.13
    - Z-score range: [8.13, 8.13]

... [more recommendations] ...

SUMMARY:
✓ Test 1 (BH ArBr+RNH2)     | Match: exact    | Precedents: 112 | Recs: 3
✓ Test 2 (Suzuki ArCl+ArB)  | Match: exact    | Precedents: 107 | Recs: 3
✓ Test 4 (C-O ArBr+ROH)     | Match: category | Precedents: 110 | Recs: 3
```

**Interpretation:**
- ✅ System is working correctly
- ✅ Exact matches found for common reactions
- ✅ Fallback to category works for less common substrates
- ✅ High z-scores (>3.0) in recommendations
- ✅ Multiple recommendations available

---

## Programmatic Testing (For Your Own Code)

### Test a Specific Reaction

```python
from simple_condition_recommender import SimpleConditionRecommender
from pathlib import Path

# Initialize
csv_path = "data-processor/other_data/z-Score Peaks with FG_STANDARDIZED.csv"
recommender = SimpleConditionRecommender(csv_path)

# Test
result = recommender.recommend(
    reaction_type="Buchwald-Hartwig",
    electrophile="ArBr",
    nucleophile="RNH2"
)

# Check results
print(f"Match level: {result['match_level']}")
print(f"Precedents: {result['total_precedents']}")

if result['recommendations']:
    top = result['recommendations'][0]
    print(f"\nTop recommendation:")
    print(f"  Catalyst: {top['catalyst']}")
    print(f"  Z-score: {top['evidence']['avg_zscore']}")
    print(f"  Confidence: {top['confidence_score']}")
else:
    print("No recommendations found")
```

### Test Multiple Reactions at Once

```python
test_cases = [
    ("Buchwald-Hartwig", "ArBr", "RNH2"),
    ("Suzuki-Miyaura", "ArCl", "ArB(OH)2"),
    ("CN-Coupling", "ArBr", "ArNH2"),
]

for rxn, elec, nuc in test_cases:
    result = recommender.recommend(rxn, elec, nuc)
    status = "✓" if result['recommendations'] else "✗"
    print(f"{status} {rxn:20s} {result['match_level']:15s} {result['total_precedents']:4d} precedents")
```

---

## Troubleshooting

### "No recommendations found"

**Possible causes:**
1. Reaction type not in database → Check spelling, use standardized names
2. Substrate combination too rare → Try category-level search
3. Z-score threshold too high → Lower from 1.0 to 0.5

**How to fix:**
```python
# Check available reaction types
import pandas as pd
df = pd.read_csv("z-Score Peaks with FG_STANDARDIZED.csv")
print(df['Reaction_Type_Standardized'].value_counts().head(20))

# Check available electrophiles
print(df['Reactant_Type_Electrophile'].value_counts().head(20))

# Check available nucleophiles
print(df['Reactant_Type_Nucleophile'].value_counts().head(20))
```

### "Only getting reaction_type matches"

**Possible causes:**
1. Substrate types not matching exactly → Check exact string matches
2. Reactant types not standardized → Verify against reactant_types.json
3. Very rare combination → Expected behavior

**How to check:**
```python
# Check if exact combination exists
filtered = df[
    (df['Reaction_Type_Standardized'] == 'Buchwald-Hartwig') &
    (df['Reactant_Type_Electrophile'] == 'ArBr') &
    (df['Reactant_Type_Nucleophile'] == 'RNH2')
]
print(f"Found {len(filtered)} reactions with this combination")
print(f"High performers (z>1): {len(filtered[filtered['z-Score']>1.0])}")
```

### "Low confidence scores"

**This is normal when:**
- Many different conditions work equally well
- No single dominant catalyst/ligand combination
- Dataset has high diversity for this substrate

**What to do:**
- Look at top 3-5 recommendations (not just top 1)
- Consider running a small screen with top options
- Check if functional groups might be affecting selectivity

---

## Next Steps After Testing

### If Results Look Good:
1. ✅ Test with your own specific reactions
2. ✅ Validate recommendations match your domain knowledge
3. ✅ Consider integrating into your workflow
4. ✅ Document any issues or edge cases

### If Results Need Improvement:
1. 🔧 Adjust z-score threshold (try 0.5 or 1.5 instead of 1.0)
2. 🔧 Modify scoring weights in `_aggregate_conditions()`
3. 🔧 Add functional group filtering
4. 🔧 Improve category matching logic

### Before Production Use:
1. ⚠️ Test on hold-out data (reactions not used in development)
2. ⚠️ Validate with expert chemists
3. ⚠️ Document known limitations
4. ⚠️ Set up monitoring for recommendation quality

---

## Summary

**To test the recommender quickly:**
```powershell
python data-processor/other_data/quick_test.py
```

**To test your specific reactions:**
```powershell
python data-processor/other_data/test_recommender.py
# Select option 3 (Interactive)
```

**To run comprehensive tests:**
```powershell
python data-processor/other_data/test_recommender.py
# Select option 6 (Run all tests)
```

**Good results mean:**
- Exact matches for common reactions
- Reasonable z-scores (>1.5)
- Multiple recommendations available
- Chemically sensible suggestions

**The system is production-ready when:**
- Test results meet your quality standards
- Recommendations align with expert knowledge
- Edge cases handled gracefully
- Performance is acceptable (<100ms per query)

Then you can proceed with integration! 🚀
