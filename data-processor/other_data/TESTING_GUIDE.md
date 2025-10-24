# Testing Guide for Simple Condition Recommender

## Quick Start

You have a complete testing framework in `test_recommender.py`. Here's how to use it:

### Run the Interactive Test Menu

```powershell
python data-processor/other_data/test_recommender.py
```

This gives you 7 testing options:

---

## Test Options Explained

### **Test 1: Known Combinations**
Tests if the recommender can reproduce conditions that actually worked in the dataset.

**What it does:**
- Samples random successful reactions from the dataset
- Gets recommendations for those substrate combinations
- Checks if the actual successful conditions appear in top 3/5 recommendations

**How to interpret:**
- High "Found in top 3/5" rate (>50%) = Good recommendation accuracy
- Low rate might mean we need to adjust the scoring algorithm

**Run manually:**
```python
from test_recommender import RecommenderTester
from pathlib import Path

tester = RecommenderTester("data-processor/other_data/z-Score Peaks with FG_STANDARDIZED.csv")
tester.test_known_combinations(n_samples=20)
```

---

### **Test 2: Coverage by Reaction Type**
Shows which reaction types have good/poor data coverage.

**What it shows:**
- How many high-performing examples exist for each reaction type
- Number of unique substrate combinations
- Average precedents per combination
- Top substrate combinations for each reaction

**Why it matters:**
- Reactions with many precedents = more reliable recommendations
- Reactions with few precedents = might need category or general fallback

**Example output:**
```
Buchwald-Hartwig:
  High performers: 3,077
  Unique substrate combos: 145
  Avg precedents per combo: 21.2
  Top combinations:
    ArBr + RNH2-a-branch: 112 precedents
    ArCl + RNH2: 98 precedents
```

---

### **Test 3: Interactive Testing**
Manually test any substrate combination you want.

**How to use:**
1. Shows you the top 10 available reaction types, electrophiles, nucleophiles
2. You input your desired combination
3. Get instant recommendations with full details

**Example session:**
```
Reaction Type: Buchwald-Hartwig
Electrophile Type: ArBr
Nucleophile Type: RNH2

→ Returns top 3-5 recommendations with z-scores, precedent counts, ELN IDs
```

**When to use:**
- Testing specific chemistry you're interested in
- Exploring what's available in the dataset
- Validating recommendations make chemical sense

---

### **Test 4: Edge Cases & Error Handling**
Tests how the system handles bad inputs.

**Tests:**
- Non-existent reaction types
- Invalid substrate types
- Rare combinations (should fall back gracefully)
- Empty strings

**Why it matters:**
- Ensures the system doesn't crash on bad input
- Validates fallback logic works correctly

---

### **Test 5: Recommendation Quality**
Analyzes if we're recommending the BEST conditions available.

**What it measures:**
- Average z-score of top-1 recommendation
- Average z-score of top-3 recommendations
- Gap between our recommendation and the best possible z-score

**Good results:**
- Top-1 recommendation z-score > 2.0
- Gap (max - top1) < 1.0

**If results are poor:**
- May need to adjust scoring weights
- Consider adding more features (functional groups, etc.)

---

### **Test 6: Run All Tests**
Runs tests 1, 2, 4, 5 in sequence (non-interactive).

**When to use:**
- Quick validation after code changes
- Before committing changes
- Generating comprehensive test results

---

### **Test 7: Generate Test Report**
Creates a markdown report with:
- Dataset coverage statistics
- Top reaction types with details
- Unique combinations available

**Output:** `test_report.md`

**When to use:**
- Documentation
- Sharing results with team
- Tracking changes over time

---

## Simple Testing Examples

### Example 1: Quick Validation

```python
from simple_condition_recommender import SimpleConditionRecommender
from pathlib import Path

# Initialize
csv_path = Path("data-processor/other_data/z-Score Peaks with FG_STANDARDIZED.csv")
recommender = SimpleConditionRecommender(csv_path)

# Test a known good combination
result = recommender.recommend(
    reaction_type="Buchwald-Hartwig",
    electrophile="ArBr",
    nucleophile="RNH2"
)

print(f"Match level: {result['match_level']}")
print(f"Precedents found: {result['total_precedents']}")
print(f"Top recommendation: {result['recommendations'][0]['catalyst']}")
print(f"Confidence: {result['recommendations'][0]['confidence_score']}")
```

---

### Example 2: Test Multiple Reactions

```python
test_cases = [
    ("Buchwald-Hartwig", "ArBr", "RNH2"),
    ("Suzuki-Miyaura", "ArCl", "ArB(OH)2"),
    ("CN-Coupling", "ArBr", "ArNH2"),
    ("Amide-coupling", "RCO2H", "RNH2"),
]

for rxn, elec, nuc in test_cases:
    result = recommender.recommend(rxn, elec, nuc)
    print(f"\n{rxn} ({elec} + {nuc}):")
    print(f"  Match: {result['match_level']}")
    print(f"  Precedents: {result['total_precedents']}")
    if result['recommendations']:
        top = result['recommendations'][0]
        print(f"  Top: {top['catalyst']}/{top['base']}/{top['solvent']}")
        print(f"  Z-score: {top['evidence']['avg_zscore']}")
```

---

### Example 3: Validate Fallback Logic

```python
# Test exact match
result1 = recommender.recommend("Buchwald-Hartwig", "ArBr", "RNH2")
print(f"Common combo: {result1['match_level']}")  # Should be 'exact'

# Test category match
result2 = recommender.recommend("Buchwald-Hartwig", "ArF", "R2NH-a-branch")
print(f"Less common: {result2['match_level']}")  # Might be 'category'

# Test reaction-only match
result3 = recommender.recommend("Stille", "FakeElec", "FakeNuc")
print(f"Novel combo: {result3['match_level']}")  # Should be 'reaction_type'
```

---

## Interpreting Test Results

### Good Signs ✅
- **Match level "exact"** for common reactions (Buchwald-Hartwig, Suzuki, etc.)
- **Total precedents > 20** for exact matches
- **Confidence scores > 0.3** for top recommendations
- **Average z-scores > 1.5** for recommended conditions
- **Multiple recommendations** available (diversity)

### Warning Signs ⚠️
- **Match level "reaction_type"** for supposedly common reactions
  - → May need more data or better reactant type standardization
- **Very few precedents** (<5) for exact matches
  - → Recommendations less reliable, consider screening
- **Low confidence scores** (<0.2)
  - → Many competing options, unclear winner
- **Very high z-score gaps** (>2.0 between max and top1)
  - → We're missing the best conditions, check scoring algorithm

### Red Flags 🚩
- **No recommendations** returned
  - → Reaction type not in database
- **All recommendations have z-score < 1.5**
  - → No clearly successful conditions found
- **Recommended conditions don't make chemical sense**
  - → Data quality issue or aggregation bug

---

## Advanced Testing

### Test with Your Own Data

If you have experimental results:

```python
# Your experiment
your_rxn = "Buchwald-Hartwig"
your_elec = "ArBr"
your_nuc = "RNH2"
your_actual_catalyst = "XantPhosPdCl2"
your_actual_ligand = "XantPhos"
your_actual_base = "Cs2CO3"
your_actual_solvent = "Dioxane"

# Get recommendations
result = recommender.recommend(your_rxn, your_elec, your_nuc, top_n=10)

# Check if your conditions are recommended
for i, rec in enumerate(result['recommendations'], 1):
    if (rec['catalyst'] == your_actual_catalyst and
        rec['ligand'] == your_actual_ligand and
        rec['base'] == your_actual_base and
        rec['solvent'] == your_actual_solvent):
        print(f"✓ Your conditions found at rank {i}")
        print(f"  Z-score: {rec['evidence']['avg_zscore']}")
        print(f"  Success rate: {rec['evidence']['success_rate']}")
        break
else:
    print("✗ Your conditions not in top 10 recommendations")
    print("  This might indicate:")
    print("  - Your conditions are novel (not in database)")
    print("  - Your conditions underperformed historically")
    print("  - Need more data for this substrate combination")
```

---

### Batch Testing

Test many reactions from a file:

```python
import pandas as pd

# Your reactions to test (CSV with columns: reaction_type, electrophile, nucleophile)
test_df = pd.read_csv("my_reactions_to_test.csv")

results = []
for idx, row in test_df.iterrows():
    result = recommender.recommend(
        row['reaction_type'],
        row['electrophile'],
        row['nucleophile']
    )
    
    results.append({
        'reaction': f"{row['reaction_type']} ({row['electrophile']} + {row['nucleophile']})",
        'match_level': result['match_level'],
        'precedents': result['total_precedents'],
        'top_catalyst': result['recommendations'][0]['catalyst'] if result['recommendations'] else None,
        'top_zscore': result['recommendations'][0]['evidence']['avg_zscore'] if result['recommendations'] else None
    })

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv("recommendation_test_results.csv", index=False)
print(results_df)
```

---

## Next Steps After Testing

### If Results Look Good (>60% exact matches, reasonable z-scores):
1. ✅ Document successful test cases
2. ✅ Add the recommender to your main application
3. ✅ Create API endpoint
4. ✅ Build simple UI

### If Results Need Improvement:
1. 🔧 Adjust scoring weights (currently 50% frequency, 30% z-score, 20% consistency)
2. 🔧 Try different z-score thresholds (currently >1.0 for "high performer")
3. 🔧 Add functional group filtering
4. 🔧 Refine category matching logic
5. 🔧 Consider combining multiple conditions (e.g., catalyst+ligand as a pair)

---

## Troubleshooting

**"No precedents found"**
- Check spelling of reaction type, electrophile, nucleophile
- Use exact IDs from the standardized dataset
- Try interactive test (option 3) to see available options

**"All recommendations have low z-scores"**
- The substrate combination might be inherently difficult
- Consider screening multiple conditions experimentally
- Check if functional groups are problematic

**"Test 1 shows 0% found in top 3"**
- The scoring algorithm may need adjustment
- Current grouping by all 4 conditions might be too strict
- Consider ranking catalyst+ligand combinations separately

**"Memory errors"**
- Dataset is large (66K reactions)
- Close other applications
- Use smaller test samples (n_samples=10 instead of 50)

---

## Summary

**To test the recommender:**
```powershell
python data-processor/other_data/test_recommender.py
```

**Quick tests to run:**
1. Option 2 (Coverage) - See what's in the dataset
2. Option 3 (Interactive) - Test your specific reactions
3. Option 6 (All tests) - Comprehensive validation

**Good results indicate:**
- Exact matches for common reactions
- High z-scores (>1.5) in recommendations
- Clear confidence differentiation
- Chemically sensible suggestions

**The system is ready when:**
- You're confident in the recommendations
- Test results are reproducible
- Edge cases are handled gracefully
- Documentation is complete

Then you can integrate into your main application! 🚀
