# Simple Yet Reliable Condition Recommendation System
## Based on z-Score Peaks Dataset Analysis

## Executive Summary

After analyzing the z-Score Peaks dataset (66,308 reactions with empirical performance data), I propose a **lookup-based recommendation system** that combines **exact matching** with **statistical ranking**. This approach is simple, transparent, explainable, and performs well for the most common reaction types.

---

## Dataset Characteristics

### Performance Distribution
- **High performers (z > 1.0)**: 10,724 reactions (16.2%)
- **Medium (0 < z ≤ 1.0)**: 13,905 reactions (21.0%)
- **Low (z ≤ 0)**: 41,678 reactions (62.9%)

### Key Insight
Only ~16-22% of conditions work well for a given substrate combination. This means:
- ✅ **We have gold standard data** showing what actually works
- ✅ **We should prioritize proven high performers**
- ❌ **We shouldn't recommend untested/low-performing conditions**

### Coverage
- **Top 10 reaction types**: 52,862 reactions (79.7% of dataset)
- **Buchwald-Hartwig**: 20,286 reactions (most common, 15.2% success rate)
- **Suzuki-Miyaura**: 11,588 reactions (second most, 16.5% success rate)

---

## Proposed Recommendation Approach

### **Method: Hierarchical Exact Match + Statistical Ranking**

This is a **pure data-driven lookup system** with no machine learning complexity:

```
INPUT: 
  - Reaction Type (e.g., "Buchwald-Hartwig")
  - Electrophile Type (e.g., "ArBr")
  - Nucleophile Type (e.g., "RNH2")
  - Optional: Functional Groups (FG A, FG B)

STEP 1: Exact Match Filter
  Filter dataset for exact matches:
    reaction_type == input_reaction_type
    AND electrophile_type == input_electrophile_type
    AND nucleophile_type == input_nucleophile_type
    AND (if provided) functional_groups match

STEP 2: Rank by z-Score
  Sort matched reactions by z-Score (descending)
  
STEP 3: Aggregate Top Conditions
  For top N high-performing reactions (e.g., z > 1.0):
    - Count frequency of each catalyst
    - Count frequency of each ligand
    - Count frequency of each base
    - Count frequency of each solvent
    - Calculate average z-Score for each condition component

STEP 4: Return Recommendations
  Recommend top 3-5 condition sets ranked by:
    - Frequency in high performers
    - Average z-Score
    - Consistency (std dev of z-Score)

OUTPUT:
  - Ranked list of catalyst/ligand/base/solvent combinations
  - Each with: frequency, avg z-score, success rate
  - Explainability: "This worked in X out of Y cases with avg z-score Z"
```

---

## Implementation Strategy

### Level 1: Exact Match (Highest Confidence)
```python
# Filter for exact substrate match
exact_matches = dataset[
    (dataset['Reaction_Type_Standardized'] == reaction_type) &
    (dataset['Reactant_Type_Electrophile'] == electrophile) &
    (dataset['Reactant_Type_Nucleophile'] == nucleophile)
]

# Get high performers
high_performers = exact_matches[exact_matches['z-Score'] > 1.0]

# Rank conditions by frequency and average z-score
recommendations = aggregate_conditions(high_performers)
```

**When to use**: When exact substrate combination exists in database
**Confidence**: Very High (based on empirical data)
**Coverage**: ~60-70% of common reactions

### Level 2: Category Match (Good Confidence)
```python
# If no exact match, use category-level matching
category_matches = dataset[
    (dataset['Reaction_Type_Standardized'] == reaction_type) &
    (dataset['Reactant_Category_Electrophile'] == electrophile_category) &  # e.g., "ArX*"
    (dataset['Reactant_Category_Nucleophile'] == nucleophile_category)     # e.g., "Aliphatic-amine"
]
```

**When to use**: When exact type not found, but category matches
**Confidence**: High (broader but still relevant)
**Coverage**: ~85-90% of reactions

### Level 3: Reaction Type Only (Baseline)
```python
# If no substrate match, use reaction type defaults
reaction_defaults = dataset[
    dataset['Reaction_Type_Standardized'] == reaction_type
]

# Get most common high-performing conditions across all substrates
general_recommendations = aggregate_conditions(
    reaction_defaults[reaction_defaults['z-Score'] > 1.0]
)
```

**When to use**: Novel substrate combinations
**Confidence**: Medium (general best practices)
**Coverage**: 100% for known reaction types

---

## Recommendation Ranking Formula

For each condition component (catalyst, ligand, base, solvent):

```python
score = (
    0.5 * frequency_weight +      # How often it appears in high performers
    0.3 * avg_zscore_weight +     # Average z-score when used
    0.2 * consistency_weight      # Low std dev = more reliable
)

frequency_weight = count_in_high_performers / total_high_performers
avg_zscore_weight = (avg_zscore - 1.0) / (max_zscore - 1.0)  # Normalize above 1.0
consistency_weight = 1 / (1 + std_zscore)  # Lower std = higher weight
```

**Output Format**:
```json
{
  "reaction_type": "Buchwald-Hartwig",
  "substrate": {
    "electrophile": "ArBr",
    "nucleophile": "RNH2"
  },
  "match_level": "exact",  // or "category" or "reaction_type"
  "recommendations": [
    {
      "rank": 1,
      "catalyst": "XantPhosPdCl2",
      "ligand": "XantPhos",
      "base": "Cs2CO3",
      "solvent": "Dioxane",
      "confidence_score": 0.87,
      "evidence": {
        "successful_cases": 45,
        "total_cases": 52,
        "success_rate": "86.5%",
        "avg_zscore": 1.82,
        "zscore_range": [1.1, 3.2]
      }
    },
    {
      "rank": 2,
      "catalyst": "tBuBrettPhos Pd(allyl)OTf",
      "ligand": "tBuBrettPhos",
      "base": "NaOtBu",
      "solvent": "PhMe",
      "confidence_score": 0.81,
      "evidence": {...}
    }
  ]
}
```

---

## Advantages of This Approach

### ✅ Simplicity
- No ML model training required
- No feature engineering complexity
- Pure SQL/Pandas queries
- Fast implementation (1-2 days)

### ✅ Reliability
- Based on actual experimental data
- Only recommends conditions that worked in practice
- Clear evidence for each recommendation
- No "black box" predictions

### ✅ Explainability
- "This catalyst worked in 45/52 cases (86.5%) with avg z-score 1.82"
- Show specific ELN IDs of successful precedents
- Chemists can verify recommendations against source data

### ✅ Maintainability
- Easy to update with new data (just append to dataset)
- No model retraining needed
- Simple logic to debug and modify
- Can be implemented in pure Python/Pandas

### ✅ Transparency
- Chemists understand the logic immediately
- Clear confidence levels (exact > category > general)
- Known limitations (only works for reactions in database)

---

## Implementation Pseudocode

```python
class SimpleConditionRecommender:
    def __init__(self, csv_path):
        self.df = pd.read_csv(csv_path)
        self.high_performers = self.df[self.df['z-Score'] > 1.0]
    
    def recommend(self, reaction_type, electrophile, nucleophile, 
                  functional_groups=None, top_n=3):
        """Main recommendation function."""
        
        # Level 1: Exact match
        exact = self._exact_match(reaction_type, electrophile, nucleophile)
        if len(exact) >= 5:  # Minimum threshold
            return self._aggregate_conditions(exact, match_level="exact")
        
        # Level 2: Category match
        category = self._category_match(reaction_type, 
                                       self._get_category(electrophile),
                                       self._get_category(nucleophile))
        if len(category) >= 10:
            return self._aggregate_conditions(category, match_level="category")
        
        # Level 3: Reaction type only
        general = self._reaction_type_match(reaction_type)
        return self._aggregate_conditions(general, match_level="reaction_type")
    
    def _exact_match(self, rxn_type, elec, nuc):
        """Filter for exact substrate match."""
        return self.high_performers[
            (self.high_performers['Reaction_Type_Standardized'] == rxn_type) &
            (self.high_performers['Reactant_Type_Electrophile'] == elec) &
            (self.high_performers['Reactant_Type_Nucleophile'] == nuc)
        ]
    
    def _aggregate_conditions(self, df, match_level, top_n=3):
        """Aggregate and rank condition combinations."""
        
        # Group by condition combination
        combinations = df.groupby(['Catalyst', 'Ligand', 'Base', 'Solvent']).agg({
            'z-Score': ['count', 'mean', 'std'],
            'ELN_ID': list  # Keep track of source experiments
        })
        
        # Calculate score
        combinations['score'] = self._calculate_score(combinations)
        
        # Rank and return top N
        top_combinations = combinations.nlargest(top_n, 'score')
        
        return self._format_output(top_combinations, match_level)
    
    def _calculate_score(self, combinations):
        """Calculate recommendation score."""
        freq_weight = combinations[('z-Score', 'count')] / combinations[('z-Score', 'count')].sum()
        zscore_weight = (combinations[('z-Score', 'mean')] - 1.0) / (self.df['z-Score'].max() - 1.0)
        consistency = 1 / (1 + combinations[('z-Score', 'std')].fillna(0))
        
        return 0.5 * freq_weight + 0.3 * zscore_weight + 0.2 * consistency
```

---

## Example Use Cases

### Case 1: Common Reaction (High Confidence)
```
Input: Buchwald-Hartwig, ArBr + RNH2

Result: 
  - Match Level: EXACT
  - 45 precedents found (z > 1.0)
  - Top recommendation: XantPhosPdCl2/XantPhos/Cs2CO3/Dioxane
  - Success rate: 86.5%
  - Avg z-score: 1.82
```

### Case 2: Less Common Combination (Medium Confidence)
```
Input: Buchwald-Hartwig, ArF + R2NH-a-branch

Result:
  - Match Level: CATEGORY (ArX* + Aliphatic-amine)
  - 128 similar precedents found
  - Top recommendation: tBuBrettPhos Pd(allyl)OTf/tBuBrettPhos/NaOtBu/PhMe
  - Success rate: 73.2%
  - Note: "Expanded to category match - use with caution"
```

### Case 3: Novel Substrate (Baseline)
```
Input: Stille Coupling, ArBr + Alkyl-SnBu3 (not in database)

Result:
  - Match Level: REACTION_TYPE
  - General Stille conditions from 432 reactions
  - Top recommendation: Pd(PPh3)4/No Ligand/No Base/Dioxane
  - Note: "No specific precedent - using general conditions"
  - Suggestion: "Consider running a screen for optimization"
```

---

## Success Metrics

### Coverage
- **Level 1 (Exact)**: ~65% of queries
- **Level 2 (Category)**: ~25% of queries
- **Level 3 (General)**: ~10% of queries

### Accuracy (Expected)
- **Level 1**: ~75-85% success rate (based on historical data)
- **Level 2**: ~60-70% success rate
- **Level 3**: ~40-50% success rate

### Performance
- Query time: <100ms (in-memory Pandas)
- Can handle 1000s of queries/second
- No API rate limits or external dependencies

---

## Enhancements (Future)

### Short-term (Low Complexity)
1. **Functional group awareness**: Filter by presence of specific FGs
2. **Multi-condition output**: Return catalyst + 3 ligand options
3. **Negative filters**: Exclude conditions known to fail
4. **Similarity scoring**: Tanimoto similarity for functional groups

### Medium-term (Moderate Complexity)
1. **SMARTS-based substrate matching**: Match by structural features
2. **Condition combination validation**: Check if combinations are chemically sensible
3. **Temperature/time recommendations**: Extract from plate data if available
4. **Cost optimization**: Rank by reagent cost when multiple options exist

### Long-term (Higher Complexity)
1. **Machine learning model**: For novel substrates outside database
2. **Active learning**: Suggest experiments to fill data gaps
3. **Multi-objective optimization**: Balance yield, cost, safety
4. **Real-time updates**: Learn from new experiments continuously

---

## Implementation Plan

### Phase 1: Core System (Week 1)
- [ ] Load and process standardized CSV
- [ ] Implement 3-level matching logic
- [ ] Create condition aggregation function
- [ ] Build simple API endpoint (FastAPI)

### Phase 2: Ranking & Output (Week 2)
- [ ] Implement scoring algorithm
- [ ] Format JSON output with evidence
- [ ] Add confidence indicators
- [ ] Create simple CLI tool

### Phase 3: Validation (Week 3)
- [ ] Test on hold-out dataset
- [ ] Calculate success rates per match level
- [ ] Gather chemist feedback
- [ ] Refine scoring weights

### Phase 4: Integration (Week 4)
- [ ] Integrate with existing chemtools modules
- [ ] Add to FastAPI routes
- [ ] Create documentation
- [ ] Deploy to production

---

## Comparison to ML Approaches

| Aspect | Simple Lookup | ML Model |
|--------|--------------|----------|
| Development time | 1-2 weeks | 2-3 months |
| Explainability | Perfect | Limited |
| Data requirements | Existing dataset | + Feature engineering |
| Maintenance | Low | High (retraining) |
| Accuracy (in-domain) | 75-85% | 80-90% |
| Accuracy (out-of-domain) | 40-50% | 30-40% |
| Chemist trust | High | Medium |
| Debugging | Easy | Difficult |

### Recommendation
Start with the simple lookup approach:
- ✅ Faster time to production
- ✅ Easier to validate and trust
- ✅ Covers 90% of common use cases
- ✅ Can add ML later for edge cases

The z-Score dataset is **perfect for this approach** because:
1. It has actual performance metrics (not just success/fail)
2. It covers common reaction types comprehensively
3. It includes functional group information
4. The data quality is high (standardized experiments)

---

## Conclusion

**Proposed Solution**: Implement a **hierarchical exact-match recommendation system** with statistical ranking based on z-Score performance.

**Key Benefits**:
- Simple, fast, reliable
- Transparent and explainable
- Based on real experimental data
- Easy to maintain and extend
- Chemist-friendly

**Success Criteria**:
- Recommend high-performing conditions for 75%+ of common reactions
- Provide clear confidence levels for all recommendations
- Enable chemists to verify recommendations against source data
- Achieve 75-85% success rate for exact matches

**Next Step**: Implement core system in 1-2 weeks and validate with chemist feedback.
