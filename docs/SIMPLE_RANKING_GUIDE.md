# Simple Precedent-Centric Recommendation

## Philosophy

**Precedents are the most important** - they are real experimental data. However, dataset quality issues can lead to errors (e.g., Ullmann getting Suzuki results). We need simple methods to **filter/rerank precedents** using chemical knowledge.

## The Simplified Approach

### Core Workflow

```
1. Find k similar precedents (DRFP similarity)
2. Optionally rerank precedents using:
   - Rule-based: Boost precedents matching rule patterns
   - Analytics-based: Boost precedents using popular reagents
   - None: Use similarity only
3. Extract top conditions from reranked precedents
```

### Three Strategies

#### Strategy 1: Similarity Only (Baseline)
```python
recommend_simple(
    reaction_smiles=rxn,
    family="Suzuki",
    rerank_strategy='none'
)
```
- Simple DRFP k-NN ranking
- No additional filtering
- Fast, straightforward
- **Use when**: High-quality precedents available

#### Strategy 2: Rule-Based Reranking
```python
recommend_simple(
    reaction_smiles=rxn,
    family="Suzuki",
    rerank_strategy='rule'
)
```
- Match reaction against rule database
- Boost precedents whose reagents match rule recommendations
- Scoring: `similarity × (1 + boost)`
  - Catalyst match: +30% boost
  - Base match: +20% boost
  - Solvent match: +20% boost
- **Use when**: Strong chemical rules available (e.g., Ullmann → Cu, Buchwald → Pd)

#### Strategy 3: Analytics-Based Reranking
```python
recommend_simple(
    reaction_smiles=rxn,
    family="Suzuki",
    rerank_strategy='analytics'
)
```
- Get most popular reagents from full dataset
- Boost precedents using popular reagents
- Scoring: `similarity × (1 + boost)`
  - Popular catalyst: +30% boost (scaled by rank)
  - Popular base: +20% boost
  - Popular solvent: +20% boost
- **Use when**: Large, high-quality dataset with clear trends

## Comparison with Old Fusion Approach

### Old Fusion (Complex)
```python
# ❌ Complex adaptive weighting
weights = compute_adaptive_weights(evidence)
  → α (precedents) adjusted by count, diversity, similarity
  → β (analytics) adjusted by dataset size
  → γ (rules) adjusted by confidence
  → δ (ML) adjusted by model availability

# ❌ Multiple scoring functions
score_from_precedents(candidate, evidence)
score_from_analytics(candidate, evidence)
score_from_rules(candidate, evidence)
score_from_ml(candidate, evidence)

# ❌ Fusion formula
total_score = α×precedent + β×analytics + γ×rules + δ×ML

# ❌ Multi-pass ranking
1. Initial fusion scoring
2. Rule-alignment reranking
3. Deduplication
```

### New Simple Approach (Clear)
```python
# ✅ Simple precedent ranking
precedents = find_similar_reactions(k=30)

# ✅ Optional reranking
if strategy == 'rule':
    precedents = boost_by_rule_match(precedents)
elif strategy == 'analytics':
    precedents = boost_by_popularity(precedents)

# ✅ Extract top conditions
top_conditions = Counter(precedents[:20])
```

## Test Results

For Suzuki reaction `Clc1ccc(C#N)cc1 + furan-boronic-acid`:

| Strategy | #1 Catalyst | #1 Base | Top Precedent Yield |
|----------|------------|---------|-------------------|
| Similarity only | Pd | K3PO4 (7778-53-2) | 71% |
| Rule-reranked | Pd | K3PO4 (7778-53-2) | 71% |
| Analytics-reranked | Pd | K3PO4 (7778-53-2) | 98% |

**Notice**: Analytics reranking moved the 98% yield precedent to #1 position!

## When to Use Which Strategy

### Similarity Only
- ✅ High-quality, curated datasets
- ✅ Need fast recommendations
- ✅ Query very similar to precedents
- ❌ Dataset has quality issues
- ❌ Query is edge case

### Rule-Based Reranking  
- ✅ Strong chemical rules exist (Ullmann→Cu, Buchwald→Pd)
- ✅ Dataset has systematic errors (wrong catalyst family)
- ✅ Need chemical validity guarantees
- ❌ Rules are incomplete/outdated
- ❌ Novel reaction types

### Analytics-Based Reranking
- ✅ Large dataset with clear trends
- ✅ Want to favor popular/robust conditions
- ✅ Dataset is high-quality overall
- ❌ Small dataset (<1000 reactions)
- ❌ Dataset has batch effects

## Implementation Details

### Boost Calculation

Both strategies use multiplicative boost:
```python
combined_score = similarity × (1.0 + boost)
```

**Why multiplicative?**
- Preserves relative similarity order
- Precedent with 0.9 similarity + 30% boost = 1.17
- Precedent with 0.5 similarity + 30% boost = 0.65
- High similarity still wins, but matches get lifted

**Boost ranges:**
- Catalyst match: 0.0 - 0.3 (30%)
- Base match: 0.0 - 0.2 (20%)
- Solvent match: 0.0 - 0.2 (20%)
- **Maximum total boost**: 70% (all three match)

### Rule Matching

```python
# Load rule database for family
rule_db = chem.rules.load_database("suzuki_db.json")
match_result = chem.rules.match(rule_db, reaction_smiles)

# Extract recommended reagents
rule_catalyst = "Pd(PPh3)4"
rule_base = "584-08-7"  # K2CO3
rule_solvent = "7732-18-5"  # Water

# Boost precedents matching these
```

### Analytics Ranking

```python
# Get dataset statistics
top_cores = chem.analytics.get_condition_cores("Suzuki", top_n=10)
  → [(core1, count1, avg_yield1), (core2, count2, avg_yield2), ...]

# Create popularity ranks
core_ranks = {core1: 10, core2: 9, core3: 8, ...}

# Boost calculation (scaled by rank)
if precedent.core in core_ranks:
    boost += (core_ranks[precedent.core] / 10) * 0.3
```

## Next Steps

1. **Evaluate on validation set**: Which strategy gives best results?
2. **Add to CLI**: Let users choose strategy with `--rerank rule` or `--rerank analytics`
3. **Hybrid approach**: Use rule-reranking for families with strong rules (Ullmann, Buchwald), analytics for others (Suzuki, Heck)
4. **Monitor performance**: Track which strategy works best for each family

## Code Location

- **Implementation**: `chemtools/ml/simple_precedent_ranker.py`
- **Test script**: `test_simple_ranker.py`
- **Usage**:
  ```python
  from chemtools.ml.simple_precedent_ranker import recommend_simple
  
  result = recommend_simple(
      reaction_smiles="Br...>>...",
      family="Suzuki",
      k=30,
      rerank_strategy='rule'  # or 'analytics' or 'none'
  )
  
  print(result['top_cores'])  # Most common catalysts
  print(result['reasoning'])  # Explanation
  ```
