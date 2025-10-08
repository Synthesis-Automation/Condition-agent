# Implementation Guide: Multi-Source Evidence Fusion

## Executive Summary

**Problem**: Current ML recommendations over-rely on top-k precedents, which can be biased and unrepresentative of the full dataset.

**Solution**: Multi-source evidence fusion that balances:
- Precedent search (DRFP k-NN)
- Dataset analytics (frequency & yield statistics)
- Rule-based matching (SMARTS schemes)
- ML yield prediction

**Key Innovation**: Adaptive weighting adjusts trust in each evidence source based on data quality.

## Quick Start

### Running the Demo

```bash
python demo_fusion_recommendation.py
```

### Example Output

```
⚖️  Adaptive Weights:
   α (precedents): 0.412
   β (analytics):  0.441
   γ (rules):      0.147
   δ (ML):         0.000

💡 Reasoning:
   • Low similarity (0.50) → precedents less relevant
   • No strong rule match → low rule weight
   • No ML model available → δ=0
```

## Architecture

### Evidence Collection (`collect_evidence`)

```python
evidence = {
    'precedents': {
        'reactions': [...],       # Top-k similar reactions
        'similarities': [...],     # DRFP similarity scores
        'diversity_score': 0.536,  # How diverse are precedents?
        'coverage': 10,            # Number found
        'avg_similarity': 0.500    # Average similarity
    },
    'analytics': {
        'catalytic_systems': [(system, count, yield), ...],
        'bases': [(base, count, yield), ...],
        'solvents': [(solvent, count, yield), ...],
        'dataset_size': 1343
    },
    'rules': {
        'matched_schemes': [...],
        'confidence': 0.0
    },
    'ml': {
        'model_available': False
    }
}
```

### Adaptive Weighting (`compute_adaptive_weights`)

Adjusts weights α, β, γ, δ based on:

| Condition | Weight Adjustment | Reasoning |
|-----------|------------------|-----------|
| `n_precedents < 10` | α ×0.5, β ×1.5 | Few precedents → trust analytics |
| `diversity < 0.3` | α ×0.7 | Low diversity → precedents biased |
| `similarity < 0.6` | α ×0.8 | Low similarity → precedents less relevant |
| `dataset_size > 10,000` | β ×1.3 | Large dataset → trust statistics |
| `rule_confidence > 0.9` | γ ×2.0 | Exact match → trust chemistry |

**Example**:
```python
# Evidence quality:
n_prec = 10
diversity = 0.536
avg_sim = 0.500
dataset_size = 1,343

# Resulting weights:
α = 0.412  # Precedents (moderate - low similarity)
β = 0.441  # Analytics (slightly higher - validates precedents)
γ = 0.147  # Rules (low - no strong match)
δ = 0.000  # ML (none - no model available)
```

### Component Scoring

#### 1. Precedent Score (PS)

```python
def score_from_precedents(candidate, precedent_evidence):
    """
    Score based on precedent support.
    
    Combines:
        - Core match: Does precedent use same catalyst?
        - Base match: Does precedent use same base?
        - Solvent match: Does precedent use same solvent?
        - Weighted by DRFP similarity
    
    Returns: Float in [0, 1]
    """
```

**Example**:
- 6/10 precedents use Pd/PtBu3 (core match)
- 9/10 use same base
- 9/10 use same solvent
- Average similarity: 0.7
- **PS = 0.7 × (6/10 + 9/10 + 9/10) / 3 = 0.633**

#### 2. Analytics Score (AS)

```python
def score_from_analytics(candidate, analytics_evidence):
    """
    Score based on dataset-level statistics.
    
    Combines:
        - Frequency score: How common is this condition?
        - Yield score: How successful is this condition?
    
    Formula: 0.6 × freq_score + 0.4 × yield_score
    
    Returns: Float in [0, 1]
    """
```

**Example**:
- Pd(OAc)2 + ligand: 45 reactions, 90.1% yield (rank 1/20)
- NaOtBu: 493 reactions, 75.6% yield (rank 1/15)
- Dioxane: 200 reactions, 78.0% yield (rank 2/15)
- **AS = 0.6 × 0.95 + 0.4 × 0.81 = 0.894**

#### 3. Rule Score (RS)

```python
def score_from_rules(candidate, rule_evidence):
    """
    Score based on chemistry rules.
    
    Returns:
        - 1.0 if exact SMARTS match
        - 0.5-0.8 if partial match
        - 0.0 if no match
    """
```

#### 4. ML Score (MS)

```python
def score_from_ml(candidate, ml_evidence, reaction_smiles):
    """
    Score based on ML yield prediction.
    
    Returns: predicted_yield / 100.0
    """
```

### Score Fusion

```python
total_score = α×PS + β×AS + γ×RS + δ×MS
```

**Example**:
```
Candidate: Pd/PtBu3, NaOtBu, Toluene, 100°C, 12h

Component Scores:
  PS = 0.633 (6/10 precedents match)
  AS = 0.200 (less common in dataset)
  RS = 0.000 (no rule match)
  MS = 0.500 (neutral - no model)

Weights:
  α = 0.412
  β = 0.441
  γ = 0.147
  δ = 0.000

Total Score:
  = 0.412 × 0.633 + 0.441 × 0.200 + 0.147 × 0.000 + 0.000 × 0.500
  = 0.261 + 0.088 + 0.000 + 0.000
  = 0.349
  
Confidence: LOW (score < 0.5)
```

## Integration Steps

### Step 1: Replace Current Recommendation Logic

**Current** (`chemtools/recommend/core.py`):
```python
# Simple k-NN voting
core_counts = Counter([p.get("core") for p in precs])
chosen_core = max(core_counts, key=core_counts.get)
```

**New** (with fusion):
```python
from chemtools.ml.fusion_recommender import recommend_with_fusion

results = recommend_with_fusion(
    reaction_smiles=rxn_smiles_norm,
    family=fam,
    k=k,
    top_n=max_variants,
    relax=relax
)
```

### Step 2: Add Analytics Integration

Already implemented! ✅
- `chem.analytics.get_common_catalytic_systems()`
- `chem.analytics.get_common_bases()`
- `chem.analytics.get_common_solvents()`

### Step 3: Add Rule Integration (TODO)

```python
# In collect_evidence()
if family in RULE_DATABASES:
    from chemtools import rules
    db = rules.load_database(f"{family}_db.json")
    match_result = rules.match(db, reaction_smiles)
    
    evidence['rules'] = {
        'matched_schemes': match_result.schemes,
        'confidence': match_result.confidence,
        'recommended': match_result.conditions
    }
```

### Step 4: Add ML Integration (TODO)

```python
# In collect_evidence()
try:
    from chemtools.ml.drfp_predictor import DRFPYieldPredictor
    
    model_path = "models/drfp_yield_v1.pkl"
    predictor = DRFPYieldPredictor.load(model_path)
    
    evidence['ml'] = {
        'model_available': True,
        'model_confidence': 0.85,  # From validation
        'predictor': predictor
    }
except:
    evidence['ml'] = {'model_available': False}
```

## Testing Strategy

### Unit Tests

```python
def test_diversity_score():
    """Test precedent diversity calculation."""
    # Same conditions → low diversity
    precs_same = [
        {'core': 'Pd', 'base_uid': 'K2CO3', 'solvent_uid': 'Water', 'T_C': 80}
    ] * 10
    assert compute_diversity_score(precs_same) < 0.1
    
    # Different conditions → high diversity
    precs_diff = [
        {'core': f'Cat{i}', 'base_uid': f'Base{i}', 
         'solvent_uid': f'Solv{i}', 'T_C': 80 + i*30}
        for i in range(10)
    ]
    assert compute_diversity_score(precs_diff) > 0.8


def test_adaptive_weights_sparse_precedents():
    """Test weight adaptation for sparse precedents."""
    evidence = {
        'precedents': {
            'coverage': 5,  # Few precedents
            'diversity_score': 0.8,
            'avg_similarity': 0.5
        },
        'analytics': {'dataset_size': 50000},  # Large dataset
        'rules': {'confidence': 0.0},
        'ml': {'model_available': False}
    }
    
    weights = compute_adaptive_weights(evidence)['weights']
    
    # Should favor analytics over sparse precedents
    assert weights['β'] > weights['α']
    assert weights['β'] > 0.5  # Analytics should dominate


def test_adaptive_weights_strong_precedents():
    """Test weight adaptation for strong precedents."""
    evidence = {
        'precedents': {
            'coverage': 100,  # Many precedents
            'diversity_score': 0.7,  # Diverse
            'avg_similarity': 0.85  # Very similar
        },
        'analytics': {'dataset_size': 1000},
        'rules': {'confidence': 0.0},
        'ml': {'model_available': False}
    }
    
    weights = compute_adaptive_weights(evidence)['weights']
    
    # Should favor strong precedents
    assert weights['α'] > weights['β']
    assert weights['α'] > 0.45
```

### Integration Tests

```python
def test_fusion_vs_baseline():
    """Compare fusion recommendation with baseline."""
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    family = "C_N_Coupling_Pd"
    
    # Baseline (current system)
    from chemtools.recommend import recommend_from_reaction
    baseline = recommend_from_reaction(reaction, k=25)
    
    # Fusion (new system)
    fusion = recommend_with_fusion(reaction, family, k=25, top_n=5)
    
    # Both should return recommendations
    assert len(baseline['formatted']['recommended_conditions']) > 0
    assert len(fusion['recommended_conditions']) > 0
    
    # Fusion should provide component scores
    top = fusion['recommended_conditions'][0]
    assert 'component_scores' in top
    assert 'weights' in top
    assert 'reasoning' in top
```

### Benchmark Tests

```python
def test_benchmark_reactions():
    """Test on known reactions with expected outcomes."""
    test_cases = [
        {
            'reaction': "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            'family': "C_N_Coupling_Pd",
            'expected_core': "Pd",
            'expected_base': "NaOtBu",
            'expected_confidence': "HIGH"
        },
        # Add more test cases
    ]
    
    for case in test_cases:
        results = recommend_with_fusion(
            case['reaction'],
            case['family'],
            k=50,
            top_n=5
        )
        
        top = results['recommended_conditions'][0]
        
        # Check if expected conditions appear in top 5
        top5_cores = [r.candidate.core for r in results['recommended_conditions']]
        assert case['expected_core'] in top5_cores
```

## Performance Considerations

### Timing Analysis

```
collect_evidence():     ~0.5s
  - precedent.knn():    ~0.3s (DRFP search)
  - analytics calls:    ~0.2s (dataset loading)

compute_weights():      ~0.001s

score_candidates():     ~0.01s per candidate
  - 30 candidates:      ~0.3s

Total: ~0.8s (acceptable for interactive use)
```

### Optimization Strategies

1. **Cache Analytics Data**
   ```python
   # Cache analytics for each family
   _analytics_cache = {}
   
   def get_analytics_evidence(family):
       if family not in _analytics_cache:
           _analytics_cache[family] = {
               'systems': chem.analytics.get_common_catalytic_systems(...),
               ...
           }
       return _analytics_cache[family]
   ```

2. **Limit Candidate Generation**
   - Max 30-40 candidates total
   - From precedents: Top 10 unique conditions
   - From analytics: Top 3 systems × Top 2 bases × Top 2 solvents = 12
   - From rules: 1-2 exact matches

3. **Parallel Scoring**
   ```python
   from concurrent.futures import ThreadPoolExecutor
   
   with ThreadPoolExecutor(max_workers=4) as executor:
       scored = list(executor.map(
           lambda c: score_candidate(c, evidence, weights, rxn),
           candidates
       ))
   ```

## Deployment Plan

### Phase 1: Prototype (✅ Complete)
- [x] Implement fusion_recommender.py
- [x] Create demo script
- [x] Document architecture

### Phase 2: Integration (Next)
- [ ] Integrate with recommend/core.py
- [ ] Add rule-based evidence collection
- [ ] Add ML predictor integration
- [ ] Add analytics caching

### Phase 3: Testing (Future)
- [ ] Unit tests for all components
- [ ] Integration tests
- [ ] Benchmark on test reactions
- [ ] Compare with baseline system

### Phase 4: Validation (Future)
- [ ] Run on 100 test reactions
- [ ] Measure yield prediction accuracy
- [ ] Measure user satisfaction
- [ ] A/B testing if possible

### Phase 5: Production (Future)
- [ ] Make fusion the default method
- [ ] Add configuration options
- [ ] Update documentation
- [ ] Monitor performance

## Configuration Options

```python
# Future: Allow users to configure fusion behavior

FUSION_CONFIG = {
    'precedent_threshold': 10,   # Min precedents to trust them
    'diversity_threshold': 0.3,   # Min diversity for healthy precedents
    'analytics_weight_boost': 1.3,  # Boost for large datasets
    'rule_weight_boost': 2.0,     # Boost for exact matches
    'ml_confidence_threshold': 0.8,  # Min ML confidence to use
    'max_candidates': 30,         # Max candidates to score
}
```

## Monitoring & Logging

```python
# Add logging to track fusion decisions

import logging

logger = logging.getLogger('fusion_recommender')

def recommend_with_fusion(...):
    logger.info(f"Fusion recommendation for {family}")
    
    evidence = collect_evidence(...)
    logger.info(f"Evidence: {evidence['precedents']['coverage']} prec, "
                f"diversity={evidence['precedents']['diversity_score']:.2f}")
    
    weights = compute_adaptive_weights(evidence)
    logger.info(f"Adaptive weights: α={weights['weights']['α']:.3f}, "
                f"β={weights['weights']['β']:.3f}")
    
    # ... rest of function
```

## Conclusion

The multi-source evidence fusion system:

✅ **Solves the top-k precedent problem** by balancing with dataset analytics

✅ **Adapts to data quality** with dynamic weighting

✅ **Provides interpretability** through component scores and reasoning

✅ **Prevents bad recommendations** by validating precedents against dataset patterns

✅ **Ready for integration** with existing pipeline

**Next Action**: Integrate `recommend_with_fusion()` into `chemtools/recommend/core.py` as an optional recommendation mode.
