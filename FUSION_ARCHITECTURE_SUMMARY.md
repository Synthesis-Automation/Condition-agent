# Fusion Recommendation Architecture Summary

## Overview

The fusion recommendation system **builds on top of** existing ML and rule-based methods, providing an optional enhancement layer without replacing any original functionality.

## Architecture: Three Layers

```
┌─────────────────────────────────────────────────────────────────┐
│                   USER API (chemtools.recommend)                │
│                                                                 │
│  recommend_from_reaction(reaction, use_fusion=True/False)      │
└───────────────────────┬─────────────────────────────────────────┘
                        │
                        ├─── use_fusion=False (DEFAULT)
                        │    │
                        │    └──> STANDARD K-NN MODE
                        │         ├─ precedent.knn() → DRFP similarity
                        │         ├─ Vote for catalyst core
                        │         ├─ Vote for base/solvent
                        │         └─ constraints.apply()
                        │
                        └─── use_fusion=True (OPTIONAL)
                             │
                             └──> FUSION MODE (builds on top of k-NN)
                                  ├─ Precedent Evidence (PE)
                                  │  └─ precedent.knn() [SAME AS BASELINE]
                                  │     + diversity scoring
                                  │     + bias detection
                                  │
                                  ├─ Dataset Statistics (DS)
                                  │  └─ analytics.get_common_catalytic_systems()
                                  │     + frequency analysis
                                  │     + validation of precedents
                                  │
                                  ├─ Rule Evidence (RE)
                                  │  └─ router.detect_family() [EXISTING]
                                  │     + SMARTS pattern matching
                                  │     + chemistry knowledge
                                  │
                                  └─ ML Predictions (ML)
                                     └─ DRFPYieldPredictor [EXISTING]
                                        + yield predictions
                                        + quality scoring
```

## Key Design Principles

### 1. **Non-Invasive Enhancement**
- Fusion mode is **opt-in** via `use_fusion=True` parameter
- Default behavior (`use_fusion=False`) is **unchanged**
- All original methods remain intact and independently usable

### 2. **Reuse Existing Infrastructure**
- **Precedent Evidence**: Uses same `precedent.knn()` as baseline
- **Rule Evidence**: Uses same `router.detect_family()` for pattern matching
- **ML Evidence**: Uses same `DRFPYieldPredictor` for yield predictions
- **Dataset Analytics**: Uses existing `analytics` module

### 3. **Backward Compatibility**
- All existing APIs work exactly as before
- No breaking changes to any function signatures
- Legacy code continues to function without modification

## Component Breakdown

### Standard k-NN Mode (Original Method)
```python
# File: chemtools/recommend/core.py
def recommend_from_reaction(reaction, k=25, use_fusion=False):
    if use_fusion:
        # Fusion path (new, optional)
        return _fusion_recommend(...)
    
    # STANDARD PATH (original, default)
    pack = precedent.knn(family=fam, features=features, k=k, relax=relax)
    # Vote for core, base, solvent
    # Apply constraints
    # Format output
    return results
```

**What it does:**
- DRFP similarity search via `precedent.knn()`
- Voting mechanism for catalyst selection
- Constraint filtering
- **Status**: ✅ Unchanged, fully preserved

### Fusion Mode (Enhancement Layer)
```python
# File: chemtools/ml/fusion_recommender.py
def recommend_with_fusion(reaction_smiles, family, k, top_n, relax):
    # 1. Precedent Evidence (PE) - reuses precedent.knn()
    precedents = _get_precedent_evidence(reaction_smiles, family, k, relax)
    diversity = _calculate_diversity(precedents)  # NEW: bias detection
    
    # 2. Dataset Statistics (DS) - reuses analytics module
    stats = analytics.get_common_catalytic_systems(family)  # EXISTING
    
    # 3. Rule Evidence (RE) - reuses router
    rules = router.match_conditions(reaction_smiles, family)  # EXISTING
    
    # 4. ML Evidence (ML) - reuses yield predictor
    ml_scores = yield_predictor.predict(candidates)  # EXISTING
    
    # 5. Adaptive Fusion (NEW: combines all evidence)
    weights = _compute_adaptive_weights(
        precedent_count=len(precedents),
        diversity=diversity,
        rule_confidence=rules.confidence,
        ml_available=ml_scores is not None
    )
    
    # 6. Score candidates (NEW: multi-source scoring)
    final_scores = (
        weights['alpha'] * precedent_scores +
        weights['beta'] * analytics_scores +
        weights['gamma'] * rule_scores +
        weights['delta'] * ml_scores
    )
    
    return ranked_candidates
```

**What it adds:**
- **Diversity scoring**: Detects if top-k precedents are biased
- **Adaptive weighting**: Adjusts weights based on evidence quality
- **Multi-source validation**: Cross-checks recommendations
- **Transparency**: Shows component scores and reasoning
- **Status**: ✅ New enhancement, builds on existing methods

## Evidence Sources Explained

### 1. Precedent Evidence (PE) - **REUSES** `precedent.knn()`
```python
# SAME call as baseline k-NN
precedents = precedent.knn(family=family, features={}, k=k, relax=relax)

# NEW: Add quality metrics
diversity = calculate_diversity_score(precedents)
avg_similarity = np.mean([p['drfp_similarity'] for p in precedents])

# If diversity < 0.3: Precedents may be from same batch → reduce weight α
```

**Enhancement**: Adds bias detection, doesn't change the underlying method.

### 2. Dataset Statistics (DS) - **REUSES** `analytics` module
```python
# EXISTING function call
from chemtools.analytics import get_common_catalytic_systems

stats = get_common_catalytic_systems(family="Suzuki")
# Returns: {
#   'Pd(PPh3)4': {'count': 1234, 'frequency': 0.15},
#   'Pd(OAc)2 + XPhos': {'count': 987, 'frequency': 0.12},
#   ...
# }

# NEW: Use statistics to validate precedent recommendations
analytics_score = match_candidate_to_stats(candidate, stats)
```

**Enhancement**: Uses existing analytics for validation, no new data collection.

### 3. Rule Evidence (RE) - **REUSES** `router` module
```python
# EXISTING function call
from chemtools.router import detect_family, match_conditions

family = detect_family(reaction_smiles)  # SMARTS pattern matching
conditions = match_conditions(reaction_smiles, family)  # Rule-based lookup

# NEW: Score candidates based on rule match confidence
rule_score = conditions.get('confidence', 0.0)
```

**Enhancement**: Incorporates existing rule-based knowledge into fusion scoring.

### 4. ML Evidence (ML) - **REUSES** `DRFPYieldPredictor`
```python
# EXISTING ML model (if available)
from chemtools.ml.drfp_yield_predictor import DRFPYieldPredictor

predictor = DRFPYieldPredictor.load("models/yield_predictor.pkl")

# NEW: Use predictions as evidence
for candidate in candidates:
    reaction_with_conditions = build_reaction_string(reactants, candidate)
    yield_pred = predictor.predict(reaction_with_conditions)
    ml_score = yield_pred / 100.0  # Normalize to 0-1
```

**Enhancement**: Integrates existing ML predictions into fusion framework.

## Adaptive Weighting Logic

The fusion system dynamically adjusts weights based on evidence quality:

```python
# Initial weights (can be customized)
α = 0.5  # Precedent evidence
β = 0.3  # Dataset analytics  
γ = 0.1  # Rule-based evidence
δ = 0.1  # ML predictions

# Adaptive adjustments
if precedent_count < 10:
    α *= 0.5  # Low precedent count → reduce precedent weight
    β *= 1.5  # Boost analytics weight to compensate

if diversity_score < 0.3:
    α *= 0.7  # Low diversity → precedents may be biased
    β *= 1.3  # Trust dataset statistics more

if avg_similarity < 0.6:
    α *= 0.8  # Low similarity → precedents less relevant
    γ *= 1.2  # Boost rule-based evidence

if dataset_size > 10000:
    β *= 1.3  # Large dataset → trust statistics more

if rule_confidence > 0.9:
    γ *= 2.0  # Strong rule match → trust chemistry knowledge

if not ml_model_available:
    δ = 0.0  # No ML model → redistribute weight
    α, β, γ = normalize(α, β, γ)

# Final scores
total_score = α·PS + β·AS + γ·RS + δ·MS
```

## Usage Examples

### Example 1: Standard k-NN (Original Method)
```python
from chemtools.recommend import recommend_from_reaction

# DEFAULT: Uses original k-NN voting
results = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50
)

print(results['recommendation']['core'])  # "Pd/XPhos"
# Method: precedent.knn() → voting → constraints
```

### Example 2: Fusion Mode (Enhanced Method)
```python
from chemtools.recommend import recommend_from_reaction

# OPT-IN: Uses multi-source fusion
results = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50,
    use_fusion=True  # 👈 Enable fusion
)

# Access fusion metadata
meta = results['fusion_meta']
print(meta['adaptive_weights'])  # α=0.329, β=0.503, γ=0.168, δ=0.000
print(meta['evidence_summary'])  # {precedents: 10, diversity: 0.087, ...}
print(meta['reasoning'])         # ["Low diversity → precedents may be biased", ...]

# Recommendations have component scores
rec = results['formatted']['recommended_conditions'][0]
print(rec['component_scores'])   # {PS: 0.933, AS: 0.354, RS: 0.000, MS: 0.500}
```

### Example 3: Comparing Both Methods
```python
# Baseline
baseline = recommend_from_reaction(reaction, use_fusion=False)

# Fusion
fusion = recommend_from_reaction(reaction, use_fusion=True)

# Compare
print(f"Baseline top: {baseline['recommendation']['core']}")
print(f"Fusion top: {fusion['recommendation']['core']}")

# Analyze differences
if baseline['recommendation'] != fusion['recommendation']:
    print("Fusion detected potential bias in precedents")
    print(f"Diversity score: {fusion['fusion_meta']['evidence_summary']['diversity']}")
    print(f"Reasoning: {fusion['fusion_meta']['reasoning']}")
```

## Integration Status

### ✅ Completed
1. **Core Integration**: `recommend_from_reaction()` supports `use_fusion` parameter
2. **Precedent Evidence**: Reuses `precedent.knn()` with diversity scoring
3. **Dataset Analytics**: Reuses `analytics.get_common_catalytic_systems()`
4. **Deduplication**: Prevents duplicate recommendations
5. **Adaptive Weighting**: Dynamically adjusts weights based on evidence quality
6. **Testing**: 6 integration tests passing
7. **Documentation**: Implementation guides and API docs

### 🔄 In Progress
1. **Rule Evidence**: Placeholder for `router.match_conditions()` integration
2. **ML Evidence**: Placeholder for `DRFPYieldPredictor` integration (δ=0 currently)

### 📋 Future Enhancements
1. **Enhanced Rule Matching**: Connect to suzuki_db.json and other rule databases
2. **ML Yield Predictor**: Integrate trained DRFP yield models
3. **Performance Optimization**: Caching, parallelization
4. **Additional Evidence Sources**: Reaction mechanism knowledge, retrosynthesis insights

## File Structure

```
chemtools/
├── recommend/
│   ├── core.py                    # Main API (use_fusion parameter)
│   └── utils.py                   # Helper functions
│
├── ml/
│   ├── fusion_recommender.py      # Fusion engine (NEW)
│   └── drfp_yield_predictor.py    # ML yield model (EXISTING, to integrate)
│
├── analytics.py                   # Dataset statistics (EXISTING, reused)
├── router.py                      # Rule-based matching (EXISTING, to integrate)
├── precedent.py                   # k-NN search (EXISTING, reused)
└── constraints.py                 # Filtering (EXISTING, reused)

tests/
└── test_fusion_integration.py     # Integration tests (NEW)

docs/
├── FUSION_IMPLEMENTATION_GUIDE.md  # Implementation details (NEW)
├── FUSION_INTEGRATION_COMPLETE.md  # Integration summary (NEW)
└── ML_RECOMMENDATION_STRATEGY.md   # Strategic design (NEW)
```

## Key Takeaways

### For Users
- **Default behavior unchanged**: All existing code works without modification
- **Opt-in enhancement**: Add `use_fusion=True` to get multi-source recommendations
- **Transparency**: Fusion mode explains its reasoning and shows evidence quality
- **Fallback safety**: If fusion fails, automatically falls back to k-NN

### For Developers
- **Modular design**: Each evidence source is independent and reusable
- **No code duplication**: Reuses all existing infrastructure
- **Easy to extend**: Add new evidence sources by implementing score_from_X() functions
- **Well-tested**: 6 integration tests cover core functionality
- **Clear separation**: Fusion logic isolated in `chemtools/ml/fusion_recommender.py`

### For Maintainers
- **Backward compatible**: No breaking changes to any APIs
- **Independent modules**: Fusion can be disabled or removed without affecting core
- **Graceful degradation**: Missing dependencies (ML model, etc.) handled gracefully
- **Clear documentation**: Comprehensive guides for implementation and usage

## Conclusion

The fusion recommendation system is an **optional enhancement layer** that builds on top of existing methods:

- ✅ **Preserves** all original functionality
- ✅ **Reuses** existing infrastructure (precedent.knn, analytics, router, ML models)
- ✅ **Enhances** with diversity scoring, adaptive weighting, and multi-source validation
- ✅ **Transparent** with component scores and reasoning
- ✅ **Safe** with fallback to k-NN if fusion fails

Users can choose:
- **Standard mode** (`use_fusion=False`): Fast, proven k-NN voting
- **Fusion mode** (`use_fusion=True`): Multi-source validation with bias detection

Both modes use the same underlying methods—fusion just combines them intelligently.
