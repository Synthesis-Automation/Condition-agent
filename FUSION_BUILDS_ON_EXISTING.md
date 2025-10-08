# Fusion System: Building on Top of Existing Methods

## Visual Overview

This document clarifies that the fusion system **BUILDS ON TOP** of existing ML and rule-based methods, rather than replacing them.

## The Three Modes of Operation

```
┌──────────────────────────────────────────────────────────────────┐
│                        USER CHOOSES MODE                         │
└──────────────────────────────────────────────────────────────────┘
                                 │
                 ┌───────────────┼───────────────┐
                 │               │               │
                 ▼               ▼               ▼
         ┌───────────┐   ┌───────────┐   ┌───────────┐
         │  k-NN     │   │   Rules   │   │  Fusion   │
         │  Only     │   │   Only    │   │ (Combined)│
         └───────────┘   └───────────┘   └───────────┘
              │               │               │
              │               │               │
    use_fusion=False    (separate)    use_fusion=True
         (DEFAULT)        module       (ENHANCEMENT)
```

## Current Implementation Status

### ✅ What Exists and Works Independently

```python
# 1. k-NN Precedent Search (EXISTING - chemtools/precedent.py)
def knn(family, features, k, relax):
    """
    DRFP-based k-nearest neighbors search.
    Returns k most similar precedents.
    """
    # Load DRFP fingerprints
    # Compute similarity
    # Return top-k precedents
    pass

# 2. Rule-Based Matching (EXISTING - chemtools/router.py)
def detect_family(reaction_smiles):
    """
    SMARTS pattern-based family detection.
    Returns matched reaction family and confidence.
    """
    # Match SMARTS patterns
    # Return family + confidence
    pass

# 3. Dataset Analytics (EXISTING - chemtools/analytics.py)
def get_common_catalytic_systems(family):
    """
    Statistical analysis of dataset.
    Returns frequency distribution of conditions.
    """
    # Analyze dataset
    # Return statistics
    pass

# 4. ML Yield Predictor (EXISTING - chemtools/ml/drfp_yield_predictor.py)
class DRFPYieldPredictor:
    """
    Machine learning model for yield prediction.
    Uses DRFP features to predict reaction yield.
    """
    def predict(self, reaction_smiles):
        # DRFP encoding
        # ML model inference
        return predicted_yield
```

### 🆕 What Fusion Adds (NEW Layer)

```python
# Fusion Recommender (NEW - chemtools/ml/fusion_recommender.py)
def recommend_with_fusion(reaction_smiles, family, k, top_n, relax):
    """
    Multi-source evidence fusion that BUILDS ON existing methods.
    
    Does NOT replace anything—just combines them intelligently.
    """
    
    # ============================================================
    # STEP 1: Gather Evidence from EXISTING Methods
    # ============================================================
    
    # 1a. Call EXISTING precedent search
    from chemtools import precedent
    precedent_data = precedent.knn(
        family=family,
        features={},
        k=k,
        relax=relax
    )  # ← SAME function as baseline mode!
    
    # 1b. Call EXISTING analytics
    from chemtools.analytics import get_common_catalytic_systems
    analytics_data = get_common_catalytic_systems(family)
    # ← SAME function used elsewhere!
    
    # 1c. Call EXISTING rule matcher
    from chemtools.router import detect_family
    rule_data = detect_family(reaction_smiles)
    # ← SAME function as baseline!
    
    # 1d. Call EXISTING ML predictor (if available)
    try:
        from chemtools.ml.drfp_yield_predictor import DRFPYieldPredictor
        predictor = DRFPYieldPredictor.load()
        ml_data = predictor.predict(candidates)
        # ← SAME model used elsewhere!
    except:
        ml_data = None  # Graceful fallback
    
    # ============================================================
    # STEP 2: NEW Analysis - Assess Evidence Quality
    # ============================================================
    
    # NEW: Calculate diversity of precedents
    diversity = calculate_diversity_score(precedent_data['precedents'])
    
    # NEW: Detect if precedents are biased (e.g., all from same batch)
    if diversity < 0.3:
        warning = "Low diversity - precedents may be biased"
    
    # NEW: Assess dataset size
    dataset_size = analytics_data.get('total_reactions', 0)
    
    # NEW: Check rule confidence
    rule_confidence = rule_data.get('confidence', 0.0)
    
    # ============================================================
    # STEP 3: NEW Scoring - Adaptive Weight Adjustment
    # ============================================================
    
    # NEW: Dynamic weight computation based on evidence quality
    weights = compute_adaptive_weights(
        precedent_count=len(precedent_data['precedents']),
        diversity=diversity,
        avg_similarity=precedent_data['avg_similarity'],
        dataset_size=dataset_size,
        rule_confidence=rule_confidence,
        ml_available=(ml_data is not None)
    )
    
    # Example output:
    # weights = {
    #     'alpha': 0.329,  # Precedents (reduced due to low diversity)
    #     'beta': 0.503,   # Analytics (increased to compensate)
    #     'gamma': 0.168,  # Rules (moderate confidence)
    #     'delta': 0.000   # ML (not available)
    # }
    
    # ============================================================
    # STEP 4: NEW Fusion - Combine Evidence Sources
    # ============================================================
    
    candidates = build_candidate_list(
        precedent_data,
        analytics_data,
        rule_data
    )
    
    for candidate in candidates:
        # NEW: Multi-source scoring
        PS = score_from_precedents(candidate, precedent_data)
        AS = score_from_analytics(candidate, analytics_data)
        RS = score_from_rules(candidate, rule_data)
        MS = score_from_ml(candidate, ml_data) if ml_data else 0.5
        
        # NEW: Weighted combination
        total_score = (
            weights['alpha'] * PS +
            weights['beta'] * AS +
            weights['gamma'] * RS +
            weights['delta'] * MS
        )
        
        # NEW: Store component scores for transparency
        candidate['score'] = total_score
        candidate['component_scores'] = {
            'PS': PS, 'AS': AS, 'RS': RS, 'MS': MS
        }
    
    # ============================================================
    # STEP 5: NEW Output - Enhanced Metadata
    # ============================================================
    
    return {
        'recommended_conditions': sorted(candidates, key=lambda c: c['score'], reverse=True),
        'evidence': {
            'precedents': precedent_data,  # From EXISTING method
            'analytics': analytics_data,    # From EXISTING method
            'rules': rule_data,            # From EXISTING method
            'ml': ml_data                  # From EXISTING method (if available)
        },
        'evidence_summary': {
            'precedent_count': len(precedent_data['precedents']),
            'diversity_score': diversity,  # NEW metric
            'dataset_size': dataset_size,
            'rule_confidence': rule_confidence
        },
        'adaptive_weights': weights,  # NEW: shows how evidence was weighted
        'reasoning': [  # NEW: explains weight adjustments
            f"Low diversity ({diversity:.2f}) → precedents may be biased",
            f"Low similarity ({avg_sim:.2f}) → precedents less relevant",
            "No ML model available → δ=0"
        ]
    }
```

## How It All Fits Together

### Mode 1: Standard k-NN (Original, DEFAULT)

```python
# User code
from chemtools.recommend import recommend_from_reaction

results = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50
    # use_fusion=False is the DEFAULT
)

# What happens internally:
# 1. precedent.knn() → get top-50 similar reactions
# 2. Vote for catalyst core
# 3. Vote for base/solvent
# 4. Apply constraints
# 5. Return recommendation

# Result:
# {
#   'recommendation': {'core': 'Pd/XPhos', 'base': 'K3PO4', ...},
#   'alternatives': {...},
#   'precedent_pack': {...},
#   'formatted': {...}
# }
```

**Uses**: `precedent.knn()` only  
**Fast**: Single method call  
**Proven**: Existing, well-tested approach

### Mode 2: Fusion Enhancement (NEW, OPT-IN)

```python
# User code  
from chemtools.recommend import recommend_from_reaction

results = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50,
    use_fusion=True  # 👈 OPT-IN to fusion
)

# What happens internally:
# 1. precedent.knn() → get top-50 similar reactions [SAME as Mode 1]
# 2. analytics.get_common_catalytic_systems() → get dataset stats [EXISTING]
# 3. router.detect_family() → get rule matches [EXISTING]
# 4. DRFPYieldPredictor.predict() → get ML predictions [EXISTING, if available]
# 5. Calculate diversity of precedents [NEW]
# 6. Compute adaptive weights [NEW]
# 7. Score candidates with multi-source evidence [NEW]
# 8. Return enhanced results

# Result:
# {
#   'recommendation': {'core': 'Pd/XPhos', 'base': 'K3PO4', ...},
#   'fusion_meta': {  # 👈 NEW metadata
#       'adaptive_weights': {'alpha': 0.329, 'beta': 0.503, ...},
#       'evidence_summary': {'precedents': 10, 'diversity': 0.087, ...},
#       'reasoning': ["Low diversity → precedents may be biased", ...]
#   },
#   'formatted': {
#       'recommended_conditions': [
#           {
#               'summary': {...},
#               'component_scores': {  # 👈 NEW transparency
#                   'PS': 0.933,  # Precedent score
#                   'AS': 0.354,  # Analytics score
#                   'RS': 0.000,  # Rule score
#                   'MS': 0.500   # ML score
#               }
#           }
#       ]
#   }
# }
```

**Uses**: All existing methods + fusion logic  
**Comprehensive**: Multi-source validation  
**Transparent**: Shows component scores and reasoning  
**Safe**: Falls back to k-NN if fusion fails

### Mode 3: Individual Methods (Can Still Use Directly)

```python
# Users can STILL call individual methods directly if needed

# Just precedent search
from chemtools import precedent
precs = precedent.knn(family="Suzuki", features={}, k=50)

# Just analytics
from chemtools.analytics import get_common_catalytic_systems
stats = get_common_catalytic_systems("Suzuki")

# Just rules
from chemtools.router import detect_family
family = detect_family("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")

# Just ML
from chemtools.ml.drfp_yield_predictor import DRFPYieldPredictor
predictor = DRFPYieldPredictor.load()
yield_pred = predictor.predict("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")
```

**All original methods remain accessible and unchanged!**

## Key Takeaway

```
┌─────────────────────────────────────────────────────────────────┐
│                     FUSION DOES NOT REPLACE                     │
│                                                                 │
│  Instead, it BUILDS ON TOP by:                                 │
│  1. Calling all existing methods                               │
│  2. Assessing evidence quality (NEW)                           │
│  3. Adapting weights dynamically (NEW)                         │
│  4. Combining scores intelligently (NEW)                       │
│  5. Providing transparency (NEW)                               │
│                                                                 │
│  Original methods are:                                         │
│  ✅ Unchanged                                                  │
│  ✅ Still independently usable                                 │
│  ✅ Still the default (k-NN mode)                              │
│                                                                 │
│  Fusion is:                                                    │
│  🆕 Optional enhancement (opt-in)                              │
│  🆕 Combines existing methods                                  │
│  🆕 Adds quality assessment                                    │
│  🆕 Provides transparency                                      │
└─────────────────────────────────────────────────────────────────┘
```

## Recommended Usage

### When to use Standard k-NN (use_fusion=False)
- Fast recommendations needed
- Well-established reaction types with many precedents
- High-throughput screening
- When you trust the precedent database

### When to use Fusion (use_fusion=True)
- Novel or challenging substrates
- When precedent count is low (<10)
- When you suspect bias in precedents
- When you want transparency and reasoning
- When you want to validate recommendations with multiple evidence sources

### Example Decision Tree

```
Is this a novel/challenging reaction?
│
├─ No → Use standard k-NN (fast, proven)
│
└─ Yes → Are precedents available?
    │
    ├─ Many (>50) with good diversity → k-NN is fine
    │
    └─ Few (<10) or low diversity → Use fusion
        │
        └─ Fusion will:
            ├─ Detect bias (diversity scoring)
            ├─ Validate with dataset stats
            ├─ Incorporate chemistry rules
            └─ Show you the reasoning
```

## Summary Table

| Feature | k-NN Mode | Fusion Mode |
|---------|-----------|-------------|
| **Speed** | Fast (single method) | Moderate (multiple methods) |
| **Precedent Search** | ✅ precedent.knn() | ✅ precedent.knn() (SAME) |
| **Dataset Analytics** | ❌ Not used | ✅ analytics.get_common_*() |
| **Rule Matching** | ❌ Not used | ✅ router.detect_family() |
| **ML Predictions** | ❌ Not used | ✅ DRFPYieldPredictor (if available) |
| **Bias Detection** | ❌ No | ✅ Diversity scoring (NEW) |
| **Adaptive Weights** | ❌ No | ✅ Quality-based (NEW) |
| **Transparency** | ⚠️ Limited | ✅ Component scores + reasoning |
| **Default** | ✅ Yes (use_fusion=False) | ❌ No (opt-in) |
| **Backward Compatible** | ✅ Original behavior | ✅ Non-invasive addition |
| **Fallback** | N/A | ✅ Falls back to k-NN if error |

## Conclusion

**The fusion system is an enhancement layer that preserves all original methods:**

1. ✅ **k-NN remains the default** (fast, proven, unchanged)
2. ✅ **All existing methods work independently** (precedent, analytics, router, ML)
3. ✅ **Fusion is opt-in** (add `use_fusion=True` to try it)
4. ✅ **No breaking changes** (backward compatible)
5. ✅ **Safe fallback** (if fusion fails, reverts to k-NN)

Users have **three options**:
- Use k-NN alone (default, fast)
- Use individual methods directly (precedent, analytics, rules, ML)
- Use fusion (combines all methods intelligently)

**All methods coexist peacefully—fusion just orchestrates them together.**
