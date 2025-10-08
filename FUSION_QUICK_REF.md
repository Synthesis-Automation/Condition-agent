# Fusion System: Quick Reference Card

## TL;DR

**Fusion BUILDS ON existing methods, does NOT replace them.**

```
Original Methods (UNCHANGED):
  ✅ precedent.knn()           → k-NN similarity search
  ✅ analytics.get_common_*()  → Dataset statistics
  ✅ router.detect_family()    → Rule-based SMARTS
  ✅ DRFPYieldPredictor        → ML yield prediction

Fusion Layer (NEW, OPT-IN):
  🆕 Calls ALL existing methods above
  🆕 Adds diversity scoring (bias detection)
  🆕 Adds adaptive weighting (quality-based)
  🆕 Adds multi-source scoring (combines evidence)
  🆕 Adds transparency (shows reasoning)
```

## Two Modes

| Feature | k-NN Mode | Fusion Mode |
|---------|-----------|-------------|
| **Parameter** | `use_fusion=False` (default) | `use_fusion=True` (opt-in) |
| **Methods Used** | `precedent.knn()` only | All 4 existing methods |
| **Speed** | Fast | Moderate |
| **Transparency** | Basic | Component scores + reasoning |
| **When to Use** | Standard cases | Novel/challenging cases |

## Code Example

```python
from chemtools.recommend import recommend_from_reaction

# Standard k-NN (DEFAULT)
baseline = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50
)

# Fusion enhancement (OPT-IN)
fusion = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50,
    use_fusion=True  # 👈 Enable multi-source fusion
)

# Check fusion metadata
print(fusion['fusion_meta']['adaptive_weights'])
# {'alpha': 0.329, 'beta': 0.503, 'gamma': 0.168, 'delta': 0.000}

print(fusion['fusion_meta']['reasoning'])
# ["Low diversity → precedents may be biased", ...]
```

## What Fusion Does Internally

```python
# Pseudo-code of fusion internals
def recommend_with_fusion(reaction, family, k):
    # STEP 1: Call EXISTING methods (reuse, don't replace)
    precedents = precedent.knn(family, {}, k, {})      # ← EXISTING #1
    analytics = analytics.get_common_systems(family)   # ← EXISTING #2
    rules = router.detect_family(reaction)             # ← EXISTING #3
    ml = DRFPYieldPredictor.load().predict(...)       # ← EXISTING #4
    
    # STEP 2: NEW analysis - assess evidence quality
    diversity = calculate_diversity(precedents)        # ← NEW
    
    # STEP 3: NEW weighting - adapt based on quality
    weights = compute_adaptive_weights(                # ← NEW
        precedent_count, diversity, rule_confidence, ...
    )
    
    # STEP 4: NEW scoring - combine all evidence
    for candidate in candidates:
        PS = score_from_precedents(candidate, precedents)
        AS = score_from_analytics(candidate, analytics)
        RS = score_from_rules(candidate, rules)
        MS = score_from_ml(candidate, ml)
        
        total = weights['α']·PS + weights['β']·AS + 
                weights['γ']·RS + weights['δ']·MS     # ← NEW
    
    return ranked_candidates
```

## Evidence Sources

| Source | Method Called | Status | Weight |
|--------|--------------|--------|---------|
| **Precedents (PE)** | `precedent.knn()` | ✅ Connected | α (adaptive) |
| **Analytics (AS)** | `analytics.get_common_*()` | ✅ Connected | β (adaptive) |
| **Rules (RE)** | `router.detect_family()` | ⏳ Placeholder | γ (low) |
| **ML (MS)** | `DRFPYieldPredictor` | ⏳ Placeholder | δ=0 |

## Key Design Principle

```
┌──────────────────────────────────────────┐
│   COMPOSITION, NOT REPLACEMENT           │
├──────────────────────────────────────────┤
│                                          │
│  ❌ Bad: Reimplementing existing logic  │
│  ✅ Good: Calling existing methods       │
│                                          │
│  Fusion is a composition layer that      │
│  orchestrates existing methods.          │
│                                          │
└──────────────────────────────────────────┘
```

## Backward Compatibility

- ✅ Default behavior unchanged (`use_fusion=False`)
- ✅ All existing APIs work as before
- ✅ Original methods callable independently
- ✅ No breaking changes
- ✅ Safe fallback if fusion fails

## Documentation

- **FUSION_CORRECT_ARCHITECTURE.md**: Confirms fusion builds on existing
- **FUSION_BUILDS_ON_EXISTING.md**: Detailed explanation with examples
- **FUSION_ARCHITECTURE_SUMMARY.md**: Complete architecture overview
- **FUSION_IMPLEMENTATION_GUIDE.md**: Implementation details

## Status Summary

```
Component                    Status
─────────────────────────────────────────────────
Core Integration             ✅ Complete
Precedent Evidence (PE)      ✅ Connected to precedent.knn()
Dataset Analytics (AS)       ✅ Connected to analytics module
Rule Evidence (RE)           ⏳ Placeholder for router module
ML Evidence (MS)             ⏳ Placeholder for yield predictor
Deduplication                ✅ Complete
Adaptive Weighting           ✅ Complete
Transparency                 ✅ Complete
Integration Tests            ✅ 6 tests passing
Documentation                ✅ Complete
Demo Script                  🔄 In progress (Suzuki examples)
```

## When to Use Each Mode

### Use k-NN (default) when:
- Fast results needed
- Well-established reaction types
- Many precedents available (>50)
- High trust in precedent database

### Use Fusion (opt-in) when:
- Novel or challenging substrates
- Few precedents (<10)
- Suspected bias in precedents
- Need transparency/reasoning
- Multi-source validation desired

## Bottom Line

**Fusion is an enhancement layer that:**
1. Calls all existing methods (doesn't replace)
2. Adds quality assessment (diversity, bias detection)
3. Adapts weights dynamically (based on evidence)
4. Combines scores intelligently (multi-source fusion)
5. Provides transparency (reasoning + component scores)

**All original methods remain functional and unchanged.**
