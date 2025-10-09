# All Recommendation Functions - Complete Overview

## Summary: You Have 5 Different Recommendation Functions

Here's what each one does and when to use it:

---

## 1. **`recommend_from_reaction()`** ⭐ MAIN ENTRY POINT
**Location**: `chemtools/recommend/core.py`

**What it is**: The **primary public API** - routes to other methods based on settings

```python
from chemtools.recommend import recommend_from_reaction

result = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50,
    use_fusion=False  # ← Controls which method is used
)
```

**Workflow**:
```
recommend_from_reaction()
    ↓
    if use_fusion == True:
        → calls recommend_with_fusion()  (Complex multi-source)
    else:
        → uses simple precedent-based logic (Default)
```

**When to use**: This is the **main function** you should call from your code.

---

## 2. **`recommend_with_fusion()`** 🔬 COMPLEX FUSION
**Location**: `chemtools/ml/fusion_recommender.py`

**What it is**: Multi-source evidence fusion with adaptive weighting

```python
from chemtools.ml.fusion_recommender import recommend_with_fusion

result = recommend_with_fusion(
    reaction_smiles=rxn,
    family="Suzuki",
    k=50,
    use_rule_alignment=True  # Default: True
)
```

**Workflow**:
```
1. Collect evidence from 4 sources:
   - Precedents (DRFP k-NN)
   - Analytics (dataset statistics)
   - Rules (SMARTS matching)
   - ML (yield prediction models)

2. Compute adaptive weights (α, β, γ, δ)

3. Score candidates: α×precedent + β×analytics + γ×rules + δ×ML

4. Rerank by rule-alignment (if enabled)
```

**Complexity**: ❌ HIGH (10+ functions, adaptive weighting, multi-pass ranking)

**When to use**: 
- Called internally by `recommend_from_reaction(use_fusion=True)`
- Don't call directly unless you need full fusion complexity

---

## 3. **`recommend_simple()`** ✅ NEW SIMPLE RANKER
**Location**: `chemtools/ml/simple_precedent_ranker.py`

**What it is**: Simple precedent-centric ranking with optional reranking

```python
from chemtools.ml.simple_precedent_ranker import recommend_simple

result = recommend_simple(
    reaction_smiles=rxn,
    family="Suzuki",
    k=30,
    rerank_strategy='rule'  # 'rule', 'analytics', or 'none'
)
```

**Workflow**:
```
1. Find k similar precedents (DRFP k-NN)

2. Optionally rerank:
   - 'rule': Boost by rule match
   - 'analytics': Boost by popularity
   - 'none': Use similarity only

3. Extract top conditions from reranked precedents
```

**Complexity**: ✅ LOW (3 simple functions, transparent)

**When to use**:
- **Recommended!** When you want simple, transparent ranking
- When precedents are primary but need quality filtering
- When you want to choose rule vs analytics reranking

---

## 4. **`recommend_with_yield_prediction()`** 📊 ML PREDICTOR
**Location**: `chemtools/ml/recommender.py`

**What it is**: Uses ML models to predict yields for candidates

```python
from chemtools.ml.recommender import recommend_with_yield_prediction

result = recommend_with_yield_prediction(
    reaction_smiles=rxn,
    family="Suzuki",
    k=50
)
```

**Workflow**:
```
1. Get precedents
2. Load ML yield prediction model
3. Predict yields for all candidate conditions
4. Rank by predicted yield
```

**When to use**:
- When you have trained ML models available
- When you want yield-based ranking
- Not commonly used (most families don't have models)

---

## 5. **`recommend_conditions_structured()`** 📋 FORMATTER
**Location**: `chemtools/recommend/core.py`

**What it is**: Wrapper that formats output into structured JSON

```python
from chemtools.recommend import recommend_conditions_structured

result = recommend_conditions_structured(
    reaction=rxn,
    k=50,
    use_fusion=False
)
```

**What it does**: 
- Calls `recommend_from_reaction()` internally
- Formats output with chemicals list, conditions dict, summary
- Used by API endpoints

**When to use**: When you need nicely formatted JSON output (API, UI)

---

## Quick Decision Guide

```
┌─────────────────────────────────────────────────┐
│ What do you want to do?                         │
└─────────────────────────────────────────────────┘
         ↓
    ┌────────────────────────────────────┐
    │ Need structured JSON for API/UI?   │
    └────────────────────────────────────┘
         Yes ↓                  No ↓
    ┌──────────────┐      ┌──────────────┐
    │ Use:         │      │ Need simple  │
    │ recommend_   │      │ transparent  │
    │ conditions_  │      │ ranking?     │
    │ structured() │      └──────────────┘
    └──────────────┘         Yes ↓      No ↓
                        ┌─────────┐  ┌─────────┐
                        │ Use:    │  │ Use:    │
                        │ recommend│ │ recommend│
                        │ _simple()│ │ _from_  │
                        │         │  │ reaction│
                        │ ✅ NEW! │  │ ()      │
                        └─────────┘  └─────────┘
```

---

## Function Comparison Table

| Function | Complexity | Speed | Transparency | Use Case |
|----------|-----------|-------|--------------|----------|
| `recommend_from_reaction()` | Medium | Fast | Medium | **Main entry point** |
| `recommend_with_fusion()` | ❌ Very High | Slow | ❌ Low | Multi-source fusion (complex) |
| `recommend_simple()` | ✅ Low | ✅ Fast | ✅ High | **Simple precedent ranking** ⭐ |
| `recommend_with_yield_prediction()` | Medium | Medium | Medium | ML yield prediction |
| `recommend_conditions_structured()` | Medium | Fast | Medium | JSON formatting wrapper |

---

## Current Default Behavior

### What happens when you call the main function:

```python
from chemtools.recommend import recommend_from_reaction

# Default: use_fusion=False
result = recommend_from_reaction(reaction=rxn, k=50)
```

**Current flow**:
1. `recommend_from_reaction(use_fusion=False)` called
2. Uses simple precedent-based logic (NOT fusion, NOT the new simple ranker)
3. Returns recommendations based on frequency in top-k precedents

### To use fusion:
```python
result = recommend_from_reaction(reaction=rxn, use_fusion=True)
# → calls recommend_with_fusion() internally
# → complex multi-source fusion with rule-alignment
```

### To use new simple ranker:
```python
from chemtools.ml.simple_precedent_ranker import recommend_simple

result = recommend_simple(
    reaction_smiles=rxn,
    family="Suzuki",
    rerank_strategy='rule'  # or 'analytics'
)
# → simple, transparent precedent ranking
```

---

## Recommendation: Which Should You Use?

### For Most Users ✅
```python
from chemtools.ml.simple_precedent_ranker import recommend_simple

result = recommend_simple(
    reaction_smiles=rxn,
    family=family,
    rerank_strategy='rule'  # Simple rule-based boost
)
```

**Why**: Simple, fast, transparent, achieves your goal of precedent-centric with quality filtering

### For API/Production
```python
from chemtools.recommend import recommend_conditions_structured

result = recommend_conditions_structured(
    reaction=rxn,
    k=50,
    use_fusion=False  # Keep it simple
)
```

**Why**: Structured output, stable API, works with existing infrastructure

### Avoid
```python
# ❌ Don't use directly unless you understand the complexity
from chemtools.ml.fusion_recommender import recommend_with_fusion
```

**Why**: Overly complex, hard to debug, opaque reasoning

---

## Next Steps

**Option 1**: Integrate `recommend_simple()` into `recommend_from_reaction()`
```python
def recommend_from_reaction(..., rerank_strategy='rule'):
    if rerank_strategy in ['rule', 'analytics']:
        return recommend_simple(..., rerank_strategy=rerank_strategy)
    else:
        # existing logic
```

**Option 2**: Add new CLI argument
```bash
python scripts/local_recommendation_cli.py \
    --rxn "..." \
    --family Suzuki \
    --rerank rule  # or analytics
```

**Option 3**: Keep them separate
- Use `recommend_simple()` for new features
- Keep `recommend_from_reaction()` for backward compatibility

---

## Summary

**You have 5 functions, but realistically:**

1. **`recommend_simple()`** ← ✅ **Use this for new code** (simple, transparent)
2. **`recommend_from_reaction()`** ← Main entry point (backward compatible)
3. **`recommend_conditions_structured()`** ← API wrapper (formatted output)
4. **`recommend_with_fusion()`** ← ❌ Avoid (overly complex)
5. **`recommend_with_yield_prediction()`** ← Specialized (ML models)

**The new `recommend_simple()` achieves your goal**: precedent-centric with simple rule/analytics filtering, without the fusion complexity!
