# Quick Answer: How Many Recommendation Functions?

## You Have 5 Functions (But Only Need 2)

### The 5 Functions:

```
1. recommend_from_reaction()          ← Main entry point (routes to others)
2. recommend_with_fusion()            ← Complex fusion (avoid)
3. recommend_simple()                 ← NEW! Simple ranker (use this!)
4. recommend_with_yield_prediction()  ← ML predictor (specialized)
5. recommend_conditions_structured()  ← JSON formatter (API wrapper)
```

---

## Visual Map

```
YOUR CODE
   ↓
┌──────────────────────────────────────────────────────┐
│ recommend_from_reaction(use_fusion=?)                │  ← Main Entry
└──────────────────────────────────────────────────────┘
   ↓                                    ↓
use_fusion=True                    use_fusion=False
   ↓                                    ↓
┌─────────────────────┐         ┌──────────────────┐
│ recommend_with_     │         │ Simple precedent │
│ fusion()            │         │ frequency logic  │
│                     │         │ (default)        │
│ ❌ Complex!         │         └──────────────────┘
│ - Adaptive weights  │
│ - 4 evidence sources│
│ - Rule alignment    │
└─────────────────────┘


NEW SIMPLE RANKER (not connected to main entry yet):

┌─────────────────────────────────────────────────────┐
│ recommend_simple(rerank_strategy='rule')            │
│                                                     │
│ ✅ Simple precedent ranking                         │
│ ✅ Optional rule/analytics reranking                │
│ ✅ Transparent reasoning                            │
└─────────────────────────────────────────────────────┘
```

---

## What You Should Use

### Option 1: New Simple Ranker (Recommended) ✅

```python
from chemtools.ml.simple_precedent_ranker import recommend_simple

result = recommend_simple(
    reaction_smiles="Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1",
    family="Ullmann_CN",
    k=30,
    rerank_strategy='rule'  # 'rule', 'analytics', or 'none'
)

# Clear, simple, transparent!
```

### Option 2: Main Entry Point (Existing) 📍

```python
from chemtools.recommend import recommend_from_reaction

result = recommend_from_reaction(
    reaction="Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1",
    k=50,
    use_fusion=False  # Keep it simple
)

# Works, but less control over ranking
```

### Option 3: API Wrapper (For JSON output) 📋

```python
from chemtools.recommend import recommend_conditions_structured

result = recommend_conditions_structured(
    reaction="Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1",
    k=50
)

# Returns nicely formatted JSON
```

---

## Summary

**Don't be confused!** You really only need to know about:

1. **`recommend_simple()`** ← Use for new code (simple + transparent)
2. **`recommend_from_reaction()`** ← Existing main API (backward compatible)

The other 3 functions are either:
- ❌ Too complex (`recommend_with_fusion`)
- 🔧 Specialized (`recommend_with_yield_prediction`)
- 📦 Just a wrapper (`recommend_conditions_structured`)

**For your use case** (precedents + simple reranking), use `recommend_simple()`!
