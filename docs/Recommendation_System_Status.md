# Recommendation System Analysis

## Current State (No Obsolete Code Found)

After comprehensive analysis of the codebase, **there is NO obsolete recommendation system to remove**. Here's what exists:

---

## 📁 Active Recommendation Systems

### 1. **`chemtools/recommend.py`** - k-NN Precedent-Based System ✅ ACTIVE

**Status:** **Core foundation**, actively used by API and ML system

**Key Functions:**
- `recommend_from_reaction()` - Main k-NN recommendation function
- `recommend_conditions_structured()` - API-friendly wrapper
- `design_plate_from_reaction()` - Plate design feature

**Used By:**
- ✅ API endpoints: `/api/v1/recommend`, `/api/v1/recommend/conditions`, `/api/v1/design_plate`
- ✅ `chemtools/recommend_ml.py` (imports and wraps it)
- ✅ UI (Gradio): `app/ui_gradio.py`, `scripts/ui_gradio.py`
- ✅ Scripts: `scripts/recommend_from_rxn.py`, `scripts/design_plate.py`
- ✅ Tests: `tests/test_recommend.py`, `tests/test_selector_payloads_amide.py`

**Why It's Not Obsolete:**
- The ML system (`recommend_ml.py`) **wraps this** - it doesn't replace it
- It's the **fallback** when ML conditions aren't met (n_precedents < threshold)
- It provides **precedent retrieval** which ML needs for interpretability

---

### 2. **`chemtools/recommend_ml.py`** - Hybrid ML + k-NN System ✅ ACTIVE

**Status:** **New addition**, uses k-NN as foundation

**Key Function:**
- `hybrid_recommend()` - Combines ML yield prediction with k-NN precedents

**Strategy:**
```python
if n_precedents >= ml_threshold and model_available:
    # Use ML to predict yields, re-rank k-NN results
    use_ml = True
else:
    # Fall back to pure k-NN
    use_ml = False
```

**Dependencies:**
- Imports `recommend_from_reaction` from `chemtools/recommend.py`
- Uses k-NN results as base, adds ML yield predictions on top

---

### 3. **`chemtools/rules/api.py`** - Rule-Based CRL System ✅ ACTIVE

**Status:** **Complementary system** for structured rule matching

**Key Function:**
- `recommend()` - Selects playbooks from Condition Rule Library (CRL)

**Used By:**
- `chemtools/condition_rules.py:recommend_rule_based()`
- Amide formation selectors and rule-based workflows

**Purpose:** Different from k-NN/ML - uses expert-defined rules, not precedent matching

---

## 🔍 Verification Results

### Functions Checked:
| Function | Location | Status | Usage Count |
|----------|----------|--------|-------------|
| `recommend_from_reaction` | `chemtools/recommend.py` | ✅ ACTIVE | 15+ references |
| `recommend_conditions_structured` | `chemtools/recommend.py` | ✅ ACTIVE | API endpoint |
| `design_plate_from_reaction` | `chemtools/recommend.py` | ✅ ACTIVE | API endpoint |
| `hybrid_recommend` | `chemtools/recommend_ml.py` | ✅ ACTIVE | New ML system |
| `rules.recommend` | `chemtools/rules/api.py` | ✅ ACTIVE | CRL system |

### No Obsolete Functions Found ❌
- No functions named `get_recommendation` or similar deprecated variants
- No unused recommendation code detected
- All recommendation systems serve distinct purposes

---

## 📊 System Architecture

```
User Request
    │
    ├─→ Pure k-NN (recommend.py)
    │   └─→ Precedent matching → Vote-based selection
    │
    ├─→ Hybrid ML (recommend_ml.py)
    │   ├─→ Calls recommend.py for base precedents
    │   └─→ Adds ML yield predictions → Re-ranks
    │
    └─→ Rule-Based (rules/api.py)
        └─→ Playbook matching → Expert rules
```

---

## ✅ Conclusion

**NO CODE TO REMOVE**

All recommendation systems are:
1. **Actively maintained**
2. **Used by API/UI**
3. **Serve complementary purposes**

The systems work together:
- **k-NN** (recommend.py) = Foundation
- **ML** (recommend_ml.py) = Enhancement layer
- **Rules** (rules/api.py) = Expert knowledge

---

## 💡 What You Might Be Thinking Of

If you're looking to clean up code, these are candidates (but NOT recommendation systems):

### Potentially Unused Code (Unrelated to Recommendations):
1. **Test fixtures** that are no longer referenced
2. **Old data processing scripts** in `data-processor/` (check if outdated)
3. **Temporary scripts** like `scripts/tmp_*.py` files
4. **Old documentation** in `docs/old/` directory

### Next Steps:
If you want to clean up:
1. ✅ Remove `scripts/tmp_*.py` temporary scripts
2. ✅ Archive `docs/old/` deprecated documentation
3. ✅ Check `data-processor/` for obsolete processing scripts
4. ❌ **DO NOT** remove `chemtools/recommend.py` - it's core infrastructure

---

## 📝 Recommendation

**Keep all three recommendation systems intact:**
- `chemtools/recommend.py` - Core k-NN engine
- `chemtools/recommend_ml.py` - ML enhancement
- `chemtools/rules/api.py` - Expert rules

They form a **complementary ecosystem**, not competing alternatives.
