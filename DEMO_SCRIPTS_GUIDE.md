# ChemTools Demo Scripts - Usage Guide

## 📚 Available Demo Scripts

We've created comprehensive demonstration scripts to showcase all ChemTools features:

### 1. **`demo_chemtools_quick.py`** ⭐ **START HERE**

**Quick, working examples of core features**

```bash
python demo_chemtools_quick.py
```

**What it demonstrates:**
- ✅ SMILES normalization (molecules and reactions)
- ✅ Reaction family detection
- ✅ Molecular featurization
- ✅ Precedent search (k-NN)
- ✅ Condition recommendations (NEW modular API)
- ✅ ML-enhanced recommendations
- ✅ DRFP similarity calculation
- ✅ Property lookup

**Includes:**
- NEW import patterns (refactored modules)
- Function signatures (correct parameter names)
- Performance tips
- Links to documentation

**Runtime:** ~5 seconds (works without dataset loading)

---

### 2. **`demo_chemtools_complete.py`** 

**Comprehensive showcase with 17 demonstrations**

```bash
python demo_chemtools_complete.py
```

**What it demonstrates:**
- All features from quick demo PLUS:
- Structured recommendation API
- Plate design generation  
- Resource management (ChemTools class)
- Backward compatibility
- FastAPI endpoints reference
- Gradio UI tabs overview
- Performance optimization
- Testing approaches
- Troubleshooting guide

**Includes:**
- Detailed code examples
- Error handling demonstrations
- Migration guide references
- Complete API documentation

**Runtime:** ~10 seconds

---

## 🚀 Quick Start Examples

### Example 1: Recommend Conditions

```python
from chemtools.recommend import recommend_from_reaction

result = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=25,
    relax={"use_drfp": True}
)

print(f"Core: {result['recommendation']['core']}")
print(f"Base: {result['recommendation']['base']}")
print(f"Confidence: {result['recommendation']['confidence']:.2f}")
```

### Example 2: ML-Enhanced Recommendations

```python
from chemtools.ml.recommender import hybrid_recommend

result = hybrid_recommend(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50
)

for cond in result['recommended_conditions'][:3]:
    print(f"{cond['core']}: {cond.get('predicted_yield_pct', 'N/A')}%")
```

### Example 3: Detect Reaction Family

```python
from chemtools.router import detect_family

result = detect_family(["Brc1ccccc1", "Nc1ccccc1"])
print(f"Family: {result['family']}")
print(f"Confidence: {result['confidence']:.2f}")
```

### Example 4: Selective Loading (Fast)

```python
from chemtools import ChemTools

# Only load what you need - 50-100x faster!
fast = ChemTools(
    datasets=["C_N_Coupling_Pd"],
    reagent_dbs=["ligand", "base"]
)

precedents = fast.precedent.knn(
    family="C_N_Coupling_Pd",
    features={"bin": "LG:Br|NUC:aniline"},
    k=5
)
```

---

## 🔑 Key API Changes (October 2025 Refactoring)

### NEW Module Structure

```python
# ✅ Recommendations - NEW modular package
from chemtools.recommend import recommend_from_reaction
from chemtools.recommend import recommend_conditions_structured
from chemtools.recommend import design_plate_from_reaction

# ✅ ML - NEW location
from chemtools.ml.recommender import hybrid_recommend

# ✅ DRFP - Correct function name
from chemtools.reaction_similarity import drfp_tanimoto
```

### Important Parameter Names

**⚠️ Parameter name varies by function:**

```python
# recommend_from_reaction uses "reaction"
recommend_from_reaction(
    reaction="...",  # ← "reaction" (not "reaction_smiles")
    k=25
)

# hybrid_recommend uses "reaction_smiles"
hybrid_recommend(
    reaction_smiles="...",  # ← "reaction_smiles"
    k=50
)

# design_plate_from_reaction uses "reaction"
design_plate_from_reaction(
    reaction="...",  # ← "reaction"
    top_cores=5
)
```

---

## 📊 Demo Output Examples

### SMILES Normalization
```
✅ normalize('c1ccccc1O')
   → Oc1ccccc1

✅ normalize_reaction('Brc1ccccc1.Nc1ccccc1>>...')
   → Reactants: 2
   → Products: 1
```

### Reaction Family Detection
```
✅ detect_family(['Brc1ccccc1', 'Nc1ccccc1'])
   → Family: Ullmann_CN
   → Confidence: 0.90
```

### Featurization
```
✅ featurize('Brc1ccccc1', 'Nc1ccccc1')
   → LG: Br
   → Nuc class: amine_primary
   → Bin: LG:Br|NUC:amine_primary
```

### Recommendation
```
✅ recommend_from_reaction('...')
   → Family: C_N_Coupling_Cu
   → Core: Cu
   → Base: K3PO4
   → Confidence: 0.38
```

---

## 🛠️ Troubleshooting

### Issue: Dataset not loaded

**Solution:**
```python
from chemtools import ChemTools

# Use selective loading
chem = ChemTools(datasets=["C_N_Coupling_Pd"])
```

### Issue: DRFP not available

**Solution:**
```bash
pip install drfp
```

Or disable DRFP:
```python
relax = {"use_drfp": False}
```

### Issue: Slow first request

**Expected!** Dataset loads on first use (60-120s).

**Solutions:**
- Use `preload=True` for API servers
- Keep server running to reuse cache
- Use selective dataset loading

---

## 📖 Documentation References

### Getting Started
1. **`CHEMTOOLS_QUICKSTART.md`** - 15 min quick start
2. **`CHEMTOOLS_CLASS_GUIDE.md`** - Comprehensive API guide (30 min)
3. **`README.md`** - Project overview and endpoints

### Migration
4. **`MIGRATION_GUIDE.md`** - Migrating from old API
5. **`RECOMMEND_REFACTORING_SUCCESS.md`** - What changed in refactoring

### Advanced
6. **API Docs:** `http://127.0.0.1:8000/docs` (start server first)
7. **Test Suite:** `pytest -q` or `pytest tests/chemtools/ -v`

---

## 🎯 Next Steps

### 1. Run Quick Demo
```bash
python demo_chemtools_quick.py
```

### 2. Try Gradio UI
```bash
python app/ui_gradio.py
# Open: http://127.0.0.1:7860
```

### 3. Start API Server
```bash
uvicorn app.main:app --reload --port 8000
# Open: http://127.0.0.1:8000/docs
```

### 4. Run Tests
```bash
pytest -q                        # All tests
pytest test_new_recommend.py -v  # Recommendation tests
pytest tests/chemtools/ -v       # Package tests
```

### 5. Read Documentation
- Start with `CHEMTOOLS_QUICKSTART.md`
- Deep dive with `CHEMTOOLS_CLASS_GUIDE.md`
- Check migration guide if updating old code

---

## ✅ Demo Script Checklist

- ✅ `demo_chemtools_quick.py` - Quick working examples
- ✅ `demo_chemtools_complete.py` - Comprehensive showcase
- ✅ Both scripts test actual refactored APIs
- ✅ Include correct import patterns
- ✅ Show proper function signatures
- ✅ Demonstrate performance tips
- ✅ Link to documentation

---

## 💡 Pro Tips

1. **Start with quick demo** - Get working examples fast
2. **Use selective loading** - 50-100x faster for specific reactions
3. **Keep server running** - Reuse cached data
4. **Read the signatures** - Parameter names vary (`reaction` vs `reaction_smiles`)
5. **Check docs** - CHEMTOOLS_QUICKSTART.md has everything you need

---

**Created:** October 7, 2025  
**ChemTools Version:** v2 (Refactored Modular Architecture)  
**Status:** ✅ Production Ready
