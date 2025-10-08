# ChemTools Demo Scripts - Quick Reference

## 📚 Two Focused Demos

### 1. **`demo_basic_tools.py`** - Core Utilities
**Runtime:** ~2 seconds  
**Dependencies:** Core only (no ML models needed)

```powershell
python demo_basic_tools.py
```

**Demonstrates:**
- ✅ SMILES normalization (`normalize`, `normalize_reaction`)
- ✅ Family detection (`detect_family`, `detect_family_from_reaction`)
- ✅ Molecular featurization (`featurize`)
- ✅ Property lookup (`lookup`)
- ✅ DRFP similarity (`drfp_tanimoto`, optional)

**Perfect for:**
- Learning the basics
- Testing core functionality
- Quick validation
- No ML dependencies needed

---

### 2. **`demo_recommendations.py`** - Condition Recommendations
**Runtime:** ~5-10 seconds  
**Dependencies:** Core + optional ML models

```powershell
python demo_recommendations.py
```

**Demonstrates:**
- ✅ Precedent search (`knn`)
- ✅ Condition recommendations (`recommend_from_reaction`)
- ✅ Structured outputs (`recommend_conditions_structured`)
- ✅ ML-enhanced predictions (`hybrid_recommend`)
- ✅ Plate design (`design_plate_from_reaction`)
- ✅ Advanced options (DRFP, constraints, variants)

**Perfect for:**
- Condition optimization
- HTS plate design
- ML-enhanced predictions
- Full feature exploration

---

## 🎯 Quick Comparison

| Feature | Basic Tools | Recommendations |
|---------|-------------|-----------------|
| **SMILES normalization** | ✅ | - |
| **Family detection** | ✅ (both methods) | - |
| **Featurization** | ✅ | - |
| **Property lookup** | ✅ | - |
| **DRFP similarity** | ✅ | - |
| **Precedent search** | - | ✅ |
| **Recommendations** | - | ✅ |
| **ML predictions** | - | ✅ |
| **Plate design** | - | ✅ |
| **Runtime** | ~2s | ~5-10s |
| **ML models required** | ❌ | Optional |

---

## 📖 Key Learnings from Demos

### Family Detection - TWO Methods

**Method 1:** `detect_family(reactants)` - Rule-based only
```python
from chemtools.router import detect_family

result = detect_family(['Brc1ccccc1', 'Nc1ccccc1'])
# → {'family': 'Ullmann_CN', 'confidence': 0.90}
```

**Method 2:** `detect_family_from_reaction(reaction_smiles)` - ⭐ Better!
```python
from chemtools.router import detect_family_from_reaction

result = detect_family_from_reaction('Brc1ccccc1.Nc1ccccc1>Pd>...')
# → Detects Pd catalyst → Buchwald_CN (higher accuracy)
# → Detects Cu catalyst → Ullmann_CN
```

**Why Method 2 is better:**
- ✅ Accepts full reaction SMILES (reactants>>agents>>products)
- ✅ Detects catalysts from agents (Pd, Cu, Ni, Co)
- ✅ Applies catalyst override logic
- ✅ Optional ML-based detection via `use_rxn_insight=True`

### API Parameter Inconsistency ⚠️

**Different parameter names across functions:**

```python
# These use 'reaction' parameter:
recommend_from_reaction(reaction=rxn, k=5)
design_plate_from_reaction(reaction=rxn, plate_size=24)

# This uses 'reaction_smiles' parameter:
hybrid_recommend(reaction_smiles=rxn, k=5)
```

**Always check the function signature!**

### normalize_reaction() Output

```python
result = normalize_reaction('Brc1ccccc1.Nc1ccccc1>>...')

# ✅ Correct:
normalized = result.get('normalized')

# ❌ Wrong (returns None):
normalized = result.get('reaction_smiles')
```

### knn() Returns Dict, Not List

```python
result = knn(family='Ullmann C-N coupling', features={...}, k=5)

# ✅ Correct:
precedents = result.get('precedents', [])
support = result.get('support', 0)

# ❌ Wrong (KeyError):
top = result[0]
```

### Plate Design Parameters

```python
# ✅ Correct (October 2025):
design_plate_from_reaction(
    reaction=rxn,
    plate_size=12  # Number of wells
)

# ❌ Old API (deprecated):
design_plate_from_reaction(
    reaction=rxn,
    top_cores=3,
    variants_per_core=2
)
```

---

## 🚀 Performance Tips

### Selective Dataset Loading (50-100x faster)
```python
from chemtools import ChemTools

fast = ChemTools(
    datasets=["C_N_Coupling_Pd"],  # Only this reaction type
    reagent_dbs=["ligand", "base"] # Only these reagents
)

precedents = fast.precedent.knn(...)
```

### DRFP Optimization
```python
relax = {
    "use_drfp": True,
    "drfp_threshold": 0.5,
    "precompute_drfp": "candidates",  # Not "all"
    "drfp_n_bits": 2048               # Smaller = faster
}
```

### Environment Variables
```powershell
# Windows (faster startup, testing only)
$env:CHEMTOOLS_DISABLE_RDKIT='1'
$env:CHEMTOOLS_LOAD_DATASET='0'

# Linux/Mac
export CHEMTOOLS_DISABLE_RDKIT=1
export CHEMTOOLS_LOAD_DATASET=0
```

---

## 🔧 Import Patterns

### Basic Tools
```python
from chemtools.smiles import normalize, normalize_reaction
from chemtools.router import detect_family, detect_family_from_reaction
from chemtools.featurizers.molecular import featurize
from chemtools.properties import lookup
from chemtools.reaction_similarity import drfp_tanimoto  # Optional
```

### Recommendations
```python
from chemtools.precedent import knn
from chemtools.recommend import (
    recommend_from_reaction,
    recommend_conditions_structured,
    design_plate_from_reaction
)
from chemtools.ml.recommender import hybrid_recommend  # Optional
```

---

## 📚 Next Steps

### After Running Demos:

1. **Interactive UI**
   ```powershell
   python app/ui_gradio.py
   ```

2. **Documentation**
   - `CHEMTOOLS_QUICKSTART.md` - 15 min start
   - `CHEMTOOLS_CLASS_GUIDE.md` - 30 min deep dive
   - `DEMO_SPLIT_GUIDE.md` - Detailed demo guide

3. **Tests**
   ```powershell
   pytest -q
   ```

4. **API Documentation**
   ```powershell
   make run  # Start server
   # Visit: http://127.0.0.1:8000/docs
   ```

---

## ✅ Summary

**Before (Single Demo):**
- ❌ 17 demonstrations in one file
- ❌ 600+ lines
- ❌ ~10-15 seconds runtime
- ❌ Information overload

**After (Split Demos):**
- ✅ Focused demos (5 + 6 demonstrations)
- ✅ Clearer organization
- ✅ Faster basic demo (~2s)
- ✅ Easier to find examples
- ✅ Optional ML dependencies

**Both demos work perfectly!** 🎉

---

## 🐛 Common Issues

### "DRFP not installed"
```powershell
pip install drfp
```

### "No precedents found"
- Dataset may be empty (sample data)
- Try different family name
- Enable DRFP relaxation

### Windows encoding errors
- Both demos include UTF-8 fix
- Should work on Windows/Mac/Linux

### "ML module not available"
- ML models optional
- Demos handle gracefully

---

**Created:** October 2025  
**Last Updated:** Demo split (basic tools + recommendations)
