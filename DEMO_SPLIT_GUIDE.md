# Demo Scripts Guide - Split Version

## Overview

The ChemTools demonstration scripts have been split into two focused files:

1. **`demo_basic_tools.py`** - Core utilities (SMILES, router, featurization)
2. **`demo_recommendations.py`** - Condition recommendations and ML

This separation makes it easier to:
- Learn specific features without information overload
- Test individual components independently
- Find relevant examples faster
- Run lighter demos when ML models aren't needed

---

## Quick Start

### Basic Tools Demo
```powershell
# Fast (~2 seconds) - No ML dependencies required
python demo_basic_tools.py
```

**What it demonstrates:**
- ✅ SMILES normalization (molecules & reactions)
- ✅ Family detection (with & without catalysts)
- ✅ Molecular featurization (leaving groups, nucleophiles)
- ✅ Property lookup
- ✅ DRFP similarity (optional)

### Recommendations Demo
```powershell
# Slower (~5-10 seconds) - May require ML models
python demo_recommendations.py
```

**What it demonstrates:**
- ✅ Precedent search (k-NN)
- ✅ Condition recommendations (NEW modular API)
- ✅ Structured outputs
- ✅ ML-enhanced predictions
- ✅ Plate design for HTS
- ✅ Advanced options (DRFP, constraints, variants)

---

## Demo 1: Basic Tools (`demo_basic_tools.py`)

### Features Covered

#### 1. SMILES Normalization
```python
from chemtools.smiles import normalize, normalize_reaction

# Canonicalize molecule
result = normalize('c1ccccc1O')  # → 'Oc1ccccc1'

# Normalize reaction
result = normalize_reaction('Brc1ccccc1.Nc1ccccc1>>...')
# → {'normalized': '...', 'reactants': [...], 'products': [...]}
```

#### 2. Family Detection (TWO METHODS)
```python
from chemtools.router import detect_family, detect_family_from_reaction

# Method 1: From reactant list (rule-based)
result = detect_family(['Brc1ccccc1', 'Nc1ccccc1'])

# Method 2: From reaction SMILES (detects catalysts!)
result = detect_family_from_reaction('Brc1ccccc1.Nc1ccccc1>Pd>...')
# ⭐ Pd catalyst → Buchwald_CN, Cu catalyst → Ullmann_CN
```

#### 3. Molecular Featurization
```python
from chemtools.featurizers.molecular import featurize

result = featurize('Brc1ccccc1', 'Nc1ccccc1')
# → {'lg': 'Br', 'nuc': 'amine_primary', 'bin': 'LG:Br|NUC:amine_primary'}
```

#### 4. Property Lookup
```python
from chemtools.properties import lookup

result = lookup('K3PO4')
# → {'name': '...', 'cas': '...'}
```

#### 5. DRFP Similarity
```python
from chemtools.reaction_similarity import drfp_tanimoto

similarity = drfp_tanimoto(rxn1, rxn2)
# → 0.850 (high similarity)
```

### When to Use
- ✅ Learning the basics
- ✅ Testing SMILES/featurization only
- ✅ No ML dependencies available
- ✅ Quick validation (~2 seconds runtime)

---

## Demo 2: Recommendations (`demo_recommendations.py`)

### Features Covered

#### 1. Precedent Search
```python
from chemtools.precedent import knn

precedents = knn(
    family="Ullmann C-N coupling",
    features={"lg": "Br", "nuc": "amine_primary"},
    k=5
)
```

#### 2. Basic Recommendations
```python
from chemtools.recommend import recommend_from_reaction

result = recommend_from_reaction(
    reaction='Brc1ccccc1.Nc1ccccc1>>...',  # ⚠️ 'reaction' not 'reaction_smiles'
    k=5
)
# → {'family': '...', 'conditions': [...], 'confidence': 0.90}
```

#### 3. Structured Outputs
```python
from chemtools.recommend import recommend_conditions_structured

result = recommend_conditions_structured(reaction=rxn, k=3)
# → Nested dicts with full reagent details
```

#### 4. ML-Enhanced Recommendations
```python
from chemtools.ml.recommender import hybrid_recommend

result = hybrid_recommend(
    reaction_smiles=rxn,  # ⚠️ Different param name!
    k=5
)
# → Includes predicted yields
```

#### 5. Plate Design
```python
from chemtools.recommend import design_plate_from_reaction

result = design_plate_from_reaction(
    reaction=rxn,
    top_cores=3,
    variants_per_core=2
)
# → HTS plate layout
```

#### 6. Advanced Options
```python
# DRFP relaxation
relax = {
    "use_drfp": True,
    "drfp_threshold": 0.5,
    "precompute_drfp": "candidates"
}

# Constraint rules
constraints = {
    "core": ["CuI", "Cu2O"],
    "solvent": ["DMF", "DMSO"]
}

result = recommend_from_reaction(
    reaction=rxn,
    k=5,
    relax=relax,
    constraint_rules=constraints
)
```

### When to Use
- ✅ Learning condition recommendations
- ✅ Testing ML models
- ✅ Designing HTS experiments
- ✅ Full feature exploration

---

## Comparison: Old vs New Demos

### Before (Single Demo)
```
demo_chemtools_complete.py
├── 17 demonstrations
├── 600+ lines
├── ~10-15 seconds runtime
└── Information overload for beginners
```

### After (Split Demos)
```
demo_basic_tools.py
├── 5 demonstrations
├── 250 lines
├── ~2 seconds runtime
└── Perfect for learning basics

demo_recommendations.py
├── 6 demonstrations
├── 400 lines
├── ~5-10 seconds runtime
└── Focus on recommendations
```

**Benefits:**
- ✅ Faster startup for basic tools
- ✅ Clearer organization
- ✅ Easier to find specific examples
- ✅ Run without ML dependencies (basic tools)

---

## API Inconsistencies to Remember

### Parameter Naming (October 2025)
```python
# ⚠️ Inconsistent parameter names across functions

# These use 'reaction' parameter:
recommend_from_reaction(reaction=rxn, ...)
design_plate_from_reaction(reaction=rxn, ...)

# This uses 'reaction_smiles' parameter:
hybrid_recommend(reaction_smiles=rxn, ...)

# Always check function signature!
```

### normalize_reaction() Output
```python
# ⚠️ Key is 'normalized', not 'reaction_smiles'
result = normalize_reaction(rxn)
normalized = result.get('normalized')  # ✅ Correct
normalized = result.get('reaction_smiles')  # ❌ Wrong (returns None)
```

---

## Running the Demos

### Basic Tools (Fast)
```powershell
# Full demo
python demo_basic_tools.py

# Expected output:
# ✅ SMILES normalization examples
# ✅ Family detection (both methods)
# ✅ Featurization results
# ✅ Property lookups
# ✅ DRFP similarity (if installed)
# Runtime: ~2 seconds
```

### Recommendations (Comprehensive)
```powershell
# Full demo
python demo_recommendations.py

# Expected output:
# ✅ Precedent search results
# ✅ Condition recommendations
# ✅ Structured outputs
# ✅ ML predictions (if models available)
# ✅ Plate designs
# ✅ Advanced options
# Runtime: ~5-10 seconds
```

---

## Troubleshooting

### "DRFP not installed"
```powershell
pip install drfp
```

### "ML module not available"
- Ensure LightGBM models are in `models/` directory
- Check `chemtools/ml/recommender.py` imports

### "No precedents found"
- Check dataset loading: `ChemTools(datasets=[...])`
- Verify family name matches dataset
- Try DRFP relaxation for broader search

### Windows Encoding Errors
Both demos now include UTF-8 encoding fix:
```python
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
```

---

## Next Steps

### After Running Demos
1. **Try Interactive UI**: `python app/ui_gradio.py`
2. **Read Documentation**:
   - `CHEMTOOLS_QUICKSTART.md` (15 min)
   - `CHEMTOOLS_CLASS_GUIDE.md` (30 min)
3. **Run Tests**: `pytest -q`
4. **Explore API**: Start server and visit `http://127.0.0.1:8000/docs`

### For Development
- See `MIGRATION_GUIDE.md` for old → new API migration
- See `RECOMMEND_REFACTORING_SUCCESS.md` for architecture details
- Check `tests/` for comprehensive test examples

---

## Summary

| Demo | Focus | Runtime | Dependencies |
|------|-------|---------|--------------|
| **demo_basic_tools.py** | SMILES, router, featurization | ~2s | Core only |
| **demo_recommendations.py** | Precedents, conditions, ML | ~5-10s | ML models |

**Key Improvements:**
- ✅ Clearer separation of concerns
- ✅ Faster basic demo (no ML overhead)
- ✅ Easier to find specific examples
- ✅ Better learning progression

**Both demos include:**
- ✅ Working examples
- ✅ Import patterns
- ✅ Function signatures
- ✅ Tips & best practices
- ✅ Error handling
- ✅ Windows UTF-8 support
