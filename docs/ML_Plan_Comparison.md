# ML Implementation Plan: Comparison Summary

## What Changed? ⚡

### Original Plan → Updated Plan

| Aspect | ❌ Original (Hand-Crafted) | ✅ Updated (End-to-End Learning) |
|--------|---------------------------|----------------------------------|
| **Features** | Categorical (LG, nuc_class, ortho_count, para_EWG, etc.) | DRFP → Neural embeddings → GNN |
| **Generalization** | Ullmann-specific rules | Universal (all reaction families) |
| **Maintenance** | High (new rules for each family) | Low (data-driven) |
| **Timeline** | 12 weeks | **8-10 weeks** |
| **Phase 1** | XGBoost on categorical features | DRFP + LightGBM |
| **Phase 2** | Neural embeddings | Neural network + learned condition embeddings |
| **Phase 3** | Multi-task learning | Graph Neural Network (GNN) |
| **Expected MAE** | 10-12% → 8-10% | 12-15% → 8-10% (same target) |

---

## Why Avoid Hand-Crafted Features?

### Problems with Current Categorical System

```python
# Current (Ullmann-specific):
features = {
    "LG": "Br",              # ❌ What about OTf? OMs? Complex leaving groups?
    "elec_class": "aryl",    # ❌ Loses all structural information
    "nuc_class": "aniline",  # ❌ All anilines treated the same
    "ortho_count": 1,        # ❌ Doesn't capture which ortho position
    "para_EWG": False,       # ❌ What about meta? Mixed substituents?
}

# Issues:
# 1. Ullmann rules don't work for Suzuki/Amide/etc.
# 2. Need expert to define rules for every new reaction type
# 3. Loses 90% of structural information
# 4. Can't transfer learn from other datasets
```

### Benefits of Learned Representations

```python
# Updated (Universal):
representations = {
    "DRFP": [0, 1, 0, ..., 1],     # 4096 bits - captures all substructures
    # OR
    "GNN": molecular_graph,         # Learns optimal features from data
    # OR  
    "Transformer": pretrained_emb,  # ChemBERTa, MolCLR, etc.
}

# Benefits:
# ✅ Works for ALL reaction families
# ✅ No expert rules needed
# ✅ Captures full molecular structure
# ✅ Can use pretrained models
# ✅ Learns what's important from data
```

---

## Updated 3-Phase Approach

### Phase 1: DRFP + LightGBM (2-3 weeks)
**Simple and effective baseline**

```python
Input:
  Reaction SMILES → DRFP (4096 bits)
  Conditions → One-hot encoding

Model:
  LightGBM Gradient Boosting

Output:
  Predicted yield ± uncertainty

Performance:
  MAE: 12-15%
  Training: 5-10 minutes (CPU)
  Inference: <100ms
```

**Advantages:**
- ✅ No hand-crafted features needed
- ✅ Works across all reaction families
- ✅ Fast training on CPU
- ✅ Easy to deploy

---

### Phase 2: Neural Network + Learned Embeddings (3-4 weeks)
**Better representations through learning**

```python
Input:
  Reaction SMILES → DRFP (4096 bits)
  Core/Base/Solvent → Learned embeddings (64 dims each)
  T/time → Normalized numeric

Model:
  Neural network (3 hidden layers + residual connections)

Output:
  Predicted yield + log variance (calibrated uncertainty)

Performance:
  MAE: 10-12%
  Training: 1-2 hours (GPU)
  Inference: <100ms
```

**Key improvement:**
- Core/Base/Solvent embeddings learn chemical similarity
- Example: Pd(PPh3)4 and Pd(dppf)Cl2 have similar embeddings
- Better uncertainty via negative log-likelihood loss

---

### Phase 3: Graph Neural Network (3-4 weeks, optional)
**State-of-the-art with end-to-end learning**

```python
Input:
  Reactants → Molecular graphs (atoms + bonds)
  Products → Molecular graphs
  Conditions → Learned embeddings

Model:
  GNN encoder (4 graph conv layers)
  Reaction = Product_embedding - Reactant_embedding
  MLP prediction head

Output:
  Predicted yield + uncertainty

Performance:
  MAE: 8-10%
  Training: 3-5 hours (GPU)
  Inference: 200-500ms
```

**Best performance:**
- Learns optimal molecular representations
- Can capture 3D geometry
- Transfer learning from ChemBERTa/MolCLR
- State-of-the-art accuracy

---

## Implementation Comparison

### Original Plan (Hand-Crafted Features)

```python
# Phase 1: Extract categorical features
def featurize_ullmann(electrophile, nucleophile):
    # ❌ Expert rules for Ullmann only
    lg = detect_leaving_group(electrophile)  # Br, Cl, I, OTf
    elec_class = classify_electrophile(electrophile)  # aryl, vinyl, alkyl
    nuc_class = classify_nucleophile(nucleophile)  # aniline, amine_primary, etc.
    # ... more rules ...
    return {"LG": lg, "elec_class": elec_class, "nuc_class": nuc_class}

# Phase 2: Train on categorical features
model = XGBRegressor()
model.fit(categorical_features, yields)
```

**Problems:**
- Need to write `featurize_suzuki()`, `featurize_amide()`, etc.
- Each new reaction type needs expert rules
- Loses molecular structure information

---

### Updated Plan (Learned Representations)

```python
# Phase 1: Use DRFP (no rules needed!)
from drfp import DrfpEncoder

encoder = DrfpEncoder()

def encode_reaction(reaction_smiles):
    # ✅ Works for ANY reaction
    fp = encoder.encode(reaction_smiles, n_folded_length=4096)
    return np.array(fp)

# Train directly on fingerprints
model = LGBMRegressor()
model.fit(reaction_fingerprints, yields)
```

**Advantages:**
- No expert rules needed
- Works for all reaction types
- Captures full structural information
- Easy to extend

---

## Quick Start: Updated Steps

### Week 1: Data Preparation
```bash
# Extract and split dataset
python scripts/ml/prepare_dataset.py \
  --source data/reaction_dataset \
  --output data/ml_splits

# Output:
# data/ml_train.jsonl (70%)
# data/ml_val.jsonl (15%)
# data/ml_test.jsonl (15%)
```

### Week 2: Train DRFP Baseline
```bash
# Train baseline model (no hand-crafted features!)
python scripts/ml/train_drfp_model.py \
  --train data/ml_train.jsonl \
  --val data/ml_val.jsonl \
  --output models/drfp_yield_v1.pkl

# Expected output:
# Validation MAE: 12-15%
# Training time: 5-10 minutes
```

### Week 3: Integrate with API
```python
# chemtools/recommend_ml.py
from chemtools.ml.drfp_predictor import DRFPYieldPredictor

model = DRFPYieldPredictor.load('models/drfp_yield_v1.pkl')

# Predict yield for any reaction (no family-specific code!)
yield_pred, uncertainty = model.predict(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    conditions=[{
        'core': 'Cu/1,10-phenanthroline',
        'base_uid': 'K2CO3',
        'solvent_uid': 'DMSO',
        'T_C': 110,
        'time_h': 24
    }]
)[0]

print(f"Predicted yield: {yield_pred:.1f}% ± {uncertainty:.1f}%")
```

---

## What You Get

### Immediate (Phase 1, Week 1-3):
✅ Yield prediction for **all reaction families**  
✅ No hand-crafted features to maintain  
✅ Works on CPU (no GPU needed)  
✅ Easy to deploy (<100MB model)

### Short-term (Phase 2, Week 4-7):
✅ Improved accuracy (10-12% MAE)  
✅ Calibrated uncertainty estimates  
✅ Learned condition embeddings (chemical similarity)  
✅ Active learning for experiment design

### Long-term (Phase 3, Week 8-10):
✅ State-of-the-art performance (8-10% MAE)  
✅ End-to-end learning from molecular graphs  
✅ Transfer learning from pretrained models  
✅ Best-in-class recommendation system

---

## Decision Tree: Which Phase to Implement?

```
Start with Phase 1 (DRFP + LightGBM)
│
├─ Good enough? (MAE < 15%)
│  └─ YES → Deploy to production ✅
│  └─ NO → Continue to Phase 2
│
├─ Phase 2 (Neural + Learned Embeddings)
│  │
│  ├─ Good enough? (MAE < 12%)
│  │  └─ YES → Deploy to production ✅
│  │  └─ NO → Continue to Phase 3
│  │
│  └─ Phase 3 (Graph Neural Network)
│     └─ State-of-the-art (MAE < 10%) → Deploy ✅
```

**Recommendation:** Start with Phase 1, evaluate, then decide whether Phase 2/3 is worth the effort.

---

## Files to Create (Updated)

```
chemtools/ml/
├── __init__.py
├── drfp_predictor.py           # Phase 1: DRFP + LightGBM
├── neural_yield_predictor.py   # Phase 2: Neural network
├── gnn_yield_predictor.py      # Phase 3: GNN
├── data_loader.py              # PyTorch Dataset
└── evaluation.py               # Metrics, plots

scripts/ml/
├── prepare_dataset.py          # Extract and split data
├── train_drfp_model.py         # Phase 1 training
├── train_neural_model.py       # Phase 2 training
├── train_gnn_model.py          # Phase 3 training
└── evaluate_models.py          # Benchmark all models

chemtools/
└── recommend_ml.py             # Hybrid recommender (ML + k-NN fallback)

models/
├── drfp_yield_v1.pkl           # Phase 1 model (100MB)
├── neural_yield_v1.pt          # Phase 2 model (50MB)
└── gnn_yield_v1.pt             # Phase 3 model (80MB)
```

---

## Timeline Summary

| Phase | Weeks | Deliverables | Expected MAE |
|-------|-------|--------------|--------------|
| **Phase 1** | 2-3 | DRFP baseline | 12-15% |
| **Phase 2** | 3-4 | Neural network | 10-12% |
| **Phase 3** | 3-4 | GNN (optional) | 8-10% |
| **Total** | **8-10** | Production system | **8-15%** |

**vs Original Plan:** 12 weeks → 8-10 weeks (20% faster!)

---

## Key Takeaways

### ❌ What We're Removing:
- Hand-crafted categorical features (LG, nuc_class, ortho_count, etc.)
- Family-specific featurization functions
- Expert rules for each reaction type
- Maintenance burden of feature engineering

### ✅ What We're Adding:
- Universal DRFP fingerprints (works for all reactions)
- Learned condition embeddings (chemical similarity)
- Optional GNN for end-to-end learning
- Faster path to production (8-10 weeks)

### 🎯 What Stays the Same:
- Yield prediction as primary goal
- Hybrid system (ML + k-NN fallback)
- Calibrated uncertainty estimates
- Integration with existing API

---

## Next Steps

1. **Review** this updated plan with your team
2. **Prepare** dataset with `scripts/ml/prepare_dataset.py`
3. **Start Phase 1** (DRFP baseline) - 2 weeks
4. **Evaluate** and decide: Deploy or continue to Phase 2?

**Ready to implement!** 🚀

Read full details: `docs/ML_Updated_Plan_End_to_End_Learning.md`
