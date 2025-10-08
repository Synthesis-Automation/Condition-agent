# ML Recommendation System: Quick Summary

## Current System Analysis

### What You Have Now ✅
```
k-NN Precedent System (chemtools/precedent.py + recommend.py)
├─ Categorical feature matching (LG, nuc_class, ortho, etc.)
├─ Optional DRFP fingerprint similarity (whole-reaction Tanimoto)
├─ Vote-based aggregation (Laplace smoothing)
├─ Simple confidence scores (heuristic from vote share)
└─ Median T/time from precedents
```

**Strengths:**
- ✅ Interpretable and explainable
- ✅ Works with small datasets
- ✅ Fast inference (<100ms)
- ✅ No training required

**Limitations:**
- ❌ No yield prediction
- ❌ Hand-crafted features don't generalize
- ❌ Simple similarity metrics
- ❌ Cold start with new reaction types
- ❌ No uncertainty quantification

---

## Proposed ML Architecture 🚀

### 5-Phase Implementation Plan

**Phase 1: Gradient Boosting (2-3 weeks)**
```python
XGBoost/LightGBM Models
├─ Yield prediction (MAE target: <10%)
├─ Core/base/solvent classification
└─ Learned feature interactions
```

**Phase 2: Neural Embeddings (3-4 weeks)**
```python
Graph Neural Network / Transformer
├─ Learn reaction representations
├─ Contrastive learning (triplet loss)
└─ Better similarity than DRFP
```

**Phase 3: Multi-Task Learning (2-3 weeks)**
```python
Single Model Predicts All:
├─ Yield (regression)
├─ Core (classification)
├─ Base (classification)
├─ Solvent (classification)
├─ Temperature (regression)
└─ Time (regression)
```

**Phase 4: Uncertainty Quantification (1-2 weeks)**
```python
Calibrated Confidence
├─ MC Dropout / Ensembles
├─ Conformal prediction intervals
└─ Active learning suggestions
```

**Phase 5: Production Integration (1-2 weeks)**
```python
Hybrid System
├─ ML for high-data scenarios (>50 precedents)
├─ k-NN fallback for sparse data
├─ FastAPI endpoints
└─ Gradio UI updates
```

---

## Expected Impact 📊

| Metric | Current | With ML | Improvement |
|--------|---------|---------|-------------|
| Yield Prediction | None | MAE ~10% | New capability |
| Top-5 Accuracy | ~60% | ~75-80% | +15-20% |
| Confidence | Heuristic | Calibrated | Trustworthy |
| Cold Start | Poor | Transfer learning | Much better |
| Inference Time | 50ms | 80ms | Acceptable |

**Business Impact:**
- 📉 Reduce screening plate sizes by 30-40%
- ⏱️ Faster time to first hit (less trial & error)
- 💰 Cost savings from fewer failed experiments
- 🎯 Active learning identifies high-value experiments

---

## Hybrid System Design

```python
class HybridRecommender:
    def recommend(self, reaction):
        precedents = knn_search(reaction)  # Always retrieve
        
        if len(precedents) >= 50:
            # Use ML predictions
            ml_output = self.ml_model.predict(reaction)
            return {
                'method': 'ml',
                'yield_pred': ml_output.yield,  # NEW
                'uncertainty': ml_output.std,    # NEW
                'recommendations': ml_ranked_conditions,
                'precedents': precedents,        # Still available
            }
        else:
            # Fallback to k-NN
            return knn_based_recommendation(precedents)
```

**Best of both worlds:**
- ML predictions when data is rich
- k-NN fallback for sparse scenarios
- Precedents always shown for interpretability

---

## Quick Start: Phase 1 Implementation

### Step 1: Data Preparation
```bash
# Extract training data
python scripts/ml/prepare_dataset.py \
  --source data/reaction_dataset \
  --output data/ml_train.jsonl \
  --split 0.7/0.15/0.15

# Output: train.jsonl, val.jsonl, test.jsonl
```

### Step 2: Train Baseline Model
```bash
# Train XGBoost yield predictor
python scripts/ml/train_yield_model.py \
  --train data/ml_train.jsonl \
  --val data/ml_val.jsonl \
  --output models/yield_v1.pkl

# Expected: MAE 12-15% on validation
```

### Step 3: Evaluate
```bash
# Benchmark against k-NN
python scripts/ml/evaluate_models.py \
  --test data/ml_test.jsonl \
  --knn-baseline \
  --ml-model models/yield_v1.pkl

# Compare metrics: MAE, top-5 accuracy, inference time
```

### Step 4: Deploy
```python
# Update recommend.py
from chemtools.ml import YieldPredictor

predictor = YieldPredictor.load('models/yield_v1.pkl')

def recommend_with_ml(reaction, k=50):
    # ... existing k-NN code ...
    
    # Add ML prediction
    if n_precedents >= 50:
        yield_pred, uncertainty = predictor.predict(reaction)
        recommendation['yield_prediction'] = yield_pred
        recommendation['confidence'] = 1 - (uncertainty / 100)
    
    return recommendation
```

---

## Data Requirements

**Current:** ~1,000 reactions (estimated)  
**Minimum for ML:** 1,000 per family  
**Recommended:** 5,000+ per family

**Data sources:**
- ✅ Current dataset: `data/reaction_dataset/*.jsonl`
- 📖 Literature mining: USPTO, PubMed (tools: RxnScribe, RoboRXN)
- 🧪 ELN exports: Benchling, SciNote
- 🔬 Vendor data: SciFinder, Reaxys (if licensed)
- 🤖 Synthetic generation: Template-based enumeration

---

## Timeline & Resources

**Total Time:** 12 weeks (3 months)

**Week 1-2:** Baseline XGBoost (Phase 1a)  
**Week 3-4:** Condition ranking (Phase 1b)  
**Week 5-6:** Neural embeddings (Phase 2)  
**Week 7-8:** Multi-task model (Phase 3)  
**Week 9-10:** Uncertainty + active learning (Phase 4)  
**Week 11-12:** Production deployment (Phase 5)

**Team:**
- 1 ML engineer (full-time)
- 1 chemist (part-time, for validation)
- Optional: GPU instance (AWS p3.2xlarge ~$3/hour)

**Budget:**
- GPU compute: ~$500-1000/month
- Experiment tracking: Free (Weights & Biases community tier)
- Total: ~$1,500 for 3 months

---

## Success Metrics

### Technical Metrics (Phase 1)
- [ ] Yield prediction MAE < 15%
- [ ] Top-5 condition accuracy > 60%
- [ ] Inference time < 100ms
- [ ] Model size < 100MB (for deployment)

### Business Metrics (After Deployment)
- [ ] Lab hit rate improves by >10%
- [ ] Average screening plate size reduces 30%
- [ ] Time to first hit reduces by 20%
- [ ] 80%+ chemists prefer ML recommendations (user survey)

---

## Risk Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Insufficient data | ML underperforms | Transfer learning, data augmentation |
| Overfitting | Poor generalization | Cross-validation, regularization |
| Slow inference | Production bottleneck | Model distillation, caching |
| User distrust | Low adoption | SHAP explanations, show precedents |

---

## Next Actions

### Immediate (This Week)
1. ✅ Review this analysis with team
2. ⬜ Assess current dataset size per family
3. ⬜ Set up ML infrastructure (GPU, experiment tracking)
4. ⬜ Start Phase 1: baseline XGBoost

### Short-term (Month 1)
1. ⬜ Train yield predictor
2. ⬜ A/B test against k-NN
3. ⬜ Gather user feedback
4. ⬜ Prioritize Phase 2 vs Phase 3

### Long-term (Months 2-3)
1. ⬜ Neural embeddings or multi-task (based on Phase 1 results)
2. ⬜ Production deployment
3. ⬜ Documentation & training materials
4. ⬜ Monitor performance & iterate

---

## Files to Create

```
chemtools/ml/
├── __init__.py
├── feature_engineering.py  # Feature encoding pipeline
├── yield_predictor.py      # XGBoost/LightGBM wrapper
├── data_loader.py          # PyTorch Dataset
└── evaluation.py           # Metrics, plots

scripts/ml/
├── prepare_dataset.py      # Train/val/test split
├── train_yield_model.py    # Training script
└── evaluate_models.py      # Benchmark suite

models/
└── yield_v1.pkl            # Trained model (gitignored)

tests/ml/
└── test_yield_predictor.py # Unit tests
```

---

## Questions?

**Technical:** See full analysis in `docs/ML_Recommendation_Analysis_and_Plan.md`  
**Implementation:** Phase 1 code examples included  
**Timeline:** 12 weeks total, can start with 2-week Phase 1 pilot

**Ready to start?** Begin with Phase 1 baseline (2 weeks, ~$200 budget)
