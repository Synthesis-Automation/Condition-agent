# ML-Based Recommendation: Analysis & Improvement Plan

**Generated:** 2025-10-06  
**Status:** Analysis Complete, Implementation Roadmap Proposed

---

## Executive Summary

**Current State:** The system uses a **k-NN precedent-based approach** with categorical feature matching + optional DRFP fingerprint similarity. This is **deterministic and interpretable** but has several limitations.

**Recommendation:** Implement a **hybrid ML architecture** that combines:
1. **Gradient Boosting Models** (XGBoost/LightGBM) for yield prediction and condition ranking
2. **Neural Fingerprint Embeddings** for reaction similarity learning
3. **Multi-task Learning** for simultaneous prediction of core/base/solvent/temperature/time
4. **Keep existing k-NN as fallback** for interpretability and cold-start scenarios

**Expected Impact:**
- **15-25% improvement** in condition recommendation accuracy
- **Yield prediction** within ±10% (vs current no prediction)
- **Confidence calibration** based on model uncertainty
- **Active learning** to identify high-value experiments

---

## Current Architecture Analysis

### 1. Data Flow (`chemtools/recommend.py`)

```python
# Current workflow:
reaction_smiles 
  → normalize_reaction()         # chemtools.smiles
  → detect_family()              # chemtools.router (rule-based)
  → featurize(elec, nuc)         # chemtools.featurizers.molecular
  → precedent.knn()              # k-NN retrieval
  → vote on core/base/solvent    # Laplace smoothing
  → median T/time                # Simple statistics
  → confidence from vote_share   # Heuristic
```

**Key components:**
- `recommend_from_reaction()`: Main entry (1000+ lines)
- `recommend_conditions_structured()`: API-friendly wrapper
- `design_plate_from_reaction()`: Diversified plate design

### 2. Precedent Retrieval (`chemtools/precedent.py`)

**Current approach:**
```python
def knn(family, features, k=50, relax=None):
    # 1. Candidate pool: exact bin match → fallbacks (nuc_class, LG, any)
    cands = _candidate_pool(rows, family_txt, features, k, relax)
    
    # 2. Categorical similarity scoring
    for candidate in cands:
        sim_cat = _similarity(target_feat, candidate.features)  # 0-1 score
        
        # Optional: DRFP whole-reaction similarity
        if use_drfp:
            sim_fp = tanimoto(query_drfp, candidate_drfp)
            sim_total = blend(sim_cat, sim_fp, drfp_weight)
        
        # Yield-weighted neighbor score
        neighbor_score = sim_total * (0.5 + 0.5 * yield_norm)
    
    # 3. Return top-k precedents
    return sorted(scored)[:k]
```

**Categorical similarity weights:**
```python
weights = {
    "LG": 0.35,           # Leaving group
    "nuc_class": 0.35,    # Nucleophile class
    "ortho_count": 0.10,  # Ortho substituents
    "para_EWG": 0.10,     # Para electron-withdrawing
    "heteroaryl": 0.10,   # Heteroaromatic flag
}
```

### 3. Feature Extraction (`chemtools/featurizers/molecular.py`)

**Current features (Ullmann-focused):**
```python
{
    "LG": "Br",                    # Categorical
    "elec_class": "aryl",          # Categorical
    "nuc_class": "aniline",        # Categorical
    "ortho_count": 1,              # Numeric
    "para_EWG": False,             # Boolean
    "heteroaryl": False,           # Boolean
    "n_basicity": "aromatic_primary",  # Categorical
    "steric_alpha": "low",         # Categorical
    "bin": "LG:Br|NUC:aniline"     # Composite key
}
```

**Optional role-aware features** (`chemtools.features.role`):
- Global descriptors (MW, logP, HBA, HBD, TPSA, etc.)
- Role-specific patterns (amine/alcohol/aryl_halide SMARTS)
- Centered ECFP fingerprints (512-2048 bits)

**DRFP fingerprints** (`chemtools/reaction_similarity.py`):
- Whole-reaction difference fingerprints (4096 bits default)
- Tanimoto similarity for reaction comparison
- LRU-cached for speed (200K cache size)

### 4. Integrations

**Optional ML components:**
- `chemtools.integrations.molpipeline`: Role-based feature aggregation (sklearn pipelines)
- `chemtools.reaction_type_detector`: rxn-insight integration for classification
- Precomputed DRFP index (`artifacts/all_drfp_4096.npz`)

---

## Limitations of Current Approach

### 1. **No Yield Prediction**
- Current system only retrieves precedents but doesn't predict expected yield
- Confidence score is heuristic (vote share), not calibrated
- No uncertainty quantification

### 2. **Hand-Crafted Features**
- Ullmann-specific categorical features don't generalize to other reaction types
- Missing important structural information (e.g., steric bulk beyond categorical bins)
- DRFP is global; doesn't capture local reactivity patterns

### 3. **Simple Similarity Metrics**
- Categorical exact matching is brittle
- Fixed weights (LG: 0.35, nuc: 0.35, etc.) not learned from data
- DRFP Tanimoto doesn't account for reaction context

### 4. **Vote-Based Aggregation**
- Laplace smoothing is generic; doesn't consider substrate-specific preferences
- Median T/time ignores substrate reactivity
- No multi-objective optimization (yield vs cost vs time)

### 5. **Cold Start Problem**
- New reaction types with sparse precedents give low confidence
- No transfer learning from related reactions

### 6. **Scalability**
- O(N) similarity computation for all candidates
- DRFP precomputation helps but requires manual regeneration
- No indexing structure (e.g., LSH, FAISS)

---

## Proposed ML Architecture

### Phase 1: Gradient Boosting Models (2-3 weeks)

**Goal:** Predict yield and rank conditions using structured features

#### 1.1 Yield Prediction Model

```python
# Model: XGBoost or LightGBM regression
# Input features:
features = {
    # Substrate features (from featurizers)
    "LG": one_hot(["Cl", "Br", "I", "OTf"]),
    "elec_class": one_hot(["aryl", "vinyl", "alkyl"]),
    "nuc_class": one_hot(["aniline", "amine_primary", "amine_secondary"]),
    "ortho_count": numeric(0-2),
    "para_EWG": binary,
    "heteroaryl": binary,
    
    # Condition features
    "core": label_encode(cores),      # 50-100 unique cores
    "base": label_encode(bases),      # 20-50 bases
    "solvent": label_encode(solvents),# 30-60 solvents
    "T_C": numeric(0-150),
    "time_h": log_transform(0.1-48),
    
    # Interaction features (learned)
    "core_x_elec_class": interaction,
    "base_x_nuc_class": interaction,
    
    # DRFP embeddings (optional)
    "drfp_pca": dimensionality_reduce(4096 → 64),
}

# Target: yield_pct (0-100)
# Loss: MAE or Huber (robust to outliers)
```

**Implementation:**
```python
# chemtools/ml/yield_predictor.py
from lightgbm import LGBMRegressor
import numpy as np

class YieldPredictor:
    def __init__(self, n_estimators=500, max_depth=8):
        self.model = LGBMRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=0.05,
            objective='regression',
            metric='mae',
            random_state=42
        )
        self.feature_encoder = FeatureEncoder()
    
    def fit(self, reactions: List[Dict], yields: List[float]):
        X = self.feature_encoder.transform(reactions)
        self.model.fit(X, yields)
    
    def predict(self, reaction: Dict) -> Tuple[float, float]:
        X = self.feature_encoder.transform([reaction])
        pred = self.model.predict(X)[0]
        
        # Uncertainty from tree variance
        uncertainty = self._estimate_uncertainty(X)
        return pred, uncertainty
    
    def _estimate_uncertainty(self, X):
        # Use quantile regression or ensemble variance
        leaves = self.model.predict(X, pred_leaf=True)
        # Compute variance across trees
        return np.std(leaves)
```

#### 1.2 Condition Ranking Model

```python
# Multi-class classification for core/base/solvent
# Or: Learning-to-rank model for candidate reranking

class ConditionRanker:
    def __init__(self):
        self.core_model = LGBMClassifier(...)
        self.base_model = LGBMClassifier(...)
        self.solvent_model = LGBMClassifier(...)
    
    def rank_candidates(self, reaction: Dict, candidates: List[Dict]):
        # Score each candidate combination
        scores = []
        for cand in candidates:
            feat = self._combine_features(reaction, cand)
            
            # Predict probabilities
            p_core = self.core_model.predict_proba(feat)
            p_base = self.base_model.predict_proba(feat)
            p_solv = self.solvent_model.predict_proba(feat)
            
            # Combined score (log-likelihood)
            score = np.log(p_core + 1e-9) + \
                    np.log(p_base + 1e-9) + \
                    np.log(p_solv + 1e-9)
            
            scores.append((score, cand))
        
        return sorted(scores, reverse=True)
```

**Expected performance:**
- Yield prediction MAE: **8-12%** (literature: 10-15%)
- Top-5 accuracy for core: **70-80%**
- Inference time: **<50ms** per reaction

---

### Phase 2: Neural Fingerprint Embeddings (3-4 weeks)

**Goal:** Learn reaction representations that capture reactivity patterns

#### 2.1 Reaction Encoder

```python
# Model: Graph Neural Network or Transformer
# Architecture options:
# 1. MPNN (Message Passing Neural Network) - good for small molecules
# 2. Transformer (Molecular Transformer) - better for complex reactions
# 3. DRFP + MLP (simple baseline)

class ReactionEncoder(nn.Module):
    def __init__(self, hidden_dim=256, num_layers=4):
        super().__init__()
        # Option 1: Use pretrained molecular embeddings
        self.mol_encoder = ChemBERTa.from_pretrained()
        
        # Option 2: Learn from scratch with GNN
        # self.gnn = GraphConv(layers=num_layers, hidden=hidden_dim)
        
        # Reaction-level aggregation
        self.reaction_head = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, 128)  # Final embedding
        )
    
    def forward(self, reactants, products):
        # Encode reactants and products separately
        r_emb = self.mol_encoder(reactants)
        p_emb = self.mol_encoder(products)
        
        # Difference representation
        diff = p_emb - r_emb
        
        # Concatenate
        rxn_emb = torch.cat([r_emb.mean(0), diff.mean(0)])
        
        return self.reaction_head(rxn_emb)
```

#### 2.2 Contrastive Learning

```python
# Train embeddings to cluster similar reactions
# Use triplet loss or NT-Xent (SimCLR)

class ContrastiveLearner:
    def __init__(self, encoder):
        self.encoder = encoder
        self.temperature = 0.07
    
    def compute_loss(self, batch):
        # Anchor: original reaction
        # Positive: same core, different substrates
        # Negative: different core
        
        anchor_emb = self.encoder(batch['anchor'])
        pos_emb = self.encoder(batch['positive'])
        neg_emb = self.encoder(batch['negative'])
        
        # Cosine similarity
        sim_pos = F.cosine_similarity(anchor_emb, pos_emb)
        sim_neg = F.cosine_similarity(anchor_emb, neg_emb)
        
        # Contrastive loss
        loss = -torch.log(
            torch.exp(sim_pos / self.temperature) /
            (torch.exp(sim_pos / self.temperature) + 
             torch.exp(sim_neg / self.temperature))
        )
        
        return loss.mean()
```

**Data augmentation strategies:**
- **Chemical:** Add/remove protecting groups, vary stereochemistry
- **Condition:** Perturb T/time within ±10%, swap similar bases/solvents
- **Synthetic:** Generate reactions from retrosynthesis tools

---

### Phase 3: Multi-Task Learning (2-3 weeks)

**Goal:** Jointly optimize yield + core + base + solvent + T + time

```python
class MultiTaskRecommender(nn.Module):
    def __init__(self, reaction_encoder):
        super().__init__()
        self.encoder = reaction_encoder  # Shared backbone
        
        # Task-specific heads
        self.yield_head = nn.Linear(128, 1)
        self.core_head = nn.Linear(128, num_cores)
        self.base_head = nn.Linear(128, num_bases)
        self.solvent_head = nn.Linear(128, num_solvents)
        self.temp_head = nn.Linear(128, 1)
        self.time_head = nn.Linear(128, 1)
    
    def forward(self, reaction):
        emb = self.encoder(reaction)
        
        return {
            'yield': self.yield_head(emb),
            'core': self.core_head(emb),
            'base': self.base_head(emb),
            'solvent': self.solvent_head(emb),
            'temp': self.temp_head(emb),
            'time': self.time_head(emb),
        }
    
    def loss(self, outputs, targets, weights=None):
        # Weighted multi-task loss
        if weights is None:
            weights = {
                'yield': 1.0,
                'core': 0.5,
                'base': 0.3,
                'solvent': 0.3,
                'temp': 0.2,
                'time': 0.2,
            }
        
        loss = 0
        loss += weights['yield'] * F.mse_loss(outputs['yield'], targets['yield'])
        loss += weights['core'] * F.cross_entropy(outputs['core'], targets['core'])
        loss += weights['base'] * F.cross_entropy(outputs['base'], targets['base'])
        loss += weights['solvent'] * F.cross_entropy(outputs['solvent'], targets['solvent'])
        loss += weights['temp'] * F.mse_loss(outputs['temp'], targets['temp'])
        loss += weights['time'] * F.mse_loss(outputs['time'], targets['time'])
        
        return loss
```

**Advantages:**
- **Shared representations** improve sample efficiency
- **Regularization** from correlated tasks
- **Faster inference** (single forward pass)

---

### Phase 4: Uncertainty Quantification (1-2 weeks)

**Goal:** Calibrated confidence scores for recommendations

#### 4.1 Bayesian Neural Networks

```python
# Use dropout at inference time (MC Dropout)
class BayesianYieldPredictor:
    def predict_with_uncertainty(self, reaction, n_samples=100):
        self.model.train()  # Enable dropout
        
        predictions = []
        for _ in range(n_samples):
            pred = self.model(reaction)
            predictions.append(pred)
        
        predictions = torch.stack(predictions)
        
        mean = predictions.mean(dim=0)
        std = predictions.std(dim=0)
        
        return mean, std  # Epistemic uncertainty
```

#### 4.2 Ensemble Methods

```python
# Train multiple models with different initializations
class EnsembleRecommender:
    def __init__(self, n_models=5):
        self.models = [
            MultiTaskRecommender(...) 
            for _ in range(n_models)
        ]
    
    def predict_with_uncertainty(self, reaction):
        outputs = [model(reaction) for model in self.models]
        
        # Mean prediction
        mean_yield = torch.stack([o['yield'] for o in outputs]).mean()
        
        # Uncertainty from model disagreement
        std_yield = torch.stack([o['yield'] for o in outputs]).std()
        
        return mean_yield, std_yield
```

#### 4.3 Conformal Prediction

```python
# Calibrate prediction intervals on holdout set
class ConformalPredictor:
    def calibrate(self, val_reactions, val_yields):
        # Compute residuals on validation set
        residuals = []
        for rxn, y_true in zip(val_reactions, val_yields):
            y_pred, _ = self.model.predict(rxn)
            residuals.append(abs(y_true - y_pred))
        
        # Store quantiles
        self.quantiles = np.quantile(residuals, [0.05, 0.5, 0.95])
    
    def predict_interval(self, reaction, confidence=0.9):
        y_pred, _ = self.model.predict(reaction)
        
        # Symmetric interval
        alpha = 1 - confidence
        lower = y_pred - self.quantiles[int((1-alpha/2) * 100)]
        upper = y_pred + self.quantiles[int((1-alpha/2) * 100)]
        
        return y_pred, (lower, upper)
```

---

### Phase 5: Active Learning & Experiment Design (1-2 weeks)

**Goal:** Identify most informative experiments to run next

```python
class ActiveLearner:
    def __init__(self, recommender):
        self.recommender = recommender
    
    def suggest_next_experiments(self, candidate_reactions, n=10):
        # Acquisition functions:
        # 1. Uncertainty sampling: max(uncertainty)
        # 2. Expected improvement: max(E[max(0, y - y_best)])
        # 3. Thompson sampling: sample from posterior
        
        scores = []
        for rxn in candidate_reactions:
            yield_pred, uncertainty = self.recommender.predict_with_uncertainty(rxn)
            
            # Acquisition score
            # Option 1: Pure exploration (high uncertainty)
            acq = uncertainty
            
            # Option 2: Exploitation + exploration (UCB)
            # acq = yield_pred + 2 * uncertainty
            
            # Option 3: Expected improvement
            # acq = self._expected_improvement(yield_pred, uncertainty, best_yield)
            
            scores.append((acq, rxn))
        
        # Return top-n
        return sorted(scores, reverse=True)[:n]
    
    def _expected_improvement(self, mu, sigma, best_y, xi=0.01):
        # Probability of improvement
        from scipy.stats import norm
        
        z = (mu - best_y - xi) / (sigma + 1e-9)
        ei = (mu - best_y - xi) * norm.cdf(z) + sigma * norm.pdf(z)
        
        return ei
```

---

## Implementation Roadmap

### Week 1-2: Data Preparation & Baseline
- [ ] Extract training data from `data/reaction_dataset/*.jsonl`
- [ ] Create train/val/test splits (70/15/15) stratified by family
- [ ] Implement feature engineering pipeline
- [ ] Train baseline XGBoost yield predictor
- [ ] Evaluate on holdout set (target: MAE < 15%)

### Week 3-4: Condition Ranking
- [ ] Implement multi-class classifiers for core/base/solvent
- [ ] Add interaction features (substrate × condition)
- [ ] Hyperparameter tuning with Optuna
- [ ] A/B test against current k-NN system

### Week 5-6: Neural Embeddings
- [ ] Implement reaction encoder (MPNN or Transformer)
- [ ] Set up contrastive learning with triplet loss
- [ ] Pretrain on full dataset (unsupervised)
- [ ] Fine-tune on yield prediction task

### Week 7-8: Multi-Task Learning
- [ ] Build multi-task architecture
- [ ] Implement task-specific heads
- [ ] Train with weighted loss
- [ ] Compare single-task vs multi-task performance

### Week 9-10: Uncertainty & Active Learning
- [ ] Add MC Dropout / ensemble
- [ ] Calibrate with conformal prediction
- [ ] Implement acquisition functions
- [ ] Test on simulated experiment loop

### Week 11-12: Production Integration
- [ ] Wrap models in FastAPI endpoints
- [ ] Add caching layer (Redis)
- [ ] Update Gradio UI with ML predictions
- [ ] Write documentation & tests

---

## Evaluation Metrics

### 1. Yield Prediction
- **MAE** (Mean Absolute Error): target < 10%
- **R²**: target > 0.65
- **Calibration**: prediction intervals at 90% coverage

### 2. Condition Recommendation
- **Top-1 accuracy**: % times best condition is rank 1
- **Top-5 accuracy**: % times best condition in top 5
- **MRR** (Mean Reciprocal Rank): average 1/rank of best condition

### 3. Business Metrics
- **Lab hit rate**: % reactions with yield > 50%
- **Cost savings**: reduced screening plate sizes
- **Time to first hit**: days to find working conditions

---

## Hybrid System Architecture

```python
# chemtools/recommend_ml.py

class HybridRecommender:
    """Combines ML predictions with k-NN fallback."""
    
    def __init__(self):
        self.ml_model = MultiTaskRecommender.load('models/v1.pt')
        self.knn_engine = precedent.knn  # Existing system
        self.use_ml_threshold = 50  # Min precedents to use ML
    
    def recommend(self, reaction, k=50):
        # 1. Get precedents (always)
        prec_pack = self.knn_engine(
            family=detect_family(reaction),
            features=featurize(reaction),
            k=k
        )
        
        n_precedents = prec_pack.get('support', 0)
        
        # 2. Decide: ML or k-NN
        if n_precedents >= self.use_ml_threshold:
            # Use ML predictions
            ml_output = self.ml_model(reaction)
            
            # Rerank precedents by predicted yield
            ranked = self._rerank_by_ml(prec_pack, ml_output)
            
            return {
                'method': 'ml',
                'yield_prediction': ml_output['yield'],
                'uncertainty': ml_output['uncertainty'],
                'recommendations': ranked[:5],
                'precedents': prec_pack,
            }
        else:
            # Fall back to k-NN
            knn_rec = self._vote_based_recommendation(prec_pack)
            
            return {
                'method': 'knn_fallback',
                'reason': f'insufficient_precedents (n={n_precedents})',
                'recommendations': knn_rec,
                'precedents': prec_pack,
            }
    
    def _rerank_by_ml(self, prec_pack, ml_output):
        # Score each precedent condition by predicted yield
        precedents = prec_pack.get('precedents', [])
        scored = []
        
        for p in precedents:
            # Predict yield if we used this condition
            pred_yield = self.ml_model.predict_yield(
                reaction=...,
                core=p['core'],
                base=p['base_uid'],
                solvent=p['solvent_uid'],
                T=p['T_C'],
                time=p['time_h']
            )
            
            scored.append((pred_yield, p))
        
        return [p for _, p in sorted(scored, reverse=True)]
```

---

## Data Requirements

### Current Dataset Size
- **Ullmann C-N**: ~500 reactions (based on sample files)
- **Suzuki**: ~200 reactions (estimated)
- **Amide Formation**: ~300 reactions (estimated)

### ML Training Requirements
**Minimum:** 1,000 reactions per family
**Recommended:** 5,000+ reactions per family

### Data Augmentation Strategy
1. **Literature mining:** PubMed, USPTO patents (using tools like RxnScribe, IBM RoboRXN)
2. **ELN integration:** Export from Benchling/SciNote
3. **Vendor data:** SciFinder, Reaxys (license permitting)
4. **Synthetic generation:** Template-based enumeration

---

## Risk Mitigation

### Risk 1: Insufficient Data
**Impact:** ML models underperform k-NN baseline  
**Mitigation:**
- Use transfer learning from pretrained models (ChemBERTa, MoleculeNet)
- Data augmentation (chemical & condition perturbations)
- Semi-supervised learning with unlabeled reactions

### Risk 2: Overfitting to Conditions
**Impact:** Model memorizes dataset-specific trends  
**Mitigation:**
- Cross-validation by publication/batch
- Regularization (L2, dropout, early stopping)
- Ensemble multiple models with different architectures

### Risk 3: Computational Cost
**Impact:** Slow inference in production  
**Mitigation:**
- Model distillation (compress large model → small model)
- ONNX export for CPU inference
- Batch prediction with caching

### Risk 4: User Trust
**Impact:** Chemists don't trust "black box" predictions  
**Mitigation:**
- SHAP/LIME explanations for predictions
- Show similar precedents alongside ML predictions
- Confidence intervals and uncertainty visualization

---

## Success Criteria

### Phase 1 (Baseline ML): 
✅ Yield prediction MAE < 15%  
✅ Top-5 condition accuracy > 60%  
✅ Inference time < 100ms

### Phase 2 (Neural Embeddings):
✅ Reaction similarity captures chemical intuition (validated by chemists)  
✅ Embedding clustering aligns with reaction families  
✅ Transfer learning improves cold-start scenarios

### Phase 3 (Multi-Task):
✅ Multi-task outperforms single-task by ≥5%  
✅ All condition components predicted with >50% top-5 accuracy

### Phase 4 (Production):
✅ Hybrid system deployed to API  
✅ Gradio UI shows ML predictions + explanations  
✅ 90% of users prefer ML recommendations over k-NN (user survey)

---

## Code Structure

```
chemtools/
├── ml/
│   ├── __init__.py
│   ├── feature_engineering.py   # Encode categorical + numeric features
│   ├── yield_predictor.py       # XGBoost/LightGBM regressor
│   ├── condition_ranker.py      # Multi-class classifiers
│   ├── reaction_encoder.py      # GNN/Transformer for embeddings
│   ├── multitask_model.py       # Joint prediction architecture
│   ├── uncertainty.py           # MC Dropout, ensembles, conformal
│   ├── active_learning.py       # Acquisition functions
│   ├── data_loader.py           # PyTorch Dataset for reactions
│   ├── training.py              # Training loops, checkpointing
│   └── evaluation.py            # Metrics, visualization
├── recommend_ml.py              # HybridRecommender (ML + k-NN)
└── recommend.py                 # Existing k-NN system (keep as fallback)

models/
├── yield_predictor_v1.pkl       # Trained XGBoost
├── reaction_encoder_v1.pt       # Trained PyTorch encoder
└── multitask_v1.pt              # Multi-task model

scripts/
├── train_yield_model.py         # CLI for training
├── train_encoder.py             # Contrastive learning script
└── evaluate_models.py           # Benchmark suite

tests/
└── ml/
    ├── test_yield_predictor.py
    ├── test_reaction_encoder.py
    └── test_hybrid_recommender.py
```

---

## Next Steps

1. **Review this plan** with domain experts (chemists + ML engineers)
2. **Prioritize phases** based on business impact
3. **Gather additional data** if needed (target: 5K+ reactions)
4. **Set up ML infrastructure:**
   - GPU instance (AWS p3.2xlarge or GCP T4)
   - Experiment tracking (Weights & Biases, MLflow)
   - Model registry (DVC, MLflow)
5. **Start with Phase 1** (baseline XGBoost) - quickest ROI

---

## References

- **Yield Prediction:** Schwaller et al. (2021) "Prediction of chemical reaction yields using deep learning"
- **Reaction Embeddings:** Coley et al. (2019) "A graph-convolutional neural network model for the prediction of chemical reactivity"
- **Multi-Task Learning:** Varnek et al. (2020) "Multi-task learning for chemical reaction prediction"
- **Active Learning:** Häse et al. (2021) "Chimera: enabling hierarchy based multi-objective optimization for self-driving laboratories"

---

**Questions? Contact:** [Your team]  
**Status:** Ready for implementation  
**Estimated Timeline:** 12 weeks (3 months)  
**Expected ROI:** 15-25% improvement in recommendation accuracy
