# ML Recommendation: End-to-End Learning Plan (Updated)

**Generated:** 2025-10-06  
**Update:** Prioritize yield prediction with learned representations (avoid hand-crafted features)  
**Timeline:** 8-10 weeks

---

## Executive Summary

**Key Changes from Original Plan:**
- ❌ **Skip hand-crafted categorical features** (LG, nuc_class, ortho_count, etc.)
- ✅ **Use learned molecular representations** from the start (DRFP, GNN, or pretrained embeddings)
- ✅ **Focus on yield prediction** as primary objective
- ✅ **Faster path to production** (8-10 weeks vs 12 weeks)

**New Architecture:**
```
Reaction SMILES
  → Learned Representations (DRFP/GNN/Transformer)
  → Neural Network (yield + conditions)
  → Predictions with uncertainty
```

---

## Why Avoid Hand-Crafted Features?

### Problems with Current Categorical Features

**Current system:**
```python
features = {
    "LG": "Br",                    # ❌ Doesn't generalize (Suzuki ≠ Ullmann)
    "elec_class": "aryl",          # ❌ Too coarse (loses structural info)
    "nuc_class": "aniline",        # ❌ Binary categories miss subtleties
    "ortho_count": 1,              # ❌ Ignores 3D geometry
    "para_EWG": False,             # ❌ Hard to define for all substrates
    "bin": "LG:Br|NUC:aniline"     # ❌ Family-specific
}
```

**Issues:**
1. **Family-specific:** Ullmann features don't work for Suzuki, Amide, etc.
2. **Information loss:** Reduces complex molecules to 5-10 categories
3. **Hard to maintain:** Requires expert rules for every new reaction type
4. **Limited expressiveness:** Can't capture steric/electronic effects accurately
5. **No transfer learning:** Can't leverage pretrained models

### Benefits of Learned Representations

**Molecular fingerprints (DRFP, ECFP, Morgan):**
- ✅ Universal across all reaction types
- ✅ Captures substructure patterns automatically
- ✅ No expert rules needed
- ✅ Easy to extend to new reactions

**Graph Neural Networks (GNN):**
- ✅ Learns from atom connectivity directly
- ✅ Captures 3D geometry (if using 3D GNN)
- ✅ Transfer learning from pretrained models (ChemBERTa, MolCLR)
- ✅ End-to-end differentiable

**Transformer embeddings:**
- ✅ State-of-the-art on benchmarks
- ✅ Pretrained on millions of molecules
- ✅ Handles reaction context naturally
- ✅ Easy fine-tuning

---

## Updated 3-Phase Architecture

### Phase 1: DRFP + Gradient Boosting (2-3 weeks)

**Simplest path to production - No hand-crafted features!**

```python
# chemtools/ml/drfp_predictor.py

from drfp import DrfpEncoder
from lightgbm import LGBMRegressor
import numpy as np

class DRFPYieldPredictor:
    """Yield prediction using DRFP fingerprints only."""
    
    def __init__(self, n_bits=4096, radius=3):
        self.encoder = DrfpEncoder()
        self.n_bits = n_bits
        self.radius = radius
        
        # Gradient boosting on DRFP features
        self.model = LGBMRegressor(
            n_estimators=1000,
            max_depth=8,
            learning_rate=0.03,
            num_leaves=64,
            min_child_samples=20,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_alpha=0.1,
            reg_lambda=0.1,
            random_state=42,
        )
        
        # Optional: condition embeddings (learned)
        self.condition_encoder = None  # Will be learned
    
    def _encode_reaction(self, reaction_smiles):
        """Encode reaction to DRFP fingerprint."""
        fp = self.encoder.encode(
            reaction_smiles,
            n_folded_length=self.n_bits,
            radius=self.radius,
            rings=True,
            mapping=False,
        )
        return np.array(fp[0] if isinstance(fp, list) else fp, dtype=np.float32)
    
    def _encode_conditions(self, core, base, solvent, temp, time):
        """Encode conditions using learned embeddings or one-hot."""
        # Option 1: One-hot (simple baseline)
        core_vec = self._one_hot(core, self.core_vocab_)
        base_vec = self._one_hot(base, self.base_vocab_)
        solv_vec = self._one_hot(solvent, self.solv_vocab_)
        
        # Normalize numerical conditions
        temp_norm = (temp - 80) / 50  # Typical range: 20-150°C
        time_norm = np.log1p(time) / 5  # Log scale, typical: 0.1-48h
        
        # Concatenate all
        return np.concatenate([
            core_vec, base_vec, solv_vec, 
            [temp_norm, time_norm]
        ])
    
    def fit(self, reactions_df):
        """
        Train on reaction dataset.
        
        reactions_df columns:
        - reaction_smiles: str (reactants>>products)
        - core: str
        - base_uid: str
        - solvent_uid: str
        - T_C: float
        - time_h: float
        - yield_pct: float
        """
        # Build vocabularies
        self.core_vocab_ = sorted(reactions_df['core'].unique())
        self.base_vocab_ = sorted(reactions_df['base_uid'].unique())
        self.solv_vocab_ = sorted(reactions_df['solvent_uid'].unique())
        
        # Encode all reactions
        X_list = []
        for _, row in reactions_df.iterrows():
            # Reaction fingerprint
            rxn_fp = self._encode_reaction(row['reaction_smiles'])
            
            # Condition features
            cond_vec = self._encode_conditions(
                row['core'],
                row['base_uid'],
                row['solvent_uid'],
                row['T_C'],
                row['time_h']
            )
            
            # Concatenate
            X_list.append(np.concatenate([rxn_fp, cond_vec]))
        
        X = np.vstack(X_list)
        y = reactions_df['yield_pct'].values
        
        # Train model
        self.model.fit(X, y)
        
        return self
    
    def predict(self, reaction_smiles, conditions_list):
        """
        Predict yield for reaction with different conditions.
        
        Args:
            reaction_smiles: str
            conditions_list: List[Dict] with keys: core, base_uid, solvent_uid, T_C, time_h
        
        Returns:
            predictions: List[Tuple[float, float]] - (yield, uncertainty)
        """
        rxn_fp = self._encode_reaction(reaction_smiles)
        
        X_list = []
        for cond in conditions_list:
            cond_vec = self._encode_conditions(
                cond['core'],
                cond['base_uid'],
                cond['solvent_uid'],
                cond['T_C'],
                cond['time_h']
            )
            X_list.append(np.concatenate([rxn_fp, cond_vec]))
        
        X = np.vstack(X_list)
        
        # Predict yields
        y_pred = self.model.predict(X)
        
        # Estimate uncertainty from tree ensemble
        # (LightGBM doesn't have built-in uncertainty, use quantile regression)
        # For now, use simple heuristic: uncertainty = 1.0 / num_training_samples_in_leaf
        
        # Better: train quantile regression models for 5th and 95th percentiles
        uncertainty = np.full_like(y_pred, 10.0)  # Placeholder: ±10%
        
        return list(zip(y_pred, uncertainty))
    
    def _one_hot(self, value, vocab):
        vec = np.zeros(len(vocab), dtype=np.float32)
        if value in vocab:
            vec[vocab.index(value)] = 1.0
        return vec
```

**Training script:**
```python
# scripts/ml/train_drfp_model.py

import pandas as pd
from chemtools.ml.drfp_predictor import DRFPYieldPredictor
import joblib

def main():
    # Load dataset
    df = pd.read_json('data/ml_train.jsonl', lines=True)
    
    # Filter out low-quality data
    df = df[df['yield_pct'] > 0]  # Remove failed reactions
    df = df[df['yield_pct'] <= 100]  # Remove invalid yields
    
    # Train model
    model = DRFPYieldPredictor(n_bits=4096, radius=3)
    model.fit(df)
    
    # Save
    joblib.dump(model, 'models/drfp_yield_v1.pkl')
    
    # Evaluate on validation set
    val_df = pd.read_json('data/ml_val.jsonl', lines=True)
    evaluate(model, val_df)

def evaluate(model, val_df):
    y_true = []
    y_pred = []
    
    for _, row in val_df.iterrows():
        pred, _ = model.predict(
            row['reaction_smiles'],
            [row[['core', 'base_uid', 'solvent_uid', 'T_C', 'time_h']].to_dict()]
        )[0]
        
        y_true.append(row['yield_pct'])
        y_pred.append(pred)
    
    from sklearn.metrics import mean_absolute_error, r2_score
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    
    print(f"Validation MAE: {mae:.2f}%")
    print(f"Validation R²: {r2:.3f}")

if __name__ == '__main__':
    main()
```

**Expected performance:**
- **MAE: 12-15%** (good baseline without hand-crafted features)
- **R²: 0.50-0.65** (reasonable correlation)
- **Training time: 5-10 minutes** (on CPU)
- **Inference: <100ms** per reaction

**Advantages:**
- ✅ No hand-crafted features
- ✅ Works across all reaction families
- ✅ Fast training and inference
- ✅ Easy to interpret (SHAP on DRFP bits)

**Limitations:**
- DRFP is a fixed representation (not learned end-to-end)
- One-hot conditions don't capture chemical similarity

---

### Phase 2: Neural Network with Learned Embeddings (3-4 weeks)

**Replace one-hot conditions with learned embeddings**

```python
# chemtools/ml/neural_yield_predictor.py

import torch
import torch.nn as nn
from drfp import DrfpEncoder

class NeuralYieldPredictor(nn.Module):
    """
    Neural network for yield prediction.
    
    Architecture:
    - Reaction: DRFP fingerprint (4096 bits)
    - Core/Base/Solvent: Learned embeddings (64 dims each)
    - T/time: Numeric inputs
    - MLP: 3 hidden layers with residual connections
    """
    
    def __init__(self, n_cores, n_bases, n_solvents, 
                 drfp_dim=4096, embed_dim=64, hidden_dim=512):
        super().__init__()
        
        # Condition embeddings (learned representations)
        self.core_embed = nn.Embedding(n_cores, embed_dim)
        self.base_embed = nn.Embedding(n_bases, embed_dim)
        self.solv_embed = nn.Embedding(n_solvents, embed_dim)
        
        # Input projection
        input_dim = drfp_dim + 3 * embed_dim + 2  # DRFP + 3 embeddings + T + time
        self.input_proj = nn.Linear(input_dim, hidden_dim)
        
        # MLP with residual connections
        self.layers = nn.ModuleList([
            ResidualBlock(hidden_dim, hidden_dim)
            for _ in range(3)
        ])
        
        # Output heads
        self.yield_head = nn.Linear(hidden_dim, 1)
        self.uncertainty_head = nn.Linear(hidden_dim, 1)  # Log variance
        
        # Batch normalization
        self.bn = nn.BatchNorm1d(hidden_dim)
        
        # Dropout for uncertainty
        self.dropout = nn.Dropout(0.2)
    
    def forward(self, drfp, core_idx, base_idx, solv_idx, temp, time):
        """
        Args:
            drfp: (batch, 4096) - DRFP fingerprints
            core_idx: (batch,) - Core indices
            base_idx: (batch,) - Base indices
            solv_idx: (batch,) - Solvent indices
            temp: (batch,) - Temperature (normalized)
            time: (batch,) - Time (log-normalized)
        
        Returns:
            yield_pred: (batch,) - Predicted yield
            log_var: (batch,) - Log variance (uncertainty)
        """
        # Embed conditions
        core_emb = self.core_embed(core_idx)
        base_emb = self.base_embed(base_idx)
        solv_emb = self.solv_embed(solv_idx)
        
        # Concatenate all features
        x = torch.cat([
            drfp,
            core_emb,
            base_emb,
            solv_emb,
            temp.unsqueeze(1),
            time.unsqueeze(1)
        ], dim=1)
        
        # Forward through network
        x = self.input_proj(x)
        x = self.bn(x)
        x = torch.relu(x)
        
        for layer in self.layers:
            x = layer(x)
        
        x = self.dropout(x)
        
        # Predict yield and uncertainty
        yield_pred = self.yield_head(x).squeeze(1)
        log_var = self.uncertainty_head(x).squeeze(1)
        
        return yield_pred, log_var
    
    def loss(self, yield_pred, log_var, yield_true):
        """
        Negative log-likelihood loss with uncertainty.
        
        Assumes Gaussian likelihood: p(y | x) = N(y; μ, σ²)
        Loss = 0.5 * (log(2π) + log(σ²) + (y - μ)² / σ²)
        """
        # Clamp log_var to avoid numerical issues
        log_var = torch.clamp(log_var, min=-5, max=5)
        
        # NLL loss
        mse = (yield_pred - yield_true) ** 2
        loss = 0.5 * (log_var + mse / torch.exp(log_var))
        
        return loss.mean()


class ResidualBlock(nn.Module):
    def __init__(self, dim, hidden_dim=None):
        super().__init__()
        if hidden_dim is None:
            hidden_dim = dim
        
        self.layers = nn.Sequential(
            nn.Linear(dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, dim),
        )
        self.ln = nn.LayerNorm(dim)
    
    def forward(self, x):
        return self.ln(x + self.layers(x))
```

**Training:**
```python
# scripts/ml/train_neural_model.py

import torch
from torch.utils.data import DataLoader
from chemtools.ml.neural_yield_predictor import NeuralYieldPredictor

def train_epoch(model, dataloader, optimizer):
    model.train()
    total_loss = 0
    
    for batch in dataloader:
        drfp, core, base, solv, temp, time, yield_true = batch
        
        # Forward
        yield_pred, log_var = model(drfp, core, base, solv, temp, time)
        
        # Loss
        loss = model.loss(yield_pred, log_var, yield_true)
        
        # Backward
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        
        total_loss += loss.item()
    
    return total_loss / len(dataloader)

def main():
    # Load data
    train_dataset = ReactionDataset('data/ml_train.jsonl')
    val_dataset = ReactionDataset('data/ml_val.jsonl')
    
    train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=128)
    
    # Create model
    model = NeuralYieldPredictor(
        n_cores=len(train_dataset.core_vocab),
        n_bases=len(train_dataset.base_vocab),
        n_solvents=len(train_dataset.solv_vocab),
    )
    
    # Optimizer with weight decay
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=1e-3,
        weight_decay=1e-4
    )
    
    # Learning rate scheduler
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=5
    )
    
    # Training loop
    best_val_loss = float('inf')
    for epoch in range(100):
        train_loss = train_epoch(model, train_loader, optimizer)
        val_loss = evaluate(model, val_loader)
        
        scheduler.step(val_loss)
        
        print(f"Epoch {epoch}: train_loss={train_loss:.4f}, val_loss={val_loss:.4f}")
        
        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), 'models/neural_yield_v1.pt')
            print(f"  → Saved new best model")
        
        # Early stopping
        if optimizer.param_groups[0]['lr'] < 1e-6:
            print("Learning rate too small, stopping")
            break

if __name__ == '__main__':
    main()
```

**Expected performance:**
- **MAE: 10-12%** (better than DRFP + LightGBM)
- **R²: 0.60-0.70** (learned embeddings capture chemical similarity)
- **Calibrated uncertainty** (from NLL loss)
- **Training time: 1-2 hours** (on GPU)

**Advantages:**
- ✅ Learned condition embeddings (e.g., similar cores have similar embeddings)
- ✅ Calibrated uncertainty from likelihood
- ✅ End-to-end differentiable
- ✅ Can add more components (catalysts, additives) easily

---

### Phase 3: Molecular Graph Neural Network (3-4 weeks)

**Replace DRFP with GNN for even better representations**

```python
# chemtools/ml/gnn_yield_predictor.py

import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool
from rdkit import Chem

class MolecularGNN(nn.Module):
    """
    Graph Neural Network for molecular encoding.
    
    Uses atom and bond features directly (no DRFP).
    """
    
    def __init__(self, hidden_dim=256, num_layers=4):
        super().__init__()
        
        # Atom feature embedding
        self.atom_embed = nn.Linear(9, hidden_dim)  # 9 atom features
        
        # Graph convolution layers
        self.convs = nn.ModuleList([
            GCNConv(hidden_dim, hidden_dim)
            for _ in range(num_layers)
        ])
        
        self.batch_norms = nn.ModuleList([
            nn.BatchNorm1d(hidden_dim)
            for _ in range(num_layers)
        ])
    
    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        
        # Embed atom features
        x = self.atom_embed(x)
        
        # Graph convolutions with residual connections
        for conv, bn in zip(self.convs, self.batch_norms):
            x_new = conv(x, edge_index)
            x_new = bn(x_new)
            x_new = torch.relu(x_new)
            x = x + x_new  # Residual
        
        # Global pooling
        x = global_mean_pool(x, batch)
        
        return x


class ReactionGNN(nn.Module):
    """
    Encode reaction as difference of product and reactant graphs.
    """
    
    def __init__(self, hidden_dim=256):
        super().__init__()
        self.mol_encoder = MolecularGNN(hidden_dim)
    
    def forward(self, reactants_data, products_data):
        # Encode reactants and products
        r_emb = self.mol_encoder(reactants_data)
        p_emb = self.mol_encoder(products_data)
        
        # Reaction embedding as difference
        rxn_emb = p_emb - r_emb
        
        return rxn_emb


class GNNYieldPredictor(nn.Module):
    """Full model with GNN reaction encoder + condition embeddings."""
    
    def __init__(self, n_cores, n_bases, n_solvents):
        super().__init__()
        
        self.rxn_encoder = ReactionGNN(hidden_dim=256)
        
        # Condition embeddings
        self.core_embed = nn.Embedding(n_cores, 64)
        self.base_embed = nn.Embedding(n_bases, 64)
        self.solv_embed = nn.Embedding(n_solvents, 64)
        
        # Prediction head
        input_dim = 256 + 3*64 + 2  # Rxn + conditions + T/time
        self.mlp = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, 2),  # [yield, log_var]
        )
    
    def forward(self, reactants, products, core, base, solv, temp, time):
        # Encode reaction
        rxn_emb = self.rxn_encoder(reactants, products)
        
        # Encode conditions
        core_emb = self.core_embed(core)
        base_emb = self.base_embed(base)
        solv_emb = self.solv_embed(solv)
        
        # Concatenate
        x = torch.cat([
            rxn_emb,
            core_emb,
            base_emb,
            solv_emb,
            temp.unsqueeze(1),
            time.unsqueeze(1)
        ], dim=1)
        
        # Predict
        out = self.mlp(x)
        yield_pred = out[:, 0]
        log_var = out[:, 1]
        
        return yield_pred, log_var
```

**Expected performance:**
- **MAE: 8-10%** (state-of-the-art)
- **R²: 0.70-0.80** (excellent correlation)
- **Training time: 3-5 hours** (on GPU)

**Advantages:**
- ✅ No fixed fingerprints - learns representations end-to-end
- ✅ Can capture 3D geometry (with 3D GNN variants)
- ✅ Transfer learning from pretrained molecular models (ChemBERTa, MolCLR)
- ✅ State-of-the-art performance

**Challenges:**
- More complex implementation (needs PyTorch Geometric)
- Requires more data for stable training (3K+ reactions recommended)
- Slower inference (~200-500ms per reaction)

---

## Updated Implementation Roadmap

### Week 1-2: Phase 1 - DRFP + LightGBM Baseline

**Goal:** Working yield predictor without hand-crafted features

- [ ] Extract training data from `data/reaction_dataset/*.jsonl`
- [ ] Create train/val/test splits (70/15/15)
- [ ] Implement `DRFPYieldPredictor` class
- [ ] Train baseline model
- [ ] **Target: MAE < 15%, inference < 100ms**

**Deliverables:**
```
models/drfp_yield_v1.pkl          # Trained model
scripts/ml/train_drfp_model.py    # Training script
chemtools/ml/drfp_predictor.py    # Model class
```

### Week 3-4: Condition Recommendation Integration

**Goal:** Use yield predictions to rank conditions

- [ ] Integrate with existing `precedent.knn()` system
- [ ] Create `HybridRecommender` that uses ML predictions
- [ ] A/B test: ML vs k-NN baseline
- [ ] Update API endpoints
- [ ] **Target: Top-5 accuracy > 65%**

**Deliverables:**
```
chemtools/recommend_ml.py         # Hybrid recommender
app/main.py                       # Updated API routes
```

### Week 5-7: Phase 2 - Neural Network with Learned Embeddings

**Goal:** Improve accuracy with learned condition representations

- [ ] Implement `NeuralYieldPredictor` (PyTorch)
- [ ] Create `ReactionDataset` class
- [ ] Train with NLL loss (uncertainty calibration)
- [ ] Hyperparameter tuning (hidden dims, learning rate, etc.)
- [ ] **Target: MAE < 12%, calibrated uncertainty**

**Deliverables:**
```
models/neural_yield_v1.pt
chemtools/ml/neural_yield_predictor.py
chemtools/ml/data_loader.py
scripts/ml/train_neural_model.py
```

### Week 8-10: Phase 3 - Graph Neural Network (Optional)

**Goal:** State-of-the-art performance with GNN

- [ ] Install PyTorch Geometric
- [ ] Implement `MolecularGNN` and `ReactionGNN`
- [ ] Preprocess molecules to graph format
- [ ] Train GNN model
- [ ] **Target: MAE < 10%, R² > 0.70**

**Optional enhancements:**
- [ ] Pretrain on public datasets (USPTO, Reaxys)
- [ ] 3D conformer generation for geometry
- [ ] Attention mechanisms for interpretability

**Deliverables:**
```
models/gnn_yield_v1.pt
chemtools/ml/gnn_yield_predictor.py
chemtools/ml/mol_to_graph.py
```

---

## Data Preparation Script

```python
# scripts/ml/prepare_dataset.py

import json
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split

def load_reaction_dataset(dataset_dir='data/reaction_dataset'):
    """Load all JSONL files from dataset directory."""
    rows = []
    
    for jsonl_file in Path(dataset_dir).glob('*.jsonl'):
        with open(jsonl_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    try:
                        row = json.loads(line)
                        rows.append(row)
                    except json.JSONDecodeError:
                        continue
    
    return rows

def extract_features(row):
    """Extract features from dataset row."""
    # Get reaction SMILES
    smiles_block = row.get('smiles', {})
    reactants = smiles_block.get('reactants', '')
    products = smiles_block.get('products', '')
    reaction_smiles = f"{reactants}>>{products}"
    
    # Get conditions
    conditions = row.get('conditions', {})
    core = row.get('condition_core', 'Unknown')
    
    # Get base (first base reagent)
    base_uid = None
    for reagent in row.get('reagents', []):
        if reagent.get('role', '').upper() == 'BASE':
            base_uid = reagent.get('cas') or reagent.get('uid') or reagent.get('name')
            break
    
    # Get solvent (first solvent)
    solvents = row.get('solvents', [])
    solvent_uid = solvents[0].get('cas') or solvents[0].get('name') if solvents else None
    
    # Get yield
    yield_pct = conditions.get('yield_pct')
    
    # Get T and time
    temp_c = conditions.get('temperature_c')
    time_h = conditions.get('time_h')
    
    return {
        'reaction_id': row.get('reaction_id'),
        'reaction_smiles': reaction_smiles,
        'reaction_type': row.get('reaction_type'),
        'core': core,
        'base_uid': base_uid,
        'solvent_uid': solvent_uid,
        'T_C': temp_c,
        'time_h': time_h,
        'yield_pct': yield_pct,
    }

def main():
    # Load data
    print("Loading dataset...")
    rows = load_reaction_dataset()
    print(f"Loaded {len(rows)} reactions")
    
    # Extract features
    print("Extracting features...")
    df = pd.DataFrame([extract_features(row) for row in rows])
    
    # Filter valid data
    df = df.dropna(subset=['reaction_smiles', 'yield_pct'])
    df = df[df['yield_pct'] > 0]
    df = df[df['yield_pct'] <= 100]
    print(f"After filtering: {len(df)} reactions")
    
    # Statistics
    print("\nDataset statistics:")
    print(f"  Reaction types: {df['reaction_type'].nunique()}")
    print(f"  Cores: {df['core'].nunique()}")
    print(f"  Bases: {df['base_uid'].nunique()}")
    print(f"  Solvents: {df['solvent_uid'].nunique()}")
    print(f"  Mean yield: {df['yield_pct'].mean():.1f}%")
    print(f"  Yield std: {df['yield_pct'].std():.1f}%")
    
    # Train/val/test split (stratified by reaction_type)
    train_df, temp_df = train_test_split(
        df, test_size=0.3, random_state=42, 
        stratify=df['reaction_type']
    )
    val_df, test_df = train_test_split(
        temp_df, test_size=0.5, random_state=42,
        stratify=temp_df['reaction_type']
    )
    
    print(f"\nSplit sizes:")
    print(f"  Train: {len(train_df)} ({len(train_df)/len(df)*100:.1f}%)")
    print(f"  Val: {len(val_df)} ({len(val_df)/len(df)*100:.1f}%)")
    print(f"  Test: {len(test_df)} ({len(test_df)/len(df)*100:.1f}%)")
    
    # Save splits
    Path('data').mkdir(exist_ok=True)
    train_df.to_json('data/ml_train.jsonl', orient='records', lines=True)
    val_df.to_json('data/ml_val.jsonl', orient='records', lines=True)
    test_df.to_json('data/ml_test.jsonl', orient='records', lines=True)
    
    print("\nSaved splits to data/ml_*.jsonl")

if __name__ == '__main__':
    main()
```

---

## Evaluation Framework

```python
# chemtools/ml/evaluation.py

import numpy as np
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
import matplotlib.pyplot as plt

def evaluate_yield_predictor(model, test_df):
    """Comprehensive evaluation of yield predictor."""
    
    y_true = []
    y_pred = []
    uncertainties = []
    
    for _, row in test_df.iterrows():
        pred, unc = model.predict(
            row['reaction_smiles'],
            [row[['core', 'base_uid', 'solvent_uid', 'T_C', 'time_h']].to_dict()]
        )[0]
        
        y_true.append(row['yield_pct'])
        y_pred.append(pred)
        uncertainties.append(unc)
    
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    uncertainties = np.array(uncertainties)
    
    # Metrics
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)
    
    # Within ±10% accuracy
    within_10 = np.mean(np.abs(y_true - y_pred) < 10) * 100
    
    # Calibration: % of true yields within predicted intervals
    lower = y_pred - 1.96 * uncertainties
    upper = y_pred + 1.96 * uncertainties
    coverage = np.mean((y_true >= lower) & (y_true <= upper)) * 100
    
    print("=" * 60)
    print("EVALUATION RESULTS")
    print("=" * 60)
    print(f"MAE:              {mae:.2f}%")
    print(f"RMSE:             {rmse:.2f}%")
    print(f"R²:               {r2:.3f}")
    print(f"Within ±10%:      {within_10:.1f}%")
    print(f"95% CI Coverage:  {coverage:.1f}% (target: 95%)")
    print("=" * 60)
    
    # Plots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # 1. Predicted vs True
    axes[0].scatter(y_true, y_pred, alpha=0.5)
    axes[0].plot([0, 100], [0, 100], 'r--', label='Perfect prediction')
    axes[0].set_xlabel('True Yield (%)')
    axes[0].set_ylabel('Predicted Yield (%)')
    axes[0].set_title(f'Predicted vs True (R²={r2:.3f})')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # 2. Residuals
    residuals = y_true - y_pred
    axes[1].scatter(y_pred, residuals, alpha=0.5)
    axes[1].axhline(0, color='r', linestyle='--')
    axes[1].set_xlabel('Predicted Yield (%)')
    axes[1].set_ylabel('Residual (True - Pred)')
    axes[1].set_title(f'Residual Plot (MAE={mae:.2f}%)')
    axes[1].grid(True, alpha=0.3)
    
    # 3. Uncertainty calibration
    axes[2].errorbar(y_pred, y_true, yerr=1.96*uncertainties, 
                     fmt='o', alpha=0.3, capsize=2)
    axes[2].plot([0, 100], [0, 100], 'r--')
    axes[2].set_xlabel('Predicted Yield (%)')
    axes[2].set_ylabel('True Yield (%)')
    axes[2].set_title(f'Uncertainty (Coverage={coverage:.1f}%)')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('evaluation_results.png', dpi=150)
    print("\nSaved plots to evaluation_results.png")
    
    return {
        'mae': mae,
        'rmse': rmse,
        'r2': r2,
        'within_10': within_10,
        'coverage': coverage,
    }
```

---

## Integration with Existing System

```python
# chemtools/recommend_ml.py

from chemtools import precedent, recommend
from chemtools.ml.drfp_predictor import DRFPYieldPredictor
import joblib

class HybridRecommender:
    """
    Hybrid recommender combining ML predictions with k-NN fallback.
    
    Strategy:
    - Use ML when sufficient precedents (>50)
    - Fall back to k-NN for sparse data
    - Always show precedents for interpretability
    """
    
    def __init__(self, model_path='models/drfp_yield_v1.pkl'):
        self.ml_model = joblib.load(model_path)
        self.use_ml_threshold = 50
    
    def recommend(self, reaction, k=50, constraint_rules=None):
        """
        Recommend conditions for reaction.
        
        Returns same format as chemtools.recommend.recommend_from_reaction
        but with additional ML predictions.
        """
        # Get base recommendation from k-NN system
        base_rec = recommend.recommend_from_reaction(
            reaction,
            k=k,
            constraint_rules=constraint_rules,
            max_variants=5
        )
        
        prec_pack = base_rec.get('precedent_pack', {})
        n_precedents = prec_pack.get('support', 0)
        
        # Decide whether to use ML
        if n_precedents >= self.use_ml_threshold:
            # Use ML predictions
            variants = base_rec.get('formatted', {}).get('recommended_conditions', [])
            
            # Predict yield for each variant
            for variant in variants:
                combo = variant.get('combo', {})
                conditions = variant.get('conditions', {})
                
                # Extract T and time
                temp_str = conditions.get('temperature', '')
                time_str = conditions.get('time', '')
                
                try:
                    temp = float(temp_str.split()[0]) if temp_str else 80.0
                except:
                    temp = 80.0
                
                try:
                    time = float(time_str.split()[0]) if time_str else 12.0
                except:
                    time = 12.0
                
                # Predict yield
                pred_yield, uncertainty = self.ml_model.predict(
                    reaction,
                    [{
                        'core': variant['summary']['core'],
                        'base_uid': combo.get('base_uid'),
                        'solvent_uid': combo.get('solvent_uid'),
                        'T_C': temp,
                        'time_h': time,
                    }]
                )[0]
                
                # Add ML predictions to variant
                variant['ml_prediction'] = {
                    'predicted_yield': round(float(pred_yield), 1),
                    'uncertainty': round(float(uncertainty), 1),
                    'confidence': round(max(0, min(1, 1 - uncertainty/50)), 3),
                }
            
            # Re-rank by predicted yield
            variants_sorted = sorted(
                variants,
                key=lambda v: v.get('ml_prediction', {}).get('predicted_yield', 0),
                reverse=True
            )
            
            # Update ranks
            for i, v in enumerate(variants_sorted, 1):
                v['rank'] = i
                v['summary']['rank'] = i
            
            base_rec['formatted']['recommended_conditions'] = variants_sorted
            base_rec['ml_method'] = 'neural_network'
            base_rec['ml_enabled'] = True
        
        else:
            # k-NN fallback
            base_rec['ml_method'] = 'knn_fallback'
            base_rec['ml_enabled'] = False
            base_rec['ml_reason'] = f'insufficient_precedents (n={n_precedents} < {self.use_ml_threshold})'
        
        return base_rec
```

**Usage:**
```python
from chemtools.recommend_ml import HybridRecommender

recommender = HybridRecommender('models/drfp_yield_v1.pkl')

result = recommender.recommend("Brc1ccccc1.Nc1ccccc1>>")

# Access ML predictions
for variant in result['formatted']['recommended_conditions']:
    print(f"Rank {variant['rank']}: "
          f"Predicted yield = {variant['ml_prediction']['predicted_yield']}% "
          f"± {variant['ml_prediction']['uncertainty']}%")
```

---

## Success Criteria (Updated)

### Phase 1 (Week 1-2): DRFP Baseline
- ✅ MAE < 15% on test set
- ✅ Inference time < 100ms per reaction
- ✅ Model size < 100MB
- ✅ Works across all reaction families (no hand-crafted features)

### Phase 2 (Week 3-7): Neural Network
- ✅ MAE < 12% on test set
- ✅ Calibrated uncertainty (95% CI coverage > 90%)
- ✅ Top-5 condition accuracy > 70%
- ✅ Learned embeddings cluster similar conditions

### Phase 3 (Week 8-10): GNN (Optional)
- ✅ MAE < 10% on test set
- ✅ R² > 0.70
- ✅ Interpretable attention weights
- ✅ Transfer learning improves cold-start

---

## Summary: What Changed?

| Aspect | Original Plan | Updated Plan |
|--------|---------------|--------------|
| **Features** | Hand-crafted (LG, nuc_class, etc.) | ❌ **Removed** - Use DRFP/GNN |
| **Phase 1** | XGBoost on categorical features | DRFP + LightGBM (universal) |
| **Phase 2** | Neural embeddings (after Phase 1) | Neural network with learned condition embeddings |
| **Phase 3** | Multi-task learning | GNN for end-to-end learning |
| **Timeline** | 12 weeks | **8-10 weeks** (faster) |
| **Generalization** | Family-specific | **Universal** (all reactions) |
| **Maintenance** | High (rules for each family) | **Low** (data-driven) |

---

## Next Steps

1. **Week 1:** Prepare dataset with `scripts/ml/prepare_dataset.py`
2. **Week 1-2:** Implement and train DRFP baseline
3. **Week 2:** Evaluate baseline, compare to k-NN
4. **Week 3:** Integrate ML predictions into API
5. **Week 4:** User testing and feedback
6. **Week 5+:** Decide: Continue to Phase 2 (Neural) or Phase 3 (GNN)?

**Ready to start!** 🚀
