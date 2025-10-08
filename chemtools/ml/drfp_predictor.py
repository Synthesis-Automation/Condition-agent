"""
DRFP-based yield predictor using LightGBM (Phase 1).

This predictor uses differential reaction fingerprints (DRFP) to encode
the reaction transformation and one-hot encoded condition features to
predict reaction yield. It avoids hand-crafted features entirely.

Example:
    from chemtools.ml.drfp_predictor import DRFPYieldPredictor
    
    predictor = DRFPYieldPredictor()
    predictor.fit(train_df)
    predictions = predictor.predict(test_df)
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Dict, List, Tuple, Any
import sys
import warnings

import numpy as np

try:
    from drfp import DrfpEncoder
except ImportError:
    print("ERROR: drfp not installed. Install with: pip install drfp")
    sys.exit(1)

try:
    import lightgbm as lgb
except ImportError:
    print("ERROR: lightgbm not installed. Install with: pip install lightgbm")
    sys.exit(1)


class DRFPYieldPredictor:
    """
    Yield predictor using DRFP fingerprints and LightGBM.
    
    Features:
        - DRFP (4096 bits): Differential reaction fingerprint
        - One-hot encoded core, base, solvent
        - Numerical T_C, time_h
    
    Model: LightGBM gradient boosting regressor
    
    Target: MAE < 15%, inference < 100ms
    """
    
    def __init__(
        self,
        n_bits: int = 4096,
        radius: int = 3,
        n_estimators: int = 1000,
        max_depth: int = 8,
        learning_rate: float = 0.05,
        random_state: int = 42,
    ):
        """
        Initialize DRFP yield predictor.
        
        Args:
            n_bits: DRFP fingerprint size (default: 4096)
            radius: DRFP radius (default: 3)
            n_estimators: Number of boosting rounds (default: 1000)
            max_depth: Maximum tree depth (default: 8)
            learning_rate: Learning rate (default: 0.05)
            random_state: Random seed (default: 42)
        """
        self.n_bits = n_bits
        self.radius = radius
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.random_state = random_state
        
        # DRFP encoder (Note: DRFP 0.3.x uses fixed 2048 bits, ignoring n_bits parameter)
        # We'll post-process to target n_bits if needed
        self.drfp_encoder = DrfpEncoder()
        
        # LightGBM model
        self.model = None
        
        # Vocabularies for one-hot encoding (learned from training data)
        self.core_vocab: List[str] = []
        self.base_vocab: List[str] = []
        self.solvent_vocab: List[str] = []
        
        # Training statistics
        self.train_stats: Dict[str, Any] = {}
    
    def _encode_reaction(self, reaction_smiles: str) -> np.ndarray:
        """
        Encode reaction SMILES to DRFP fingerprint.
        
        Args:
            reaction_smiles: Reaction SMILES (reactants>>products)
        
        Returns:
            DRFP fingerprint as binary array (n_bits,)
        """
        try:
            fp = self.drfp_encoder.encode([reaction_smiles])[0]
            fp_array = np.array(fp, dtype=np.float32)
            
            # DRFP 0.3.x returns 2048 bits by default
            # If user requested different n_bits, we'll use actual size
            actual_size = len(fp_array)
            if actual_size != self.n_bits:
                warnings.warn(
                    f"DRFP returned {actual_size} bits, using actual size (requested {self.n_bits} ignored)"
                )
                self.n_bits = actual_size  # Update to actual size
            
            return fp_array
        except Exception as e:
            warnings.warn(f"DRFP encoding failed for {reaction_smiles}: {e}")
            return np.zeros(self.n_bits, dtype=np.float32)
    
    def _one_hot(self, value: str, vocab: List[str]) -> np.ndarray:
        """
        One-hot encode categorical value.
        
        Args:
            value: Categorical value
            vocab: Vocabulary list
        
        Returns:
            One-hot vector (len(vocab),)
        """
        vec = np.zeros(len(vocab), dtype=np.float32)
        if value in vocab:
            idx = vocab.index(value)
            vec[idx] = 1.0
        return vec
    
    def _encode_conditions(self, row: Dict[str, Any]) -> np.ndarray:
        """
        Encode reaction conditions to feature vector.
        
        Args:
            row: Reaction data with core, base_uid, solvent_uid, T_C, time_h
        
        Returns:
            Condition feature vector
        """
        features = []
        
        # One-hot encoded categorical features
        features.append(self._one_hot(row.get('core', 'Unknown'), self.core_vocab))
        features.append(self._one_hot(row.get('base_uid', 'Unknown'), self.base_vocab))
        features.append(self._one_hot(row.get('solvent_uid', 'Unknown'), self.solvent_vocab))
        
        # Numerical features (normalized)
        temp = (row.get('T_C', 80.0) - self.train_stats.get('temp_mean', 80.0)) / max(
            self.train_stats.get('temp_std', 20.0), 1e-6
        )
        time = (row.get('time_h', 12.0) - self.train_stats.get('time_mean', 12.0)) / max(
            self.train_stats.get('time_std', 5.0), 1e-6
        )
        
        features.append(np.array([temp, time], dtype=np.float32))
        
        return np.concatenate(features)
    
    def _build_vocabularies(self, df) -> None:
        """
        Build vocabularies from training data.
        
        Args:
            df: Training dataframe with core, base_uid, solvent_uid columns
        """
        # Get unique values, sorted for reproducibility
        self.core_vocab = sorted(df['core'].unique().tolist())
        self.base_vocab = sorted(df['base_uid'].unique().tolist())
        self.solvent_vocab = sorted(df['solvent_uid'].unique().tolist())
        
        print(f"Vocabularies built:")
        print(f"  Cores:    {len(self.core_vocab)}")
        print(f"  Bases:    {len(self.base_vocab)}")
        print(f"  Solvents: {len(self.solvent_vocab)}")
    
    def _compute_train_stats(self, df) -> None:
        """
        Compute training set statistics for normalization.
        
        Args:
            df: Training dataframe with T_C, time_h columns
        """
        self.train_stats = {
            'temp_mean': df['T_C'].mean(),
            'temp_std': df['T_C'].std(),
            'time_mean': df['time_h'].mean(),
            'time_std': df['time_h'].std(),
            'yield_mean': df['yield_pct'].mean(),
            'yield_std': df['yield_pct'].std(),
        }
        
        print(f"Training statistics:")
        print(f"  T_C:   {self.train_stats['temp_mean']:.1f} ± {self.train_stats['temp_std']:.1f}")
        print(f"  time_h: {self.train_stats['time_mean']:.1f} ± {self.train_stats['time_std']:.1f}")
        print(f"  yield:  {self.train_stats['yield_mean']:.1f} ± {self.train_stats['yield_std']:.1f}")
    
    def _featurize(self, df) -> np.ndarray:
        """
        Featurize dataframe to model inputs.
        
        Args:
            df: Dataframe with reaction_smiles and condition columns
        
        Returns:
            Feature matrix (n_samples, n_features)
        """
        features = []
        
        print(f"Featurizing {len(df)} reactions...")
        for idx, row in df.iterrows():
            # DRFP fingerprint
            drfp = self._encode_reaction(row['reaction_smiles'])
            
            # Condition features
            cond = self._encode_conditions(row)
            
            # Concatenate
            features.append(np.concatenate([drfp, cond]))
        
        return np.array(features, dtype=np.float32)
    
    def fit(
        self,
        train_df,
        val_df=None,
        verbose: bool = True,
    ) -> Dict[str, float]:
        """
        Train yield predictor on training data.
        
        Args:
            train_df: Training dataframe with reaction_smiles, conditions, yield_pct
            val_df: Optional validation dataframe for early stopping
            verbose: Print training progress
        
        Returns:
            Training metrics
        """
        print("=" * 70)
        print("Training DRFP Yield Predictor")
        print("=" * 70)
        
        # Build vocabularies and stats
        self._build_vocabularies(train_df)
        self._compute_train_stats(train_df)
        print()
        
        # Featurize
        X_train = self._featurize(train_df)
        y_train = train_df['yield_pct'].values
        
        print(f"Feature matrix: {X_train.shape}")
        print()
        
        # Validation set
        eval_set = None
        if val_df is not None:
            X_val = self._featurize(val_df)
            y_val = val_df['yield_pct'].values
            eval_set = [(X_val, y_val)]
            print(f"Validation set: {X_val.shape}")
            print()
        
        # LightGBM parameters
        params = {
            'objective': 'regression',
            'metric': 'mae',
            'boosting_type': 'gbdt',
            'n_estimators': self.n_estimators,
            'max_depth': self.max_depth,
            'learning_rate': self.learning_rate,
            'num_leaves': 2 ** self.max_depth - 1,
            'min_child_samples': 20,
            'subsample': 0.8,
            'colsample_bytree': 0.8,
            'random_state': self.random_state,
            'verbosity': 0 if not verbose else 1,
        }
        
        print("Training LightGBM model...")
        print(f"Parameters: {params}")
        print()
        
        # Train
        self.model = lgb.LGBMRegressor(**params)
        
        self.model.fit(
            X_train,
            y_train,
            eval_set=eval_set,
            eval_metric='mae',
            callbacks=[
                lgb.early_stopping(stopping_rounds=50, verbose=verbose),
                lgb.log_evaluation(period=100 if verbose else 0),
            ] if eval_set else None,
        )
        
        # Training metrics
        y_pred_train = self.model.predict(X_train)
        train_mae = np.mean(np.abs(y_train - y_pred_train))
        train_rmse = np.sqrt(np.mean((y_train - y_pred_train) ** 2))
        
        metrics = {
            'train_mae': train_mae,
            'train_rmse': train_rmse,
        }
        
        if val_df is not None:
            y_pred_val = self.model.predict(X_val)
            val_mae = np.mean(np.abs(y_val - y_pred_val))
            val_rmse = np.sqrt(np.mean((y_val - y_pred_val) ** 2))
            metrics['val_mae'] = val_mae
            metrics['val_rmse'] = val_rmse
        
        print()
        print("Training complete!")
        print(f"  Train MAE:  {metrics['train_mae']:.2f}%")
        print(f"  Train RMSE: {metrics['train_rmse']:.2f}%")
        if val_df is not None:
            print(f"  Val MAE:    {metrics['val_mae']:.2f}%")
            print(f"  Val RMSE:   {metrics['val_rmse']:.2f}%")
        print()
        
        return metrics
    
    def predict(self, df) -> np.ndarray:
        """
        Predict yields for reactions.
        
        Args:
            df: Dataframe with reaction_smiles and condition columns
        
        Returns:
            Predicted yields (n_samples,)
        """
        if self.model is None:
            raise ValueError("Model not trained. Call fit() first.")
        
        X = self._featurize(df)
        return self.model.predict(X)
    
    def save(self, filepath: str) -> None:
        """
        Save model to disk.
        
        Args:
            filepath: Output path (.pkl)
        """
        path = Path(filepath)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(path, 'wb') as f:
            pickle.dump({
                'drfp_encoder': self.drfp_encoder,
                'model': self.model,
                'core_vocab': self.core_vocab,
                'base_vocab': self.base_vocab,
                'solvent_vocab': self.solvent_vocab,
                'train_stats': self.train_stats,
                'n_bits': self.n_bits,
                'radius': self.radius,
            }, f)
        
        print(f"Model saved to: {path}")
    
    @classmethod
    def load(cls, filepath: str) -> DRFPYieldPredictor:
        """
        Load model from disk.
        
        Args:
            filepath: Model path (.pkl)
        
        Returns:
            Loaded predictor
        """
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        
        predictor = cls(n_bits=data['n_bits'], radius=data['radius'])
        predictor.drfp_encoder = data['drfp_encoder']
        predictor.model = data['model']
        predictor.core_vocab = data['core_vocab']
        predictor.base_vocab = data['base_vocab']
        predictor.solvent_vocab = data['solvent_vocab']
        predictor.train_stats = data['train_stats']
        
        return predictor
