#!/usr/bin/env python
"""
Train DRFP-based yield predictor (Phase 1).

Trains a LightGBM model using DRFP fingerprints and condition features
to predict reaction yields without hand-crafted features.

Usage:
    python scripts/ml/train_drfp_model.py
    python scripts/ml/train_drfp_model.py --train data/ml_train.jsonl --val data/ml_val.jsonl --output models/drfp_yield_v1.pkl
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas not installed. Install with: pip install pandas")
    sys.exit(1)

try:
    from chemtools.ml.drfp_predictor import DRFPYieldPredictor
except ImportError:
    print("ERROR: chemtools.ml not found. Check PYTHONPATH or install dependencies.")
    sys.exit(1)


def load_jsonl(filepath: str) -> pd.DataFrame:
    """Load JSONL file to DataFrame."""
    rows = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Train DRFP yield predictor"
    )
    parser.add_argument(
        '--train',
        default='data/ml_train.jsonl',
        help='Training data (JSONL)'
    )
    parser.add_argument(
        '--val',
        default='data/ml_val.jsonl',
        help='Validation data (JSONL)'
    )
    parser.add_argument(
        '--output',
        default='models/drfp_yield_v1.pkl',
        help='Output model path'
    )
    parser.add_argument(
        '--n-bits',
        type=int,
        default=4096,
        help='DRFP fingerprint size (default: 4096)'
    )
    parser.add_argument(
        '--radius',
        type=int,
        default=3,
        help='DRFP radius (default: 3)'
    )
    parser.add_argument(
        '--n-estimators',
        type=int,
        default=1000,
        help='Number of boosting rounds (default: 1000)'
    )
    parser.add_argument(
        '--max-depth',
        type=int,
        default=8,
        help='Maximum tree depth (default: 8)'
    )
    parser.add_argument(
        '--learning-rate',
        type=float,
        default=0.05,
        help='Learning rate (default: 0.05)'
    )
    parser.add_argument(
        '--random-state',
        type=int,
        default=42,
        help='Random seed (default: 42)'
    )
    
    args = parser.parse_args(argv)
    
    print("=" * 70)
    print("DRFP Yield Predictor Training")
    print("=" * 70)
    print(f"Train: {args.train}")
    print(f"Val:   {args.val}")
    print(f"Output: {args.output}")
    print()
    
    # Load data
    print("Loading data...")
    train_df = load_jsonl(args.train)
    val_df = load_jsonl(args.val) if Path(args.val).exists() else None
    
    print(f"  Train: {len(train_df)} reactions")
    if val_df is not None:
        print(f"  Val:   {len(val_df)} reactions")
    print()
    
    # Initialize predictor
    predictor = DRFPYieldPredictor(
        n_bits=args.n_bits,
        radius=args.radius,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        random_state=args.random_state,
    )
    
    # Train
    start_time = time.time()
    metrics = predictor.fit(train_df, val_df=val_df, verbose=True)
    train_time = time.time() - start_time
    
    print("=" * 70)
    print("Training Results")
    print("=" * 70)
    print(f"Training time: {train_time:.1f}s")
    print()
    print("Metrics:")
    for key, value in metrics.items():
        print(f"  {key:15s} {value:.2f}%")
    print()
    
    # Save model
    predictor.save(args.output)
    print()
    
    # Quick inference test
    print("Inference speed test...")
    start_time = time.time()
    _ = predictor.predict(train_df.head(100))
    inference_time = (time.time() - start_time) / 100 * 1000  # ms per sample
    
    print(f"  Inference time: {inference_time:.1f}ms per reaction")
    print()
    
    print("=" * 70)
    print("✓ Training complete!")
    print("=" * 70)
    print()
    print("Next steps:")
    print("  1. Evaluate on test set: python scripts/ml/evaluate_models.py")
    print("  2. Integrate with recommender: See chemtools/recommend_ml.py")
    print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
