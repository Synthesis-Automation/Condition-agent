#!/usr/bin/env python
"""
Evaluate trained yield predictor on test set.

Loads a trained model and evaluates it on test data with comprehensive
metrics and visualizations.

Usage:
    python scripts/ml/evaluate_models.py
    python scripts/ml/evaluate_models.py --model models/drfp_yield_v1.pkl --test data/ml_test.jsonl --output results/
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas not installed. Install with: pip install pandas")
    sys.exit(1)

try:
    from chemtools.ml.drfp_predictor import DRFPYieldPredictor
    from chemtools.ml.evaluation import evaluate_yield_predictor
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
        description="Evaluate trained yield predictor on test set"
    )
    parser.add_argument(
        '--model',
        default='models/drfp_yield_v1.pkl',
        help='Trained model path (.pkl)'
    )
    parser.add_argument(
        '--test',
        default='data/ml_test.jsonl',
        help='Test data (JSONL)'
    )
    parser.add_argument(
        '--output',
        default='results',
        help='Output directory for plots and metrics'
    )
    
    args = parser.parse_args(argv)
    
    print("=" * 70)
    print("Model Evaluation")
    print("=" * 70)
    print(f"Model: {args.model}")
    print(f"Test:  {args.test}")
    print(f"Output: {args.output}")
    print()
    
    # Check files exist
    if not Path(args.model).exists():
        print(f"ERROR: Model not found: {args.model}")
        print()
        print("Train a model first:")
        print("  python scripts/ml/train_drfp_model.py")
        return 1
    
    if not Path(args.test).exists():
        print(f"ERROR: Test data not found: {args.test}")
        print()
        print("Prepare dataset first:")
        print("  python scripts/ml/prepare_dataset.py")
        return 1
    
    # Load model
    print("Loading model...")
    predictor = DRFPYieldPredictor.load(args.model)
    print(f"  Model loaded successfully")
    print()
    
    # Load test data
    print("Loading test data...")
    test_df = load_jsonl(args.test)
    print(f"  Test set: {len(test_df)} reactions")
    print()
    
    # Evaluate
    results = evaluate_yield_predictor(
        predictor,
        test_df,
        output_dir=args.output,
        verbose=True,
    )
    
    # Save metrics to JSON
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    metrics_file = output_path / 'test_metrics.json'
    with open(metrics_file, 'w') as f:
        json.dump(results['metrics'], f, indent=2)
    
    print(f"Metrics saved to: {metrics_file}")
    print()
    
    # Check if target MAE achieved
    target_mae = 15.0
    if results['metrics']['mae'] <= target_mae:
        print(f"✓ Target MAE achieved: {results['metrics']['mae']:.2f}% <= {target_mae}%")
    else:
        print(f"⚠ Target MAE not achieved: {results['metrics']['mae']:.2f}% > {target_mae}%")
    print()
    
    print("=" * 70)
    print("✓ Evaluation complete!")
    print("=" * 70)
    print()
    print(f"Results saved in: {args.output}/")
    print("  - evaluation_results.png (plots)")
    print("  - test_metrics.json (metrics)")
    print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
