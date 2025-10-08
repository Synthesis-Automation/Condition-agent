"""
Evaluation tools for yield predictors.

Provides metrics and visualization for assessing yield prediction performance.

Example:
    from chemtools.ml.evaluation import evaluate_yield_predictor
    
    results = evaluate_yield_predictor(predictor, test_df)
    print(f"MAE: {results['mae']:.2f}%")
"""

from __future__ import annotations

from typing import Dict, Any, Tuple
import sys

import numpy as np

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas not installed. Install with: pip install pandas")
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
except ImportError:
    print("WARNING: matplotlib not installed. Plotting disabled.")
    plt = None


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """
    Compute regression metrics for yield prediction.
    
    Args:
        y_true: True yields (n_samples,)
        y_pred: Predicted yields (n_samples,)
    
    Returns:
        Dictionary of metrics:
            - mae: Mean absolute error
            - rmse: Root mean squared error
            - r2: R-squared score
            - within_10pct: Fraction within ±10 percentage points
            - within_15pct: Fraction within ±15 percentage points
            - mean_error: Mean signed error (bias)
    """
    errors = y_pred - y_true
    abs_errors = np.abs(errors)
    
    mae = np.mean(abs_errors)
    rmse = np.sqrt(np.mean(errors ** 2))
    
    # R-squared
    ss_res = np.sum(errors ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r2 = 1 - (ss_res / (ss_tot + 1e-10))
    
    # Within ±X% thresholds
    within_10 = np.mean(abs_errors <= 10.0)
    within_15 = np.mean(abs_errors <= 15.0)
    
    # Bias
    mean_error = np.mean(errors)
    
    return {
        'mae': mae,
        'rmse': rmse,
        'r2': r2,
        'within_10pct': within_10 * 100,  # Convert to percentage
        'within_15pct': within_15 * 100,
        'mean_error': mean_error,
    }


def plot_predictions(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    title: str = "Yield Predictions",
    output_path: str | None = None,
) -> None:
    """
    Plot predicted vs true yields with residuals.
    
    Args:
        y_true: True yields
        y_pred: Predicted yields
        title: Plot title
        output_path: Save path (if None, display plot)
    """
    if plt is None:
        print("WARNING: matplotlib not available, skipping plot")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Predicted vs True
    ax = axes[0]
    ax.scatter(y_true, y_pred, alpha=0.5, s=20, edgecolors='k', linewidths=0.5)
    ax.plot([0, 100], [0, 100], 'r--', linewidth=2, label='Perfect prediction')
    ax.plot([0, 100], [10, 110], 'orange', linestyle='--', alpha=0.5, label='±10%')
    ax.plot([0, 100], [-10, 90], 'orange', linestyle='--', alpha=0.5)
    ax.set_xlabel('True Yield (%)', fontsize=12)
    ax.set_ylabel('Predicted Yield (%)', fontsize=12)
    ax.set_title('Predicted vs True Yields', fontsize=14)
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    
    # Residuals
    ax = axes[1]
    residuals = y_pred - y_true
    ax.scatter(y_true, residuals, alpha=0.5, s=20, edgecolors='k', linewidths=0.5)
    ax.axhline(0, color='r', linestyle='--', linewidth=2, label='Zero error')
    ax.axhline(10, color='orange', linestyle='--', alpha=0.5, label='±10%')
    ax.axhline(-10, color='orange', linestyle='--', alpha=0.5)
    ax.set_xlabel('True Yield (%)', fontsize=12)
    ax.set_ylabel('Residual (Predicted - True) (%)', fontsize=12)
    ax.set_title('Residual Plot', fontsize=14)
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 100)
    
    fig.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {output_path}")
    else:
        plt.show()
    
    plt.close()


def evaluate_yield_predictor(
    predictor,
    test_df: pd.DataFrame,
    output_dir: str | None = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Comprehensive evaluation of yield predictor.
    
    Args:
        predictor: Trained yield predictor with predict() method
        test_df: Test dataframe with yield_pct column
        output_dir: Directory for saving plots (if None, no plots saved)
        verbose: Print results
    
    Returns:
        Dictionary with:
            - metrics: Computed metrics dict
            - predictions: Predicted yields
            - true_yields: True yields
    """
    if verbose:
        print("=" * 70)
        print("Yield Predictor Evaluation")
        print("=" * 70)
    
    # Predict
    y_true = test_df['yield_pct'].values
    y_pred = predictor.predict(test_df)
    
    # Compute metrics
    metrics = compute_metrics(y_true, y_pred)
    
    if verbose:
        print(f"Test set: {len(test_df)} reactions")
        print()
        print("Metrics:")
        print(f"  MAE:             {metrics['mae']:.2f}%")
        print(f"  RMSE:            {metrics['rmse']:.2f}%")
        print(f"  R²:              {metrics['r2']:.3f}")
        print(f"  Within ±10%:     {metrics['within_10pct']:.1f}%")
        print(f"  Within ±15%:     {metrics['within_15pct']:.1f}%")
        print(f"  Mean Error:      {metrics['mean_error']:.2f}% (bias)")
        print()
    
    # Plot
    if output_dir:
        from pathlib import Path
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        
        plot_path = out_path / 'evaluation_results.png'
        plot_predictions(
            y_true,
            y_pred,
            title=f"Yield Prediction (MAE={metrics['mae']:.1f}%, R²={metrics['r2']:.2f})",
            output_path=str(plot_path),
        )
    
    if verbose:
        print("=" * 70)
    
    return {
        'metrics': metrics,
        'predictions': y_pred,
        'true_yields': y_true,
    }


def compare_models(
    results_dict: Dict[str, Dict[str, Any]],
    output_path: str | None = None,
) -> None:
    """
    Compare multiple models side-by-side.
    
    Args:
        results_dict: Dictionary mapping model name -> evaluation results
        output_path: Save path for comparison plot
    """
    if plt is None:
        print("WARNING: matplotlib not available, skipping comparison plot")
        return
    
    print("=" * 70)
    print("Model Comparison")
    print("=" * 70)
    
    # Extract metrics
    model_names = list(results_dict.keys())
    maes = [results_dict[name]['metrics']['mae'] for name in model_names]
    rmses = [results_dict[name]['metrics']['rmse'] for name in model_names]
    r2s = [results_dict[name]['metrics']['r2'] for name in model_names]
    
    # Print table
    print(f"{'Model':<30s} {'MAE':>8s} {'RMSE':>8s} {'R²':>8s}")
    print("-" * 70)
    for name, mae, rmse, r2 in zip(model_names, maes, rmses, r2s):
        print(f"{name:<30s} {mae:>7.2f}% {rmse:>7.2f}% {r2:>8.3f}")
    print()
    
    # Bar plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # MAE
    axes[0].bar(model_names, maes, color='steelblue', edgecolor='k')
    axes[0].set_ylabel('MAE (%)', fontsize=12)
    axes[0].set_title('Mean Absolute Error', fontsize=14)
    axes[0].grid(axis='y', alpha=0.3)
    
    # RMSE
    axes[1].bar(model_names, rmses, color='coral', edgecolor='k')
    axes[1].set_ylabel('RMSE (%)', fontsize=12)
    axes[1].set_title('Root Mean Squared Error', fontsize=14)
    axes[1].grid(axis='y', alpha=0.3)
    
    # R²
    axes[2].bar(model_names, r2s, color='seagreen', edgecolor='k')
    axes[2].set_ylabel('R²', fontsize=12)
    axes[2].set_title('R² Score', fontsize=14)
    axes[2].grid(axis='y', alpha=0.3)
    axes[2].set_ylim(0, 1)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Comparison plot saved to: {output_path}")
    else:
        plt.show()
    
    plt.close()
    
    print("=" * 70)
