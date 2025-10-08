"""
Machine learning models for reaction condition prediction.

This module provides ML-based yield predictors and condition recommenders
using learned representations (DRFP, neural networks, GNNs) instead of
hand-crafted features.
"""

from __future__ import annotations

__all__ = []

# Phase 1: DRFP baseline (available when dependencies installed)
try:
    from .drfp_predictor import DRFPYieldPredictor
    __all__.append('DRFPYieldPredictor')
except ImportError:
    pass

# Phase 2: Neural network (available when PyTorch installed)
try:
    from .neural_yield_predictor import NeuralYieldPredictor
    __all__.append('NeuralYieldPredictor')
except ImportError:
    pass

# Always available
try:
    from .evaluation import evaluate_yield_predictor
    __all__.append('evaluate_yield_predictor')
except ImportError:
    pass
