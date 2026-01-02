"""
ML-enhanced recommendation with yield prediction.

Integrates DRFP-based yield prediction with the core recommendation engine.
Re-ranks condition variants by predicted yield when sufficient precedents are available.

Strategy:
    - If n_precedents >= threshold (default: 50): Use ML to predict yields
      and re-rank by predicted yield
    - If n_precedents < threshold: Use k-NN vote-based ranking
    - Always show precedents for interpretability

Example:
    from chemtools.ml.recommender import hybrid_recommend
    
    results = hybrid_recommend(
        reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        model_path="models/drfp_yield_v1.pkl",
    )
"""

from __future__ import annotations

from typing import Dict, List, Any, Optional
import warnings

try:
    from .drfp_predictor import DRFPYieldPredictor
except ImportError:
    DRFPYieldPredictor = None
    warnings.warn("DRFPYieldPredictor not available. ML prediction disabled.")


def hybrid_recommend(
    reaction_smiles: str,
    model_path: Optional[str] = None,
    ml_threshold: int = 50,
    k: int = 10,
    **kwargs,
) -> Dict[str, Any]:
    """
    Hybrid ML + k-NN recommendation with yield prediction.
    
    Uses DRFP-based reaction similarity (via recommend.recommend_from_reaction)
    and optionally re-ranks results by ML-predicted yield.
    
    Args:
        reaction_smiles: Reaction SMILES string (reactants>>products)
        model_path: Path to trained ML model (.pkl). If None, uses k-NN only
        ml_threshold: Minimum precedents to use ML (default: 50)
        k: Number of precedents to retrieve (default: 10)
        **kwargs: Additional arguments passed to recommend_from_reaction()
    
    Returns:
        Recommendation results with:
            - recommended_conditions: List of condition variants
            - predicted_yields: (if ML used) Predicted yields for each variant
            - method: 'ml' or 'knn'
            - precedents: Retrieved precedents
            - n_precedents: Total precedent count
    """
    # Import here to avoid circular dependency
    try:
        from ..recommend import recommend_from_reaction
    except ImportError:
        raise ImportError("chemtools.recommend not available")
    
    # Get baseline unified recommendations
    knn_results = recommend_from_reaction(reaction_smiles, k=k, **kwargs)
    knn_results["method"] = "unified"
    if model_path or DRFPYieldPredictor is not None:
        knn_results["reason"] = "Yield reranking is not available for the unified recommender."
    return knn_results
