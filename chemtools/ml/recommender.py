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
    
    # Get baseline DRFP-based recommendations
    knn_results = recommend_from_reaction(reaction_smiles, k=k, **kwargs)
    
    # Extract from formatted output structure
    formatted = knn_results.get('formatted', {})
    precedent_pack = knn_results.get('precedent_pack', {})
    recommendation = knn_results.get('recommendation', {})
    
    n_precedents = precedent_pack.get('support', 0)
    recommended = formatted.get('recommended_conditions', [])
    precedents = precedent_pack.get('precedents', [])
    
    # Get default temperature and time from recommendation
    default_T = recommendation.get('T_C', 80.0)
    default_time = recommendation.get('time_h', 12.0)
    
    # Decision: Use ML or stick with k-NN?
    use_ml = (
        model_path is not None
        and DRFPYieldPredictor is not None
        and n_precedents >= ml_threshold
        and len(recommended) > 0
    )
    
    if not use_ml:
        # Return k-NN results as-is
        knn_results['method'] = 'knn'
        knn_results['reason'] = (
            f'Using k-NN (n_precedents={n_precedents} < threshold={ml_threshold})'
            if n_precedents < ml_threshold
            else 'ML model not available'
        )
        return knn_results
    
    # Load ML model
    try:
        predictor = DRFPYieldPredictor.load(model_path)
    except Exception as e:
        warnings.warn(f"Failed to load ML model: {e}. Falling back to k-NN.")
        knn_results['method'] = 'knn'
        knn_results['reason'] = f'Model load failed: {e}'
        return knn_results
    
    # Predict yields for each recommended variant
    predicted_yields = []
    
    for variant in recommended:
        # Extract from formatted structure
        summary = variant.get('summary', {})
        combo = variant.get('combo', {})
        
        # Get base and solvent from summary
        base_info = summary.get('base', {})
        solvent_info = summary.get('solvent', {})
        
        # Build test dataframe row
        test_row = {
            'reaction_smiles': reaction_smiles,
            'core': summary.get('core', 'Unknown'),
            'base_uid': base_info.get('cas') or base_info.get('name', 'Unknown') if base_info else 'Unknown',
            'solvent_uid': solvent_info.get('cas') or solvent_info.get('name', 'Unknown') if solvent_info else 'Unknown',
            'T_C': default_T,
            'time_h': default_time,
        }
        
        try:
            import pandas as pd
            test_df = pd.DataFrame([test_row])
            yield_pred = predictor.predict(test_df)[0]
            predicted_yields.append(float(yield_pred))
        except Exception as e:
            warnings.warn(f"Yield prediction failed for variant: {e}")
            predicted_yields.append(None)
    
    # Re-rank by predicted yield (descending)
    valid_indices = [i for i, y in enumerate(predicted_yields) if y is not None]
    
    if not valid_indices:
        # Prediction failed for all variants, fall back to k-NN
        knn_results['method'] = 'knn'
        knn_results['reason'] = 'Yield prediction failed for all variants'
        return knn_results
    
    # Sort by predicted yield
    sorted_indices = sorted(
        valid_indices,
        key=lambda i: predicted_yields[i],
        reverse=True,
    )
    
    # Reorder recommendations
    reranked_conditions = [recommended[i] for i in sorted_indices]
    reranked_yields = [predicted_yields[i] for i in sorted_indices]
    
    # Attach predicted yields to each variant
    for variant, yield_pred in zip(reranked_conditions, reranked_yields):
        variant['predicted_yield_pct'] = yield_pred
    
    # Update the formatted structure with reranked conditions
    formatted_updated = dict(formatted)
    formatted_updated['recommended_conditions'] = reranked_conditions
    formatted_updated['meta'] = formatted_updated.get('meta', {})
    formatted_updated['meta']['ml_reranked'] = True
    formatted_updated['meta']['predicted_yields'] = reranked_yields
    
    # Return full knn_results structure with updated formatted section
    ml_results = dict(knn_results)
    ml_results['formatted'] = formatted_updated
    ml_results['method'] = 'ml'
    ml_results['reason'] = f'ML re-ranking (n_precedents={n_precedents})'
    ml_results['predicted_yields'] = reranked_yields
    
    return ml_results


def recommend_with_yield_prediction(
    reaction_smiles: str,
    model: Optional[Any] = None,  # DRFPYieldPredictor instance
    model_path: Optional[str] = None,
    top_n: int = 5,
    **knn_kwargs,
) -> List[Dict[str, Any]]:
    """
    Simplified API for getting top-N conditions with yield predictions.
    
    Args:
        reaction_smiles: Reaction SMILES
        model: Pre-loaded DRFPYieldPredictor instance (optional)
        model_path: Path to model file (used if model not provided)
        top_n: Number of top conditions to return (default: 5)
        **knn_kwargs: Additional k-NN arguments
    
    Returns:
        List of top-N condition dictionaries with predicted_yield_pct field
    """
    # Load model if not provided
    if model is None and model_path is not None:
        if DRFPYieldPredictor is None:
            raise ImportError("DRFPYieldPredictor not available")
        model = DRFPYieldPredictor.load(model_path)
    
    # Get hybrid recommendations
    results = hybrid_recommend(
        reaction_smiles,
        model_path=model_path,
        **knn_kwargs,
    )
    
    # Extract from formatted structure
    formatted = results.get('formatted', {})
    conditions = formatted.get('recommended_conditions', [])
    
    # Return top-N
    return conditions[:top_n]
