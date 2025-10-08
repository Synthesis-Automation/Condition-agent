"""Similarity calculations for precedent ranking.

Provides functions to compute similarity between query and precedent reactions
based on categorical features and numeric parameters.
"""
from typing import Dict, Any, List, Tuple


def _similarity(a: Dict[str, Any], b: Dict[str, Any]) -> float:
    """Calculate similarity between two feature dictionaries.
    
    Uses weighted categorical matching for discrete features and exponential
    decay for continuous features. Exact bin match returns 1.0.
    
    Args:
        a: First feature dict (typically query)
        b: Second feature dict (typically precedent)
        
    Returns:
        Similarity score between 0.0 and 1.0
    """
    # Exact bin match gets perfect similarity
    if (a.get("bin") or "") == (b.get("bin") or "") and a.get("bin"):
        return 1.0

    # Weighted categorical matching
    weights = {
        "LG": 0.35,
        "nuc_class": 0.35,
        "ortho_count": 0.10,
        "para_EWG": 0.10,
        "heteroaryl": 0.10,
    }
    score = 0.0
    total = sum(weights.values())
    for k, w in weights.items():
        av = a.get(k)
        bv = b.get(k)
        # Normalize bools to exact equality
        if isinstance(av, bool) or isinstance(bv, bool):
            if bool(av) == bool(bv):
                score += w
        else:
            if av is not None and bv is not None and str(av).lower() == str(bv).lower():
                score += w

    # Optional small numeric distances if present in feature dicts
    # Use an exponential decay mapped to <= 0.15 extra credit total
    numeric_keys: List[Tuple[str, float, float]] = [
        ("T_C", 50.0, 0.10),  # (scale, weight)
        ("time_h", 8.0, 0.05),
    ]
    for key, scale, w in numeric_keys:
        if key in a and key in b:
            try:
                da = float(a[key])
                db = float(b[key])
                import math
                sim_num = math.exp(-abs(da - db) / max(1e-9, scale))
                score += w * sim_num
                total += w
            except Exception:
                # ignore numeric similarity if non-numeric
                pass

    if total <= 0:
        return 0.0
    return max(0.0, min(1.0, score / total))
