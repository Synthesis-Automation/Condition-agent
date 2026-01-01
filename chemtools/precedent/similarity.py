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

    # Weighted categorical matching (legacy features)
    weights = {
        "LG": 0.35,
        "nuc_class": 0.35,
        "ortho_count": 0.10,
        "para_EWG": 0.10,
        "heteroaryl": 0.10,
    }
    score = 0.0
    total = 0.0
    for k, w in weights.items():
        av = a.get(k)
        bv = b.get(k)
        if av is None or bv is None:
            continue
        total += w
        # Normalize bools to exact equality
        if isinstance(av, bool) or isinstance(bv, bool):
            if bool(av) == bool(bv):
                score += w
        else:
            if av is not None and bv is not None and str(av).lower() == str(bv).lower():
                score += w

    # Taxonomy-aware categorical matching
    taxonomy_keys = {
        "reaction_type": 0.20,
        "reaction_category": 0.10,
    }
    for key, weight in taxonomy_keys.items():
        av = a.get(key)
        bv = b.get(key)
        if av is None or bv is None:
            continue
        total += weight
        if str(av).lower() == str(bv).lower():
            score += weight

    # Optional small numeric distances if present in feature dicts
    # Use an exponential decay mapped to <= 0.25 extra credit total
    numeric_keys: List[Tuple[str, float, float]] = [
        ("T_C", 50.0, 0.10),  # (scale, weight)
        ("time_h", 8.0, 0.05),
        ("max_aryl_steric", 6.0, 0.05),
        ("avg_aryl_electronic", 5.0, 0.05),
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

    # Boolean feature overlap (new taxonomy payloads)
    bool_keys = [
        key for key in a.keys() & b.keys()
        if isinstance(a.get(key), bool) and isinstance(b.get(key), bool)
    ]
    if bool_keys:
        matches = sum(1 for key in bool_keys if bool(a.get(key)) == bool(b.get(key)))
        weight = 0.25
        score += weight * (matches / max(1, len(bool_keys)))
        total += weight

    if total <= 0:
        return 0.0
    return max(0.0, min(1.0, score / total))
