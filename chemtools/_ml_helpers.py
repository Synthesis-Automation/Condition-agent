"""
Internal helpers for ML-based reaction detection using rxn-insight.

This module provides the low-level rxn-insight integration that was previously
in reaction_type_detector.py. It's marked internal (_ml_helpers) to discourage
direct use - all detection should go through detect_reaction().
"""

from typing import Any, Dict, Optional, Set
import logging

logger = logging.getLogger(__name__)


def is_available() -> bool:
    """Check if rxn-insight package is available."""
    try:
        import rxn_insight  # type: ignore
        _ = rxn_insight
        return True
    except Exception:
        return False


def _call_insight(reaction_smiles: str) -> Any:
    """
    Attempt to call rxn_insight with the reaction SMILES.
    
    Since rxn_insight API may evolve, try several likely function names.
    
    Returns:
        Raw rxn-insight result object, or None on failure
    """
    try:
        import rxn_insight  # type: ignore
    except Exception:
        return None

    # Try common entrypoints
    candidates = [
        "analyze", "analyze_reaction", "classify", "classify_reaction",
        "analyze_smiles", "classify_smiles",
    ]
    for name in candidates:
        fn = getattr(rxn_insight, name, None)
        if callable(fn):
            try:
                return fn(reaction_smiles)
            except Exception:
                continue

    # Try class-based API (e.g., rxn_insight.Reaction(...).analyze())
    try:
        cls = getattr(rxn_insight, "Reaction", None)
        if cls is not None:
            obj = cls(reaction_smiles)  # type: ignore
            for m in ("analyze", "classify", "analyze_reaction"):
                meth = getattr(obj, m, None)
                if callable(meth):
                    try:
                        return meth()
                    except Exception:
                        pass
    except Exception:
        pass

    return None


def _extract_fields(res: Any) -> tuple[Optional[str], Optional[str], Optional[float]]:
    """
    Extract (rxn_class, rxn_name, confidence) from rxn-insight result.
    
    Handles both dict and object responses.
    
    Returns:
        (rxn_class, rxn_name, confidence) tuple
    """
    rxn_class: Optional[str] = None
    rxn_name: Optional[str] = None
    conf: Optional[float] = None

    try:
        # Dict-style response
        if isinstance(res, dict):
            rxn_class = (
                res.get("class")
                or res.get("rxn_class")
                or res.get("classification")
                or res.get("class_label")
            )
            rxn_name = (
                res.get("name")
                or res.get("rxn_name")
                or res.get("specific")
                or res.get("label")
                or res.get("class_name")
            )
            c = res.get("confidence") or res.get("score") or res.get("prob")
            try:
                conf = float(c) if c is not None else None
            except Exception:
                conf = None
            return rxn_class, rxn_name, conf

        # Object-style response
        if hasattr(res, "classification") or hasattr(res, "class_label"):
            rxn_class = getattr(res, "classification", None) or getattr(res, "class_label", None)
        if hasattr(res, "name") or hasattr(res, "label"):
            rxn_name = getattr(res, "name", None) or getattr(res, "label", None)
        if hasattr(res, "confidence") or hasattr(res, "score"):
            try:
                conf = float(getattr(res, "confidence", None) or getattr(res, "score", None))
            except Exception:
                conf = None
    except Exception:
        pass

    return rxn_class, rxn_name, conf


def call_rxn_insight(reaction_smiles: str) -> Dict[str, Any]:
    """
    Call rxn-insight and return standardized result.
    
    This is the main entry point for ML detection, consolidating the logic
    from the old _detect_reaction_type_impl().
    
    Returns:
        {
            "available": bool,       # Is rxn-insight installed?
            "success": bool,         # Did we get a result?
            "rxn_class": str | None,
            "rxn_name": str | None,
            "confidence": float | None,
            "raw": any,              # Raw response for debugging
        }
    """
    if not is_available():
        return {
            "available": False,
            "success": False,
            "rxn_class": None,
            "rxn_name": None,
            "confidence": None,
            "raw": None,
        }

    # Call rxn-insight
    raw = _call_insight(reaction_smiles)
    if raw is None:
        return {
            "available": True,
            "success": False,
            "rxn_class": None,
            "rxn_name": None,
            "confidence": None,
            "raw": None,
        }

    # Extract fields
    rxn_class, rxn_name, conf = _extract_fields(raw)
    
    # Success if we got at least a class or name
    success = bool(rxn_class or rxn_name)

    return {
        "available": True,
        "success": success,
        "rxn_class": rxn_class,
        "rxn_name": rxn_name,
        "confidence": conf,
        "raw": raw,
    }
