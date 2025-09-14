from __future__ import annotations

"""
rxn-insight integration for reaction type detection.

This module attempts to use the optional `rxn_insight` package to classify
reaction SMILES into a coarse reaction family used by chemtools, with
graceful fallback when the package is unavailable or its API changes.

Public API:
 - is_available() -> bool
 - detect_reaction_type(reaction_smiles: str) -> dict

Returned dict (best-effort):
{
  "available": bool,             # Whether rxn_insight is importable
  "success": bool,               # Whether a classification was obtained
  "rxn_class": str | None,       # Broad class (e.g., "C-C Coupling")
  "rxn_name": str | None,        # Specific name (e.g., "Suzuki coupling with boronic acids")
  "mapped_family": str | None,   # Mapped to chemtools family (e.g., "Suzuki_CC")
  "confidence": float | None,    # If provided or inferred
  "raw": any                     # Raw object from rxn_insight for debugging
}
"""

from typing import Any, Dict, Optional

try:
    # Optional classifier path available in newer rxn_insight builds
    from rxn_insight.classification import ReactionClassifier  # type: ignore
    _HAS_RC = True
except Exception:
    _HAS_RC = False
try:
    # Use our own helpers to parse reactants and apply simple functional-group hits
    from .smiles import normalize_reaction as _norm_rxn  # type: ignore
    from .router import _rule_hits as _hits  # type: ignore
except Exception:
    _norm_rxn = None  # type: ignore
    _hits = None  # type: ignore


def is_available() -> bool:
    try:
        import rxn_insight  # type: ignore
        _ = rxn_insight
        return True
    except Exception:
        return False


def _call_insight(reaction_smiles: str) -> Any:
    """Attempt to call common entrypoints in rxn_insight.

    Since rxn_insight API may evolve, try a few likely function names.
    Returns whatever the library returns, or None on failure.
    """
    try:
        import rxn_insight  # type: ignore
    except Exception:
        return None

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

    # Some packages expose a class-based API
    # e.g., rxn_insight.Reaction(...).analyze()
    try:
        cls = getattr(rxn_insight, "Reaction", None)
        if cls is not None:
            obj = cls(reaction_smiles)  # type: ignore[call-arg]
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
    """Best-effort extraction of (rxn_class, rxn_name, confidence) from result."""
    rxn_class: Optional[str] = None
    rxn_name: Optional[str] = None
    conf: Optional[float] = None

    try:
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

        # Some libs return objects with attributes
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


def _map_to_family(rxn_class: Optional[str], rxn_name: Optional[str]) -> Optional[str]:
    """Map rxn-insight's (class, name) to chemtools family labels."""
    c = (rxn_class or "").lower()
    n = (rxn_name or "").lower()

    # C-N couplings (Buchwald–Hartwig, Ullmann–Goldberg, Chan–Lam etc.)
    if (
        "heteroatom alkylation" in c
        or "heteroatom arylation" in c
        or ("n-arylation" in n)
        or ("buchwald" in n or "hartwig" in n or "ullmann" in n or "goldberg" in n or "chan-lam" in n)
    ):
        return "Ullmann_CN"

    # Suzuki C–C coupling
    if ("c-c coupling" in c and "suzuki" in n) or ("suzuki" in n):
        return "Suzuki_CC"

    # Other C–C couplings could be added here as families become supported

    # Amide formation (acylation)
    if ("acylation" in c or "amide" in n) and (
        "carboxylic acid" in n or "amine" in n or "amide" in n
    ):
        return "Amide_Coupling"

    return None


def detect_reaction_type(reaction_smiles: str) -> Dict[str, Any]:
    """Detect reaction type using rxn_insight, mapping to chemtools family.

    Always returns a dict; falls back to {available: False, success: False}
    when rxn_insight is not installed or classification fails.
    """
    avail = is_available()
    if not avail:
        return {"available": False, "success": False, "rxn_class": None, "rxn_name": None, "mapped_family": None, "confidence": None, "raw": None}

    raw = _call_insight(reaction_smiles)
    rxn_class, rxn_name, conf = _extract_fields(raw)
    mapped = _map_to_family(rxn_class, rxn_name)

    # If no mapping yet, try ReactionClassifier + lightweight heuristics
    if not mapped and _HAS_RC:
        try:
            rc = ReactionClassifier(reaction_smiles)  # type: ignore[call-arg]
            # Populate internal flags
            try:
                rc.classify_reaction()
            except Exception:
                pass
            # Derive broad class
            if getattr(rc, 'is_cc_coupling')() if hasattr(rc, 'is_cc_coupling') else False:  # type: ignore[misc]
                rxn_class = rxn_class or 'C-C Coupling'
                # Use boron presence to identify Suzuki when possible
                boron = False
                try:
                    if _norm_rxn and _hits:
                        norm = _norm_rxn(reaction_smiles)
                        reactants = [
                            (r.get('smiles_norm') or r.get('largest_smiles') or r.get('input') or '')
                            for r in (norm.get('reactants') or [])
                        ]
                        reactants = [s for s in reactants if s]
                        h = _hits(reactants)
                        boron = bool(h.get('boron'))
                except Exception:
                    boron = False
                if boron:
                    rxn_name = rxn_name or 'Suzuki coupling with boronic acids'
                    mapped = mapped or 'Suzuki_CC'
            elif getattr(rc, 'is_heteroatom_alkylation')() if hasattr(rc, 'is_heteroatom_alkylation') else False:  # type: ignore[misc]
                rxn_class = rxn_class or 'Heteroatom Alkylation and Arylation'
                mapped = mapped or 'Ullmann_CN'
            elif getattr(rc, 'is_acylation')() if hasattr(rc, 'is_acylation') else False:  # type: ignore[misc]
                rxn_class = rxn_class or 'Acylation'
                mapped = mapped or 'Amide_Coupling'
        except Exception:
            pass

    success = bool(rxn_class or rxn_name or mapped)
    return {
        "available": True,
        "success": success,
        "rxn_class": rxn_class,
        "rxn_name": rxn_name,
        "mapped_family": mapped,
        "confidence": conf,
        "raw": raw,
    }
