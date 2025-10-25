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
  "mapped_family": str | None,   # Canonical taxonomy reaction type ID (e.g., "suzuki_miyaura")
  "confidence": float | None,    # If provided or inferred
    "raw": any,                    # Raw object from rxn_insight for debugging
    "catalysts": list[str]         # Metals inferred from agents (if any)
}
"""

from typing import Any, Dict, Optional, Set

try:
    # Optional classifier path available in newer rxn_insight builds
    from rxn_insight.classification import ReactionClassifier  # type: ignore
    _HAS_RC = True
except Exception:
    _HAS_RC = False
try:
    # Use our own helpers to parse reactants and apply simple functional-group hits
    from .smiles import normalize_reaction as _norm_rxn  # type: ignore
    from .router import (  # type: ignore
        _rule_hits as _hits,
        _detect_agent_metals as _agent_metals,
        resolve_reaction_family as _resolve_family,
    )
except Exception:
    _norm_rxn = None  # type: ignore
    _hits = None  # type: ignore
    _agent_metals = None  # type: ignore
    def _resolve_family(value: Optional[str]):  # type: ignore
        return value


def _resolve_family_safe(label: Optional[str]) -> Optional[str]:
    try:
        resolved = _resolve_family(label)
    except Exception:
        resolved = None
    if resolved:
        return resolved
    if label and label not in {"Unknown", "unknown"}:
        return label
    return None


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

    alias: Optional[str] = None

    # C-N couplings (Buchwald–Hartwig, Ullmann–Goldberg, Chan–Lam etc.)
    if (
        "heteroatom alkylation" in c
        or "heteroatom arylation" in c
        or ("n-arylation" in n)
        or ("buchwald" in n or "hartwig" in n or "ullmann" in n or "goldberg" in n or "chan-lam" in n)
    ):
        alias = "Ullmann_CN"

    # Suzuki C–C coupling
    elif ("c-c coupling" in c and "suzuki" in n) or ("suzuki" in n):
        alias = "Suzuki_CC"

    # Other C–C couplings could be added here as families become supported

    # Amide formation (acylation)
    elif ("acylation" in c or "amide" in n) and (
        "carboxylic acid" in n or "amine" in n or "amide" in n
    ):
        alias = "Amide_Coupling"

    if alias:
        return _resolve_family_safe(alias)
    return None


def _refine_cn_family(
    mapped: Optional[str],
    rxn_name: Optional[str],
    rxn_class: Optional[str],
    catalysts: Set[str],
) -> tuple[Optional[str], Optional[str]]:
    if not catalysts:
        return mapped, rxn_name

    cn_related = bool(
        (mapped in {"cn_coupling", "ullmann_cn", "buchwald_hartwig_c_n"})
        or (rxn_class and "heteroatom" in rxn_class.lower())
        or (rxn_name and any(token in rxn_name.lower() for token in ("buchwald", "hartwig", "ullmann", "amination")))
    )
    if not cn_related:
        return mapped, rxn_name

    if "Pd" in catalysts:
        override = _resolve_family_safe("Buchwald_CN")
        if override:
            mapped = override
            if not rxn_name or "buchwald" not in rxn_name.lower():
                rxn_name = rxn_name or "Buchwald-Hartwig C–N coupling"
    elif "Cu" in catalysts and (mapped is None or mapped in {"Unknown", "cn_coupling"}):
        override = _resolve_family_safe("Ullmann_CN")
        if override:
            mapped = override
            if not rxn_name or "ullmann" not in rxn_name.lower():
                rxn_name = rxn_name or "Ullmann/Golberg C–N coupling"

    return mapped, rxn_name


def detect_reaction_type(reaction_smiles: str) -> Dict[str, Any]:
    """Detect reaction type using rxn_insight, returning canonical taxonomy families."""
    if not is_available():
        return {
            "available": False,
            "success": False,
            "rxn_class": None,
            "rxn_name": None,
            "mapped_family": None,
            "confidence": None,
            "raw": None,
        }

    norm = _norm_rxn(reaction_smiles) if _norm_rxn else None
    catalysts: Set[str] = set()
    if norm is not None and _agent_metals is not None:
        try:
            catalysts = _agent_metals(norm.get("agents") or [])  # type: ignore[arg-type]
        except Exception:
            catalysts = set()

    raw = _call_insight(reaction_smiles)
    rxn_class, rxn_name, conf = _extract_fields(raw)
    mapped_family = _map_to_family(rxn_class, rxn_name)

    # If no mapping yet, try ReactionClassifier + lightweight heuristics
    if not mapped_family and _HAS_RC:
        try:
            rc = ReactionClassifier(reaction_smiles)  # type: ignore[call-arg]
            try:
                rc.classify_reaction()
            except Exception:
                pass

            if getattr(rc, "is_cc_coupling")() if hasattr(rc, "is_cc_coupling") else False:  # type: ignore[misc]
                rxn_class = rxn_class or "C-C Coupling"
                boron = False
                try:
                    if norm is None and _norm_rxn:
                        norm = _norm_rxn(reaction_smiles)
                    if norm and _hits:
                        reactants = [
                            (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
                            for r in (norm.get("reactants") or [])
                        ]
                        reactants = [s for s in reactants if s]
                        boron = bool(_hits(reactants).get("boron"))
                except Exception:
                    boron = False
                if boron:
                    rxn_name = rxn_name or "Suzuki coupling with boronic acids"
                    mapped_family = _resolve_family_safe("Suzuki_CC")
            elif getattr(rc, "is_heteroatom_alkylation")() if hasattr(rc, "is_heteroatom_alkylation") else False:  # type: ignore[misc]
                rxn_class = rxn_class or "Heteroatom Alkylation and Arylation"
                mapped_family = mapped_family or _resolve_family_safe("Ullmann_CN")
            elif getattr(rc, "is_acylation")() if hasattr(rc, "is_acylation") else False:  # type: ignore[misc]
                rxn_class = rxn_class or "Acylation"
                mapped_family = mapped_family or _resolve_family_safe("Amide_Coupling")
        except Exception:
            pass

    if catalysts:
        mapped_family, rxn_name = _refine_cn_family(
            mapped_family,
            rxn_name,
            rxn_class,
            catalysts,
        )

    success = bool(rxn_class or rxn_name or mapped_family)
    return {
        "available": True,
        "success": success,
        "rxn_class": rxn_class,
        "rxn_name": rxn_name,
        "mapped_family": mapped_family,
        "confidence": conf,
        "raw": raw,
        "catalysts": sorted(catalysts) if catalysts else [],
    }
