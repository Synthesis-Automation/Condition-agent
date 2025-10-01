from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

from chemtools.rules import api as crl_api


@lru_cache(maxsize=1)
def _load_default_crl() -> Dict[str, Any]:
    return crl_api.load_default_crl()


def _load_custom_crl(path: Optional[str | Path]) -> Dict[str, Any]:
    if path is None:
        return _load_default_crl()
    return crl_api.load_crl(path)


def recommend_rule_based(
    family: str,
    features: Dict[str, Any],
    *,
    job_ctx: Optional[Dict[str, Any]] = None,
    top_n: int = 3,
    crl_path: Optional[str | Path] = None,
) -> List[Dict[str, Any]]:
    """Return ranked rule-based condition recommendations for a reaction family.

    Args:
        family: Rule family token, e.g. ``"Buchwald_CN"``.
        features: Structured substrate/context features (nested dicts allowed).
        job_ctx: Optional job-level context (mode, base, solvent, outcome, etc.).
        top_n: Number of ranked playbooks to return.
        crl_path: Optional path to an alternate CRL JSON.

    Returns:
        List of recommendation dicts ordered by descending score.
    """
    crl = _load_custom_crl(crl_path)
    return crl_api.recommend(
        crl,
        family,
        features,
        job_ctx=job_ctx,
        top_n=top_n,
    )


def list_playbooks(
    family: str,
    *,
    crl_path: Optional[str | Path] = None,
) -> List[Dict[str, Any]]:
    """Return all playbooks for a family (unordered)."""
    crl = _load_custom_crl(crl_path)
    fam = (crl.get("families") or {}).get(family, {})
    return list(fam.get("playbooks") or [])


def get_family_defaults(
    family: str,
    *,
    crl_path: Optional[str | Path] = None,
) -> Dict[str, Any]:
    crl = _load_custom_crl(crl_path)
    shared_defaults = crl.get("defaults") or {}
    fam_defaults = (crl.get("families") or {}).get(family, {}).get("defaults") or {}
    merged = dict(shared_defaults)
    merged.update(fam_defaults)
    return merged
