"""
Unified feature bundles for molecules and reactions.

This module forwards to ``chemtools.featurizers.formatters`` while keeping a
small compatibility shim for legacy callers that expect ``result["reaction"]``
to contain the reaction payload.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from .formatters import featurize_molecule as _featurize_molecule
from .formatters import featurize_reaction as _featurize_reaction


def featurize_molecule(smiles: str, *args: Any, **kwargs: Any) -> Dict[str, Any]:
    return _featurize_molecule(smiles, *args, **kwargs)


def featurize_reaction(
    reaction_smiles: str,
    *args: Any,
    **kwargs: Any,
) -> Dict[str, Any]:
    payload: Dict[str, Any] = _featurize_reaction(reaction_smiles, *args, **kwargs)
    if isinstance(payload, dict) and "reaction" not in payload:
        payload = dict(payload)
        payload["reaction"] = dict(payload)
    return payload


__all__ = ["featurize_molecule", "featurize_reaction"]
