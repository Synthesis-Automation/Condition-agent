"""
Unified feature bundles for molecules and reactions.

This module forwards to ``chemtools.featurizers.formatters`` while keeping a
small compatibility shim for legacy callers that expect ``result["reaction"]``
to contain the reaction payload.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from chemtools.reaction import featurize_reaction as _featurize_reaction
from .reaction_assist import apply_llm_reaction_assist


def featurize_reaction(
    reaction_smiles: str,
    *args: Any,
    **kwargs: Any,
) -> Dict[str, Any]:
    payload: Dict[str, Any] = _featurize_reaction(reaction_smiles, *args, **kwargs)
    options = kwargs.get("options")
    if options is None and len(args) >= 2 and isinstance(args[1], dict):
        options = args[1]
    elif options is None and len(args) >= 1 and isinstance(args[0], dict):
        options = args[0]
    if isinstance(payload, dict):
        payload = apply_llm_reaction_assist(
            payload,
            reaction_smiles=reaction_smiles,
            options=options if isinstance(options, dict) else None,
        )
    if isinstance(payload, dict) and "reaction" not in payload:
        payload = dict(payload)
        payload["reaction"] = dict(payload)
    return payload


__all__ = ["featurize_reaction"]
