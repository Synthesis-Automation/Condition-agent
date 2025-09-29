"""Substrate featurisation wrapper for Condition MCP."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Tuple

from pydantic import BaseModel, Field

from chemtools.featurizers import molecular as feat_molecular

from .base import SchemaStamped, pick_first, validate_payload


def _select_substrates(reactants: Iterable[str]) -> Tuple[str, str]:
    """Split reactants into electrophile and nucleophile heuristically."""

    clean = [r.strip() for r in reactants if isinstance(r, str) and r.strip()]
    if not clean:
        return "", ""
    if len(clean) == 1:
        return clean[0], ""

    def is_electrophile(token: str) -> bool:
        t = token.lower()
        return any(sig in t for sig in ("br", "cl", " i", "os(=o)(=o)c(f)(f)f", "otf"))

    first, second, *rest = clean + [""]
    if is_electrophile(first):
        return first, second
    if is_electrophile(second):
        return second, first
    return first, second


class FeaturizeInput(BaseModel):
    """Input payload for the ``featurize_substrates`` tool."""

    reactants: List[str] = Field(..., description="Canonical reactant SMILES strings")


class FeaturizeOutput(SchemaStamped):
    """Descriptors derived from the electrophile/nucleophile pair."""

    electrophile: str
    nucleophile: str
    descriptors: Dict[str, Any]


def featurize_substrates(data: Dict[str, Any]) -> Dict[str, Any]:
    """Run the Ullmann featuriser on the supplied reactants."""

    payload = validate_payload(FeaturizeInput, data)
    electrophile, nucleophile = _select_substrates(payload.reactants)
    descriptors = feat_molecular.featurize(electrophile, nucleophile)
    # Drop role-aware block to keep output deterministic
    if isinstance(descriptors, dict) and "role_aware" in descriptors:
        descriptors = {k: v for k, v in descriptors.items() if k != "role_aware"}

    output = FeaturizeOutput(
        electrophile=electrophile,
        nucleophile=nucleophile,
        descriptors=descriptors,
    )
    return output.model_dump()
