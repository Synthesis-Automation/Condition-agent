"""Feature mapping helpers bridging substrate analysis to MCP server schema."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional


@dataclass(frozen=True)
class ClassifiedReactants:
    electrophiles: List[str]
    nucleophiles: List[str]


AROMATIC_MARKERS = ("c1", "c2", "c3")


def _classify_electrophile(smiles: str) -> Dict[str, str]:
    smi = (smiles or "").strip()
    lowered = smi.lower()
    if "otf" in lowered or "osf" in lowered:
        return {"class": "aryl sulfonate"}
    if "cl" in lowered:
        if any(marker in lowered for marker in AROMATIC_MARKERS):
            return {"class": "aryl chloride"}
        return {"class": "alkyl chloride"}
    if "br" in lowered:
        if any(marker in lowered for marker in AROMATIC_MARKERS):
            return {"class": "aryl bromide"}
        return {"class": "alkyl bromide"}
    if "i" in lowered:
        if any(marker in lowered for marker in AROMATIC_MARKERS):
            return {"class": "aryl iodide"}
        return {"class": "alkyl iodide"}
    if "b(oh)2" in lowered or "bpin" in lowered:
        return {"class": "boron reagent"}
    return {"class": "unknown"}


def _classify_nucleophile(smiles: str) -> Dict[str, str]:
    smi = (smiles or "").strip()
    lowered = smi.lower()
    if "n" not in lowered:
        return {"class": "unknown"}
    if "s(=o)" in lowered or "so2" in lowered:
        return {"class": "sulfonamide"}
    if "c(=o)n" in lowered:
        return {"class": "amide"}
    if "nc1" in lowered:
        return {"class": "primary aniline"}
    if lowered.count("n") > 1 or "n(" in lowered:
        return {"class": "secondary amine"}
    return {"class": "primary amine"}


def map_electrophile(reactants: Iterable[str]) -> Dict[str, str]:
    for smi in reactants:
        info = _classify_electrophile(smi)
        if info.get("class") != "unknown":
            return info
    return {"class": "unknown"}


def map_nucleophile(reactants: Iterable[str]) -> Dict[str, str]:
    for smi in reactants:
        info = _classify_nucleophile(smi)
        if info.get("class") != "unknown":
            return info
    return {"class": "unknown"}


def build_features(
    reactants: Iterable[str],
    context: Optional[Dict[str, Any]] = None,
    overrides: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Compose the feature payload for the MCP rules server.

    Parameters
    ----------
    reactants:
        Iterable of reactant SMILES.
    context:
        Optional experiment context (mode, base, solvent, temperature, etc.).
    overrides:
        Optional explicit overrides for the payload.
    """

    reactant_list = [r for r in reactants if r]
    features: Dict[str, Any] = {
        "electrophile": map_electrophile(reactant_list),
        "nucleophile": map_nucleophile(reactant_list),
    }
    if context:
        features.update({k: v for k, v in context.items() if v is not None})
    if overrides:
        features.update(overrides)
    return features


__all__ = ["build_features", "map_electrophile", "map_nucleophile"]
