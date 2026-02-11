"""
Dataset-ingestion routing for taxonomy benchmarking.

This module classifies reaction rows before reaction-family benchmarking so
low-information or out-of-scope records can be excluded consistently.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Set, Tuple

from . import rdkit_helpers


_TRANSITION_OR_HEAVY_METALS: Set[str] = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Ga", "In", "Sn", "Sb", "Bi",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
}
_METAL_TOKEN_RE = re.compile(r"\[([A-Z][a-z]?)(?:[^\]]*)\]")


def _split_reaction_sides(reaction_smiles: str) -> Tuple[List[str], List[str]]:
    text = str(reaction_smiles or "").strip()
    if not text:
        return [], []
    if ">>" in text:
        left, right = text.split(">>", 1)
        return [x for x in left.split(".") if x], [x for x in right.split(".") if x]
    parts = text.split(">")
    if len(parts) == 3:
        return [x for x in parts[0].split(".") if x], [x for x in parts[2].split(".") if x]
    return [x for x in text.split(".") if x], []


def _metal_tokens(text: str) -> List[str]:
    hits = sorted({m for m in _METAL_TOKEN_RE.findall(str(text or "")) if m})
    return hits


def _count_unparseable_components(components: List[str]) -> int:
    if not rdkit_helpers.rdkit_available():
        return 0
    failures = 0
    for smi in components:
        try:
            mol = rdkit_helpers.parse_smiles(str(smi))
        except Exception:
            mol = None
        if mol is None:
            failures += 1
    return failures


def classify_reaction_for_taxonomy_benchmark(
    reaction_smiles: str,
) -> Dict[str, Any]:
    """
    Classify a row for taxonomy benchmarking inclusion.

    Route values:
      - eligible_taxonomy_benchmark
      - exclude_missing_side_or_malformed_reaction_smiles
      - exclude_organometallic_or_coordination_complex
      - exclude_unparseable_components
    """
    text = str(reaction_smiles or "").strip()
    reactants, products = _split_reaction_sides(text)

    route = "eligible_taxonomy_benchmark"
    reason = "eligible"
    excluded = False

    if not text or not reactants or not products:
        route = "exclude_missing_side_or_malformed_reaction_smiles"
        reason = "missing_side_or_malformed_reaction_smiles"
        excluded = True

    metal_tokens = _metal_tokens(text)
    has_coordination_arrow = "->" in text or "<-" in text
    has_heavy_metal = bool(set(metal_tokens) & _TRANSITION_OR_HEAVY_METALS)
    if not excluded and (has_coordination_arrow and has_heavy_metal):
        route = "exclude_organometallic_or_coordination_complex"
        reason = "coordination_arrow_and_metal_token"
        excluded = True

    components = [str(x) for x in (reactants + products)]
    unparseable_components = _count_unparseable_components(components)
    if not excluded and unparseable_components > 0:
        route = "exclude_unparseable_components"
        reason = "rdkit_component_parse_failures"
        excluded = True

    return {
        "route": route,
        "excluded": excluded,
        "reason": reason,
        "reactant_count": len(reactants),
        "product_count": len(products),
        "component_count": len(components),
        "unparseable_components": unparseable_components,
        "metal_tokens": metal_tokens,
        "has_coordination_arrow": has_coordination_arrow,
    }


__all__ = ["classify_reaction_for_taxonomy_benchmark"]
