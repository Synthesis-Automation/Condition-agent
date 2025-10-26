"""
Unified analysis helpers that combine SMILES normalisation, taxonomy-aware
reactant classification, and reaction family resolution.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any, Dict, Iterable, List, Optional

from . import reactions as _reactions
from . import reactants as _reactants
from . import smiles as _smiles
from ._registry import get_registry

__all__ = [
    "normalize",
    "normalize_reaction",
    "ReactantMatch",
    "classify_reactant_smiles",
    "classify_reactant_category",
    "classify_reactant_group",
    "classify_reactant_batch",
    "get_reactant_category_matches",
    "get_all_reactant_matches",
    "normalize_reactant_identifier",
    "normalize_reaction_type",
    "resolve_reaction_family",
    "canonical_family_label",
    "analyze_reaction",
]

# Re-export frequently used helpers for convenience.
normalize = _smiles.normalize
normalize_reaction = _smiles.normalize_reaction
ReactantMatch = _reactants.ReactantMatch
classify_reactant_smiles = _reactants.classify_reactant_smiles
classify_reactant_category = _reactants.classify_reactant_category
classify_reactant_group = _reactants.classify_reactant_group
classify_reactant_batch = _reactants.classify_reactant_batch
get_reactant_category_matches = _reactants.get_reactant_category_matches
get_all_reactant_matches = _reactants.get_all_reactant_matches
normalize_reactant_identifier = _reactants.normalize_reactant_identifier
normalize_reaction_type = _reactants.normalize_reaction_type
resolve_reaction_family = _reactions.resolve_reaction_family
canonical_family_label = _reactions.canonical_family_label


def _match_to_dict(match: Optional[_reactants.ReactantMatch]) -> Optional[Dict[str, Any]]:
    if match is None:
        return None
    if is_dataclass(match):
        return asdict(match)
    return {
        "category": match.category,
        "member_type": match.member_type,
        "name": match.name,
        "group": match.group,
        "smarts": match.smarts,
        "category_smarts": match.category_smarts,
        "description": match.description,
        "specificity": match.specificity,
        "is_general": match.is_general,
    }


def _matches_to_dict(matches: Iterable[_reactants.ReactantMatch]) -> List[Dict[str, Any]]:
    return [m for m in (_match_to_dict(match) for match in matches) if m is not None]


def analyze_reaction(
    reaction_smiles: str,
    *,
    use_rxn_insight: bool = True,
) -> Dict[str, Any]:
    """
    Analyse a reaction SMILES string, returning normalised components,
    per-reactant taxonomy matches, and canonical reaction family metadata.
    """
    norm = normalize_reaction(reaction_smiles)

    reactant_results: List[Dict[str, Any]] = []
    for item in norm.get("reactants") or []:
        smiles_value = (
            item.get("smiles_norm")
            or item.get("largest_smiles")
            or item.get("input")
            or ""
        )
        best = classify_reactant_smiles(smiles_value)
        all_matches = _reactants.get_all_reactant_matches(smiles_value)
        categories = _reactants.get_reactant_category_matches(smiles_value)
        reactant_results.append(
            {
                "normalized": item,
                "taxonomy": {
                    "best_match": _match_to_dict(best),
                    "category_matches": categories,
                    "all_matches": _matches_to_dict(all_matches),
                },
            }
        )

    # Detect reaction family using existing router heuristics.
    from .. import router as _router  # Local import to avoid circular dependency.

    detection = _router.detect_family_from_reaction(
        reaction_smiles, use_rxn_insight=use_rxn_insight
    )
    canonical_family = resolve_reaction_family(detection.get("family"))

    registry = get_registry()
    reaction_type_meta: Optional[Dict[str, Any]] = None
    required_roles: List[Dict[str, Any]] = []
    reactant_requirements: Optional[List[Dict[str, Any]]] = None

    if registry and canonical_family:
        rt = registry.reaction_types.get(canonical_family)
        if rt:
            reaction_type_meta = {
                "id": rt.id,
                "name": rt.name,
                "category_id": rt.category_id,
                "description": rt.description,
                "aliases": list(rt.aliases),
                "metadata": dict(rt.metadata),
                "source_ids": list(rt.source_ids),
            }
            required_roles = [
                {
                    "role_id": req.role_id,
                    "required": bool(req.required),
                    "default_family_id": req.default_family_id,
                    "notes": req.notes,
                }
                for req in rt.required_roles
            ]
            reactant_requirements = [
                {
                    "reactant_type_id": req.reactant_type_id,
                    "stoichiometry": req.stoichiometry,
                    "notes": req.notes,
                    "original_tokens": list(req.original_tokens),
                }
                for req in rt.reactants
            ]

    return {
        "input": reaction_smiles,
        "normalized": norm,
        "reactants": reactant_results,
        "agents": norm.get("agents") or [],
        "products": norm.get("products") or [],
        "family": {
            "detected": detection,
            "canonical_id": canonical_family,
            "reaction_type": reaction_type_meta,
            "required_roles": required_roles,
            "reactant_requirements": reactant_requirements,
        },
    }

