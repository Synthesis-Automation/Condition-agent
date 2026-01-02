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
from ...taxonomy import reaction_catalog as _reaction_catalog
from ..reaction_detection import detect_reaction_types_from_smiles

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
    reactant_smiles_for_detection: List[str] = []
    for item in norm.get("reactants") or []:
        smiles_value = (
            item.get("smiles_norm")
            or item.get("largest_smiles")
            or item.get("input")
            or ""
        )
        if smiles_value:
            reactant_smiles_for_detection.append(smiles_value)
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

    detection_result = detect_reaction_types_from_smiles(reactant_smiles_for_detection)
    detection_payload = detection_result.to_dict()
    matches = detection_result.matches
    canonical_family = matches[0].reaction_type if matches else None

    reaction_type_meta: Optional[Dict[str, Any]] = None
    slot_evidence: Optional[Dict[str, List[str]]] = None

    if canonical_family:
        rt = _reaction_catalog.get_reaction_type(canonical_family)
        if rt:
            reaction_type_meta = {
                "id": rt.id,
                "name": rt.name,
                "category": rt.category,
                "description": rt.description,
                "aliases": list(rt.aliases),
                "metadata": dict(rt.metadata),
                "catalysts": list(rt.catalysts),
                "conditions": rt.conditions,
                "reactants": {slot: list(values) for slot, values in rt.reactants.items()},
                "reference_reactions": list(rt.reference_reactions),
                "notes": rt.notes,
            }
        slot_evidence = {slot: list(values) for slot, values in matches[0].slot_evidence.items()}

    return {
        "input": reaction_smiles,
        "normalized": norm,
        "reactants": reactant_results,
        "agents": norm.get("agents") or [],
        "products": norm.get("products") or [],
        "family": {
            "detected": detection_payload,
            "canonical_id": canonical_family,
            "reaction_type": reaction_type_meta,
            "slot_evidence": slot_evidence,
        },
    }

