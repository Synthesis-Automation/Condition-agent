"""Reaction parsing, typing, atom mapping, events, CRK, and featurization."""

from .atom_mapping import (
    add_atom_mapping,
    analyze_bond_changes,
    analyze_bond_changes_hybrid,
    analyze_with_mcs,
    is_available as rxnmapper_available,
)
from .featurize import (
    CrkResult,
    build_crk,
    classify_agent_roles,
    featurize_reaction,
    format_bond_change_key,
    get_crk_options,
)
from .inference import (
    GeneralReactionAnalysis,
    ReactionDecision,
    ReactionValidation,
    analyze_reaction_general,
    classify_reaction,
)
from .models import ReactionEvents, ReactionFeatures
from .typing import (
    DetectionResult,
    ReactionMatch,
    detect_reaction_type,
    detect_reaction_types,
    extract_reaction_key,
)


def detect_reaction(reaction_smiles: str, use_ml: bool = True) -> dict:
    """Return a simple family-style detection payload for app services."""
    result = detect_reaction_type(reaction_smiles)
    best = result.best_match
    family = best.reaction_type if best else "Unknown"
    confidence = best.confidence if best else 0.0
    return {
        "family": family,
        "confidence": confidence,
        "method": "taxonomy_rule",
        "details": {
            "use_ml_requested": bool(use_ml),
            "functional_groups": {},
            "matches": [
                {
                    "reaction_type": match.reaction_type,
                    "confidence": match.confidence,
                    "electrophile": match.electrophile,
                    "nucleophile": match.nucleophile,
                    "product": match.product,
                    "slot_sources": match.slot_sources,
                }
                for match in result.matches
            ],
        },
    }

__all__ = [
    "CrkResult",
    "DetectionResult",
    "GeneralReactionAnalysis",
    "ReactionDecision",
    "ReactionEvents",
    "ReactionFeatures",
    "ReactionMatch",
    "ReactionValidation",
    "add_atom_mapping",
    "analyze_bond_changes",
    "analyze_bond_changes_hybrid",
    "analyze_reaction_general",
    "analyze_with_mcs",
    "build_crk",
    "classify_agent_roles",
    "classify_reaction",
    "detect_reaction_type",
    "detect_reaction",
    "detect_reaction_types",
    "extract_reaction_key",
    "featurize_reaction",
    "format_bond_change_key",
    "get_crk_options",
    "rxnmapper_available",
]
