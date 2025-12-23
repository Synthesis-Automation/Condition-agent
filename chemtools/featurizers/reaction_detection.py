"""
Motif-based reaction type detection for taxonomy v2.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set

from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available

from .motif_detect import detect_motifs
from .motif_registry import build_compound_registry
from ..taxonomy.reaction_catalog import ReactionTypeDefinition, load_reaction_catalog


@dataclass(frozen=True)
class ReactionMatch:
    reaction_type: str
    name: str
    category: Optional[str]
    slot_evidence: Dict[str, List[str]]
    matched_slots: int
    required_slots: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "reaction_type": self.reaction_type,
            "name": self.name,
            "category": self.category,
            "slot_evidence": {slot: list(values) for slot, values in self.slot_evidence.items()},
        }


@dataclass(frozen=True)
class ReactionDetectionResult:
    matches: List[ReactionMatch]
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        payload = {"matches": [match.to_dict() for match in self.matches]}
        if self.error:
            payload["error"] = self.error
        return payload


@lru_cache(maxsize=1)
def _load_compound_registry() -> Dict[str, Any]:
    base = Path(__file__).resolve().parent.parent / "taxonomy" / "v2_data"
    registry_paths = {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
        "templates": base / "smarts_templates.v1.json",
    }
    return build_compound_registry(registry_paths)


def detect_reaction_types(
    reaction_smiles: str,
    *,
    max_hits_per_compound: Optional[int] = None,
) -> ReactionDetectionResult:
    """
    Detect reaction types from a reaction SMILES string.
    """
    if not reaction_smiles:
        return ReactionDetectionResult(matches=[], error="empty_reaction")
    from chemtools.smiles import normalize_reaction
    normalized = normalize_reaction(reaction_smiles)
    reactants = normalized.get("reactants") or []
    reactant_smiles = [
        item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
        for item in reactants
    ]
    reactant_smiles = [s for s in reactant_smiles if s]
    if not reactant_smiles:
        return ReactionDetectionResult(matches=[], error="no_reactants")
    return detect_reaction_types_from_smiles(
        reactant_smiles, max_hits_per_compound=max_hits_per_compound
    )


def detect_reaction_types_from_smiles(
    reactant_smiles: Iterable[str],
    *,
    max_hits_per_compound: Optional[int] = None,
) -> ReactionDetectionResult:
    """
    Detect reaction types from a list of reactant SMILES strings.
    """
    if not rdkit_available():
        return ReactionDetectionResult(matches=[], error="rdkit_unavailable")

    detected_motifs = _detect_motif_ids(
        reactant_smiles, max_hits_per_compound=max_hits_per_compound
    )
    if not detected_motifs:
        return ReactionDetectionResult(matches=[])

    definitions, _ = load_reaction_catalog()
    matches: List[ReactionMatch] = []

    for definition in definitions.values():
        match = _match_reaction_definition(definition, detected_motifs)
        if match is not None:
            matches.append(match)

    matches.sort(
        key=lambda m: (
            -m.matched_slots,
            -sum(len(v) for v in m.slot_evidence.values()),
            m.reaction_type,
        )
    )
    return ReactionDetectionResult(matches=matches)


def detect_motif_ids_from_smiles(
    reactant_smiles: Iterable[str],
    *,
    max_hits_per_compound: Optional[int] = None,
) -> Set[str]:
    """Return the set of organic-compound motif IDs detected in the reactants."""
    if not rdkit_available():
        return set()
    return _detect_motif_ids(reactant_smiles, max_hits_per_compound=max_hits_per_compound)


def _detect_motif_ids(
    reactant_smiles: Iterable[str],
    *,
    max_hits_per_compound: Optional[int] = None,
) -> Set[str]:
    registry = _load_compound_registry()
    compiled = registry["compiled_compounds"]
    detected: Set[str] = set()

    for smiles in reactant_smiles:
        mol = parse_smiles(smiles)
        if mol is None:
            continue
        hits = detect_motifs(
            mol,
            compiled,
            max_hits_per_compound=max_hits_per_compound,
        )
        for hit in hits:
            compound_id = hit.get("compound_id")
            if compound_id:
                detected.add(compound_id)

    return detected


def _match_reaction_definition(
    definition: ReactionTypeDefinition,
    detected_motifs: Set[str],
) -> Optional[ReactionMatch]:
    slot_evidence: Dict[str, List[str]] = {}
    required_slots = 0
    matched_slots = 0

    for slot, allowed in definition.reactants.items():
        if not allowed:
            continue
        required_slots += 1
        hits = [motif for motif in allowed if motif in detected_motifs]
        if not hits:
            return None
        slot_evidence[slot] = hits
        matched_slots += 1

    if required_slots == 0:
        return None

    return ReactionMatch(
        reaction_type=definition.id,
        name=definition.name,
        category=definition.category or None,
        slot_evidence=slot_evidence,
        matched_slots=matched_slots,
        required_slots=required_slots,
    )
