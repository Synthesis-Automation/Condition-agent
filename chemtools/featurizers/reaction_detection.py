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
    base = Path(__file__).resolve().parent.parent / "taxonomy" / "data"
    registry_paths = {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
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
    product_smiles = [
        item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
        for item in (normalized.get("products") or [])
    ]
    product_smiles = [s for s in product_smiles if s]
    if not reactant_smiles:
        return ReactionDetectionResult(matches=[], error="no_reactants")
    return detect_reaction_types_from_smiles(
        reactant_smiles,
        product_smiles=product_smiles,
        max_hits_per_compound=max_hits_per_compound,
    )


def detect_reaction_types_from_smiles(
    reactant_smiles: Iterable[str],
    *,
    product_smiles: Optional[Iterable[str]] = None,
    max_hits_per_compound: Optional[int] = None,
) -> ReactionDetectionResult:
    """
    Detect reaction types from a list of reactant SMILES strings.
    """
    if not rdkit_available():
        return ReactionDetectionResult(matches=[], error="rdkit_unavailable")

    detected_profile = _detect_motif_profile(
        reactant_smiles, max_hits_per_compound=max_hits_per_compound
    )
    if not detected_profile:
        return ReactionDetectionResult(matches=[])

    product_profile = (
        _detect_motif_profile(product_smiles, max_hits_per_compound=max_hits_per_compound)
        if product_smiles
        else {}
    )

    definitions, _ = load_reaction_catalog()
    matches: List[ReactionMatch] = []

    for definition in definitions.values():
        match = _match_reaction_definition(definition, detected_profile, product_profile)
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
    return set(
        _detect_motif_profile(
            reactant_smiles, max_hits_per_compound=max_hits_per_compound
        ).keys()
    )


def _detect_motif_profile(
    reactant_smiles: Iterable[str],
    *,
    max_hits_per_compound: Optional[int] = None,
) -> Dict[str, Dict[str, Any]]:
    registry = _load_compound_registry()
    compiled = registry["compiled_compounds"]
    detected: Dict[str, Dict[str, Any]] = {}

    for idx, smiles in enumerate(reactant_smiles):
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
                entry = detected.setdefault(
                    compound_id,
                    {"count": 0, "molecules": set()},
                )
                entry["count"] += 1
                entry["molecules"].add(idx)

    return detected


def _match_reaction_definition(
    definition: ReactionTypeDefinition,
    detected_motifs: Dict[str, Dict[str, Any]],
    detected_products: Dict[str, Dict[str, Any]],
) -> Optional[ReactionMatch]:
    slot_evidence: Dict[str, List[str]] = {}
    required_slots = 0
    matched_slots = 0

    if not definition.reactants and not definition.products:
        return None

    def apply_requirements(
        requirements: Dict[str, Any],
        profile: Dict[str, Dict[str, Any]],
        slot_prefix: str = "",
    ) -> bool:
        nonlocal required_slots, matched_slots
        for slot, requirement in requirements.items():
            allowed = requirement.allowed
            if not allowed:
                continue
            required_slots += 1
            hits = [
                motif
                for motif in allowed
                if profile.get(motif, {}).get("count", 0) > 0
            ]
            if not hits:
                return False
            total_hits = sum(profile.get(motif, {}).get("count", 0) for motif in allowed)
            if total_hits < requirement.min_hits:
                return False
            if requirement.min_reactants > 1:
                molecule_indices = set()
                for motif in allowed:
                    entry = profile.get(motif)
                    if entry:
                        molecule_indices.update(entry.get("molecules") or set())
                if len(molecule_indices) < requirement.min_reactants:
                    return False
            slot_name = f"{slot_prefix}{slot}" if slot_prefix else slot
            slot_evidence[slot_name] = hits
            matched_slots += 1
        return True

    if definition.reactants and not apply_requirements(definition.reactants, detected_motifs):
        return None
    if definition.products and not apply_requirements(
        definition.products, detected_products, slot_prefix="product:"
    ):
        return None

    if required_slots == 0:
        return None

    # Apply top-level constraints
    constraints = definition.constraints
    if constraints:
        total_hits = 0
        all_molecule_indices = set()
        for slot_name, hits in slot_evidence.items():
            # Determine if this was a reactant or product slot
            profile = detected_products if slot_name.startswith("product:") else detected_motifs
            for motif in hits:
                entry = profile.get(motif)
                if entry:
                    total_hits += entry.get("count", 0)
                    all_molecule_indices.update(entry.get("molecules") or set())

        if "min_hits" in constraints and total_hits < int(constraints["min_hits"]):
            return None
        if "min_reactants" in constraints and len(all_molecule_indices) < int(constraints["min_reactants"]):
            return None

    return ReactionMatch(
        reaction_type=definition.id,
        name=definition.name,
        category=definition.category or None,
        slot_evidence=slot_evidence,
        matched_slots=matched_slots,
        required_slots=required_slots,
    )
