"""
Motif-based reaction type matching against taxonomy definitions.

Detects compound motifs in reactants/products and matches against reaction type patterns.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set

from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available

from ..motifs.detection import detect_motifs
from ..motifs.registry import build_compound_registry
from ...taxonomy.reaction_catalog import ReactionTypeDefinition, load_reaction_catalog

from .models import ReactionMatch, ReactionDetectionResult


@lru_cache(maxsize=1)
def _load_compound_registry() -> Dict[str, Any]:
    """Load the compiled compound registry for motif detection."""
    base = Path(__file__).resolve().parent.parent.parent / "taxonomy" / "data"
    registry_paths = {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
    }
    return build_compound_registry(registry_paths)


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
    """Build a profile of detected motifs with counts and molecule indices."""
    registry = _load_compound_registry()
    compiled = registry["compiled_compounds"]
    detected: Dict[str, Dict[str, Any]] = {}

    for idx, smiles in enumerate(reactant_smiles):
        mol = parse_smiles(smiles)
        if mol is None:
            continue
        
        from rdkit import Chem
        mol = Chem.AddHs(mol)
        
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


def match_reaction_definition(
    definition: ReactionTypeDefinition,
    detected_motifs: Dict[str, Dict[str, Any]],
    detected_products: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Optional[ReactionMatch]:
    """Match a reaction definition against detected motif profiles."""
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
    
    # Only check products if they were provided in the input
    if definition.products and detected_products is not None:
        if not apply_requirements(definition.products, detected_products, slot_prefix="product:"):
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
            is_product_slot = slot_name.startswith("product:")
            profile = detected_products if is_product_slot else detected_motifs
            
            # If it's a product slot but no products were provided, skip counting
            if is_product_slot and profile is None:
                continue
                
            for motif in hits:
                entry = profile.get(motif) if profile else None
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
