"""Local context classification shared by all reactive-site detectors."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set

from .models import ContextClassification
from .patterns import MatchIndex


_CONTEXTS_PATH = Path(__file__).with_name("definitions") / "context_facets.v2.json"


@lru_cache(maxsize=1)
def load_context_taxonomy() -> Dict[str, Any]:
    """Load context definitions owned by the reactive taxonomy."""
    with _CONTEXTS_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


CONTEXT_PRECEDENCE = tuple(load_context_taxonomy()["precedence"])


@lru_cache(maxsize=1)
def _context_definitions() -> Dict[str, Dict[str, Any]]:
    return {item["id"]: item for item in load_context_taxonomy()["contexts"]}


def _aromatic_ring_context(mol: Any, atom: Any) -> ContextClassification:
    rings = [ring for ring in mol.GetRingInfo().AtomRings() if atom.GetIdx() in ring]
    if not rings:
        return ContextClassification(
            "Ar",
            atom.GetIdx(),
            (atom.GetIdx(),),
            "aromatic_ring_system",
            facet="scaffold",
            semantic_id="carbocyclic_aromatic",
            display_token="Ar",
        )
    ring_atoms: Set[int] = set().union(*(set(ring) for ring in rings))
    heteroatoms = sorted({mol.GetAtomWithIdx(i).GetSymbol() for i in ring_atoms if mol.GetAtomWithIdx(i).GetAtomicNum() != 6})
    token = "HeteroAr" if heteroatoms else "Ar"
    return ContextClassification(
        token=token,
        attachment_atom_index=atom.GetIdx(),
        fragment_atom_indices=tuple(sorted(ring_atoms)),
        classification_method="aromatic_ring_system",
        facet="scaffold",
        semantic_id=(
            "heteroaromatic" if heteroatoms else "carbocyclic_aromatic"
        ),
        display_token=token,
        subtype="heteroaromatic_ring" if heteroatoms else "carbocyclic_aromatic_ring",
        features={
            "heteroatoms": heteroatoms,
            "ring_sizes": sorted({len(ring) for ring in rings}),
            "fused": len(rings) > 1,
        },
    )


def _sp3_carbon_context(mol: Any, atom: Any, excluded: Set[int]) -> ContextClassification:
    """Classify high-value environments adjacent to an sp3 carbon anchor."""
    neighbors = [
        neighbor for neighbor in atom.GetNeighbors()
        if neighbor.GetIdx() not in excluded and neighbor.GetAtomicNum() > 1
    ]
    benzylic = any(neighbor.GetIsAromatic() for neighbor in neighbors)
    allylic = any(
        not neighbor.GetIsAromatic()
        and any(
            bond.GetBondTypeAsDouble() == 2.0
            and bond.GetOtherAtom(neighbor).GetIdx() != atom.GetIdx()
            for bond in neighbor.GetBonds()
        )
        for neighbor in neighbors
    )
    propargylic = any(
        any(
            bond.GetBondTypeAsDouble() == 3.0
            and bond.GetOtherAtom(neighbor).GetIdx() != atom.GetIdx()
            for bond in neighbor.GetBonds()
        )
        for neighbor in neighbors
    )
    subtype = "benzylic" if benzylic else ("allylic" if allylic else ("propargylic" if propargylic else "simple_alkyl"))
    carbon_neighbors = sum(neighbor.GetAtomicNum() == 6 for neighbor in neighbors)
    attachment_classes = {0: "methyl", 1: "primary", 2: "secondary", 3: "tertiary"}
    return ContextClassification(
        token="Alkyl",
        attachment_atom_index=atom.GetIdx(),
        fragment_atom_indices=tuple(sorted({atom.GetIdx(), *(neighbor.GetIdx() for neighbor in neighbors)})),
        classification_method="sp3_attachment_environment",
        facet="scaffold",
        semantic_id="alkyl",
        display_token="Alkyl",
        subtype=subtype,
        features={
            "element": "C",
            "hybridization": str(atom.GetHybridization()),
            "aromatic": False,
            "attachment_carbon_class": attachment_classes.get(carbon_neighbors, "quaternary_or_complex"),
            "benzylic": benzylic,
            "allylic": allylic,
            "propargylic": propargylic,
        },
    )


def classify_context(
    mol: Any,
    atom_index: int,
    excluded: Iterable[int] = (),
    *,
    match_index: MatchIndex | None = None,
) -> ContextClassification:
    """Classify the retained fragment beginning at ``atom_index``."""
    atom = mol.GetAtomWithIdx(atom_index)
    excluded_set = set(excluded)
    symbol = atom.GetSymbol()
    token = "Other"

    definitions = _context_definitions()
    index = match_index or MatchIndex(mol)
    composite_tokens = [
        candidate for candidate in CONTEXT_PRECEDENCE
        if definitions[candidate].get("classification_method") == "mapped_smarts"
    ]
    for candidate in composite_tokens:
        matched_atoms = index.context_match(definitions[candidate], atom_index, excluded_set)
        if matched_atoms is not None:
            return ContextClassification(
                token=candidate,
                attachment_atom_index=atom_index,
                fragment_atom_indices=tuple(sorted(set(matched_atoms))),
                classification_method="mapped_smarts",
                facet=str(definitions[candidate].get("facet") or "fallback"),
                semantic_id=str(
                    definitions[candidate].get("semantic_id") or "other"
                ),
                display_token=str(
                    definitions[candidate].get("display_token") or candidate
                ),
                matched_pattern=candidate,
                features={"excluded_atom_indices": sorted(excluded_set)},
            )
    else:
        if symbol == "C":
            if atom.GetIsAromatic():
                return _aromatic_ring_context(mol, atom)
            elif str(atom.GetHybridization()) == "SP":
                token = "Alkynyl"
            elif str(atom.GetHybridization()) == "SP2":
                token = "Alkenyl"
            else:
                return _sp3_carbon_context(mol, atom, excluded_set)
        elif symbol in {"N", "O", "S"}:
            token = symbol

    method = "element" if symbol in {"N", "O", "S"} else ("fallback" if token == "Other" else "atom_property")
    return ContextClassification(
        token=token,
        attachment_atom_index=atom_index,
        fragment_atom_indices=(atom_index,),
        classification_method=method,
        facet=(
            "element"
            if symbol in {"N", "O", "S"}
            else ("fallback" if token == "Other" else "scaffold")
        ),
        semantic_id={
            "N": "nitrogen",
            "O": "oxygen",
            "S": "sulfur",
            "Alkenyl": "alkenyl",
            "Alkynyl": "alkynyl",
        }.get(token, "other"),
        display_token=token,
        features={"element": symbol, "hybridization": str(atom.GetHybridization()), "aromatic": atom.GetIsAromatic()},
    )


def classify_neighbor_contexts(
    mol: Any,
    center_index: int,
    excluded: Iterable[int] = (),
    *,
    match_index: MatchIndex | None = None,
) -> List[ContextClassification]:
    """Return retained heavy-atom contexts directly attached to a center."""
    excluded_set = set(excluded) | {center_index}
    center = mol.GetAtomWithIdx(center_index)
    contexts = [
        classify_context(mol, neighbor.GetIdx(), excluded_set, match_index=match_index)
        for neighbor in center.GetNeighbors()
        if neighbor.GetAtomicNum() > 1 and neighbor.GetIdx() not in excluded_set
    ]
    rank = {token: i for i, token in enumerate(CONTEXT_PRECEDENCE)}
    return sorted(contexts, key=lambda item: (rank.get(item.token, 999), item.token, item.attachment_atom_index))


__all__ = ["CONTEXT_PRECEDENCE", "classify_context", "classify_neighbor_contexts", "load_context_taxonomy"]
