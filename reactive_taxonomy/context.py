"""Local context classification shared by all reactive-site detectors."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set, Tuple

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


def _aromatic_ring_system(
    mol: Any,
    atom_index: int,
) -> Tuple[Tuple[int, ...], ...]:
    """Return the fused aromatic ring system containing ``atom_index``."""
    aromatic_rings = [
        tuple(int(index) for index in ring)
        for ring in mol.GetRingInfo().AtomRings()
        if all(mol.GetAtomWithIdx(index).GetIsAromatic() for index in ring)
    ]
    selected = [ring for ring in aromatic_rings if atom_index in ring]
    ring_atoms = set().union(*(set(ring) for ring in selected)) if selected else set()
    changed = True
    while changed:
        changed = False
        for ring in aromatic_rings:
            if ring in selected or not ring_atoms.intersection(ring):
                continue
            selected.append(ring)
            ring_atoms.update(ring)
            changed = True
    return tuple(sorted(selected, key=lambda ring: (len(ring), ring)))


def _ring_distance(
    mol: Any,
    start: int,
    target: int,
    ring_atoms: Set[int],
) -> int | None:
    """Return shortest graph distance restricted to one aromatic ring system."""
    if start == target:
        return 0
    visited = {start}
    frontier = {start}
    distance = 0
    while frontier:
        distance += 1
        frontier = {
            neighbor.GetIdx()
            for index in frontier
            for neighbor in mol.GetAtomWithIdx(index).GetNeighbors()
            if neighbor.GetIdx() in ring_atoms
            and neighbor.GetIdx() not in visited
        }
        if target in frontier:
            return distance
        visited.update(frontier)
    return None


def _aromatic_heteroatom_type(atom: Any) -> str:
    """Classify an aromatic heteroatom by observable graph state."""
    symbol = atom.GetSymbol()
    charge = int(atom.GetFormalCharge())
    hydrogen_count = int(atom.GetTotalNumHs())
    if symbol == "N":
        if charge > 0 and any(
            neighbor.GetSymbol() == "O" and neighbor.GetFormalCharge() < 0
            for neighbor in atom.GetNeighbors()
        ):
            return "pyridine_n_oxide_like"
        if charge > 0:
            return "cationic_aromatic_nitrogen"
        if hydrogen_count:
            return "pyrrole_like"
        if charge == 0:
            return "pyridine_like"
        return "anionic_aromatic_nitrogen"
    if symbol == "O":
        return "furan_like"
    if symbol == "S":
        return "thiophene_like"
    return f"aromatic_{symbol.lower()}"


def _heteroaromatic_subtype(
    rings: Tuple[Tuple[int, ...], ...],
    heteroatom_details: List[Dict[str, Any]],
) -> str:
    """Return a conservative named annotation for a graph-defined ring system."""
    if not heteroatom_details:
        return "carbocyclic_aromatic_ring"
    if len(rings) > 1:
        kinds = {str(record["aromatic_type"]) for record in heteroatom_details}
        if len(heteroatom_details) == 1:
            kind = next(iter(kinds))
            if kind == "pyridine_like":
                return "fused_pyridine_like"
            if kind == "pyrrole_like":
                return "fused_pyrrole_like"
        return "fused_heteroaromatic"

    ring_size = len(rings[0]) if rings else 0
    element_counts: Dict[str, int] = {}
    for record in heteroatom_details:
        element = str(record["element"])
        element_counts[element] = element_counts.get(element, 0) + 1
    kinds = [str(record["aromatic_type"]) for record in heteroatom_details]
    if ring_size == 6 and set(element_counts) == {"N"}:
        if element_counts["N"] == 1:
            if kinds[0] == "cationic_aromatic_nitrogen":
                return "pyridinium_like"
            return kinds[0]
        if element_counts["N"] == 2:
            return "diazine_like"
        if element_counts["N"] == 3:
            return "triazine_like"
        return "six_membered_aza_arene"
    if ring_size == 5 and len(heteroatom_details) == 1:
        return kinds[0]
    if ring_size == 5:
        return "five_membered_heteroaromatic"
    if ring_size == 6:
        return "six_membered_heteroaromatic"
    return "heteroaromatic_ring"


def _aromatic_ring_context(mol: Any, atom: Any) -> ContextClassification:
    rings = _aromatic_ring_system(mol, atom.GetIdx())
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
    heteroatom_details = []
    for index in sorted(ring_atoms):
        ring_atom = mol.GetAtomWithIdx(index)
        if ring_atom.GetAtomicNum() == 6:
            continue
        heteroatom_details.append({
            "atom_index": index,
            "element": ring_atom.GetSymbol(),
            "formal_charge": int(ring_atom.GetFormalCharge()),
            "hydrogen_count": int(ring_atom.GetTotalNumHs()),
            "aromatic_type": _aromatic_heteroatom_type(ring_atom),
            "distance_from_attachment": _ring_distance(
                mol,
                atom.GetIdx(),
                index,
                ring_atoms,
            ),
        })
    heteroatoms = sorted({
        str(record["element"]) for record in heteroatom_details
    })
    token = "HeteroAr" if heteroatoms else "Ar"
    subtype = _heteroaromatic_subtype(rings, heteroatom_details)
    element_counts = {
        element: sum(
            record["element"] == element for record in heteroatom_details
        )
        for element in heteroatoms
    }
    profile_tokens = sorted(
        f"{record['element']}:{record['aromatic_type']}:"
        f"d{record['distance_from_attachment']}"
        for record in heteroatom_details
    )
    profile = ",".join(profile_tokens)
    ring_sizes = sorted(len(ring) for ring in rings)
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
        subtype=subtype,
        features={
            "heteroatoms": heteroatoms,
            "heteroatom_counts": element_counts,
            "heteroatom_details": heteroatom_details,
            "ring_count": len(rings),
            "ring_sizes": ring_sizes,
            "fused": len(rings) > 1,
            "ring_system_key": (
                f"{token}|rings={','.join(str(size) for size in ring_sizes)}|"
                f"hetero={profile or 'none'}|fused={int(len(rings) > 1)}"
            ),
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
