"""Deterministic port-specific profiles for omitted reaction-core fragments."""

from __future__ import annotations

import json
from collections import deque
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence, Tuple

from ..chemistry.rdkit_utils import prepare_fragment_serialization_copy
from .keys import digest
from .models import (
    ReactionCoreAromaticSubstituentRelation,
    ReactionCoreRemoteClass,
    ReactionCoreSubstituentProfile,
)


_DEFINITION_PATH = (
    Path(__file__).resolve().parents[1]
    / "definitions"
    / "substituent_profiles.v1.json"
)
_SCHEMA_VERSION = "1.0"
_DEFINITION_ID = "substituent_profiles.v1"


@lru_cache(maxsize=1)
def load_substituent_profile_definition() -> Mapping[str, Any]:
    """Load and validate the substituent-profile feature vocabulary."""
    with _DEFINITION_PATH.open("r", encoding="utf-8") as handle:
        payload = dict(json.load(handle))
    if payload.get("schema_version") != _SCHEMA_VERSION:
        raise ValueError("unsupported substituent-profile definition schema")
    if payload.get("definition_id") != _DEFINITION_ID:
        raise ValueError("unexpected substituent-profile definition ID")
    if not str(payload.get("algorithm_version") or ""):
        raise ValueError("substituent-profile definition requires an algorithm")
    if not payload.get("base_classes") or not payload.get(
        "carbon_substitution_classes"
    ):
        raise ValueError("substituent-profile vocabularies cannot be empty")
    return payload


def _aromatic_system(molecule: Any, start_index: int) -> Tuple[int, ...]:
    """Return the aromatic-bond-connected system containing ``start_index``."""
    start = molecule.GetAtomWithIdx(int(start_index))
    if not start.GetIsAromatic():
        return ()
    pending = deque((int(start_index),))
    visited: set[int] = set()
    while pending:
        current = pending.popleft()
        if current in visited:
            continue
        visited.add(current)
        atom = molecule.GetAtomWithIdx(current)
        for bond in atom.GetBonds():
            if not bond.GetIsAromatic():
                continue
            neighbor = int(bond.GetOtherAtomIdx(current))
            if neighbor not in visited:
                pending.append(neighbor)
    return tuple(sorted(visited))


def _base_class(
    molecule: Any,
    *,
    fragment: set[int],
    attachment_index: int,
) -> ReactionCoreRemoteClass:
    atom = molecule.GetAtomWithIdx(int(attachment_index))
    if atom.GetIsAromatic():
        aromatic_system = _aromatic_system(molecule, attachment_index)
        if any(
            molecule.GetAtomWithIdx(index).GetAtomicNum() != 6
            for index in aromatic_system
        ):
            return "heteroaryl"
        return "aryl"
    if atom.GetAtomicNum() == 6:
        if any(
            int(neighbor.GetIdx()) in fragment
            and neighbor.GetAtomicNum() in {7, 8, 16}
            and str(
                molecule.GetBondBetweenAtoms(
                    attachment_index, int(neighbor.GetIdx())
                ).GetBondType()
            ).upper()
            == "DOUBLE"
            for neighbor in atom.GetNeighbors()
        ):
            return "acyl"
        hybridization = str(atom.GetHybridization()).upper()
        if hybridization == "SP":
            return "alkynyl"
        if hybridization == "SP2":
            return "alkenyl"
        if atom.IsInRing():
            return "ring_aliphatic"
        return "alkyl"
    if atom.GetAtomicNum() in {7, 8, 9, 15, 16, 17, 35, 53}:
        return "heteroatom"
    return "generic_R"


def _carbon_substitution(
    molecule: Any,
    *,
    fragment: set[int],
    attachment_index: int,
) -> str:
    atom = molecule.GetAtomWithIdx(int(attachment_index))
    if atom.GetAtomicNum() != 6:
        return "not_carbon"
    if atom.GetIsAromatic() or str(atom.GetHybridization()).upper() != "SP3":
        return "not_applicable"
    carbon_neighbors = sum(
        neighbor.GetAtomicNum() == 6 and int(neighbor.GetIdx()) in fragment
        for neighbor in atom.GetNeighbors()
    )
    if carbon_neighbors == 0:
        return "methyl"
    return {
        1: "primary",
        2: "secondary",
        3: "tertiary",
        4: "quaternary",
    }.get(carbon_neighbors, "unresolved")


def _adjacent_multiple_bond(
    molecule: Any,
    *,
    fragment: set[int],
    attachment_index: int,
    order: str,
) -> bool:
    attachment = molecule.GetAtomWithIdx(int(attachment_index))
    for neighbor in attachment.GetNeighbors():
        neighbor_index = int(neighbor.GetIdx())
        if neighbor_index not in fragment:
            continue
        for bond in neighbor.GetBonds():
            other = int(bond.GetOtherAtomIdx(neighbor_index))
            if other == attachment_index or other not in fragment:
                continue
            if str(bond.GetBondType()).upper() == order:
                return True
    return False


def _heteroatoms_at_distance(
    molecule: Any,
    *,
    fragment: set[int],
    attachment_index: int,
    distance: int,
) -> Tuple[str, ...]:
    pending = deque(((int(attachment_index), 0),))
    visited = {int(attachment_index)}
    elements = []
    while pending:
        current, current_distance = pending.popleft()
        if current_distance >= distance:
            continue
        for neighbor in molecule.GetAtomWithIdx(current).GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if neighbor_index not in fragment or neighbor_index in visited:
                continue
            visited.add(neighbor_index)
            next_distance = current_distance + 1
            if next_distance == distance and neighbor.GetAtomicNum() not in {1, 6}:
                elements.append(str(neighbor.GetSymbol()))
            pending.append((neighbor_index, next_distance))
    return tuple(sorted(elements))


def _ring_sizes(molecule: Any, attachment_index: int) -> Tuple[int, ...]:
    return tuple(
        sorted(
            {
                len(ring)
                for ring in molecule.GetRingInfo().AtomRings()
                if int(attachment_index) in ring
            }
        )
    )


def _aromatic_distance(
    molecule: Any,
    *,
    aromatic_system: set[int],
    start_index: int,
    end_index: int,
) -> int:
    pending = deque(((int(start_index), 0),))
    visited = {int(start_index)}
    while pending:
        current, distance = pending.popleft()
        if current == int(end_index):
            return distance
        for bond in molecule.GetAtomWithIdx(current).GetBonds():
            if not bond.GetIsAromatic():
                continue
            neighbor = int(bond.GetOtherAtomIdx(current))
            if neighbor in aromatic_system and neighbor not in visited:
                visited.add(neighbor)
                pending.append((neighbor, distance + 1))
    raise ValueError("aromatic atoms are not connected")


def _branch_component(
    molecule: Any,
    *,
    fragment: set[int],
    aromatic_system: set[int],
    start_index: int,
) -> Tuple[int, ...]:
    pending = deque((int(start_index),))
    visited: set[int] = set()
    while pending:
        current = pending.popleft()
        if current in visited or current not in fragment or current in aromatic_system:
            continue
        visited.add(current)
        pending.extend(
            int(neighbor.GetIdx())
            for neighbor in molecule.GetAtomWithIdx(current).GetNeighbors()
            if int(neighbor.GetIdx()) not in visited
        )
    return tuple(sorted(visited))


def _fragment_smiles(molecule: Any, atom_indices: Sequence[int]) -> str:
    from rdkit import Chem

    copied = prepare_fragment_serialization_copy(molecule, atom_indices)
    for atom in copied.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        return str(
            Chem.MolFragmentToSmiles(
                copied,
                atomsToUse=list(atom_indices),
                canonical=True,
                isomericSmiles=True,
            )
        )
    except Exception:
        return ""


def _aromatic_substituent_relations(
    molecule: Any,
    *,
    fragment: set[int],
    attachment_index: int,
    core_atom_index: int,
) -> Tuple[ReactionCoreAromaticSubstituentRelation, ...]:
    core_atom = molecule.GetAtomWithIdx(int(core_atom_index))
    attachment = molecule.GetAtomWithIdx(int(attachment_index))
    if not core_atom.GetIsAromatic() or not attachment.GetIsAromatic():
        return ()
    aromatic_system = set(_aromatic_system(molecule, int(core_atom_index)))
    if int(attachment_index) not in aromatic_system:
        return ()
    rings = [set(ring) for ring in molecule.GetRingInfo().AtomRings()]
    simple_six_membered = (
        len(aromatic_system) == 6
        and sum(ring == aromatic_system for ring in rings) == 1
    )
    relations = []
    consumed: set[int] = set()
    for ring_atom_index in sorted(aromatic_system):
        ring_atom = molecule.GetAtomWithIdx(ring_atom_index)
        distance = _aromatic_distance(
            molecule,
            aromatic_system=aromatic_system,
            start_index=int(core_atom_index),
            end_index=ring_atom_index,
        )
        position = (
            ("ipso", "ortho", "meta", "para")[distance]
            if simple_six_membered and distance in {0, 1, 2, 3}
            else "other"
        )
        for neighbor in ring_atom.GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if (
                neighbor_index not in fragment
                or neighbor_index in aromatic_system
                or neighbor_index in consumed
                or neighbor.GetAtomicNum() <= 1
            ):
                continue
            branch = _branch_component(
                molecule,
                fragment=fragment,
                aromatic_system=aromatic_system,
                start_index=neighbor_index,
            )
            consumed.update(branch)
            bond = molecule.GetBondBetweenAtoms(ring_atom_index, neighbor_index)
            relations.append(
                ReactionCoreAromaticSubstituentRelation(
                    reactive_atom_index=int(core_atom_index),
                    ring_atom_index=ring_atom_index,
                    aromatic_distance=distance,
                    positional_relation=position,  # type: ignore[arg-type]
                    substituent_attachment_atom_index=neighbor_index,
                    substituent_element=str(neighbor.GetSymbol()),
                    substituent_bond_order=str(bond.GetBondType()).upper(),
                    substituent_fragment_smiles=_fragment_smiles(molecule, branch),
                )
            )
    return tuple(
        sorted(
            relations,
            key=lambda value: (
                value.aromatic_distance,
                value.ring_atom_index,
                value.substituent_fragment_smiles,
            ),
        )
    )


def _neighbor_state(molecule: Any, center_index: int, neighbor_index: int) -> tuple[Any, ...]:
    neighbor = molecule.GetAtomWithIdx(int(neighbor_index))
    bond = molecule.GetBondBetweenAtoms(int(center_index), int(neighbor_index))
    return (
        neighbor.GetSymbol(),
        int(neighbor.GetFormalCharge()),
        bool(neighbor.GetIsAromatic()),
        str(neighbor.GetHybridization()).upper(),
        str(bond.GetBondType()).upper(),
        int(neighbor.IsInRing()),
    )


def _radius_environment(
    molecule: Any,
    *,
    fragment: set[int],
    attachment_index: int,
    radius: int = 2,
) -> Tuple[tuple[Any, ...], ...]:
    pending = deque(((int(attachment_index), 0),))
    visited = {int(attachment_index)}
    values = []
    while pending:
        current, distance = pending.popleft()
        if distance >= radius:
            continue
        for neighbor in molecule.GetAtomWithIdx(current).GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if neighbor_index not in fragment or neighbor_index in visited:
                continue
            visited.add(neighbor_index)
            next_distance = distance + 1
            values.append(
                (next_distance, *_neighbor_state(molecule, current, neighbor_index))
            )
            pending.append((neighbor_index, next_distance))
    return tuple(sorted(values, key=repr))


def _beta_branch_count(
    molecule: Any,
    *,
    fragment: set[int],
    attachment_index: int,
) -> int:
    total = 0
    attachment = molecule.GetAtomWithIdx(int(attachment_index))
    for neighbor in attachment.GetNeighbors():
        neighbor_index = int(neighbor.GetIdx())
        if neighbor_index not in fragment or neighbor.GetAtomicNum() != 6:
            continue
        carbon_neighbors = sum(
            other.GetAtomicNum() == 6
            and int(other.GetIdx()) in fragment
            and int(other.GetIdx()) != attachment_index
            for other in neighbor.GetNeighbors()
        )
        total += max(0, carbon_neighbors - 1)
    return total


def _feature_tokens(
    *,
    base_class: str,
    carbon_substitution: str,
    cyclic: bool,
    ring_sizes: Sequence[int],
    benzylic: bool,
    allylic: bool,
    propargylic: bool,
    alpha_branch_count: int,
    beta_branch_count: int,
    radius_1_heteroatoms: Iterable[str],
    radius_2_heteroatoms: Iterable[str],
    aromatic_substituent_relations: Sequence[
        ReactionCoreAromaticSubstituentRelation
    ],
    environment_key: str,
) -> Tuple[str, ...]:
    tokens = {
        f"L0:class:{base_class}",
        f"L1:substitution:{carbon_substitution}",
        f"L1:cyclic:{int(cyclic)}",
        f"L2:alpha_branch:{alpha_branch_count}",
        f"L2:beta_branch:{beta_branch_count}",
        f"L3:environment:{environment_key}",
    }
    tokens.update(f"L2:ring_size:{size}" for size in ring_sizes)
    tokens.update(f"L2:hetero_r1:{value}" for value in radius_1_heteroatoms)
    tokens.update(f"L2:hetero_r2:{value}" for value in radius_2_heteroatoms)
    tokens.update(
        "L2:aromatic_substituent:"
        f"{relation.positional_relation}:"
        f"{relation.substituent_element}:"
        f"{relation.substituent_fragment_smiles or '-'}"
        for relation in aromatic_substituent_relations
    )
    for name, present in (
        ("benzylic", benzylic),
        ("allylic", allylic),
        ("propargylic", propargylic),
    ):
        if present:
            tokens.add(f"L1:{name}")
    return tuple(sorted(tokens))


def build_substituent_profile(
    molecule: Any,
    *,
    fragment_atom_indices: Sequence[int],
    attachment_atom_index: int,
    core_atom_index: int,
    attachment_bond_order: str,
) -> ReactionCoreSubstituentProfile:
    """Build one graph-derived profile at a retained/omitted attachment port."""
    definition = load_substituent_profile_definition()
    fragment = {int(value) for value in fragment_atom_indices}
    attachment_index = int(attachment_atom_index)
    if attachment_index not in fragment:
        raise ValueError("substituent attachment atom must belong to the fragment")
    if int(core_atom_index) in fragment:
        raise ValueError("substituent core atom cannot belong to the omitted fragment")
    attachment = molecule.GetAtomWithIdx(attachment_index)
    base_class = _base_class(
        molecule,
        fragment=fragment,
        attachment_index=attachment_index,
    )
    carbon_substitution = _carbon_substitution(
        molecule,
        fragment=fragment,
        attachment_index=attachment_index,
    )
    carbon_neighbor_count = sum(
        neighbor.GetAtomicNum() == 6 and int(neighbor.GetIdx()) in fragment
        for neighbor in attachment.GetNeighbors()
    )
    alpha_branch_count = max(0, carbon_neighbor_count - 1)
    ring_sizes = _ring_sizes(molecule, attachment_index)
    is_saturated_carbon = (
        attachment.GetAtomicNum() == 6
        and not attachment.GetIsAromatic()
        and str(attachment.GetHybridization()).upper() == "SP3"
    )
    benzylic = is_saturated_carbon and any(
        neighbor.GetIsAromatic() and int(neighbor.GetIdx()) in fragment
        for neighbor in attachment.GetNeighbors()
    )
    allylic = is_saturated_carbon and _adjacent_multiple_bond(
        molecule,
        fragment=fragment,
        attachment_index=attachment_index,
        order="DOUBLE",
    )
    propargylic = is_saturated_carbon and _adjacent_multiple_bond(
        molecule,
        fragment=fragment,
        attachment_index=attachment_index,
        order="TRIPLE",
    )
    radius_1_heteroatoms = _heteroatoms_at_distance(
        molecule,
        fragment=fragment,
        attachment_index=attachment_index,
        distance=1,
    )
    radius_2_heteroatoms = _heteroatoms_at_distance(
        molecule,
        fragment=fragment,
        attachment_index=attachment_index,
        distance=2,
    )
    beta_branch_count = _beta_branch_count(
        molecule,
        fragment=fragment,
        attachment_index=attachment_index,
    )
    aromatic_substituent_relations = _aromatic_substituent_relations(
        molecule,
        fragment=fragment,
        attachment_index=attachment_index,
        core_atom_index=int(core_atom_index),
    )
    environment_payload = {
        "attachment": (
            attachment.GetSymbol(),
            int(attachment.GetFormalCharge()),
            bool(attachment.GetIsAromatic()),
            str(attachment.GetHybridization()).upper(),
            str(attachment_bond_order).upper(),
            int(attachment.IsInRing()),
        ),
        "radius_environment": _radius_environment(
            molecule,
            fragment=fragment,
            attachment_index=attachment_index,
        ),
        "definition_version": str(definition["definition_id"]),
    }
    environment_key = digest("RSE1", environment_payload, length=32)
    feature_tokens = _feature_tokens(
        base_class=base_class,
        carbon_substitution=carbon_substitution,
        cyclic=bool(ring_sizes),
        ring_sizes=ring_sizes,
        benzylic=benzylic,
        allylic=allylic,
        propargylic=propargylic,
        alpha_branch_count=alpha_branch_count,
        beta_branch_count=beta_branch_count,
        radius_1_heteroatoms=radius_1_heteroatoms,
        radius_2_heteroatoms=radius_2_heteroatoms,
        aromatic_substituent_relations=aromatic_substituent_relations,
        environment_key=environment_key,
    )
    profile_payload = {
        "base_class": base_class,
        "attachment_element": attachment.GetSymbol(),
        "attachment_bond_order": str(attachment_bond_order).upper(),
        "carbon_substitution": carbon_substitution,
        "feature_tokens": feature_tokens,
        "definition_version": str(definition["definition_id"]),
        "algorithm_version": str(definition["algorithm_version"]),
    }
    return ReactionCoreSubstituentProfile(
        profile_id=digest("RSP1", profile_payload, length=32),
        base_class=base_class,
        attachment_element=str(attachment.GetSymbol()),
        attachment_bond_order=str(attachment_bond_order).upper(),
        attachment_aromatic=bool(attachment.GetIsAromatic()),
        attachment_hybridization=str(attachment.GetHybridization()).upper(),
        carbon_substitution=carbon_substitution,  # type: ignore[arg-type]
        cyclic=bool(ring_sizes),
        ring_sizes=ring_sizes,
        benzylic=benzylic,
        allylic=allylic,
        propargylic=propargylic,
        alpha_branch_count=alpha_branch_count,
        beta_branch_count=beta_branch_count,
        radius_1_heteroatoms=radius_1_heteroatoms,
        radius_2_heteroatoms=radius_2_heteroatoms,
        aromatic_substituent_relations=aromatic_substituent_relations,
        local_environment_key=environment_key,
        feature_tokens=feature_tokens,
        definition_version=str(definition["definition_id"]),
        algorithm_version=str(definition["algorithm_version"]),
    )


__all__ = [
    "build_substituent_profile",
    "load_substituent_profile_definition",
]
