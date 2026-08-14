"""Graph-derived local identities for fragments departing at broken bonds."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from typing import Any, Mapping, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_parser import parse_reaction_smiles


DEPARTING_FRAGMENT_TOKEN_VERSION = "1.0"


def _atom_key(atom: Mapping[str, Any]) -> tuple[str, int, int]:
    return (
        str(atom.get("side") or ""),
        int(atom.get("component_index") or 0),
        int(atom.get("atom_index") or 0),
    )


def _atom_token(atom: Any) -> tuple[Any, ...]:
    return (
        str(atom.GetSymbol()),
        int(atom.GetFormalCharge()),
        bool(atom.GetIsAromatic()),
        str(atom.GetHybridization()),
    )


def _rooted_paths(
    molecule: Any,
    root_index: int,
    excluded_index: int,
    *,
    radius: int = 3,
) -> tuple[tuple[Any, ...], ...]:
    """Return sorted rooted graph paths on the departing side of one cut."""

    paths = []

    def visit(path: tuple[int, ...]) -> None:
        current = path[-1]
        atom = molecule.GetAtomWithIdx(current)
        encoded: list[Any] = [_atom_token(molecule.GetAtomWithIdx(path[0]))]
        for left, right in zip(path, path[1:]):
            bond = molecule.GetBondBetweenAtoms(left, right)
            encoded.extend(
                (
                    str(bond.GetBondType()),
                    _atom_token(molecule.GetAtomWithIdx(right)),
                )
            )
        paths.append(tuple(encoded))
        if len(path) - 1 >= radius:
            return
        for neighbor in atom.GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if neighbor_index == excluded_index or neighbor_index in path:
                continue
            visit((*path, neighbor_index))

    visit((root_index,))
    return tuple(sorted(paths, key=repr))


def _rooted_cycle_rank(
    molecule: Any,
    root_index: int,
    excluded_index: int,
) -> int:
    """Return the cycle rank of the departing graph after the boundary cut."""

    selected = set()
    stack = [root_index]
    while stack:
        current = stack.pop()
        if current in selected:
            continue
        selected.add(current)
        for neighbor in molecule.GetAtomWithIdx(current).GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if neighbor_index != excluded_index and neighbor_index not in selected:
                stack.append(neighbor_index)
    edge_count = sum(
        1
        for bond in molecule.GetBonds()
        if int(bond.GetBeginAtomIdx()) in selected
        and int(bond.GetEndAtomIdx()) in selected
    )
    return edge_count - len(selected) + 1


def _token(molecule: Any, root_index: int, excluded_index: int) -> str:
    ring_count = _rooted_cycle_rank(molecule, root_index, excluded_index)
    payload = {
        "paths": _rooted_paths(molecule, root_index, excluded_index),
        "ring_count": ring_count,
        "version": DEPARTING_FRAGMENT_TOKEN_VERSION,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return (
        f"DF1:R{ring_count}:"
        + hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:24]
    )


def _neighbor_descriptor(center: Any, neighbor: Any, molecule: Any) -> tuple[Any, ...]:
    bond = molecule.GetBondBetweenAtoms(
        int(center.GetIdx()),
        int(neighbor.GetIdx()),
    )
    return (
        str(neighbor.GetSymbol()),
        str(bond.GetBondType()),
        bool(neighbor.GetIsAromatic()),
        int(neighbor.GetFormalCharge()),
    )


def departing_fragment_tokens(
    reaction_smiles: str,
    signature: Mapping[str, Any],
) -> Tuple[str, ...]:
    """Return local graph tokens for endpoints lost in substitution-like edits.

    A token follows the broken-bond endpoint that does not participate in a
    formed bond and records its rooted radius-three graph after excluding the
    retained endpoint. This distinguishes alcohol, mesylate, tosylate, and
    other leaving groups without reaction names or condition information.
    """

    parsed = parse_reaction_smiles(reaction_smiles)
    if not parsed.valid:
        return ()
    components = {
        component.component_index: component for component in parsed.reactants
    }
    product_atoms_by_map = {}
    for component in parsed.products:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            map_number = int(atom.GetAtomMapNum())
            if map_number > 0:
                product_atoms_by_map[map_number] = (molecule, atom)
    product_map_numbers = set(product_atoms_by_map)
    edits = tuple(
        edit
        for edit in signature.get("edits") or ()
        if isinstance(edit, Mapping)
    )
    formed_atoms = {
        _atom_key(atom)
        for edit in edits
        if str(edit.get("edit_type") or "") == "formed"
        for atom in (edit.get("atom_1"), edit.get("atom_2"))
        if isinstance(atom, Mapping)
        and str(atom.get("side") or "") == "reactant"
    }
    tokens = []
    for edit in edits:
        if str(edit.get("edit_type") or "") != "broken":
            continue
        endpoints = tuple(
            atom
            for atom in (edit.get("atom_1"), edit.get("atom_2"))
            if isinstance(atom, Mapping)
            and str(atom.get("side") or "") == "reactant"
        )
        if len(endpoints) != 2:
            continue
        leaving = tuple(
            atom for atom in endpoints if _atom_key(atom) not in formed_atoms
        )
        retained = tuple(
            atom for atom in endpoints if _atom_key(atom) in formed_atoms
        )
        if len(leaving) != 1 or len(retained) != 1:
            continue
        component_index = int(leaving[0].get("component_index") or 0)
        if component_index != int(retained[0].get("component_index") or 0):
            continue
        component = components.get(component_index)
        molecule = parse_smiles(component.input_smiles) if component else None
        if molecule is None:
            continue
        root_index = int(leaving[0].get("atom_index") or 0)
        excluded_index = int(retained[0].get("atom_index") or 0)
        if not (
            0 <= root_index < molecule.GetNumAtoms()
            and 0 <= excluded_index < molecule.GetNumAtoms()
            and molecule.GetBondBetweenAtoms(root_index, excluded_index) is not None
        ):
            continue
        tokens.append(_token(molecule, root_index, excluded_index))
    if tokens:
        return tuple(sorted(tokens))

    # Some externally mapped datasets retain only the formed bond in the
    # serialized edit set even though the mapped reactant still explicitly
    # contains the leaving group. A mapped neighbor must be absent from every
    # reported product; an unmapped heteroatom/halogen is retained as explicit
    # reactant-side structural evidence rather than treated as correspondence.
    for edit in edits:
        if str(edit.get("edit_type") or "") != "formed":
            continue
        for endpoint in (edit.get("atom_1"), edit.get("atom_2")):
            if not isinstance(endpoint, Mapping):
                continue
            if (
                str(endpoint.get("side") or "") != "reactant"
                or str(endpoint.get("element") or "") != "C"
            ):
                continue
            component_index = int(endpoint.get("component_index") or 0)
            component = components.get(component_index)
            molecule = parse_smiles(component.input_smiles) if component else None
            if molecule is None:
                continue
            carbon_index = int(endpoint.get("atom_index") or 0)
            if not 0 <= carbon_index < molecule.GetNumAtoms():
                continue
            carbon = molecule.GetAtomWithIdx(carbon_index)
            carbon_map_number = int(carbon.GetAtomMapNum())
            retained_unmapped_neighbors: Counter[tuple[Any, ...]] = Counter()
            product_carbon = product_atoms_by_map.get(carbon_map_number)
            if product_carbon is not None:
                product_molecule, product_atom = product_carbon
                retained_unmapped_neighbors.update(
                    _neighbor_descriptor(product_atom, neighbor, product_molecule)
                    for neighbor in product_atom.GetNeighbors()
                    if int(neighbor.GetAtomMapNum()) == 0
                )
            for neighbor in carbon.GetNeighbors():
                neighbor_index = int(neighbor.GetIdx())
                map_number = int(neighbor.GetAtomMapNum())
                if (
                    neighbor.GetSymbol() == "C"
                    or (
                        map_number > 0
                        and map_number in product_map_numbers
                    )
                ):
                    continue
                descriptor = _neighbor_descriptor(carbon, neighbor, molecule)
                if map_number == 0 and retained_unmapped_neighbors[descriptor]:
                    retained_unmapped_neighbors[descriptor] -= 1
                    continue
                tokens.append(_token(molecule, neighbor_index, carbon_index))
    return tuple(sorted(tokens))


__all__ = [
    "DEPARTING_FRAGMENT_TOKEN_VERSION",
    "departing_fragment_tokens",
]
