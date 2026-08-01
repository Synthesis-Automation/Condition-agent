"""Remote-subgraph observation and continuity for minimized reaction cores."""

from __future__ import annotations

from collections import deque
from dataclasses import replace
from typing import Any, Mapping, Optional, Sequence, Tuple

from ..chemistry.rdkit_utils import parse_smiles
from ..reaction_models import ReactionComponent
from .common import Coordinate as _Coordinate, Location as _Location
from .common import atom_map_number as _atom_map_number
from .keys import digest as _digest
from .models import (
    ReactionCoreAttachmentPort,
    ReactionCoreRemoteClass,
    ReactionCoreRemoteSubgraph,
)


def _connected_remote_components(
    molecule: Any,
    active_atom_indices: set[int],
) -> Tuple[Tuple[int, ...], ...]:
    remaining = {
        int(atom.GetIdx())
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1 and int(atom.GetIdx()) not in active_atom_indices
    }
    values = []
    while remaining:
        start = min(remaining)
        queue = deque((start,))
        component = set()
        while queue:
            atom_index = queue.popleft()
            if atom_index in component or atom_index not in remaining:
                continue
            component.add(atom_index)
            atom = molecule.GetAtomWithIdx(atom_index)
            for neighbor in atom.GetNeighbors():
                neighbor_index = int(neighbor.GetIdx())
                if (
                    neighbor.GetAtomicNum() > 1
                    and neighbor_index in remaining
                    and neighbor_index not in component
                ):
                    queue.append(neighbor_index)
        remaining.difference_update(component)
        if any(
            int(neighbor.GetIdx()) in active_atom_indices
            for atom_index in component
            for neighbor in molecule.GetAtomWithIdx(atom_index).GetNeighbors()
        ):
            values.append(tuple(sorted(component)))
    return tuple(values)


def _fragment_smiles(molecule: Any, atom_indices: Sequence[int]) -> str:
    from rdkit import Chem

    copied = Chem.Mol(molecule)
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


def _remote_class(
    molecule: Any,
    atom_indices: Sequence[int],
    attachment_indices: Sequence[int],
) -> ReactionCoreRemoteClass:
    atoms = [
        molecule.GetAtomWithIdx(int(atom_index))
        for atom_index in attachment_indices
    ]
    if any(atom.GetIsAromatic() for atom in atoms):
        aromatic_atoms = [
            molecule.GetAtomWithIdx(int(atom_index))
            for atom_index in atom_indices
            if molecule.GetAtomWithIdx(int(atom_index)).GetIsAromatic()
        ]
        if any(atom.GetAtomicNum() != 6 for atom in aromatic_atoms):
            return "heteroaryl"
        return "aryl"
    if atoms and all(atom.GetAtomicNum() == 6 for atom in atoms):
        fragment = set(int(value) for value in atom_indices)
        if any(
            neighbor.GetAtomicNum() in {7, 8, 16}
            and int(neighbor.GetIdx()) in fragment
            and str(
                molecule.GetBondBetweenAtoms(
                    int(atom.GetIdx()),
                    int(neighbor.GetIdx()),
                ).GetBondType()
            ).upper()
            == "DOUBLE"
            for atom in atoms
            for neighbor in atom.GetNeighbors()
        ):
            return "acyl"
        hybridizations = {str(atom.GetHybridization()).upper() for atom in atoms}
        if hybridizations == {"SP"}:
            return "alkynyl"
        if "SP2" in hybridizations:
            return "alkenyl"
        if any(atom.IsInRing() for atom in atoms):
            return "ring_aliphatic"
        return "alkyl"
    if atoms and all(
        atom.GetAtomicNum() in {7, 8, 9, 15, 16, 17, 35, 53}
        for atom in atoms
    ):
        return "heteroatom"
    return "generic_R"


def _build_remote_subgraphs_for_side(
    *,
    side: str,
    components: Sequence[ReactionComponent],
    active_coordinates: set[_Coordinate],
    atom_map_overrides: Optional[
        Mapping[tuple[str, int, int], int]
    ] = None,
) -> Tuple[ReactionCoreRemoteSubgraph, ...]:
    values = []
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        active = {
            atom_index
            for component_index, atom_index in active_coordinates
            if component_index == component.component_index
        }
        for atom_indices in _connected_remote_components(molecule, active):
            fragment = set(atom_indices)
            ports = []
            for attachment_index in atom_indices:
                attachment = molecule.GetAtomWithIdx(attachment_index)
                for core_atom in attachment.GetNeighbors():
                    core_index = int(core_atom.GetIdx())
                    if core_index not in active:
                        continue
                    bond = molecule.GetBondBetweenAtoms(
                        core_index,
                        attachment_index,
                    )
                    ports.append(
                        ReactionCoreAttachmentPort(
                            side=side,  # type: ignore[arg-type]
                            core_component_index=component.component_index,
                            core_atom_index=core_index,
                            core_atom_map_number=(
                                _atom_map_number(
                                    core_atom,
                                    side=side,
                                    component_index=component.component_index,
                                    atom_index=core_index,
                                    atom_map_overrides=atom_map_overrides,
                                )
                                or None
                            ),
                            attachment_atom_index=attachment_index,
                            attachment_atom_map_number=(
                                _atom_map_number(
                                    attachment,
                                    side=side,
                                    component_index=component.component_index,
                                    atom_index=attachment_index,
                                    atom_map_overrides=atom_map_overrides,
                                )
                                or None
                            ),
                            attachment_element=attachment.GetSymbol(),
                            bond_order=str(bond.GetBondType()).upper(),
                        )
                    )
            if not ports:
                continue
            ports = sorted(
                ports,
                key=lambda port: (
                    port.core_component_index,
                    port.core_atom_index,
                    port.attachment_atom_index,
                    port.bond_order,
                ),
            )
            attachment_indices = tuple(
                sorted({port.attachment_atom_index for port in ports})
            )
            fragment_smiles = _fragment_smiles(molecule, atom_indices)
            remote_class = _remote_class(
                molecule,
                atom_indices,
                attachment_indices,
            )
            map_numbers = tuple(
                sorted(
                    _atom_map_number(
                        atom,
                        side=side,
                        component_index=component.component_index,
                        atom_index=atom_index,
                        atom_map_overrides=atom_map_overrides,
                    )
                    for atom_index in atom_indices
                    if _atom_map_number(
                        (atom := molecule.GetAtomWithIdx(atom_index)),
                        side=side,
                        component_index=component.component_index,
                        atom_index=atom_index,
                        atom_map_overrides=atom_map_overrides,
                    ) > 0
                )
            )
            payload = {
                "side": side,
                "component_index": component.component_index,
                "atom_indices": atom_indices,
                "fragment_smiles": fragment_smiles,
                "remote_class": remote_class,
                "ports": tuple(
                    (
                        port.core_atom_index,
                        port.attachment_atom_index,
                        port.attachment_element,
                        port.bond_order,
                    )
                    for port in ports
                ),
            }
            values.append(
                ReactionCoreRemoteSubgraph(
                    subgraph_id=_digest("RCR2", payload, length=24),
                    side=side,  # type: ignore[arg-type]
                    component_index=component.component_index,
                    atom_indices=atom_indices,
                    atom_map_numbers=map_numbers,
                    remote_class=remote_class,
                    continuity="unresolved",
                    attachment_ports=tuple(ports),
                    fragment_smiles=fragment_smiles,
                    fragment_heavy_atom_count=len(fragment),
                    fragment_heteroatom_count=sum(
                        molecule.GetAtomWithIdx(atom_index).GetAtomicNum()
                        not in {1, 6}
                        for atom_index in fragment
                    ),
                    fragment_aromatic_atom_count=sum(
                        molecule.GetAtomWithIdx(atom_index).GetIsAromatic()
                        for atom_index in fragment
                    ),
                )
            )
    return tuple(
        sorted(
            values,
            key=lambda subgraph: (
                subgraph.side,
                subgraph.remote_class,
                subgraph.fragment_smiles,
                subgraph.component_index,
                subgraph.atom_indices,
            ),
        )
    )


def _mapped_port_tokens(
    subgraph: ReactionCoreRemoteSubgraph,
) -> Tuple[Tuple[int, int, str], ...]:
    return tuple(
        sorted(
            (
                int(port.core_atom_map_number),
                int(port.attachment_atom_map_number),
                port.bond_order,
            )
            for port in subgraph.attachment_ports
            if port.core_atom_map_number is not None
            and port.attachment_atom_map_number is not None
        )
    )


def _remote_continuity(
    subgraph: ReactionCoreRemoteSubgraph,
    opposite: Sequence[ReactionCoreRemoteSubgraph],
    *,
    opposite_by_map: Mapping[int, _Location],
) -> str:
    mapped = set(subgraph.atom_map_numbers)
    if mapped:
        exact = [
            candidate
            for candidate in opposite
            if set(candidate.atom_map_numbers) == mapped
        ]
        if exact:
            candidate = min(
                exact,
                key=lambda value: (
                    value.fragment_smiles,
                    value.component_index,
                    value.atom_indices,
                ),
            )
            if (
                candidate.fragment_smiles == subgraph.fragment_smiles
                and candidate.remote_class == subgraph.remote_class
                and _mapped_port_tokens(candidate)
                == _mapped_port_tokens(subgraph)
            ):
                return "retained"
            return "changed"
        if any(mapped.intersection(candidate.atom_map_numbers) for candidate in opposite):
            return "changed"
        return "departing" if subgraph.side == "reactant" else "appearing"
    core_maps = {
        int(port.core_atom_map_number)
        for port in subgraph.attachment_ports
        if port.core_atom_map_number is not None
    }
    if core_maps and all(map_number not in opposite_by_map for map_number in core_maps):
        return "departing" if subgraph.side == "reactant" else "appearing"
    port_shapes = {
        (
            int(port.core_atom_map_number),
            port.attachment_element,
            port.bond_order,
        )
        for port in subgraph.attachment_ports
        if port.core_atom_map_number is not None
    }
    opposite_port_shapes = {
        (
            int(port.core_atom_map_number),
            port.attachment_element,
            port.bond_order,
        )
        for candidate in opposite
        for port in candidate.attachment_ports
        if port.core_atom_map_number is not None
    }
    if port_shapes and not port_shapes.intersection(opposite_port_shapes):
        return "departing" if subgraph.side == "reactant" else "appearing"
    return "unresolved"


def _with_remote_continuity(
    reactant_subgraphs: Sequence[ReactionCoreRemoteSubgraph],
    product_subgraphs: Sequence[ReactionCoreRemoteSubgraph],
    *,
    reactant_by_map: Mapping[int, _Location],
    product_by_map: Mapping[int, _Location],
) -> Tuple[ReactionCoreRemoteSubgraph, ...]:
    values = [
        replace(
            subgraph,
            continuity=_remote_continuity(
                subgraph,
                product_subgraphs,
                opposite_by_map=product_by_map,
            ),  # type: ignore[arg-type]
        )
        for subgraph in reactant_subgraphs
    ]
    values.extend(
        replace(
            subgraph,
            continuity=_remote_continuity(
                subgraph,
                reactant_subgraphs,
                opposite_by_map=reactant_by_map,
            ),  # type: ignore[arg-type]
        )
        for subgraph in product_subgraphs
    )
    return tuple(
        sorted(
            values,
            key=lambda subgraph: (
                subgraph.side,
                subgraph.remote_class,
                subgraph.fragment_smiles,
                subgraph.component_index,
                subgraph.atom_indices,
            ),
        )
    )


__all__ = [
    "_build_remote_subgraphs_for_side",
    "_with_remote_continuity",
]
