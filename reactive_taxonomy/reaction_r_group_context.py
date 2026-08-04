"""Optional functional-group context for graph-provenanced R subgraphs."""

from __future__ import annotations

from collections import deque
from typing import Any, Iterable, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_core.keys import digest
from .reaction_models import (
    ReactionComponent,
    ReactionCoreProjection,
    ReactionRGroupFunctionalContext,
    ReactionRGroupFunctionalGroup,
    ReactionRGroupPortDistance,
    ReactionSpectatorGroup,
)


def _restricted_distance(
    molecule: Any,
    *,
    allowed_atom_indices: set[int],
    source_atom_indices: Iterable[int],
    target_atom_index: int,
) -> int | None:
    """Return the shortest bond distance constrained to one remote subgraph."""
    target = int(target_atom_index)
    sources = {
        int(value)
        for value in source_atom_indices
        if int(value) in allowed_atom_indices
    }
    if target not in allowed_atom_indices or not sources:
        return None
    pending = deque((source, 0) for source in sorted(sources))
    visited = set(sources)
    while pending:
        current, distance = pending.popleft()
        if current == target:
            return distance
        for neighbor in molecule.GetAtomWithIdx(current).GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if (
                neighbor_index in allowed_atom_indices
                and neighbor_index not in visited
            ):
                visited.add(neighbor_index)
                pending.append((neighbor_index, distance + 1))
    return None


def _functional_group(
    molecule: Any,
    *,
    remote_atom_indices: set[int],
    spectator: ReactionSpectatorGroup,
    attachment_ports: Sequence[Any],
) -> ReactionRGroupFunctionalGroup | None:
    if not set(spectator.atom_indices).issubset(remote_atom_indices):
        return None
    distances = []
    for port in attachment_ports:
        distance = _restricted_distance(
            molecule,
            allowed_atom_indices=remote_atom_indices,
            source_atom_indices=spectator.atom_indices,
            target_atom_index=port.attachment_atom_index,
        )
        if distance is None:
            continue
        distances.append(
            ReactionRGroupPortDistance(
                attachment_atom_index=int(port.attachment_atom_index),
                substituent_profile_id=str(port.substituent_profile.profile_id),
                bond_distance=distance,
            )
        )
    if not distances:
        return None
    return ReactionRGroupFunctionalGroup(
        motif_id=str(spectator.group_id),
        chemist_label=str(spectator.chemist_label),
        atom_indices=tuple(sorted(int(value) for value in spectator.atom_indices)),
        tags=tuple(sorted(set(str(value) for value in spectator.tags))),
        nearest_reactive_site_id=spectator.nearest_reactive_site_id,
        distance_to_reactive_site=spectator.graph_distance,
        port_distances=tuple(
            sorted(
                distances,
                key=lambda value: (
                    value.attachment_atom_index,
                    value.substituent_profile_id,
                ),
            )
        ),
        unchanged_evidence=str(spectator.unchanged_evidence),
    )


def build_r_group_functional_contexts(
    *,
    reactants: Tuple[ReactionComponent, ...],
    core: ReactionCoreProjection | None,
    spectator_groups: Tuple[ReactionSpectatorGroup, ...],
) -> Tuple[ReactionRGroupFunctionalContext, ...]:
    """Associate unchanged motifs with containing reactant remote subgraphs."""
    if core is None or not spectator_groups:
        return ()
    components = {component.component_index: component for component in reactants}
    spectators_by_component: dict[int, list[ReactionSpectatorGroup]] = {}
    for spectator in spectator_groups:
        spectators_by_component.setdefault(spectator.component_index, []).append(
            spectator
        )
    contexts = []
    for remote in core.remote_subgraphs:
        if remote.side != "reactant":
            continue
        component = components.get(remote.component_index)
        if component is None:
            continue
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        remote_atoms = set(int(value) for value in remote.atom_indices)
        groups = tuple(
            group
            for spectator in spectators_by_component.get(remote.component_index, ())
            for group in (
                _functional_group(
                    molecule,
                    remote_atom_indices=remote_atoms,
                    spectator=spectator,
                    attachment_ports=remote.attachment_ports,
                ),
            )
            if group is not None
        )
        if not groups:
            continue
        groups = tuple(
            sorted(
                groups,
                key=lambda value: (
                    value.distance_to_reactive_site
                    if value.distance_to_reactive_site is not None
                    else 999,
                    value.motif_id,
                    value.atom_indices,
                ),
            )
        )
        profile_ids = tuple(
            sorted(
                {
                    str(port.substituent_profile.profile_id)
                    for port in remote.attachment_ports
                }
            )
        )
        payload = {
            "remote_subgraph_id": remote.subgraph_id,
            "component_index": remote.component_index,
            "profile_ids": profile_ids,
            "functional_groups": tuple(
                (
                    group.motif_id,
                    group.atom_indices,
                    tuple(
                        (
                            distance.attachment_atom_index,
                            distance.bond_distance,
                        )
                        for distance in group.port_distances
                    ),
                )
                for group in groups
            ),
        }
        contexts.append(
            ReactionRGroupFunctionalContext(
                context_id=digest("RGFC1", payload, length=24),
                remote_subgraph_id=str(remote.subgraph_id),
                side="reactant",
                component_index=int(remote.component_index),
                remote_class=remote.remote_class,
                continuity=remote.continuity,
                remote_atom_indices=tuple(sorted(remote_atoms)),
                attachment_profile_ids=profile_ids,
                functional_groups=groups,
            )
        )
    return tuple(
        sorted(
            contexts,
            key=lambda value: (
                value.component_index,
                value.remote_class,
                value.remote_subgraph_id,
            ),
        )
    )


__all__ = ["build_r_group_functional_contexts"]
