"""Topology-preserving selection of remote subgraphs used by renderers."""

from __future__ import annotations

from typing import Literal, Optional

from ..reaction_models import ReactionCoreProjection, ReactionTopology


def formed_ring_path_subgraph_ids(
    *,
    core: Optional[ReactionCoreProjection],
    topology: Optional[ReactionTopology],
    side: Literal["reactant", "product"] | None = None,
) -> frozenset[str]:
    """Return retained subgraphs required to display observed formed rings.

    A formed ring can result from multiple bonds joining originally separate
    reactant components. Consequently, the overall intramolecular versus
    intermolecular scope cannot decide whether a retained tether is safe to
    abstract. Atom-mapped ring paths are used when available; topology with
    only a ring-size or cycle-rank delta falls back conservatively to all
    retained subgraphs spanning multiple active boundary atoms.
    """
    if topology is None or core is None:
        return frozenset()
    ring_changes = tuple(topology.ring_changes or ())
    ring_sizes = tuple(topology.formed_ring_sizes or ())
    ring_count_delta = topology.ring_count_delta
    if not ring_changes and not ring_sizes and not (
        ring_count_delta is not None and int(ring_count_delta) > 0
    ):
        return frozenset()

    protected = set()
    for subgraph in core.remote_subgraphs:
        if (
            subgraph.continuity != "retained"
            or (side is not None and subgraph.side != side)
        ):
            continue
        boundary_coordinates = {
            (int(port.core_component_index), int(port.core_atom_index))
            for port in subgraph.attachment_ports
        }
        if len(boundary_coordinates) < 2:
            continue
        boundary_maps = {
            int(port.core_atom_map_number)
            for port in subgraph.attachment_ports
            if port.core_atom_map_number is not None
            and int(port.core_atom_map_number) > 0
        }
        remote_coordinates = {
            (int(subgraph.component_index), int(atom_index))
            for atom_index in subgraph.atom_indices
        }
        remote_maps = {
            int(map_number)
            for map_number in subgraph.atom_map_numbers
            if int(map_number) > 0
        }
        comparable_change_found = False
        carries_ring_path = False
        for change in ring_changes:
            references = tuple(change.atom_references or ())
            ring_maps = {
                int(reference.atom_map_number)
                for reference in references
                if reference.atom_map_number is not None
                and int(reference.atom_map_number) > 0
            }
            if ring_maps and boundary_maps:
                comparable_change_found = True
                if (
                    len(boundary_maps.intersection(ring_maps)) >= 2
                    and bool(remote_maps.intersection(ring_maps))
                ):
                    carries_ring_path = True
                    break
            if subgraph.side == "reactant":
                ring_coordinates = {
                    (int(reference.component_index), int(reference.atom_index))
                    for reference in references
                    if reference.side == "reactant"
                }
                if ring_coordinates:
                    comparable_change_found = True
                    if (
                        len(
                            boundary_coordinates.intersection(
                                ring_coordinates
                            )
                        )
                        >= 2
                        and bool(
                            remote_coordinates.intersection(ring_coordinates)
                        )
                    ):
                        carries_ring_path = True
                        break
        if carries_ring_path or not comparable_change_found:
            protected.add(str(subgraph.subgraph_id))
    return frozenset(protected)


__all__ = ["formed_ring_path_subgraph_ids"]
