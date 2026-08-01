"""Reaction-event-aware spectator functional-group derivation."""

from __future__ import annotations

from typing import Dict, List, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import (
    ReactionComponent,
    ReactionEdit,
    ReactionSpectatorGroup,
)


def _spectator_groups_from_atoms(
    components: Tuple[ReactionComponent, ...],
    reactive_by_component: Dict[int, Dict[str, set[int]]],
    *,
    evidence: str,
) -> Tuple[ReactionSpectatorGroup, ...]:
    """Return functional groups outside structurally identified reactive atoms."""
    spectators: List[ReactionSpectatorGroup] = []
    from rdkit import Chem

    for component in components:
        sites = reactive_by_component.get(component.component_index, {})
        mol = parse_smiles(component.input_smiles)
        if mol is None:
            continue
        consumed = set().union(*sites.values()) if sites else set()
        for group in component.molecule_analysis.interpretation.motifs:
            if consumed.intersection(group.atom_indices):
                continue
            nearest_id = None
            nearest_distance = None
            for site_id, atoms in sites.items():
                for left in group.atom_indices:
                    for right in atoms:
                        if int(left) == int(right):
                            distance = 0
                        else:
                            try:
                                distance = (
                                    len(
                                        Chem.GetShortestPath(
                                            mol,
                                            int(left),
                                            int(right),
                                        )
                                    )
                                    - 1
                                )
                            except Exception:
                                continue
                        if nearest_distance is None or distance < nearest_distance:
                            nearest_id, nearest_distance = site_id, distance
            spectators.append(
                ReactionSpectatorGroup(
                    group_id=group.motif_id,
                    chemist_label=group.chemist_label,
                    component_index=component.component_index,
                    atom_indices=group.atom_indices,
                    nearest_reactive_site_id=nearest_id,
                    graph_distance=nearest_distance,
                    tags=group.tags,
                    unchanged_evidence=evidence,
                )
            )
    return tuple(
        sorted(
            spectators,
            key=lambda item: (
                item.component_index,
                item.graph_distance
                if item.graph_distance is not None
                else 999,
                item.group_id,
                item.atom_indices,
            ),
        )
    )


def derive_observed_spectator_groups(
    components: Tuple[ReactionComponent, ...],
    edits: Tuple[ReactionEdit, ...],
    evidence_quality: str,
) -> Tuple[ReactionSpectatorGroup, ...]:
    """Derive spectators only from normalized edit-participating atoms."""
    reactive_by_component: Dict[int, Dict[str, set[int]]] = {}
    for edit in edits:
        for atom in (edit.atom_1, edit.atom_2):
            if atom is None or atom.side != "reactant":
                continue
            site_id = (
                f"edit_center:c{atom.component_index}:a{atom.atom_index}"
            )
            reactive_by_component.setdefault(
                atom.component_index,
                {},
            ).setdefault(site_id, set()).add(atom.atom_index)
    return _spectator_groups_from_atoms(
        components,
        reactive_by_component,
        evidence=f"{evidence_quality}_edit_exclusion",
    )


__all__ = ["derive_observed_spectator_groups"]
