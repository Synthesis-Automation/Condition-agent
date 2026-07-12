"""Reaction-event-aware spectator functional-group derivation."""

from __future__ import annotations

from typing import Dict, List, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import ReactionCandidate, ReactionComponent, ReactionSpectatorGroup


def derive_spectator_groups(
    components: Tuple[ReactionComponent, ...],
    selected: ReactionCandidate | None,
    evidence_quality: str,
) -> Tuple[ReactionSpectatorGroup, ...]:
    """Return substrate groups that do not intersect the selected site atoms."""
    if selected is None:
        return ()
    reactive_by_component: Dict[int, Dict[str, set[int]]] = {}
    for site in selected.role_assignments.values():
        site_atoms = {index for indices in site.atom_roles.values() for index in indices}
        reactive_by_component.setdefault(site.component_index, {})[site.site_id] = site_atoms
    evidence = (
        "exact_product_reconstruction_event_exclusion"
        if evidence_quality == "exact_product_reconstruction"
        else "selected_event_exclusion"
    )
    spectators: List[ReactionSpectatorGroup] = []
    from rdkit import Chem
    for component in components:
        sites = reactive_by_component.get(component.component_index, {})
        mol = parse_smiles(component.input_smiles)
        if mol is None:
            continue
        consumed = set().union(*sites.values()) if sites else set()
        for group in component.compound_analysis.functional_groups:
            if consumed.intersection(group.atom_indices):
                continue
            nearest_id = None
            nearest_distance = None
            for site_id, atoms in sites.items():
                for left in group.atom_indices:
                    for right in atoms:
                        if int(left) == int(right):
                            distance = 0
                            if nearest_distance is None or distance < nearest_distance:
                                nearest_id, nearest_distance = site_id, distance
                            continue
                        try:
                            distance = len(Chem.GetShortestPath(mol, int(left), int(right))) - 1
                        except Exception:
                            continue
                        if nearest_distance is None or distance < nearest_distance:
                            nearest_id, nearest_distance = site_id, distance
            spectators.append(ReactionSpectatorGroup(
                group_id=group.group_id,
                chemist_label=group.chemist_label,
                component_index=component.component_index,
                atom_indices=group.atom_indices,
                nearest_reactive_site_id=nearest_id,
                graph_distance=nearest_distance,
                tags=group.tags,
                unchanged_evidence=evidence,
            ))
    return tuple(sorted(spectators, key=lambda item: (item.component_index, item.graph_distance if item.graph_distance is not None else 999, item.group_id, item.atom_indices)))


__all__ = ["derive_spectator_groups"]
