"""Deterministic partitioning of normalized edits into reaction events."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from typing import Any, Optional, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_archetypes import reconcile_edit_archetype
from .reaction_label_patterns import match_reaction_label_pattern
from .reaction_models import (
    ReactionAtomReference,
    ReactionCandidate,
    ReactionComponent,
    ReactionEdit,
    ReactionEvent,
    ReactionEventRelation,
    ReactionPartner,
    ReactionTopology,
)


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def _digest(prefix: str, value: Any, *, length: int = 24) -> str:
    encoded = _canonical_json(value).encode("utf-8")
    return f"{prefix}:" + hashlib.sha256(encoded).hexdigest()[:length]


def _atom_identity(atom: ReactionAtomReference) -> Tuple[Any, ...]:
    if atom.atom_map_number is not None:
        return ("map", int(atom.atom_map_number))
    return (
        "atom",
        atom.side,
        int(atom.component_index),
        int(atom.atom_index),
    )


def _atom_chemistry(atom: Optional[ReactionAtomReference]) -> Tuple[Any, ...]:
    if atom is None:
        return ("H",)
    return (
        atom.element,
        int(atom.formal_charge),
        bool(atom.aromatic),
        atom.hybridization,
        atom.local_environment_id,
    )


def _edit_chemistry(edit: ReactionEdit) -> Tuple[Any, ...]:
    endpoints = tuple(
        sorted((_atom_chemistry(edit.atom_1), _atom_chemistry(edit.atom_2)))
    )
    return (
        edit.edit_type,
        endpoints,
        edit.old_order or "NONE",
        edit.new_order or "NONE",
    )


def _edit_provenance(edit: ReactionEdit) -> Tuple[Any, ...]:
    endpoints = tuple(
        sorted(
            _atom_identity(atom)
            for atom in (edit.atom_1, edit.atom_2)
            if atom is not None
        )
    )
    return (
        edit.edit_type,
        endpoints,
        edit.old_order or "NONE",
        edit.new_order or "NONE",
    )


def _bond_type(edit: ReactionEdit, order: Optional[str]) -> str:
    if edit.atom_2 is None:
        return f"{edit.atom_1.element}-H:{order or 'NONE'}"
    elements = sorted((edit.atom_1.element, edit.atom_2.element))
    return f"{elements[0]}-{elements[1]}:{order or 'NONE'}"


def partition_reaction_edits(
    edits: Sequence[ReactionEdit],
) -> Tuple[Tuple[ReactionEdit, ...], ...]:
    """Return connected edit groups using shared atom provenance."""
    if not edits:
        return ()
    parents = list(range(len(edits)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    owners: dict[Tuple[Any, ...], int] = {}
    for index, edit in enumerate(edits):
        for atom in (edit.atom_1, edit.atom_2):
            if atom is None:
                continue
            key = _atom_identity(atom)
            previous = owners.setdefault(key, index)
            union(index, previous)

    groups: dict[int, list[ReactionEdit]] = defaultdict(list)
    for index, edit in enumerate(edits):
        groups[find(index)].append(edit)
    normalized = [
        tuple(sorted(group, key=_edit_provenance)) for group in groups.values()
    ]
    return tuple(
        sorted(
            normalized,
            key=lambda group: (
                tuple(_edit_chemistry(edit) for edit in group),
                tuple(_edit_provenance(edit) for edit in group),
            ),
        )
    )


def _event_topology(
    reactants: Tuple[ReactionComponent, ...],
    edits: Sequence[ReactionEdit],
    *,
    evidence: str,
    confidence: float,
) -> ReactionTopology:
    components = {component.component_index: component for component in reactants}
    participating = {
        atom.component_index
        for edit in edits
        for atom in (edit.atom_1, edit.atom_2)
        if atom is not None and atom.side == "reactant"
    }
    formed_scopes = []
    tether_distances = []
    formed_ring_sizes = []
    from rdkit import Chem

    for edit in edits:
        if edit.edit_type != "formed" or edit.atom_2 is None:
            continue
        if edit.atom_1.side != "reactant" or edit.atom_2.side != "reactant":
            continue
        if edit.atom_1.component_index != edit.atom_2.component_index:
            formed_scopes.append("intermolecular")
            continue
        formed_scopes.append("intramolecular")
        component = components.get(edit.atom_1.component_index)
        molecule = parse_smiles(component.input_smiles) if component else None
        if molecule is None:
            continue
        left = int(edit.atom_1.atom_index)
        right = int(edit.atom_2.atom_index)
        if not (
            0 <= left < molecule.GetNumAtoms()
            and 0 <= right < molecule.GetNumAtoms()
        ):
            continue
        if left == right or molecule.GetBondBetweenAtoms(left, right) is not None:
            continue
        try:
            path = Chem.GetShortestPath(molecule, left, right)
        except Exception:
            continue
        if path:
            distance = len(path) - 1
            tether_distances.append(distance)
            formed_ring_sizes.append(distance + 1)

    scopes = set(formed_scopes)
    if scopes == {"intramolecular"}:
        reaction_scope = "intramolecular"
    elif scopes == {"intermolecular"}:
        reaction_scope = "intermolecular"
    elif len(scopes) > 1:
        reaction_scope = "mixed"
    elif len(participating) == 1:
        reaction_scope = "unimolecular"
    else:
        reaction_scope = "unresolved"
    return ReactionTopology(
        reaction_scope=reaction_scope,
        participating_component_indices=tuple(sorted(participating)),
        role_component_indices={},
        same_component_role_groups=(),
        formed_bond_scopes=tuple(sorted(formed_scopes)),
        reactant_tether_distances=tuple(sorted(tether_distances)),
        formed_ring_sizes=tuple(sorted(formed_ring_sizes)),
        ring_count_delta=None,
        evidence=evidence,
        confidence=confidence,
    )


def _event_sites(
    reactants: Tuple[ReactionComponent, ...], edits: Sequence[ReactionEdit]
) -> Tuple[str, ...]:
    atoms_by_component: dict[int, set[int]] = defaultdict(set)
    for edit in edits:
        for atom in (edit.atom_1, edit.atom_2):
            if atom is not None and atom.side == "reactant":
                atoms_by_component[atom.component_index].add(atom.atom_index)
    return tuple(
        sorted(
            f"r{component.component_index}:{site.site_id}"
            for component in reactants
            for site in component.compound_analysis.sites
            if atoms_by_component.get(component.component_index, set()).intersection(
                site.atom_indices
            )
        )
    )


def _candidate_atom_keys(
    candidate: ReactionCandidate,
) -> set[Tuple[int, int]]:
    keys: set[Tuple[int, int]] = set()
    for change in candidate.predicted_bond_changes:
        for role_path in (change.atom_1_role, change.atom_2_role):
            if role_path is None or "." not in role_path:
                continue
            partner_role, atom_role = role_path.split(".", 1)
            site = candidate.role_assignments.get(partner_role)
            if site is None:
                continue
            atom_indices = site.atom_roles.get(atom_role) or ()
            keys.update(
                (site.component_index, int(atom_index))
                for atom_index in atom_indices
            )
    if keys:
        return keys
    return {
        (site.component_index, int(atom_index))
        for site in candidate.role_assignments.values()
        for atom_indices in site.atom_roles.values()
        for atom_index in atom_indices
    }


def _event_relations(
    events: Sequence[ReactionEvent],
) -> Tuple[ReactionEventRelation, ...]:
    relations = []
    for index, left in enumerate(events):
        left_components = set(left.topology.participating_component_indices)
        left_maps = {
            atom.atom_map_number
            for edit in left.edits
            for atom in (edit.atom_1, edit.atom_2)
            if atom is not None and atom.atom_map_number is not None
        }
        for right in events[index + 1 :]:
            shared_components = tuple(
                sorted(left_components.intersection(
                    right.topology.participating_component_indices
                ))
            )
            right_maps = {
                atom.atom_map_number
                for edit in right.edits
                for atom in (edit.atom_1, edit.atom_2)
                if atom is not None and atom.atom_map_number is not None
            }
            shared_maps = tuple(sorted(left_maps.intersection(right_maps)))
            relation_type = (
                "shared_atom"
                if shared_maps
                else "shared_component"
                if shared_components
                else "independent_sites"
            )
            relations.append(
                ReactionEventRelation(
                    event_id_1=left.event_id,
                    event_id_2=right.event_id,
                    relation_type=relation_type,
                    shared_component_indices=shared_components,
                    shared_atom_map_numbers=shared_maps,
                    evidence="atom_provenance",
                )
            )
    return tuple(relations)


def build_reaction_events(
    *,
    reactants: Tuple[ReactionComponent, ...],
    edits: Sequence[ReactionEdit],
    partners: Sequence[ReactionPartner],
    selected: Optional[ReactionCandidate],
    selected_events: Sequence[ReactionCandidate],
    named_family: Optional[str],
    compatible_named_families: Tuple[str, ...],
    evidence: str,
    confidence: float,
) -> Tuple[Tuple[ReactionEvent, ...], Tuple[ReactionEventRelation, ...]]:
    """Build deterministic reaction events and their structural relations."""
    groups = partition_reaction_edits(edits)
    events = []
    for ordinal, group in enumerate(groups, start=1):
        chemistry = tuple(sorted(_edit_chemistry(edit) for edit in group))
        event_key = _digest(
            "E0", {"edits": chemistry, "event_schema_version": "1.2"}
        )
        event_id = f"{_digest('RE1', event_key)}:{ordinal}"
        component_indices = {
            atom.component_index
            for edit in group
            for atom in (edit.atom_1, edit.atom_2)
            if atom is not None and atom.side == "reactant"
        }
        event_partners = tuple(
            sorted(
                partner.partner_id
                for partner in partners
                if partner.component_index in component_indices
            )
        )
        formed = tuple(
            sorted(
                _bond_type(edit, edit.new_order)
                for edit in group
                if edit.edit_type == "formed"
            )
        )
        broken = tuple(
            sorted(
                _bond_type(edit, edit.old_order)
                for edit in group
                if edit.edit_type == "broken"
            )
        )
        order_changes = tuple(
            sorted(
                f"{_bond_type(edit, edit.old_order)}>{edit.new_order or 'NONE'}"
                for edit in group
                if edit.edit_type == "order_changed"
            )
        )
        hydrogen_changes = tuple(
            sorted(
                f"{_bond_type(edit, edit.old_order)}>{edit.new_order or 'NONE'}"
                for edit in group
                if edit.edit_type == "hydrogen_change"
            )
        )
        pattern = match_reaction_label_pattern(group, style="ascii")
        single_event = len(groups) == 1
        event_site_ids = _event_sites(reactants, group)
        event_atom_keys = {
            (atom.component_index, atom.atom_index)
            for edit in group
            for atom in (edit.atom_1, edit.atom_2)
            if atom is not None and atom.side == "reactant"
        }
        interpreted = tuple(
            candidate
            for candidate in selected_events
            if _candidate_atom_keys(candidate).issubset(event_atom_keys)
        )
        event_candidate = interpreted[0] if len(interpreted) == 1 else None
        event_families = (
            tuple(sorted(set(event_candidate.compatible_named_families)))
            if event_candidate is not None
            else ()
        )
        event_family = event_families[0] if len(event_families) == 1 else None
        declared_archetype = (
            event_candidate.edit_archetype
            if event_candidate is not None
            else selected.edit_archetype
            if single_event and selected is not None
            else None
        )
        edit_archetype, archetype_warnings = reconcile_edit_archetype(
            group, declared_archetype
        )
        events.append(
            ReactionEvent(
                event_id=event_id,
                event_signature_key=event_key,
                edits=group,
                partner_ids=event_partners,
                reactive_site_ids=event_site_ids,
                formed_bond_types=formed,
                broken_bond_types=broken,
                order_changes=order_changes,
                hydrogen_changes=hydrogen_changes,
                formed_connection_labels=formed,
                topology=_event_topology(
                    reactants,
                    group,
                    evidence=evidence,
                    confidence=confidence,
                ),
                edit_archetype=edit_archetype,
                transformation_class=(
                    event_candidate.transformation_class
                    if event_candidate is not None
                    else selected.transformation_class
                    if single_event and selected is not None
                    else pattern.pattern_id if pattern is not None else None
                ),
                transformation_confidence=confidence,
                named_family=(
                    event_family
                    if event_candidate is not None
                    else named_family if single_event else None
                ),
                family_confidence=(
                    1.0
                    if event_family or (single_event and named_family)
                    else 0.0
                ),
                compatible_named_families=(
                    event_families
                    if event_candidate is not None
                    else tuple(sorted(set(compatible_named_families)))
                    if single_event
                    else ()
                ),
                evidence=evidence,
                confidence=confidence,
                warnings=archetype_warnings,
            )
        )
    event_tuple = tuple(events)
    return event_tuple, _event_relations(event_tuple)


__all__ = ["build_reaction_events", "partition_reaction_edits"]
