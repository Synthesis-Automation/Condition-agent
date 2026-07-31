"""Template-free minimization of an observed reaction edit graph.

The projection is an observation contract.  It consumes normalized edits and
parsed molecular graphs, but never reaction grammars, templates, source labels,
or named families.  Every edit-participating atom remains in the active graph;
unchanged connected components are represented once as typed remote subgraphs
with one or more attachment ports.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

from ..chemistry.rdkit_utils import parse_smiles
from ..reaction_models import (
    ReactionComponent,
    ReactionEdit,
    ReactionStereoChange,
    ReactionTopology,
)
from .annotations import (
    decarboxylative_abstraction as _decarboxylative_abstraction,
)
from .common import (
    AtomIdentity as _AtomIdentity,
    Coordinate as _Coordinate,
    EditRecord as _EditRecord,
    Location as _Location,
    atom_identity as _atom_identity,
    atom_map_number as _atom_map_number,
)
from .keys import digest as _digest
from .models import (
    REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
    REACTION_CORE_PROJECTION_SCHEMA_VERSION,
    ReactionCoreAttachmentPort,
    ReactionCoreAtomState,
    ReactionCoreAtomTransition,
    ReactionCoreEvent,
    ReactionCoreProjection,
    ReactionCoreRemoteClass,
    ReactionCoreRemoteSubgraph,
)
from .rendering import (
    build_core_presentation,
    multi_center_edit_label as _multi_center_edit_label,
    single_center_transition_label as _single_center_transition_label,
    state_label as _state_label,
)
from .quality import assess_reaction_core_quality, validate_core_edits
from .remote import (
    _build_remote_subgraphs_for_side,
    _functional_group_ids,
    _with_remote_continuity,
)
from .state_changes import build_state_changes


def _component_locations(
    components: Sequence[ReactionComponent],
    *,
    side: str,
    atom_map_overrides: Optional[Mapping[tuple[str, int, int], int]] = None,
) -> tuple[Dict[int, _Location], Dict[_Coordinate, _Location]]:
    by_map: Dict[int, _Location] = {}
    by_coordinate: Dict[_Coordinate, _Location] = {}
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            atom_index = int(atom.GetIdx())
            location = (component, molecule, atom_index)
            by_coordinate[(component.component_index, atom_index)] = location
            map_number = _atom_map_number(
                atom,
                side=side,
                component_index=component.component_index,
                atom_index=atom_index,
                atom_map_overrides=atom_map_overrides,
            )
            if map_number > 0:
                by_map[map_number] = location
    return by_map, by_coordinate


def _location_for_identity(
    identity: _AtomIdentity,
    *,
    side: str,
    by_map: Mapping[int, _Location],
    by_coordinate: Mapping[_Coordinate, _Location],
) -> Optional[_Location]:
    if identity[0] == "map":
        return by_map.get(int(identity[1]))
    if identity[0] != "coordinate" or str(identity[1]) != side:
        return None
    return by_coordinate.get((int(identity[2]), int(identity[3])))


def _active_coordinates(
    identities: Iterable[_AtomIdentity],
    *,
    side: str,
    by_map: Mapping[int, _Location],
    by_coordinate: Mapping[_Coordinate, _Location],
) -> set[_Coordinate]:
    values = set()
    for identity in identities:
        location = _location_for_identity(
            identity,
            side=side,
            by_map=by_map,
            by_coordinate=by_coordinate,
        )
        if location is not None:
            component, _, atom_index = location
            values.add((component.component_index, atom_index))
    return values


def _participant_tokens(
    components: Sequence[ReactionComponent],
    active_coordinates: set[_Coordinate],
    *,
    side: str,
) -> Tuple[str, ...]:
    """Return mapping-number-independent typed handles touching active atoms."""
    values = []
    for component in components:
        active = {
            atom_index
            for component_index, atom_index in active_coordinates
            if component_index == component.component_index
        }
        if not active:
            continue
        for group in component.compound_analysis.functional_groups:
            if active.intersection(int(value) for value in group.atom_indices):
                values.append(f"{side}:group:{group.group_id}")
        for site in component.compound_analysis.sites:
            if active.intersection(int(value) for value in site.atom_indices):
                values.append(f"{side}:site:{site.site_type}")
    return tuple(sorted(set(values)))


def _neighbor_token(molecule: Any, center_index: int, neighbor_index: int) -> str:
    neighbor = molecule.GetAtomWithIdx(int(neighbor_index))
    bond = molecule.GetBondBetweenAtoms(int(center_index), int(neighbor_index))
    external_heavy_degree = sum(
        other.GetAtomicNum() > 1 and int(other.GetIdx()) != int(center_index)
        for other in neighbor.GetNeighbors()
    )
    return "|".join(
        (
            neighbor.GetSymbol(),
            "aromatic" if neighbor.GetIsAromatic() else "aliphatic",
            str(bond.GetBondType()).upper(),
            f"Q{int(neighbor.GetFormalCharge())}",
            f"H{int(neighbor.GetTotalNumHs(includeNeighbors=True))}",
            f"X{external_heavy_degree}",
            f"R{int(neighbor.IsInRing())}",
        )
    )


def _build_atom_state(
    *,
    side: str,
    location: _Location,
    active_coordinates: set[_Coordinate],
    remote_classes: Mapping[
        tuple[str, int, int, int],
        ReactionCoreRemoteClass,
    ],
    atom_map_number: Optional[int] = None,
) -> ReactionCoreAtomState:
    component, molecule, atom_index = location
    atom = molecule.GetAtomWithIdx(atom_index)
    neighbor_tokens = tuple(
        sorted(
            _neighbor_token(molecule, atom_index, int(neighbor.GetIdx()))
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() > 1
        )
    )
    state_payload = {
        "element": atom.GetSymbol(),
        "formal_charge": int(atom.GetFormalCharge()),
        "aromatic": bool(atom.GetIsAromatic()),
        "hybridization": str(atom.GetHybridization()),
        "total_hydrogens": int(atom.GetTotalNumHs(includeNeighbors=True)),
        "heavy_atom_degree": sum(
            neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors()
        ),
        "radical_electrons": int(atom.GetNumRadicalElectrons()),
        "isotope": int(atom.GetIsotope()),
        "neighbor_tokens": neighbor_tokens,
    }
    return ReactionCoreAtomState(
        side=side,  # type: ignore[arg-type]
        component_index=component.component_index,
        atom_index=atom_index,
        atom_map_number=atom_map_number or int(atom.GetAtomMapNum()) or None,
        element=atom.GetSymbol(),
        formal_charge=int(atom.GetFormalCharge()),
        aromatic=bool(atom.GetIsAromatic()),
        hybridization=str(atom.GetHybridization()),
        total_hydrogens=state_payload["total_hydrogens"],
        heavy_atom_degree=state_payload["heavy_atom_degree"],
        radical_electrons=state_payload["radical_electrons"],
        isotope=state_payload["isotope"],
        neighbor_tokens=neighbor_tokens,
        functional_group_ids=_functional_group_ids(component, (atom_index,)),
        concise_label=_state_label(
            molecule=molecule,
            component_index=component.component_index,
            atom_index=atom_index,
            active_coordinates=active_coordinates,
            remote_classes=remote_classes,
            side=side,
        ),
        state_key=_digest("RAS2", state_payload, length=24),
    )


def _edit_graph(
    edits: Sequence[ReactionEdit],
) -> tuple[
    set[_AtomIdentity],
    Dict[_AtomIdentity, set[_AtomIdentity]],
    Counter[_AtomIdentity],
    Counter[_AtomIdentity],
    tuple[object, ...],
    Tuple[str, ...],
    Tuple[_EditRecord, ...],
]:
    identities: set[_AtomIdentity] = set()
    adjacency: Dict[_AtomIdentity, set[_AtomIdentity]] = defaultdict(set)
    incidence: Counter[_AtomIdentity] = Counter()
    heavy_incidence: Counter[_AtomIdentity] = Counter()
    atom_labels: Dict[_AtomIdentity, tuple[object, ...]] = {}
    incidents: Dict[_AtomIdentity, list[tuple[object, ...]]] = defaultdict(list)
    public_tokens = []
    records = []
    for edit in edits:
        left = _atom_identity(edit.atom_1)
        identities.add(left)
        incidence[left] += 1
        atom_labels[left] = (
            edit.atom_1.element,
            edit.atom_1.formal_charge,
            edit.atom_1.aromatic,
            edit.atom_1.hybridization,
        )
        record_identities = [left]
        if edit.atom_2 is None:
            neighbor_label: tuple[object, ...] = ("H",)
            incidents[left].append(
                (
                    edit.edit_type,
                    neighbor_label,
                    edit.old_order or "NONE",
                    edit.new_order or "NONE",
                )
            )
            pair = f"{edit.atom_1.element}-H"
        else:
            right = _atom_identity(edit.atom_2)
            identities.add(right)
            incidence[right] += 1
            heavy_incidence[left] += 1
            heavy_incidence[right] += 1
            record_identities.append(right)
            adjacency[left].add(right)
            adjacency[right].add(left)
            atom_labels[right] = (
                edit.atom_2.element,
                edit.atom_2.formal_charge,
                edit.atom_2.aromatic,
                edit.atom_2.hybridization,
            )
            left_label = atom_labels[left]
            right_label = atom_labels[right]
            incidents[left].append(
                (
                    edit.edit_type,
                    right_label,
                    edit.old_order or "NONE",
                    edit.new_order or "NONE",
                )
            )
            incidents[right].append(
                (
                    edit.edit_type,
                    left_label,
                    edit.old_order or "NONE",
                    edit.new_order or "NONE",
                )
            )
            pair = "-".join(sorted((edit.atom_1.element, edit.atom_2.element)))
        token = (
            f"{edit.edit_type}:{pair}:"
            f"{edit.old_order or 'NONE'}>{edit.new_order or 'NONE'}"
        )
        public_tokens.append(token)
        records.append((token, tuple(record_identities)))
    graph_payload = tuple(
        sorted(
            (
                atom_labels[identity],
                tuple(sorted(incidents[identity], key=repr)),
            )
            for identity in identities
        )
    )
    return (
        identities,
        adjacency,
        incidence,
        heavy_incidence,
        graph_payload,
        tuple(sorted(public_tokens)),
        tuple(records),
    )


def _event_components(
    identities: set[_AtomIdentity],
    adjacency: Mapping[_AtomIdentity, set[_AtomIdentity]],
) -> Tuple[Tuple[_AtomIdentity, ...], ...]:
    remaining = set(identities)
    values = []
    while remaining:
        start = min(remaining, key=repr)
        queue = deque((start,))
        event = set()
        while queue:
            identity = queue.popleft()
            if identity in event:
                continue
            event.add(identity)
            queue.extend(adjacency.get(identity, ()))
        remaining.difference_update(event)
        values.append(tuple(sorted(event, key=repr)))
    return tuple(values)


def _stable_remote_count(
    identity: _AtomIdentity,
    *,
    reactant_by_map: Mapping[int, _Location],
    product_by_map: Mapping[int, _Location],
    active_coordinates_by_side: Mapping[str, set[_Coordinate]],
    atom_map_overrides: Optional[
        Mapping[tuple[str, int, int], int]
    ] = None,
) -> int:
    if identity[0] != "map":
        return 0
    map_number = int(identity[1])
    before = reactant_by_map.get(map_number)
    after = product_by_map.get(map_number)
    if before is None or after is None:
        return 0

    def tokens(location: _Location, side: str) -> set[tuple[int, str]]:
        component, molecule, atom_index = location
        values = set()
        for neighbor in molecule.GetAtomWithIdx(atom_index).GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if (
                component.component_index,
                neighbor_index,
            ) in active_coordinates_by_side[side]:
                continue
            neighbor_map = _atom_map_number(
                neighbor,
                side=side,
                component_index=component.component_index,
                atom_index=neighbor_index,
                atom_map_overrides=atom_map_overrides,
            )
            if neighbor_map <= 0:
                continue
            values.add(
                (
                    neighbor_map,
                    str(
                        molecule.GetBondBetweenAtoms(
                            atom_index,
                            neighbor_index,
                        ).GetBondType()
                    ).upper(),
                )
            )
        return values

    return len(tokens(before, "reactant").intersection(tokens(after, "product")))


def _state_identity(
    state: Optional[ReactionCoreAtomState],
    *,
    generic: bool,
) -> object:
    if state is None:
        return ("ABSENT",)
    if generic:
        neighbor_tokens = tuple(
            sorted("|".join(token.split("|")[:3]) for token in state.neighbor_tokens)
        )
        return (
            state.element,
            state.formal_charge,
            state.aromatic,
            state.total_hydrogens,
            state.radical_electrons,
            state.isotope,
            neighbor_tokens,
        )
    return (
        state.element,
        state.formal_charge,
        state.aromatic,
        state.hybridization,
        state.total_hydrogens,
        state.heavy_atom_degree,
        state.radical_electrons,
        state.isotope,
        state.neighbor_tokens,
    )


def _port_identity(port: ReactionCoreAttachmentPort) -> object:
    return (
        port.side,
        port.attachment_element,
        port.bond_order,
    )


def _remote_identity(
    subgraph: ReactionCoreRemoteSubgraph,
    *,
    exact: bool,
) -> object:
    base = (
        subgraph.side,
        subgraph.remote_class,
        subgraph.continuity,
        tuple(sorted((_port_identity(port) for port in subgraph.attachment_ports))),
    )
    if not exact:
        return base
    return base + (
        subgraph.fragment_smiles,
        subgraph.fragment_heavy_atom_count,
        subgraph.fragment_heteroatom_count,
        subgraph.fragment_aromatic_atom_count,
        subgraph.functional_group_ids,
    )


def _shape_remote_identity(
    subgraph: ReactionCoreRemoteSubgraph,
) -> object:
    return (
        subgraph.remote_class,
        len(subgraph.attachment_ports),
        tuple(sorted(port.attachment_element for port in subgraph.attachment_ports)),
        tuple(sorted(port.bond_order for port in subgraph.attachment_ports)),
        subgraph.functional_group_ids,
    )


def _evidence_status(evidence: str) -> str:
    if str(evidence).startswith("external"):
        return "external"
    if evidence in {
        "validated_atom_mapping",
        "validated_mapping_and_exact_reconstruction",
        "validated_mapping_and_exact_multi_event_reconstruction",
        "exact_product_reconstruction",
        "exact_multi_event_reconstruction",
    }:
        return "verified"
    if evidence in {
        "unique_scaffold_correspondence",
        "global_atom_correspondence",
        "fragmented_scaffold_correspondence",
        "partial_product_correspondence",
    }:
        return "inferred"
    return "hypothesis"


def build_reaction_core_projection(
    *,
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    edits: Sequence[ReactionEdit],
    stereo_changes: Sequence[ReactionStereoChange] = (),
    evidence: str,
    confidence: float,
    topology: Optional[ReactionTopology] = None,
    atom_map_overrides: Optional[
        Mapping[tuple[str, int, int], int]
    ] = None,
) -> Optional[ReactionCoreProjection]:
    """Build a grammar- and template-free minimized reaction-core projection."""
    if not edits:
        return None
    (
        identities,
        adjacency,
        incidence,
        heavy_incidence,
        edit_graph_payload,
        edit_tokens,
        edit_records,
    ) = _edit_graph(edits)
    reactant_by_map, reactant_by_coordinate = _component_locations(
        reactants,
        side="reactant",
        atom_map_overrides=atom_map_overrides,
    )
    product_by_map, product_by_coordinate = _component_locations(
        products,
        side="product",
        atom_map_overrides=atom_map_overrides,
    )
    shared_identities = {
        identity
        for identity in identities
        if identity[0] == "map"
        and int(identity[1]) in reactant_by_map
        and int(identity[1]) in product_by_map
    }
    if not shared_identities:
        return None
    reactant_active = _active_coordinates(
        identities,
        side="reactant",
        by_map=reactant_by_map,
        by_coordinate=reactant_by_coordinate,
    )
    product_active = _active_coordinates(
        identities,
        side="product",
        by_map=product_by_map,
        by_coordinate=product_by_coordinate,
    )
    active_coordinates_by_side = {
        "reactant": reactant_active,
        "product": product_active,
    }
    reactant_remote = _build_remote_subgraphs_for_side(
        side="reactant",
        components=reactants,
        active_coordinates=reactant_active,
        atom_map_overrides=atom_map_overrides,
    )
    product_remote = _build_remote_subgraphs_for_side(
        side="product",
        components=products,
        active_coordinates=product_active,
        atom_map_overrides=atom_map_overrides,
    )
    remote_subgraphs = _with_remote_continuity(
        reactant_remote,
        product_remote,
        reactant_by_map=reactant_by_map,
        product_by_map=product_by_map,
    )
    remote_classes = {
        (
            port.side,
            port.core_component_index,
            port.core_atom_index,
            port.attachment_atom_index,
        ): subgraph.remote_class
        for subgraph in remote_subgraphs
        for port in subgraph.attachment_ports
    }
    event_components = _event_components(identities, adjacency)
    primary_identities = set()
    for event in event_components:
        candidates = [
            identity for identity in event if identity in shared_identities
        ]
        if not candidates:
            continue
        selection_incidence = (
            heavy_incidence
            if any(heavy_incidence[identity] for identity in candidates)
            else incidence
        )
        maximum_incidence = max(
            selection_incidence[identity] for identity in candidates
        )
        candidates = [
            identity
            for identity in candidates
            if selection_incidence[identity] == maximum_incidence
        ]
        stable_counts = {
            identity: _stable_remote_count(
                identity,
                reactant_by_map=reactant_by_map,
                product_by_map=product_by_map,
                active_coordinates_by_side=active_coordinates_by_side,
                atom_map_overrides=atom_map_overrides,
            )
            for identity in candidates
        }
        maximum_stable = max(stable_counts.values(), default=0)
        primary_identities.update(
            identity
            for identity in candidates
            if stable_counts[identity] == maximum_stable
        )

    transition_by_identity: Dict[
        _AtomIdentity,
        ReactionCoreAtomTransition,
    ] = {}
    for identity in sorted(identities, key=repr):
        before_location = _location_for_identity(
            identity,
            side="reactant",
            by_map=reactant_by_map,
            by_coordinate=reactant_by_coordinate,
        )
        after_location = _location_for_identity(
            identity,
            side="product",
            by_map=product_by_map,
            by_coordinate=product_by_coordinate,
        )
        before_state = (
            _build_atom_state(
                side="reactant",
                location=before_location,
                active_coordinates=reactant_active,
                remote_classes=remote_classes,
                atom_map_number=(
                    int(identity[1]) if identity[0] == "map" else None
                ),
            )
            if before_location is not None
            else None
        )
        after_state = (
            _build_atom_state(
                side="product",
                location=after_location,
                active_coordinates=product_active,
                remote_classes=remote_classes,
                atom_map_number=(
                    int(identity[1]) if identity[0] == "map" else None
                ),
            )
            if after_location is not None
            else None
        )
        role = (
            "primary_center"
            if identity in primary_identities
            else "participant"
        )
        stable_count = _stable_remote_count(
            identity,
            reactant_by_map=reactant_by_map,
            product_by_map=product_by_map,
            active_coordinates_by_side=active_coordinates_by_side,
            atom_map_overrides=atom_map_overrides,
        )
        transition_payload = {
            "before": _state_identity(before_state, generic=False),
            "after": _state_identity(after_state, generic=False),
            "incident_edit_count": incidence[identity],
            "stable_remote_subgraph_count": stable_count,
            "role": role,
        }
        transition_by_identity[identity] = ReactionCoreAtomTransition(
            transition_id=_digest("RCA2", transition_payload, length=24),
            atom_map_number=(
                int(identity[1]) if identity[0] == "map" else None
            ),
            before_state=before_state,
            after_state=after_state,
            incident_edit_count=incidence[identity],
            stable_remote_subgraph_count=stable_count,
            role=role,  # type: ignore[arg-type]
        )
    transitions = tuple(
        sorted(
            transition_by_identity.values(),
            key=lambda transition: (
                transition.role,
                repr(_state_identity(transition.before_state, generic=False)),
                repr(_state_identity(transition.after_state, generic=False)),
                transition.transition_id,
            ),
        )
    )
    if not transitions or not primary_identities:
        return None
    state_changes = build_state_changes(
        transitions,
        stereo_changes,
        evidence=str(evidence),
    )
    state_change_payload = tuple(
        sorted(
            (
                change.change_type,
                change.elements,
                change.before_value,
                change.after_value,
            )
            for change in state_changes
        )
    )

    events = []
    for event in event_components:
        if not any(identity in shared_identities for identity in event):
            continue
        event_set = set(event)
        event_edit_tokens = tuple(
            sorted(
                token
                for token, record_identities in edit_records
                if set(record_identities).intersection(event_set)
            )
        )
        transition_ids = tuple(
            sorted(transition_by_identity[identity].transition_id for identity in event)
        )
        events.append(
            ReactionCoreEvent(
                event_id=_digest(
                    "RCE2",
                    {
                        "transitions": tuple(
                            sorted(
                                (
                                    _state_identity(
                                        transition_by_identity[identity].before_state,
                                        generic=False,
                                    ),
                                    _state_identity(
                                        transition_by_identity[identity].after_state,
                                        generic=False,
                                    ),
                                )
                                for identity in event
                            )
                        ),
                        "edits": event_edit_tokens,
                    },
                    length=24,
                ),
                transition_ids=transition_ids,
                edit_tokens=event_edit_tokens,
            )
        )
    events = sorted(events, key=lambda event: (event.edit_tokens, event.event_id))

    primary = tuple(
        transition
        for transition in transitions
        if transition.role == "primary_center"
    )
    center_transition_payload = tuple(
        sorted(
            (
                _state_identity(transition.before_state, generic=False),
                _state_identity(transition.after_state, generic=False),
            )
            for transition in primary
        )
    )
    generic_center_payload = tuple(
        sorted(
            (
                _state_identity(transition.before_state, generic=True),
                _state_identity(transition.after_state, generic=True),
            )
            for transition in primary
        )
    )
    all_transition_payload = tuple(
        sorted(
            (
                transition.role,
                _state_identity(transition.before_state, generic=False),
                _state_identity(transition.after_state, generic=False),
                transition.incident_edit_count,
            )
            for transition in transitions
        )
    )
    participant_tokens = tuple(
        sorted(
            (
                *_participant_tokens(
                    reactants,
                    reactant_active,
                    side="reactant",
                ),
                *_participant_tokens(
                    products,
                    product_active,
                    side="product",
                ),
            )
        )
    )
    retained_shape = tuple(
        sorted(
            {
                _shape_remote_identity(subgraph)
                for subgraph in remote_subgraphs
                if subgraph.continuity == "retained"
            },
            key=repr,
        )
    )
    exact_identity_payload = {
        "edit_graph": edit_graph_payload,
        "transitions": all_transition_payload,
        "state_changes": state_change_payload,
        "remote_subgraphs": tuple(
            sorted(
                (
                    _remote_identity(subgraph, exact=True)
                    for subgraph in remote_subgraphs
                ),
                key=repr,
            )
        ),
        "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
    }
    typed_identity_payload = {
        "edit_graph": edit_graph_payload,
        "transitions": all_transition_payload,
        "state_changes": state_change_payload,
        "remote_subgraphs": tuple(
            sorted(
                (
                    _remote_identity(subgraph, exact=False)
                    for subgraph in remote_subgraphs
                ),
                key=repr,
            )
        ),
        "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
    }
    exact_key = _digest("RCX2", exact_identity_payload)
    typed_key = _digest("RCT2", typed_identity_payload)
    mapping_equivalence_key = _digest("RME1", exact_identity_payload)
    shape_key = _digest(
        "RSH2",
        {
            "primary_centers": generic_center_payload,
            "participant_tokens": participant_tokens,
            "retained_remote_subgraphs": retained_shape,
            "event_count": len(events),
            "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        },
    )
    center_transition_key = _digest(
        "RCS2",
        {
            "centers": center_transition_payload,
            "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        },
    )
    generic_label = (
        _multi_center_edit_label(edits)
        if len(primary) > 1
        else _single_center_transition_label(
            primary_identity=next(iter(primary_identities)),
            transition_by_identity=transition_by_identity,
            edits=edits,
        )
    )
    abstraction = _decarboxylative_abstraction(
        edits=edits,
        reactant_by_map=reactant_by_map,
        product_by_map=product_by_map,
        topology=topology,
    )
    checked_edit_count, consistency_issues = validate_core_edits(
        edits,
        reactant_by_map=reactant_by_map,
        product_by_map=product_by_map,
    )
    remote_continuity_unresolved = any(
        subgraph.continuity == "unresolved"
        for subgraph in remote_subgraphs
    )
    no_op_primary_center = any(
        _state_identity(transition.before_state, generic=False)
        == _state_identity(transition.after_state, generic=False)
        for transition in primary
    )
    quality = assess_reaction_core_quality(
        active_atom_count=len(identities),
        mapped_active_atom_count=sum(
            identity[0] == "map" for identity in identities
        ),
        edit_count=len(edits),
        heavy_atom_edit_count=sum(edit.atom_2 is not None for edit in edits),
        checked_edit_count=checked_edit_count,
        consistency_issues=consistency_issues,
        event_count=len(events),
        remote_continuity_unresolved=remote_continuity_unresolved,
        no_op_primary_center=no_op_primary_center,
    )
    warnings = set()
    if remote_continuity_unresolved:
        warnings.add("REACTION_CORE_REMOTE_CONTINUITY_UNRESOLVED")
    if any(identity[0] != "map" for identity in identities):
        warnings.add("REACTION_CORE_PARTIAL_ATOM_PROVENANCE")
    if str(evidence).startswith("external"):
        warnings.add("REACTION_CORE_EXTERNAL_MAPPING_PROPOSAL")
    if no_op_primary_center:
        warnings.add("REACTION_CORE_NO_OP_PRIMARY_CENTER")
    if quality.status == "review":
        warnings.add("REACTION_CORE_QUALITY_REVIEW")
    elif quality.status == "blocked":
        warnings.add("REACTION_CORE_QUALITY_BLOCKED")
    presentation = build_core_presentation(
        equation=generic_label,
        edits=edits,
        state_changes=state_changes,
        remote_subgraphs=remote_subgraphs,
        evidence_status=_evidence_status(str(evidence)),
        quality=quality,
    )
    core_id = _digest(
        "RCP2",
        {
            "keys": (
                exact_key,
                typed_key,
                shape_key,
                center_transition_key,
                mapping_equivalence_key,
            ),
            "algorithm_version": REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
            "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
            "abstraction": (
                abstraction.motif_key if abstraction is not None else None
            ),
        },
        length=64,
    )
    return ReactionCoreProjection(
        core_id=core_id,
        exact_core_key=exact_key,
        typed_core_key=typed_key,
        shape_core_key=shape_key,
        center_transition_key=center_transition_key,
        mapping_equivalence_key=mapping_equivalence_key,
        atom_transitions=transitions,
        state_changes=state_changes,
        events=tuple(events),
        remote_subgraphs=remote_subgraphs,
        edit_tokens=edit_tokens,
        participant_tokens=participant_tokens,
        generic_label=generic_label,
        presentation=presentation,
        quality=quality,
        abstraction=abstraction,
        active_atom_count=len(identities),
        event_count=len(events),
        evidence=str(evidence),
        evidence_status=_evidence_status(str(evidence)),  # type: ignore[arg-type]
        confidence=float(confidence),
        warnings=tuple(sorted(warnings)),
    )


__all__ = ["build_reaction_core_projection"]
