"""Template-free minimization of an observed reaction edit graph.

The projection is an observation contract.  It consumes normalized edits and
parsed molecular graphs, but never interpretation annotations, templates, source labels,
or named families.  Every edit-participating atom remains in the active graph;
unchanged connected components are represented once as typed remote subgraphs
with one or more attachment ports.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import replace
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

from ..chemistry.rdkit_utils import parse_smiles
from ..reaction_models import (
    ReactionComponent,
    ReactionEdit,
    ReactionStereoChange,
)
from ..reaction_edits import PROJECTED_UNMAPPED_DEPARTING_BOUNDARY_EVIDENCE
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
    ReactionCoreEventPath,
    ReactionCoreEventRelation,
    ReactionCoreProjection,
    ReactionCoreRemoteClass,
    ReactionCoreRemoteSubgraph,
)
from .quality import assess_reaction_core_quality, validate_core_edits
from .remote import (
    _build_remote_subgraphs_for_side,
    _with_remote_continuity,
)
from .state_changes import build_state_changes


_REACTION_CORE_IDENTITY_SCHEMA_VERSION = "3.0"


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
    """Return graph-only atom tokens touching the active core."""
    values = []
    for component in components:
        active = {
            atom_index
            for component_index, atom_index in active_coordinates
            if component_index == component.component_index
        }
        if not active:
            continue
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom_index in sorted(active):
            atom = molecule.GetAtomWithIdx(atom_index)
            values.append(
                "|".join(
                    (
                        side,
                        "atom",
                        atom.GetSymbol(),
                        f"Q{int(atom.GetFormalCharge())}",
                        "aromatic" if atom.GetIsAromatic() else "aliphatic",
                        str(atom.GetHybridization()).upper(),
                        f"H{int(atom.GetTotalNumHs(includeNeighbors=True))}",
                        f"D{int(atom.GetDegree())}",
                    )
                )
            )
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
        state_key=_digest("RAS2", state_payload, length=24),
    )


def _edit_graph(
    edits: Sequence[ReactionEdit],
    atom_map_overrides: Optional[
        Mapping[tuple[str, int, int], int]
    ] = None,
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
    for edit_index, edit in enumerate(edits):
        left_override = (atom_map_overrides or {}).get(
            (
                edit.atom_1.side,
                edit.atom_1.component_index,
                edit.atom_1.atom_index,
            )
        )
        left = _atom_identity(edit.atom_1)
        # Externally normalized edits already carry validated map numbers in
        # mapper coordinates.  Coordinate overrides describe the original
        # input molecules and must only supply identity for an unmapped edit;
        # otherwise an unrelated original atom at the same numeric index can
        # replace the edit's true mapped identity.
        if left[0] != "map" and left_override is not None:
            left = ("map", int(left_override))
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
            right_override = (atom_map_overrides or {}).get(
                (
                    edit.atom_2.side,
                    edit.atom_2.component_index,
                    edit.atom_2.atom_index,
                )
            )
            right = _atom_identity(edit.atom_2)
            if right[0] != "map" and right_override is not None:
                right = ("map", int(right_override))
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
        records.append((edit_index, token, tuple(record_identities)))
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


def _event_state_coordinates(
    event: ReactionCoreEvent,
    transitions_by_id: Mapping[str, ReactionCoreAtomTransition],
    *,
    side: str,
) -> Dict[int, set[int]]:
    values: Dict[int, set[int]] = defaultdict(set)
    for transition_id in event.transition_ids:
        transition = transitions_by_id[transition_id]
        state = (
            transition.before_state
            if side == "reactant"
            else transition.after_state
        )
        if state is not None:
            values[int(state.component_index)].add(int(state.atom_index))
    return values


def _shortest_event_path(
    components: Sequence[ReactionComponent],
    *,
    side: str,
    component_index: int,
    left_atom_indices: set[int],
    right_atom_indices: set[int],
) -> Optional[ReactionCoreEventPath]:
    from rdkit import Chem

    component = next(
        (
            value
            for value in components
            if int(value.component_index) == int(component_index)
        ),
        None,
    )
    molecule = parse_smiles(component.input_smiles) if component is not None else None
    if molecule is None:
        return None
    candidates = []
    for left in sorted(left_atom_indices):
        for right in sorted(right_atom_indices):
            if left == right:
                candidates.append((1, (int(left),), left, right))
                continue
            try:
                path = tuple(int(value) for value in Chem.GetShortestPath(
                    molecule, int(left), int(right)
                ))
            except Exception:
                continue
            if path:
                candidates.append((len(path), path, left, right))
    if not candidates:
        return None
    _, path, start, end = min(candidates)
    return ReactionCoreEventPath(
        side=side,  # type: ignore[arg-type]
        component_index=int(component_index),
        start_atom_index=int(start),
        end_atom_index=int(end),
        atom_indices=path,
        bond_count=len(path) - 1,
    )


def _build_event_relations(
    events: Sequence[ReactionCoreEvent],
    transitions: Sequence[ReactionCoreAtomTransition],
    *,
    reactants: Sequence[ReactionComponent],
    products: Sequence[ReactionComponent],
) -> Tuple[ReactionCoreEventRelation, ...]:
    transitions_by_id = {
        transition.transition_id: transition for transition in transitions
    }
    coordinates = {
        (event.event_id, side): _event_state_coordinates(
            event,
            transitions_by_id,
            side=side,
        )
        for event in events
        for side in ("reactant", "product")
    }
    relations = []
    for index, left in enumerate(events):
        for right in events[index + 1 :]:
            left_reactant = coordinates[(left.event_id, "reactant")]
            right_reactant = coordinates[(right.event_id, "reactant")]
            left_product = coordinates[(left.event_id, "product")]
            right_product = coordinates[(right.event_id, "product")]
            shared_reactant = tuple(
                sorted(set(left_reactant).intersection(right_reactant))
            )
            shared_product = tuple(
                sorted(set(left_product).intersection(right_product))
            )
            shared_transition_ids = set(left.transition_ids).intersection(
                right.transition_ids
            )
            shared_active_coordinates = any(
                left_reactant[component].intersection(right_reactant[component])
                for component in shared_reactant
            ) or any(
                left_product[component].intersection(right_product[component])
                for component in shared_product
            )
            paths = []
            for side, components, left_coordinates, right_coordinates, shared in (
                (
                    "reactant",
                    reactants,
                    left_reactant,
                    right_reactant,
                    shared_reactant,
                ),
                (
                    "product",
                    products,
                    left_product,
                    right_product,
                    shared_product,
                ),
            ):
                for component_index in shared:
                    path = _shortest_event_path(
                        components,
                        side=side,
                        component_index=component_index,
                        left_atom_indices=left_coordinates[component_index],
                        right_atom_indices=right_coordinates[component_index],
                    )
                    if path is not None:
                        paths.append(path)
            relation_type = (
                "shared_active_atom"
                if shared_transition_ids or shared_active_coordinates
                else "same_component"
                if shared_reactant or shared_product
                else "independent"
            )
            relations.append(
                ReactionCoreEventRelation(
                    event_id_1=left.event_id,
                    event_id_2=right.event_id,
                    relation_type=relation_type,  # type: ignore[arg-type]
                    shared_reactant_component_indices=shared_reactant,
                    shared_product_component_indices=shared_product,
                    shortest_paths=tuple(
                        sorted(
                            paths,
                            key=lambda value: (
                                value.side,
                                value.component_index,
                                value.bond_count,
                                value.atom_indices,
                            ),
                        )
                    ),
                    evidence="molecular_graph_shortest_path",
                )
            )
    return tuple(
        sorted(
            relations,
            key=lambda value: (value.event_id_1, value.event_id_2),
        )
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
    )


def _shape_remote_identity(
    subgraph: ReactionCoreRemoteSubgraph,
) -> object:
    return (
        subgraph.remote_class,
        len(subgraph.attachment_ports),
        tuple(sorted(port.attachment_element for port in subgraph.attachment_ports)),
        tuple(sorted(port.bond_order for port in subgraph.attachment_ports)),
    )


def _evidence_status(evidence: str) -> str:
    if str(evidence).startswith("external"):
        return "external"
    if evidence in {
        "validated_atom_mapping",
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
    atom_map_overrides: Optional[
        Mapping[tuple[str, int, int], int]
    ] = None,
) -> Optional[ReactionCoreProjection]:
    """Build an interpretation- and template-free reaction-core projection."""
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
    ) = _edit_graph(edits, atom_map_overrides=atom_map_overrides)
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
                for _, token, record_identities in edit_records
                if set(record_identities).issubset(event_set)
            )
        )
        event_edit_indices = tuple(
            sorted(
                edit_index
                for edit_index, _, record_identities in edit_records
                if set(record_identities).issubset(event_set)
            )
        )
        transition_ids = tuple(
            sorted(transition_by_identity[identity].transition_id for identity in event)
        )
        reactant_component_indices = tuple(
            sorted(
                {
                    int(state.component_index)
                    for identity in event
                    if (
                        state := transition_by_identity[identity].before_state
                    )
                    is not None
                }
            )
        )
        product_component_indices = tuple(
            sorted(
                {
                    int(state.component_index)
                    for identity in event
                    if (
                        state := transition_by_identity[identity].after_state
                    )
                    is not None
                }
            )
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
                edit_indices=event_edit_indices,
                reactant_component_indices=reactant_component_indices,
                product_component_indices=product_component_indices,
            )
        )
    events_by_base_id: Dict[str, list[ReactionCoreEvent]] = defaultdict(list)
    for event in events:
        events_by_base_id[event.event_id].append(event)
    events = [
        replace(event, event_id=f"{base_id}:{ordinal}")
        for base_id, group in sorted(events_by_base_id.items())
        for ordinal, event in enumerate(
            sorted(
                group,
                key=lambda value: (
                    value.reactant_component_indices,
                    value.product_component_indices,
                    value.transition_ids,
                    value.edit_indices,
                ),
            ),
            start=1,
        )
    ]
    events = sorted(events, key=lambda event: (event.edit_tokens, event.event_id))
    event_relations = _build_event_relations(
        events,
        transitions,
        reactants=reactants,
        products=products,
    )

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
        "schema_version": _REACTION_CORE_IDENTITY_SCHEMA_VERSION,
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
        "schema_version": _REACTION_CORE_IDENTITY_SCHEMA_VERSION,
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
            "schema_version": _REACTION_CORE_IDENTITY_SCHEMA_VERSION,
        },
    )
    center_transition_key = _digest(
        "RCS2",
        {
            "centers": center_transition_payload,
            "schema_version": _REACTION_CORE_IDENTITY_SCHEMA_VERSION,
        },
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
    unmapped_active_identities = {
        identity for identity in identities if identity[0] != "map"
    }
    validated_departing_identities = {
        _atom_identity(atom)
        for edit in edits
        if edit.evidence == PROJECTED_UNMAPPED_DEPARTING_BOUNDARY_EVIDENCE
        for atom in (edit.atom_1, edit.atom_2)
        if atom is not None and atom.atom_map_number is None
    }
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
        unmapped_active_atoms_are_validated_departures=bool(
            unmapped_active_identities
            and unmapped_active_identities <= validated_departing_identities
        ),
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
            "schema_version": _REACTION_CORE_IDENTITY_SCHEMA_VERSION,
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
        event_relations=event_relations,
        remote_subgraphs=remote_subgraphs,
        edit_tokens=edit_tokens,
        participant_tokens=participant_tokens,
        quality=quality,
        active_atom_count=len(identities),
        event_count=len(events),
        evidence=str(evidence),
        evidence_status=_evidence_status(str(evidence)),  # type: ignore[arg-type]
        confidence=float(confidence),
        warnings=tuple(sorted(warnings)),
    )


__all__ = ["build_reaction_core_projection"]
