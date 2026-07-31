"""Template-free minimization of an observed reaction edit graph.

The projection is an observation contract.  It consumes normalized edits and
parsed molecular graphs, but never reaction grammars, templates, source labels,
or named families.  Every edit-participating atom remains in the active graph;
unchanged connected components are represented once as typed remote subgraphs
with one or more attachment ports.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import replace
import hashlib
import json
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import (
    REACTION_CORE_PROJECTION_SCHEMA_VERSION,
    REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
    ReactionAtomReference,
    ReactionComponent,
    ReactionCoreAttachmentPort,
    ReactionCoreAtomState,
    ReactionCoreAtomTransition,
    ReactionCoreEvent,
    ReactionCoreProjection,
    ReactionCoreRemoteClass,
    ReactionCoreRemoteSubgraph,
    ReactionEdit,
)


_AtomIdentity = Tuple[object, ...]
_Coordinate = Tuple[int, int]
_Location = Tuple[ReactionComponent, Any, int]
_EditRecord = Tuple[str, Tuple[_AtomIdentity, ...]]


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    )


def _digest(prefix: str, value: Any, *, length: int = 32) -> str:
    encoded = _canonical_json(value).encode("utf-8")
    return f"{prefix}:" + hashlib.sha256(encoded).hexdigest()[:length]


def _atom_identity(reference: ReactionAtomReference) -> _AtomIdentity:
    if reference.atom_map_number is not None:
        return ("map", int(reference.atom_map_number))
    return (
        "coordinate",
        str(reference.side),
        int(reference.component_index),
        int(reference.atom_index),
        str(reference.element),
    )


def _component_locations(
    components: Sequence[ReactionComponent],
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
            map_number = int(atom.GetAtomMapNum())
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


def _functional_group_ids(
    component: ReactionComponent,
    atom_indices: Iterable[int],
) -> Tuple[str, ...]:
    selected = {int(value) for value in atom_indices}
    return tuple(
        sorted(
            {
                str(group.group_id)
                for group in component.compound_analysis.functional_groups
                if selected.intersection(int(value) for value in group.atom_indices)
            }
        )
    )


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
                                int(core_atom.GetAtomMapNum()) or None
                            ),
                            attachment_atom_index=attachment_index,
                            attachment_atom_map_number=(
                                int(attachment.GetAtomMapNum()) or None
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
                    int(atom.GetAtomMapNum())
                    for atom_index in atom_indices
                    if (atom := molecule.GetAtomWithIdx(atom_index)).GetAtomMapNum()
                    > 0
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
                    functional_group_ids=_functional_group_ids(
                        component,
                        atom_indices,
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


def _remote_display(remote_class: ReactionCoreRemoteClass) -> str:
    return {
        "aryl": "Ar",
        "heteroaryl": "HetAr",
        "alkyl": "R",
        "alkenyl": "Alkenyl",
        "alkynyl": "Alkynyl",
        "acyl": "Acyl",
        "ring_aliphatic": "Cycloalkyl",
        "heteroatom": "X",
        "generic_R": "R",
    }[remote_class]


def _bond_prefix(order: str) -> str:
    return {
        "SINGLE": "",
        "DOUBLE": "=",
        "TRIPLE": "#",
        "AROMATIC": ":",
    }.get(str(order).upper(), "~")


def _active_neighbor_display(
    molecule: Any,
    center_index: int,
    neighbor_index: int,
) -> str:
    neighbor = molecule.GetAtomWithIdx(int(neighbor_index))
    symbol = neighbor.GetSymbol()
    if neighbor.GetIsAromatic():
        return "ArC" if symbol == "C" else symbol
    if symbol not in {"N", "O", "P", "S"}:
        return "R" if symbol == "C" else symbol
    external = [
        atom
        for atom in neighbor.GetNeighbors()
        if atom.GetAtomicNum() > 1 and int(atom.GetIdx()) != int(center_index)
    ]
    if external:
        if all(atom.GetAtomicNum() == 6 for atom in external):
            return f"{symbol}-R"
        suffix = ",".join(
            sorted(
                "R" if atom.GetAtomicNum() == 6 else atom.GetSymbol()
                for atom in external
            )
        )
        return f"{symbol}-({suffix})"
    if int(neighbor.GetTotalNumHs(includeNeighbors=True)) > 0:
        return f"{symbol}-H"
    return symbol


def _state_label(
    *,
    molecule: Any,
    component_index: int,
    atom_index: int,
    active_coordinates: set[_Coordinate],
    remote_classes: Mapping[
        tuple[str, int, int, int],
        ReactionCoreRemoteClass,
    ],
    side: str,
) -> str:
    atom = molecule.GetAtomWithIdx(int(atom_index))
    tokens = ["H"] * int(atom.GetTotalNumHs(includeNeighbors=True))
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() <= 1:
            continue
        neighbor_index = int(neighbor.GetIdx())
        order = str(
            molecule.GetBondBetweenAtoms(atom_index, neighbor_index).GetBondType()
        ).upper()
        if (component_index, neighbor_index) in active_coordinates:
            value = _active_neighbor_display(molecule, atom_index, neighbor_index)
        else:
            remote_class = remote_classes.get(
                (side, component_index, atom_index, neighbor_index),
                "generic_R",
            )
            value = (
                neighbor.GetSymbol()
                if remote_class == "heteroatom"
                else _remote_display(remote_class)
            )
        tokens.append(f"{_bond_prefix(order)}{value}")

    def token_order(token: str) -> tuple[int, str]:
        plain = token.lstrip("=#:~")
        if plain == "H":
            return 0, token
        if plain in {
            "Ar",
            "HetAr",
            "R",
            "Alkenyl",
            "Alkynyl",
            "Acyl",
        }:
            return 1, token
        if token.startswith(("=", "#", ":", "~")):
            return 3, token
        return 2, token

    counts = Counter(tokens)
    rendered = []
    for token in sorted(counts, key=token_order):
        count = counts[token]
        rendered.append(f"({token}){count if count > 1 else ''}")
    return f"{atom.GetSymbol()}{''.join(rendered)}"


def _build_atom_state(
    *,
    side: str,
    location: _Location,
    active_coordinates: set[_Coordinate],
    remote_classes: Mapping[
        tuple[str, int, int, int],
        ReactionCoreRemoteClass,
    ],
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
        "neighbor_tokens": neighbor_tokens,
    }
    return ReactionCoreAtomState(
        side=side,  # type: ignore[arg-type]
        component_index=component.component_index,
        atom_index=atom_index,
        atom_map_number=int(atom.GetAtomMapNum()) or None,
        element=atom.GetSymbol(),
        formal_charge=int(atom.GetFormalCharge()),
        aromatic=bool(atom.GetIsAromatic()),
        hybridization=str(atom.GetHybridization()),
        total_hydrogens=state_payload["total_hydrogens"],
        heavy_atom_degree=state_payload["heavy_atom_degree"],
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
            neighbor_map = int(neighbor.GetAtomMapNum())
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
            neighbor_tokens,
        )
    return (
        state.element,
        state.formal_charge,
        state.aromatic,
        state.hybridization,
        state.total_hydrogens,
        state.heavy_atom_degree,
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
    evidence: str,
    confidence: float,
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
    reactant_by_map, reactant_by_coordinate = _component_locations(reactants)
    product_by_map, product_by_coordinate = _component_locations(products)
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
    )
    product_remote = _build_remote_subgraphs_for_side(
        side="product",
        components=products,
        active_coordinates=product_active,
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
    exact_key = _digest(
        "RCX2",
        {
            "edit_graph": edit_graph_payload,
            "transitions": all_transition_payload,
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
        },
    )
    typed_key = _digest(
        "RCT2",
        {
            "edit_graph": edit_graph_payload,
            "transitions": all_transition_payload,
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
        },
    )
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
    generic_label = " + ".join(
        (
            f"{transition.before_state.concise_label if transition.before_state else '∅'}"
            " → "
            f"{transition.after_state.concise_label if transition.after_state else '∅'}"
        )
        for transition in primary
    )
    warnings = set()
    if any(
        subgraph.continuity == "unresolved"
        for subgraph in remote_subgraphs
    ):
        warnings.add("REACTION_CORE_REMOTE_CONTINUITY_UNRESOLVED")
    if any(identity[0] != "map" for identity in identities):
        warnings.add("REACTION_CORE_PARTIAL_ATOM_PROVENANCE")
    if str(evidence).startswith("external"):
        warnings.add("REACTION_CORE_EXTERNAL_MAPPING_PROPOSAL")
    if any(
        _state_identity(transition.before_state, generic=False)
        == _state_identity(transition.after_state, generic=False)
        for transition in primary
    ):
        warnings.add("REACTION_CORE_NO_OP_PRIMARY_CENTER")
    core_id = _digest(
        "RCP2",
        {
            "keys": (
                exact_key,
                typed_key,
                shape_key,
                center_transition_key,
            ),
            "algorithm_version": REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
            "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        },
        length=64,
    )
    return ReactionCoreProjection(
        core_id=core_id,
        exact_core_key=exact_key,
        typed_core_key=typed_key,
        shape_core_key=shape_key,
        center_transition_key=center_transition_key,
        atom_transitions=transitions,
        events=tuple(events),
        remote_subgraphs=remote_subgraphs,
        edit_tokens=edit_tokens,
        participant_tokens=participant_tokens,
        generic_label=generic_label,
        active_atom_count=len(identities),
        event_count=len(events),
        evidence=str(evidence),
        evidence_status=_evidence_status(str(evidence)),  # type: ignore[arg-type]
        confidence=float(confidence),
        warnings=tuple(sorted(warnings)),
    )


__all__ = ["build_reaction_core_projection"]
