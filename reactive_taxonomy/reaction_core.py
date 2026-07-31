"""Template-free reaction-center minimization from normalized graph edits.

The projection in this module is an observation layer.  It does not consult
reaction grammars, reaction templates, source labels, or named families.
Remote branches are cut only after their molecular graph has been inspected,
and both exact mapped-edit identity and correspondence-robust center-state
identity are retained.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
import hashlib
import json
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import (
    REACTION_CORE_PROJECTION_SCHEMA_VERSION,
    ReactionAtomReference,
    ReactionComponent,
    ReactionCoreAtomState,
    ReactionCoreBoundary,
    ReactionCoreBoundaryClass,
    ReactionCoreCenter,
    ReactionCoreProjection,
    ReactionEdit,
)


_AtomIdentity = Tuple[object, ...]
_Coordinate = Tuple[int, int]
_Location = Tuple[ReactionComponent, Any, int]


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
) -> tuple[
    Dict[int, _Location],
    Dict[_Coordinate, _Location],
]:
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


def _fragment_indices(
    molecule: Any,
    start_index: int,
    *,
    blocked: set[int],
) -> Tuple[int, ...]:
    queue = deque((int(start_index),))
    visited = set(blocked)
    values = []
    while queue:
        atom_index = queue.popleft()
        if atom_index in visited:
            continue
        visited.add(atom_index)
        atom = molecule.GetAtomWithIdx(atom_index)
        if atom.GetAtomicNum() <= 1:
            continue
        values.append(atom_index)
        for neighbor in atom.GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if neighbor.GetAtomicNum() > 1 and neighbor_index not in visited:
                queue.append(neighbor_index)
    return tuple(sorted(values))


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


def _aromatic_system_indices(
    molecule: Any,
    start_index: int,
    allowed: set[int],
) -> set[int]:
    values = set()
    queue = deque((int(start_index),))
    while queue:
        atom_index = queue.popleft()
        if atom_index in values or atom_index not in allowed:
            continue
        atom = molecule.GetAtomWithIdx(atom_index)
        if not atom.GetIsAromatic():
            continue
        values.add(atom_index)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                queue.append(int(neighbor.GetIdx()))
    return values


def _boundary_class(
    molecule: Any,
    attachment_atom_index: int,
    fragment_indices: Sequence[int],
) -> ReactionCoreBoundaryClass:
    atom = molecule.GetAtomWithIdx(int(attachment_atom_index))
    fragment_set = set(int(value) for value in fragment_indices)
    if atom.GetIsAromatic():
        aromatic_system = _aromatic_system_indices(
            molecule,
            int(attachment_atom_index),
            fragment_set,
        )
        if any(
            molecule.GetAtomWithIdx(atom_index).GetAtomicNum() != 6
            for atom_index in aromatic_system
        ):
            return "heteroaryl"
        return "aryl"
    if atom.GetAtomicNum() == 6:
        has_heteroatom_double_bond = any(
            int(neighbor.GetIdx()) in fragment_set
            and neighbor.GetAtomicNum() in {7, 8, 16}
            and str(
                molecule.GetBondBetweenAtoms(
                    int(attachment_atom_index),
                    int(neighbor.GetIdx()),
                ).GetBondType()
            ).upper()
            == "DOUBLE"
            for neighbor in atom.GetNeighbors()
        )
        if has_heteroatom_double_bond:
            return "acyl"
        hybridization = str(atom.GetHybridization()).upper()
        if hybridization == "SP":
            return "alkynyl"
        if hybridization == "SP2":
            return "alkenyl"
        if atom.IsInRing():
            return "ring_aliphatic"
        return "alkyl"
    if atom.GetAtomicNum() in {7, 8, 9, 15, 16, 17, 35, 53}:
        return "heteroatom"
    return "generic_R"


def _boundary_continuity(
    *,
    side: str,
    core_map_number: Optional[int],
    attachment_map_number: Optional[int],
    bond_order: str,
    opposite_by_map: Mapping[int, _Location],
    opposite_active_coordinates: set[_Coordinate],
) -> str:
    if core_map_number is None or attachment_map_number is None:
        return "unresolved"
    core = opposite_by_map.get(int(core_map_number))
    attachment = opposite_by_map.get(int(attachment_map_number))
    if core is None or attachment is None or core[0] is not attachment[0]:
        return "reactant_only" if side == "reactant" else "product_only"
    component, molecule, core_atom_index = core
    _, _, attachment_atom_index = attachment
    bond = molecule.GetBondBetweenAtoms(core_atom_index, attachment_atom_index)
    if (
        bond is not None
        and str(bond.GetBondType()).upper() == bond_order
        and (
            component.component_index,
            attachment_atom_index,
        )
        not in opposite_active_coordinates
    ):
        return "retained"
    return "reactant_only" if side == "reactant" else "product_only"


def _build_boundaries_for_side(
    *,
    side: str,
    components: Sequence[ReactionComponent],
    active_coordinates: set[_Coordinate],
    opposite_by_map: Mapping[int, _Location],
    opposite_active_coordinates: set[_Coordinate],
) -> Tuple[ReactionCoreBoundary, ...]:
    values = []
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        component_active = {
            atom_index
            for component_index, atom_index in active_coordinates
            if component_index == component.component_index
        }
        for core_atom_index in sorted(component_active):
            core_atom = molecule.GetAtomWithIdx(core_atom_index)
            core_map_number = int(core_atom.GetAtomMapNum()) or None
            for attachment in core_atom.GetNeighbors():
                attachment_atom_index = int(attachment.GetIdx())
                if (
                    component.component_index,
                    attachment_atom_index,
                ) in active_coordinates:
                    continue
                if attachment.GetAtomicNum() <= 1:
                    continue
                bond = molecule.GetBondBetweenAtoms(
                    core_atom_index,
                    attachment_atom_index,
                )
                bond_order = str(bond.GetBondType()).upper()
                fragment_indices = _fragment_indices(
                    molecule,
                    attachment_atom_index,
                    blocked=component_active,
                )
                boundary_class = _boundary_class(
                    molecule,
                    attachment_atom_index,
                    fragment_indices,
                )
                attachment_map_number = (
                    int(attachment.GetAtomMapNum()) or None
                )
                continuity = _boundary_continuity(
                    side=side,
                    core_map_number=core_map_number,
                    attachment_map_number=attachment_map_number,
                    bond_order=bond_order,
                    opposite_by_map=opposite_by_map,
                    opposite_active_coordinates=opposite_active_coordinates,
                )
                fragment = _fragment_smiles(molecule, fragment_indices)
                heavy_atom_count = len(fragment_indices)
                heteroatom_count = sum(
                    molecule.GetAtomWithIdx(atom_index).GetAtomicNum()
                    not in {1, 6}
                    for atom_index in fragment_indices
                )
                aromatic_atom_count = sum(
                    molecule.GetAtomWithIdx(atom_index).GetIsAromatic()
                    for atom_index in fragment_indices
                )
                functional_group_ids = _functional_group_ids(
                    component,
                    fragment_indices,
                )
                boundary_payload = {
                    "side": side,
                    "core_element": core_atom.GetSymbol(),
                    "attachment_element": attachment.GetSymbol(),
                    "bond_order": bond_order,
                    "boundary_class": boundary_class,
                    "continuity": continuity,
                    "fragment_smiles": fragment,
                    "functional_group_ids": functional_group_ids,
                }
                values.append(
                    ReactionCoreBoundary(
                        boundary_id=_digest(
                            "RCB1",
                            boundary_payload,
                            length=24,
                        ),
                        side=side,  # type: ignore[arg-type]
                        core_component_index=component.component_index,
                        core_atom_index=core_atom_index,
                        core_atom_map_number=core_map_number,
                        attachment_atom_index=attachment_atom_index,
                        attachment_atom_map_number=attachment_map_number,
                        attachment_element=attachment.GetSymbol(),
                        bond_order=bond_order,
                        boundary_class=boundary_class,
                        continuity=continuity,  # type: ignore[arg-type]
                        fragment_smiles=fragment,
                        fragment_heavy_atom_count=heavy_atom_count,
                        fragment_heteroatom_count=heteroatom_count,
                        fragment_aromatic_atom_count=aromatic_atom_count,
                        functional_group_ids=functional_group_ids,
                    )
                )
    return tuple(
        sorted(
            values,
            key=lambda boundary: (
                boundary.side,
                boundary.boundary_class,
                boundary.fragment_smiles,
                boundary.bond_order,
                boundary.core_component_index,
                boundary.core_atom_index,
                boundary.attachment_atom_index,
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


def _boundary_display(boundary_class: ReactionCoreBoundaryClass) -> str:
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
    }[boundary_class]


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
    boundary_classes: Mapping[tuple[str, int, int, int], ReactionCoreBoundaryClass],
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
            value = _active_neighbor_display(
                molecule,
                atom_index,
                neighbor_index,
            )
        else:
            boundary_class = boundary_classes.get(
                (side, component_index, atom_index, neighbor_index),
                "generic_R",
            )
            value = (
                neighbor.GetSymbol()
                if boundary_class == "heteroatom"
                else _boundary_display(boundary_class)
            )
        tokens.append(f"{_bond_prefix(order)}{value}")

    def token_order(token: str) -> tuple[int, str]:
        plain = token.lstrip("=#:~")
        if plain == "H":
            return 0, token
        if plain in {"Ar", "HetAr", "R", "Alkenyl", "Alkynyl", "Acyl"}:
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
    boundary_classes: Mapping[tuple[str, int, int, int], ReactionCoreBoundaryClass],
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
        "total_hydrogens": int(
            atom.GetTotalNumHs(includeNeighbors=True)
        ),
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
            boundary_classes=boundary_classes,
            side=side,
        ),
        state_key=_digest("RAS1", state_payload, length=24),
    )


def _edit_graph(
    edits: Sequence[ReactionEdit],
) -> tuple[
    set[_AtomIdentity],
    Dict[_AtomIdentity, set[_AtomIdentity]],
    Counter[_AtomIdentity],
    tuple[object, ...],
    Tuple[str, ...],
]:
    identities: set[_AtomIdentity] = set()
    adjacency: Dict[_AtomIdentity, set[_AtomIdentity]] = defaultdict(set)
    heavy_incidence: Counter[_AtomIdentity] = Counter()
    atom_labels: Dict[_AtomIdentity, tuple[object, ...]] = {}
    incidents: Dict[_AtomIdentity, list[tuple[object, ...]]] = defaultdict(list)
    public_tokens = []
    for edit in edits:
        left = _atom_identity(edit.atom_1)
        identities.add(left)
        atom_labels[left] = (
            edit.atom_1.element,
            edit.atom_1.formal_charge,
            edit.atom_1.aromatic,
            edit.atom_1.hybridization,
        )
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
            adjacency[left].add(right)
            adjacency[right].add(left)
            heavy_incidence[left] += 1
            heavy_incidence[right] += 1
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
            pair = "-".join(
                sorted((edit.atom_1.element, edit.atom_2.element))
            )
        public_tokens.append(
            f"{edit.edit_type}:{pair}:"
            f"{edit.old_order or 'NONE'}>{edit.new_order or 'NONE'}"
        )
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
        heavy_incidence,
        graph_payload,
        tuple(sorted(public_tokens)),
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


def _stable_boundary_count(
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
            sorted(
                "|".join(token.split("|")[:3])
                for token in state.neighbor_tokens
            )
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


def _boundary_identity(
    boundary: ReactionCoreBoundary,
    *,
    exact: bool,
) -> object:
    base = (
        boundary.side,
        boundary.attachment_element,
        boundary.bond_order,
        boundary.boundary_class,
        boundary.continuity,
    )
    if not exact:
        return base
    return base + (
        boundary.fragment_smiles,
        boundary.fragment_heavy_atom_count,
        boundary.fragment_heteroatom_count,
        boundary.fragment_aromatic_atom_count,
        boundary.functional_group_ids,
    )


def build_reaction_core_projection(
    *,
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    edits: Sequence[ReactionEdit],
    evidence: str,
    confidence: float,
) -> Optional[ReactionCoreProjection]:
    """Build a grammar- and template-free minimized reaction-core projection.

    At least one edited mapped atom must be observed on both sides.  This POC
    deliberately abstains when only a reactant-side rewrite is available,
    because it would otherwise present a predicted product state as observed.
    """
    if not edits:
        return None
    (
        identities,
        adjacency,
        heavy_incidence,
        edit_graph_payload,
        edit_tokens,
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
    boundaries = (
        _build_boundaries_for_side(
            side="reactant",
            components=reactants,
            active_coordinates=reactant_active,
            opposite_by_map=product_by_map,
            opposite_active_coordinates=product_active,
        )
        + _build_boundaries_for_side(
            side="product",
            components=products,
            active_coordinates=product_active,
            opposite_by_map=reactant_by_map,
            opposite_active_coordinates=reactant_active,
        )
    )
    boundary_classes = {
        (
            boundary.side,
            boundary.core_component_index,
            boundary.core_atom_index,
            boundary.attachment_atom_index,
        ): boundary.boundary_class
        for boundary in boundaries
    }
    event_components = _event_components(identities, adjacency)
    selected: list[tuple[_AtomIdentity, int]] = []
    for event in event_components:
        candidates = [
            identity for identity in event if identity in shared_identities
        ]
        if not candidates:
            continue
        maximum_incidence = max(heavy_incidence[identity] for identity in candidates)
        candidates = [
            identity
            for identity in candidates
            if heavy_incidence[identity] == maximum_incidence
        ]
        stable_counts = {
            identity: _stable_boundary_count(
                identity,
                reactant_by_map=reactant_by_map,
                product_by_map=product_by_map,
                active_coordinates_by_side=active_coordinates_by_side,
            )
            for identity in candidates
        }
        maximum_stable = max(stable_counts.values(), default=0)
        selected.extend(
            (identity, stable_counts[identity])
            for identity in candidates
            if stable_counts[identity] == maximum_stable
        )
    centers = []
    for identity, stable_count in sorted(selected, key=lambda value: repr(value[0])):
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
                boundary_classes=boundary_classes,
            )
            if before_location is not None
            else None
        )
        after_state = (
            _build_atom_state(
                side="product",
                location=after_location,
                active_coordinates=product_active,
                boundary_classes=boundary_classes,
            )
            if after_location is not None
            else None
        )
        center_payload = {
            "before": _state_identity(before_state, generic=False),
            "after": _state_identity(after_state, generic=False),
            "incident_edit_count": heavy_incidence[identity],
            "stable_boundary_count": stable_count,
        }
        centers.append(
            ReactionCoreCenter(
                center_id=_digest("RCC1", center_payload, length=24),
                atom_map_number=(
                    int(identity[1]) if identity[0] == "map" else None
                ),
                before_state=before_state,
                after_state=after_state,
                incident_edit_count=heavy_incidence[identity],
                stable_boundary_count=stable_count,
            )
        )
    if not centers:
        return None
    center_transition_payload = tuple(
        sorted(
            (
                _state_identity(center.before_state, generic=False),
                _state_identity(center.after_state, generic=False),
            )
            for center in centers
        )
    )
    generic_transition_payload = tuple(
        sorted(
            (
                _state_identity(center.before_state, generic=True),
                _state_identity(center.after_state, generic=True),
            )
            for center in centers
        )
    )
    exact_key = _digest(
        "RCX1",
        {
            "edit_graph": edit_graph_payload,
            "centers": center_transition_payload,
            "boundaries": tuple(
                sorted(
                    (
                        _boundary_identity(boundary, exact=True)
                        for boundary in boundaries
                    ),
                    key=repr,
                )
            ),
            "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        },
    )
    typed_key = _digest(
        "RCT1",
        {
            "edit_graph": edit_graph_payload,
            "centers": center_transition_payload,
            "boundaries": tuple(
                sorted(
                    (
                        _boundary_identity(boundary, exact=False)
                        for boundary in boundaries
                    ),
                    key=repr,
                )
            ),
            "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        },
    )
    generic_key = _digest(
        "RCG1",
        {
            "centers": generic_transition_payload,
            "edit_tokens": edit_tokens,
            "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        },
    )
    center_transition_key = _digest(
        "RCS1",
        {
            "centers": center_transition_payload,
            "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        },
    )
    generic_label = " + ".join(
        (
            f"{center.before_state.concise_label if center.before_state else '∅'}"
            " → "
            f"{center.after_state.concise_label if center.after_state else '∅'}"
        )
        for center in centers
    )
    warnings = set()
    if any(
        boundary.continuity == "unresolved" for boundary in boundaries
    ):
        warnings.add("REACTION_CORE_BOUNDARY_CONTINUITY_UNRESOLVED")
    if any(identity[0] != "map" for identity in identities):
        warnings.add("REACTION_CORE_PARTIAL_ATOM_PROVENANCE")
    if str(evidence).startswith("external"):
        warnings.add("REACTION_CORE_EXTERNAL_MAPPING_PROPOSAL")
    core_id = _digest(
        "RCP1",
        {
            "keys": (
                exact_key,
                typed_key,
                generic_key,
                center_transition_key,
            ),
            "algorithm_version": "reaction_core_projection.v1",
            "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        },
        length=64,
    )
    return ReactionCoreProjection(
        core_id=core_id,
        exact_core_key=exact_key,
        typed_core_key=typed_key,
        generic_core_key=generic_key,
        center_transition_key=center_transition_key,
        centers=tuple(centers),
        boundaries=tuple(
            sorted(
                boundaries,
                key=lambda boundary: (
                    boundary.side,
                    boundary.boundary_class,
                    boundary.fragment_smiles,
                    boundary.core_component_index,
                    boundary.core_atom_index,
                ),
            )
        ),
        edit_tokens=edit_tokens,
        generic_label=generic_label,
        active_atom_count=len(identities),
        event_count=len(
            [
                event
                for event in event_components
                if any(identity in shared_identities for identity in event)
            ]
        ),
        evidence=str(evidence),
        confidence=float(confidence),
        warnings=tuple(sorted(warnings)),
    )


__all__ = ["build_reaction_core_projection"]
