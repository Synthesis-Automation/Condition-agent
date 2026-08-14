"""General component and ring topology for reaction transformations."""

from __future__ import annotations

from collections import Counter, defaultdict, deque
import hashlib
import json
from typing import Any, Dict, Iterable, Optional, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_edits import EditNormalizationResult, reaction_atom_reference
from .reaction_models import (
    ReactionComponent,
    ReactionEdit,
    ReactionRingChange,
    ReactionSiteReference,
    ReactionTopology,
)


_Coordinate = tuple[int, int]
_Edge = tuple[_Coordinate, _Coordinate]


def _edge(left: _Coordinate, right: _Coordinate) -> _Edge:
    return tuple(sorted((left, right)))  # type: ignore[return-value]


def _shortest_path_without_edge(
    adjacency: Dict[_Coordinate, set[_Coordinate]],
    start: _Coordinate,
    end: _Coordinate,
    excluded: _Edge,
) -> tuple[_Coordinate, ...]:
    queue = deque(((start,),))
    visited = {start}
    while queue:
        path = queue.popleft()
        current = path[-1]
        for neighbor in sorted(adjacency.get(current, ())):
            if _edge(current, neighbor) == excluded:
                continue
            if neighbor == end:
                return (*path, neighbor)
            if neighbor in visited:
                continue
            visited.add(neighbor)
            queue.append((*path, neighbor))
    return ()


def _canonical_cycle(
    cycle: tuple[_Coordinate, ...],
) -> tuple[_Coordinate, ...]:
    variants = []
    for values in (cycle, tuple(reversed(cycle))):
        variants.extend(
            values[offset:] + values[:offset]
            for offset in range(len(values))
        )
    return min(variants)


def _cycle_edges(cycle: tuple[_Coordinate, ...]) -> tuple[_Edge, ...]:
    return tuple(
        _edge(cycle[index], cycle[(index + 1) % len(cycle)])
        for index in range(len(cycle))
    )


def _ring_change_id(payload: object) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return "RRG1:" + hashlib.sha256(encoded).hexdigest()[:32]


def build_reaction_ring_changes(
    *,
    reactants: Tuple[ReactionComponent, ...],
    edit_result: EditNormalizationResult,
) -> Tuple[ReactionRingChange, ...]:
    """Return minimal graph facts for cycles closed by observed formed bonds."""
    components = {component.component_index: component for component in reactants}
    adjacency: Dict[_Coordinate, set[_Coordinate]] = defaultdict(set)
    bond_orders: Dict[_Edge, str] = {}
    valid_atoms: set[_Coordinate] = set()
    for component in reactants:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() > 1:
                valid_atoms.add((component.component_index, int(atom.GetIdx())))
        for bond in molecule.GetBonds():
            left_atom = bond.GetBeginAtom()
            right_atom = bond.GetEndAtom()
            if left_atom.GetAtomicNum() <= 1 or right_atom.GetAtomicNum() <= 1:
                continue
            left = (component.component_index, int(left_atom.GetIdx()))
            right = (component.component_index, int(right_atom.GetIdx()))
            adjacency[left].add(right)
            adjacency[right].add(left)
            bond_orders[_edge(left, right)] = str(bond.GetBondType()).upper()

    formed_edges: Dict[_Edge, ReactionEdit] = {}
    for edit in edit_result.edits:
        if edit.atom_2 is None or {
            edit.atom_1.side,
            edit.atom_2.side,
        } != {"reactant"}:
            continue
        left = (edit.atom_1.component_index, edit.atom_1.atom_index)
        right = (edit.atom_2.component_index, edit.atom_2.atom_index)
        if left not in valid_atoms or right not in valid_atoms or left == right:
            continue
        edge = _edge(left, right)
        if edit.edit_type == "broken":
            adjacency[left].discard(right)
            adjacency[right].discard(left)
            bond_orders.pop(edge, None)
        elif edit.edit_type == "formed":
            adjacency[left].add(right)
            adjacency[right].add(left)
            bond_orders[edge] = str(edit.new_order or "SINGLE").upper()
            formed_edges[edge] = edit
        elif edit.edit_type == "order_changed":
            bond_orders[edge] = str(edit.new_order or "SINGLE").upper()

    cycles: Dict[tuple[_Edge, ...], tuple[_Coordinate, ...]] = {}
    for formed_edge in sorted(formed_edges):
        left, right = formed_edge
        path = _shortest_path_without_edge(
            adjacency,
            left,
            right,
            formed_edge,
        )
        if len(path) < 3:
            continue
        cycle = _canonical_cycle(path)
        edges = tuple(sorted(_cycle_edges(cycle)))
        if any(edge in formed_edges for edge in edges):
            cycles.setdefault(edges, cycle)

    changes = []
    for edges, cycle in sorted(cycles.items()):
        references = tuple(
            reaction_atom_reference(
                components[component_index],
                atom_index,
                side="reactant",
            )
            for component_index, atom_index in cycle
        )
        orders = tuple(
            bond_orders.get(edge, "UNKNOWN") for edge in _cycle_edges(cycle)
        )
        formed_types = tuple(
            sorted(
                "-".join(
                    sorted(
                        (
                            formed_edges[edge].atom_1.element,
                            formed_edges[edge].atom_2.element,  # type: ignore[union-attr]
                        )
                    )
                )
                for edge in edges
                if edge in formed_edges
            )
        )
        payload = {
            "atoms": cycle,
            "elements": tuple(reference.element for reference in references),
            "orders": orders,
            "formed_bonds": formed_types,
            "evidence": edit_result.evidence,
        }
        changes.append(
            ReactionRingChange(
                change_id=_ring_change_id(payload),
                change_type="formed",
                atom_references=references,
                element_sequence=tuple(
                    reference.element for reference in references
                ),
                bond_orders_after=orders,
                source_component_indices=tuple(
                    sorted({reference.component_index for reference in references})
                ),
                formed_bond_types=formed_types,
                aromatic_after=bool(orders) and all(
                    order == "AROMATIC" for order in orders
                ),
                evidence=edit_result.evidence,
                confidence=edit_result.confidence,
            )
        )
    return tuple(changes)


def _cycle_rank(
    components: Tuple[ReactionComponent, ...],
    included_atoms: Optional[Dict[int, set[int]]] = None,
) -> Optional[int]:
    """Return cycle rank for the full or atom-corresponded component graphs."""
    from rdkit import Chem

    total = 0
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            return None
        if included_atoms is None:
            fragment_count = len(Chem.GetMolFrags(molecule, asMols=False))
            total += molecule.GetNumBonds() - molecule.GetNumAtoms() + fragment_count
            continue
        selected = set(included_atoms.get(component.component_index, set()))
        if not selected:
            continue
        adjacency = {atom_index: set() for atom_index in selected}
        edge_count = 0
        for bond in molecule.GetBonds():
            left = int(bond.GetBeginAtomIdx())
            right = int(bond.GetEndAtomIdx())
            if left not in selected or right not in selected:
                continue
            adjacency[left].add(right)
            adjacency[right].add(left)
            edge_count += 1
        remaining = set(selected)
        fragment_count = 0
        while remaining:
            fragment_count += 1
            stack = [remaining.pop()]
            while stack:
                current = stack.pop()
                neighbors = adjacency[current] & remaining
                remaining.difference_update(neighbors)
                stack.extend(neighbors)
        total += edge_count - len(selected) + fragment_count
    return total


def _formed_edits(edits: Iterable[ReactionEdit]) -> Tuple[ReactionEdit, ...]:
    return tuple(
        edit
        for edit in edits
        if edit.edit_type == "formed" and edit.atom_2 is not None
    )


def _neighbor_descriptor(
    center: Any,
    neighbor: Any,
    molecule: Any,
) -> tuple[Any, ...]:
    bond = molecule.GetBondBetweenAtoms(
        int(center.GetIdx()),
        int(neighbor.GetIdx()),
    )
    return (
        str(neighbor.GetSymbol()),
        str(bond.GetBondType()),
        bool(neighbor.GetIsAromatic()),
        int(neighbor.GetFormalCharge()),
    )


def _has_one_sided_departure(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    edits: Iterable[ReactionEdit],
) -> bool:
    """Return whether a substitution center has an observed departing branch."""

    edits = tuple(edits)
    formed_atoms = {
        (atom.component_index, atom.atom_index)
        for edit in edits
        if edit.edit_type == "formed"
        for atom in (edit.atom_1, edit.atom_2)
        if atom is not None and atom.side == "reactant"
    }
    for edit in edits:
        if edit.edit_type != "broken" or edit.atom_2 is None:
            continue
        endpoints = tuple(
            (atom.component_index, atom.atom_index)
            for atom in (edit.atom_1, edit.atom_2)
            if atom.side == "reactant"
        )
        if len(endpoints) == 2 and sum(
            endpoint in formed_atoms for endpoint in endpoints
        ) == 1:
            return True

    reactant_components = {
        component.component_index: component for component in reactants
    }
    product_atoms_by_map = {}
    for component in products:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            map_number = int(atom.GetAtomMapNum())
            if map_number > 0:
                product_atoms_by_map[map_number] = (molecule, atom)

    for edit in edits:
        if edit.edit_type != "formed":
            continue
        for endpoint in (edit.atom_1, edit.atom_2):
            if (
                endpoint is None
                or endpoint.side != "reactant"
                or endpoint.element != "C"
            ):
                continue
            component = reactant_components.get(endpoint.component_index)
            molecule = parse_smiles(component.input_smiles) if component else None
            if molecule is None:
                continue
            carbon = molecule.GetAtomWithIdx(endpoint.atom_index)
            product_carbon = product_atoms_by_map.get(
                int(carbon.GetAtomMapNum())
            )
            retained_unmapped: Counter[tuple[Any, ...]] = Counter()
            if product_carbon is not None:
                product_molecule, product_atom = product_carbon
                retained_unmapped.update(
                    _neighbor_descriptor(product_atom, neighbor, product_molecule)
                    for neighbor in product_atom.GetNeighbors()
                    if int(neighbor.GetAtomMapNum()) == 0
                )
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetSymbol() == "C":
                    continue
                map_number = int(neighbor.GetAtomMapNum())
                if map_number > 0 and map_number in product_atoms_by_map:
                    continue
                descriptor = _neighbor_descriptor(carbon, neighbor, molecule)
                if map_number == 0 and retained_unmapped[descriptor]:
                    retained_unmapped[descriptor] -= 1
                    continue
                return True
    return False


def build_reaction_topology(
    *,
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    edit_result: EditNormalizationResult,
) -> Optional[ReactionTopology]:
    """Derive topology from atom-provenanced edits, independent of family names."""
    if not edit_result.edits:
        return None

    components_by_index = {
        component.component_index: component for component in reactants
    }
    participating = {
        atom.component_index
        for edit in edit_result.edits
        for atom in (edit.atom_1, edit.atom_2)
        if atom is not None and atom.side == "reactant"
    }
    formed_scopes = []
    tether_distances = []
    ring_changes = build_reaction_ring_changes(
        reactants=reactants,
        edit_result=edit_result,
    )
    from rdkit import Chem

    for edit in _formed_edits(edit_result.edits):
        assert edit.atom_2 is not None
        if {
            edit.atom_1.side,
            edit.atom_2.side,
        } != {"reactant"}:
            # Product-only endpoints can occur in partial external mappings
            # when a reported product atom has no supplied reactant origin.
            # They cannot establish reactant component scope or tether length.
            continue
        if edit.atom_1.component_index != edit.atom_2.component_index:
            formed_scopes.append("intermolecular")
            continue
        formed_scopes.append("intramolecular")
        component = components_by_index.get(edit.atom_1.component_index)
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

    scope_set = set(formed_scopes)
    if scope_set == {"intramolecular"}:
        reaction_scope = "intramolecular"
    elif scope_set == {"intermolecular"}:
        reaction_scope = "intermolecular"
    elif len(scope_set) > 1:
        reaction_scope = "mixed"
    elif len(participating) == 1:
        reaction_scope = "unimolecular"
    else:
        reaction_scope = "unresolved"

    project_one_sided_fragments = any(
        warning.startswith(
            "REACTANT_ONLY_FRAGMENT_INTERNAL_BONDS_IGNORED:"
        )
        for warning in edit_result.warnings
    ) or _has_one_sided_departure(
        reactants,
        products,
        edit_result.edits,
    )
    if project_one_sided_fragments:
        reactant_included: Dict[int, set[int]] = defaultdict(set)
        product_included: Dict[int, set[int]] = defaultdict(set)
        reactant_locations: Dict[int, tuple[int, int]] = {}
        product_locations: Dict[int, tuple[int, int]] = {}
        if edit_result.atom_correspondence:
            for (
                reactant_component,
                reactant_atom,
                product_component,
                product_atom,
            ) in edit_result.atom_correspondence:
                reactant_included[reactant_component].add(reactant_atom)
                product_included[product_component].add(product_atom)
        else:
            for component in reactants:
                molecule = parse_smiles(component.input_smiles)
                if molecule is None:
                    continue
                for atom in molecule.GetAtoms():
                    map_number = int(atom.GetAtomMapNum())
                    if map_number > 0:
                        reactant_locations[map_number] = (
                            component.component_index,
                            int(atom.GetIdx()),
                        )
            for component in products:
                molecule = parse_smiles(component.input_smiles)
                if molecule is None:
                    continue
                for atom in molecule.GetAtoms():
                    map_number = int(atom.GetAtomMapNum())
                    if map_number > 0:
                        product_locations[map_number] = (
                            component.component_index,
                            int(atom.GetIdx()),
                        )
            for map_number in set(reactant_locations).intersection(
                product_locations
            ):
                reactant_component, reactant_atom = reactant_locations[
                    map_number
                ]
                product_component, product_atom = product_locations[map_number]
                reactant_included[reactant_component].add(reactant_atom)
                product_included[product_component].add(product_atom)
        reactant_rank = _cycle_rank(reactants, reactant_included)
        product_rank = _cycle_rank(products, product_included)
    else:
        reactant_rank = _cycle_rank(reactants)
        product_rank = _cycle_rank(products)
    ring_count_delta = (
        product_rank - reactant_rank
        if reactant_rank is not None and product_rank is not None
        else None
    )
    return ReactionTopology(
        reaction_scope=reaction_scope,
        participating_component_indices=tuple(sorted(participating)),
        formed_bond_scopes=tuple(sorted(formed_scopes)),
        reactant_tether_distances=tuple(sorted(tether_distances)),
        # A short path in the unedited reactant is not sufficient evidence of
        # ring formation: rearrangements can break that path while forming the
        # new bond. Only cycles present after applying the complete edit set
        # are reported as formed rings.
        formed_ring_sizes=tuple(
            sorted({change.ring_size for change in ring_changes})
        ),
        ring_count_delta=ring_count_delta,
        evidence=edit_result.evidence,
        confidence=edit_result.confidence,
        ring_changes=ring_changes,
    )


def topology_label_prefix(topology: Optional[ReactionTopology]) -> str:
    """Return a concise display-only prefix for intramolecular transformations."""
    if topology is None or topology.reaction_scope != "intramolecular":
        return ""
    if len(topology.formed_ring_sizes) == 1:
        return f"intramolecular ({topology.formed_ring_sizes[0]}-membered ring) "
    return "intramolecular "


def assignment_component_scope(
    assignment: dict[str, ReactionSiteReference],
) -> str:
    """Describe whether interpreted roles share reactant components."""
    component_indices = [site.component_index for site in assignment.values()]
    if len(component_indices) < 2:
        return "unresolved"
    unique_count = len(set(component_indices))
    if unique_count == 1:
        return "intramolecular"
    if unique_count == len(component_indices):
        return "intermolecular"
    return "mixed"


__all__ = [
    "assignment_component_scope",
    "build_reaction_ring_changes",
    "build_reaction_topology",
    "topology_label_prefix",
]
