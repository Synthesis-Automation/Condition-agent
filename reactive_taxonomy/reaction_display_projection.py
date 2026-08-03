"""Build a display-only minimum reaction graph from observed reaction edits.

The projection is deliberately separate from :class:`ReactionCoreProjection`.
It is a chemist-facing drawing aid and must not be used as reaction identity or
as an input to recommendation, admission, or retrieval.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Literal, Mapping, Optional, Tuple

from .chemistry.rdkit_utils import (
    parse_smiles,
    prepare_fragment_serialization_copy,
    rdkit_available,
)
from .reaction_core.substituents import build_substituent_profile
from .reaction_models import ReactionCoreSubstituentProfile
from .reaction_render_context import ReactionRenderContext

try:  # pragma: no cover - exercised through the public availability check
    from rdkit import Chem
except ImportError:  # pragma: no cover
    Chem = None  # type: ignore[assignment]


_DEFINITION_PATH = (
    Path(__file__).with_name("definitions")
    / "reaction_display_projection.v1.json"
)


@dataclass(frozen=True)
class ReactionDisplayAromaticRelation:
    """Position of one hidden substituent relative to an active ring atom."""

    reactive_atom_index: int
    reactive_atom_map_number: Optional[int]
    attachment_ring_atom_index: int
    aromatic_ring_distance: int
    same_ring: bool
    positional_relation: Literal[
        "ipso", "ortho", "meta", "para", "other"
    ]


@dataclass(frozen=True)
class ReactionDisplaySubstituent:
    """Exact omitted fragment behind an R group or aromatic hydrogen cap."""

    substituent_id: str
    side: Literal["reactant", "product"]
    source_component_index: int
    center_atom_index: int
    center_atom_map_number: Optional[int]
    center_element: str
    center_hybridization: str
    attachment_atom_index: int
    attachment_atom_map_number: Optional[int]
    attachment_element: str
    attachment_bond_order: str
    atom_indices: Tuple[int, ...]
    atom_map_numbers: Tuple[int, ...]
    fragment_smiles: str
    remote_class: Optional[str]
    substituent_profile: ReactionCoreSubstituentProfile
    continuity: str
    boundary_action: Literal["r_placeholder", "aromatic_hydrogen_cap"]
    placeholder_index: Optional[int]
    display_label: Optional[str]
    aromatic_relations: Tuple[ReactionDisplayAromaticRelation, ...]


@dataclass(frozen=True)
class ReactionDisplayConnector:
    """One omitted connected subgraph shared by multiple display ports."""

    connector_id: str
    display_label: str
    side: Literal["reactant", "product"]
    source_component_index: int
    port_substituent_ids: Tuple[str, ...]
    placeholder_indices: Tuple[int, ...]
    port_display_labels: Tuple[str, ...]
    attachment_atom_indices: Tuple[int, ...]
    atom_indices: Tuple[int, ...]
    atom_map_numbers: Tuple[int, ...]
    fragment_smiles: str
    hidden_atom_count: int
    shortest_path_atom_indices: Tuple[int, ...]
    shortest_path_bond_count: Optional[int]


@dataclass(frozen=True)
class ReactionDisplayComponent:
    """One source component after display-only graph abstraction."""

    side: Literal["reactant", "product"]
    source_component_index: int
    display_smiles: str
    render_smiles: str
    retained_atom_indices: Tuple[int, ...]
    active_atom_indices: Tuple[int, ...]
    retained_aromatic_system_count: int
    removed_substituent_count: int
    r_group_count: int
    substituents: Tuple[ReactionDisplaySubstituent, ...]
    connectors: Tuple[ReactionDisplayConnector, ...]


@dataclass(frozen=True)
class ReactionDisplayProjection:
    """Deterministic minimum reaction SMILES intended only for rendering."""

    minimum_reaction_smiles: str
    render_reaction_smiles: str
    reactants: Tuple[ReactionDisplayComponent, ...]
    products: Tuple[ReactionDisplayComponent, ...]
    substituents: Tuple[ReactionDisplaySubstituent, ...]
    connectors: Tuple[ReactionDisplayConnector, ...]
    reaction_scope: str
    formed_ring_sizes: Tuple[int, ...]
    annotation: str | None
    evidence_status: str
    confidence: float
    warnings: Tuple[str, ...]
    definition_id: str
    schema_version: str


@lru_cache(maxsize=1)
def load_reaction_display_projection_definition() -> Dict[str, Any]:
    """Load and validate the display-minimization policy."""
    with _DEFINITION_PATH.open("r", encoding="utf-8") as handle:
        definition = dict(json.load(handle))
    expected = {
        "schema_version": "1.0",
        "definition_id": "reaction_display_projection.v1.5",
        "aromatic_system_policy": "retain_aromatic_bond_component",
        "aromatic_valence_completion_policy": (
            "retain_exocyclic_multiple_bonds"
        ),
        "multiple_bond_policy": "retain_contiguous_multiple_bond_unit",
        "active_heteroatom_shell_policy": (
            "retain_direct_noncarbon_neighbors_and_adjacent_"
            "nonaromatic_multiple_bond_systems"
        ),
        "saturated_carbon_policy": "retain_active_atom_only",
        "aromatic_carbon_boundary_policy": "remove_and_hydrogen_cap",
        "aromatic_heteroatom_boundary_policy": "retain_R_attachment",
        "nonaromatic_boundary_policy": "retain_R_attachment",
        "placeholder_smiles": "*",
        "placeholder_label_template": "R{index}",
        "hidden_connector_policy": "collapse_shared_omitted_component",
        "hidden_connector_bond_type": "zero_order_dashed",
        "hidden_connector_label_template": "S{index}",
    }
    for key, value in expected.items():
        if definition.get(key) != value:
            raise ValueError(
                f"unexpected reaction display projection setting: {key}"
            )
    if "{ring_size}" not in str(
        definition.get("intramolecular_note_template") or ""
    ):
        raise ValueError("intramolecular note template requires {ring_size}")
    return definition


def _require_rdkit() -> None:
    if not rdkit_available() or Chem is None:
        raise RuntimeError("RDKit is required for reaction display projection.")


def _active_coordinates(
    analysis: ReactionRenderContext,
    *,
    side: Literal["reactant", "product"],
) -> tuple[Dict[int, set[int]], Dict[int, set[int]]]:
    core = analysis.reaction_core
    if core is None:
        return {}, {}
    active: Dict[int, set[int]] = {}
    inherited_aromatic: Dict[int, set[int]] = {}
    for transition in core.atom_transitions:
        state = (
            transition.before_state
            if side == "reactant"
            else transition.after_state
        )
        if state is not None:
            active.setdefault(int(state.component_index), set()).add(
                int(state.atom_index)
            )
            aromatic_before = (
                transition.before_state is not None
                and bool(transition.before_state.aromatic)
            )
            if side == "reactant":
                follows_aromatic_policy = bool(state.aromatic)
            else:
                follows_aromatic_policy = aromatic_before
            if follows_aromatic_policy:
                inherited_aromatic.setdefault(
                    int(state.component_index), set()
                ).add(int(state.atom_index))
    return active, inherited_aromatic


def _aromatic_system(molecule: Any, seed: int) -> set[int]:
    """Return the aromatic-bond-connected system containing ``seed``."""
    atom = molecule.GetAtomWithIdx(seed)
    if not atom.GetIsAromatic():
        return set()
    retained = {seed}
    pending = [seed]
    while pending:
        current = pending.pop()
        current_atom = molecule.GetAtomWithIdx(current)
        for bond in current_atom.GetBonds():
            if not bond.GetIsAromatic():
                continue
            neighbor = int(bond.GetOtherAtomIdx(current))
            if neighbor in retained:
                continue
            neighbor_atom = molecule.GetAtomWithIdx(neighbor)
            if neighbor_atom.GetIsAromatic():
                retained.add(neighbor)
                pending.append(neighbor)
    return retained


def _copy_atom(atom: Any) -> Any:
    copied = Chem.Atom(atom)
    copied.SetAtomMapNum(0)
    return copied


@dataclass
class _SubstituentDraft:
    """Mutable placeholder record until cross-side labels are assigned."""

    side: Literal["reactant", "product"]
    source_component_index: int
    center_atom_index: int
    center_atom_map_number: Optional[int]
    center_element: str
    center_hybridization: str
    attachment_atom_index: int
    attachment_atom_map_number: Optional[int]
    attachment_element: str
    attachment_bond_order: str
    atom_indices: Tuple[int, ...]
    atom_map_numbers: Tuple[int, ...]
    fragment_smiles: str
    remote_class: Optional[str]
    substituent_profile: ReactionCoreSubstituentProfile
    continuity: str
    boundary_action: Literal["r_placeholder", "aromatic_hydrogen_cap"]
    aromatic_relations: Tuple[ReactionDisplayAromaticRelation, ...]
    dummy_atom_index: Optional[int]
    placeholder_index: Optional[int] = None
    display_label: Optional[str] = None


@dataclass
class _ConnectorDraft:
    """Mutable shared-omission record until scaffold labels are assigned."""

    side: Literal["reactant", "product"]
    source_component_index: int
    ports: Tuple[_SubstituentDraft, ...]
    atom_indices: Tuple[int, ...]
    atom_map_numbers: Tuple[int, ...]
    fragment_smiles: str
    shortest_path_atom_indices: Tuple[int, ...]
    shortest_path_bond_count: Optional[int]
    connector_index: Optional[int] = None
    display_label: Optional[str] = None


@dataclass
class _ComponentDraft:
    """One drawable component before global R-group label assignment."""

    side: Literal["reactant", "product"]
    source_component_index: int
    molecule: Any
    retained_atom_indices: Tuple[int, ...]
    active_atom_indices: Tuple[int, ...]
    retained_aromatic_system_count: int
    removed_substituent_count: int
    r_group_count: int
    substituents: list[_SubstituentDraft]
    connectors: list[_ConnectorDraft]


def _reaction_interface_closure(
    molecule: Any,
    *,
    active_atom_indices: set[int],
    aromatic_seed_indices: set[int],
    explicit_remote_atom_indices: set[int],
) -> tuple[set[int], set[tuple[int, ...]], set[int]]:
    """Expand edit endpoints by deterministic local structural rules."""
    retained = set(active_atom_indices).union(explicit_remote_atom_indices)
    aromatic_systems: set[tuple[int, ...]] = set()
    aromatic_policy_atoms: set[int] = set()
    for atom_index in sorted(aromatic_seed_indices):
        system = _aromatic_system(molecule, atom_index)
        if system:
            retained.update(system)
            aromatic_policy_atoms.update(system)
            aromatic_systems.add(tuple(sorted(system)))

    # An edited heteroatom can be the mapped endpoint of cleavage while the
    # functional group it belongs to remains unchanged. Preserve an adjacent
    # carbon pi-system so ester hydrolysis is shown as C(=O)OR -> C(=O)OH,
    # rather than the misleading isolated OR -> OH. Aromatic neighbors follow
    # Aromatic substituent handling remains governed by its separate policy.
    adjacent_pi_seeds: set[int] = set()
    for atom_index in sorted(active_atom_indices):
        atom = molecule.GetAtomWithIdx(atom_index)
        if atom.GetAtomicNum() in {0, 1, 6}:
            continue
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue
            neighbor_index = int(neighbor.GetIdx())
            if neighbor.GetIsAromatic():
                continue
            if any(
                not bond.GetIsAromatic()
                and bond.GetBondTypeAsDouble() >= 2.0
                for bond in neighbor.GetBonds()
            ):
                retained.add(neighbor_index)
                adjacent_pi_seeds.add(neighbor_index)

    pending = [
        value
        for value in sorted(active_atom_indices.union(adjacent_pi_seeds))
        if not molecule.GetAtomWithIdx(value).GetIsAromatic()
    ]
    visited = set(pending)
    while pending:
        atom_index = pending.pop()
        atom = molecule.GetAtomWithIdx(atom_index)
        for bond in atom.GetBonds():
            if bond.GetIsAromatic() or bond.GetBondTypeAsDouble() < 2.0:
                continue
            neighbor_index = int(bond.GetOtherAtomIdx(atom_index))
            retained.add(neighbor_index)
            neighbor = molecule.GetAtomWithIdx(neighbor_index)
            if not neighbor.GetIsAromatic() and neighbor_index not in visited:
                visited.add(neighbor_index)
                pending.append(neighbor_index)

    # Complete retained aromatic valence before spectator abstraction. An
    # exocyclic C=O, C=N, S=O, or analogous multiple bond is part of the
    # drawable pi-system even when that bond is unchanged by the reaction.
    for atom_index in tuple(sorted(aromatic_policy_atoms)):
        atom = molecule.GetAtomWithIdx(atom_index)
        for bond in atom.GetBonds():
            if bond.GetIsAromatic() or bond.GetBondTypeAsDouble() < 2.0:
                continue
            retained.add(int(bond.GetOtherAtomIdx(atom_index)))

    # Preserve the elemental shell of a local reaction interface. This keeps
    # C(=O)OH, B(OH)2, S(=O)2, and analogous handles explicit while ordinary
    # carbon frameworks remain eligible for R-group abstraction.
    for atom_index in tuple(sorted(retained)):
        if atom_index in aromatic_policy_atoms:
            continue
        atom = molecule.GetAtomWithIdx(atom_index)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() not in {0, 1, 6}:
                retained.add(int(neighbor.GetIdx()))
    return retained, aromatic_systems, aromatic_policy_atoms


def _omitted_components(
    molecule: Any,
    retained: set[int],
) -> Tuple[Tuple[int, ...], ...]:
    omitted = {
        int(atom.GetIdx())
        for atom in molecule.GetAtoms()
        if int(atom.GetIdx()) not in retained
    }
    components = []
    while omitted:
        seed = min(omitted)
        pending = [seed]
        component = {seed}
        omitted.remove(seed)
        while pending:
            current = pending.pop()
            atom = molecule.GetAtomWithIdx(current)
            for neighbor in atom.GetNeighbors():
                neighbor_index = int(neighbor.GetIdx())
                if neighbor_index in omitted:
                    omitted.remove(neighbor_index)
                    component.add(neighbor_index)
                    pending.append(neighbor_index)
        components.append(tuple(sorted(component)))
    return tuple(components)


def _connected_component_by_atom(
    molecule: Any,
    atom_indices: set[int],
) -> Dict[int, int]:
    """Index connected components within an induced atom subgraph."""
    remaining = set(atom_indices)
    component_by_atom: Dict[int, int] = {}
    component_index = 0
    while remaining:
        seed = min(remaining)
        remaining.remove(seed)
        pending = [seed]
        component_by_atom[seed] = component_index
        while pending:
            current = pending.pop()
            for neighbor in molecule.GetAtomWithIdx(current).GetNeighbors():
                neighbor_index = int(neighbor.GetIdx())
                if neighbor_index in remaining:
                    remaining.remove(neighbor_index)
                    component_by_atom[neighbor_index] = component_index
                    pending.append(neighbor_index)
        component_index += 1
    return component_by_atom


def _shortest_path_within(
    molecule: Any,
    *,
    allowed_atom_indices: set[int],
    start: int,
    end: int,
) -> Tuple[int, ...]:
    """Return a deterministic shortest path constrained to omitted atoms."""
    pending = [start]
    previous: Dict[int, Optional[int]] = {start: None}
    while pending:
        current = pending.pop(0)
        if current == end:
            path = []
            cursor: Optional[int] = current
            while cursor is not None:
                path.append(cursor)
                cursor = previous[cursor]
            return tuple(reversed(path))
        neighbors = sorted(
            int(value.GetIdx())
            for value in molecule.GetAtomWithIdx(current).GetNeighbors()
            if int(value.GetIdx()) in allowed_atom_indices
        )
        for neighbor in neighbors:
            if neighbor not in previous:
                previous[neighbor] = current
                pending.append(neighbor)
    raise ValueError("shared omitted ports are not connected")


def _fragment_smiles(molecule: Any, atom_indices: Tuple[int, ...]) -> str:
    copy = prepare_fragment_serialization_copy(molecule, atom_indices)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    fragment_smiles = str(
        Chem.MolFragmentToSmiles(
            copy,
            atomsToUse=list(atom_indices),
            canonical=True,
            isomericSmiles=True,
        )
    )
    fragment = Chem.MolFromSmiles(fragment_smiles)
    if fragment is None:
        return fragment_smiles
    return str(
        Chem.MolToSmiles(
            fragment,
            canonical=True,
            isomericSmiles=True,
        )
    )


def _remote_metadata(
    core: Any,
    *,
    side: str,
    component_index: int,
    atom_indices: Tuple[int, ...],
) -> tuple[Optional[str], str, Tuple[int, ...]]:
    omitted = set(atom_indices)
    matches = []
    for subgraph in core.remote_subgraphs:
        if (
            str(subgraph.side) != side
            or int(subgraph.component_index) != component_index
        ):
            continue
        subgraph_atoms = set(subgraph.atom_indices)
        overlap = len(omitted.intersection(subgraph_atoms))
        if overlap:
            matches.append((overlap, omitted == subgraph_atoms, subgraph))
    if not matches:
        return None, "unresolved", ()
    matches.sort(
        key=lambda value: (
            not value[1],
            -value[0],
            str(value[2].subgraph_id),
        )
    )
    _, exact, best = matches[0]
    remote_class = str(best.remote_class) if exact else None
    correspondence_maps = (
        tuple(int(value) for value in best.atom_map_numbers) if exact else ()
    )
    return remote_class, str(best.continuity), correspondence_maps


def _ring_distance(
    molecule: Any,
    *,
    system: set[int],
    start: int,
    end: int,
) -> int:
    pending = [(start, 0)]
    visited = {start}
    while pending:
        current, distance = pending.pop(0)
        if current == end:
            return distance
        atom = molecule.GetAtomWithIdx(current)
        for bond in atom.GetBonds():
            if not bond.GetIsAromatic():
                continue
            neighbor = int(bond.GetOtherAtomIdx(current))
            if neighbor in system and neighbor not in visited:
                visited.add(neighbor)
                pending.append((neighbor, distance + 1))
    raise ValueError("aromatic atoms are not connected within their system")


def _aromatic_relations(
    molecule: Any,
    *,
    attachment_ring_atom_index: int,
    aromatic_systems: set[tuple[int, ...]],
    aromatic_seed_indices: set[int],
) -> Tuple[ReactionDisplayAromaticRelation, ...]:
    containing = next(
        (
            set(system)
            for system in aromatic_systems
            if attachment_ring_atom_index in system
        ),
        None,
    )
    if containing is None:
        return ()
    rings = [set(ring) for ring in molecule.GetRingInfo().AtomRings()]
    simple_six_membered = (
        len(containing) == 6
        and sum(1 for ring in rings if ring == containing) == 1
    )
    relations = []
    for reactive_index in sorted(aromatic_seed_indices.intersection(containing)):
        distance = _ring_distance(
            molecule,
            system=containing,
            start=reactive_index,
            end=attachment_ring_atom_index,
        )
        common_rings = [
            ring
            for ring in rings
            if reactive_index in ring and attachment_ring_atom_index in ring
        ]
        if simple_six_membered and distance in {0, 1, 2, 3}:
            position = ("ipso", "ortho", "meta", "para")[distance]
        else:
            position = "other"
        reactive_atom = molecule.GetAtomWithIdx(reactive_index)
        relations.append(
            ReactionDisplayAromaticRelation(
                reactive_atom_index=reactive_index,
                reactive_atom_map_number=(
                    int(reactive_atom.GetAtomMapNum()) or None
                ),
                attachment_ring_atom_index=attachment_ring_atom_index,
                aromatic_ring_distance=distance,
                same_ring=bool(common_rings),
                positional_relation=position,
            )
        )
    return tuple(relations)


def _substituent_id(draft: _SubstituentDraft) -> str:
    payload = {
        "side": draft.side,
        "component": draft.source_component_index,
        "center": draft.center_atom_index,
        "attachment": draft.attachment_atom_index,
        "atoms": draft.atom_indices,
        "action": draft.boundary_action,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
    ).hexdigest()[:24]
    return f"RDS1:{digest}"


def _display_component_draft(
    molecule: Any,
    *,
    side: Literal["reactant", "product"],
    component_index: int,
    active_atom_indices: set[int],
    aromatic_seed_indices: set[int],
    explicit_remote_atom_indices: set[int],
    core: Any,
) -> _ComponentDraft:
    retained, aromatic_systems, aromatic_policy_atom_indices = (
        _reaction_interface_closure(
            molecule,
            active_atom_indices=active_atom_indices,
            aromatic_seed_indices=aromatic_seed_indices,
            explicit_remote_atom_indices=explicit_remote_atom_indices,
        )
    )
    editable = Chem.RWMol()
    old_to_new: Dict[int, int] = {}
    removed_substituent_count = 0
    r_group_count = 0
    substituents: list[_SubstituentDraft] = []
    aromatic_carbons_to_cap: set[int] = set()
    for atom_index in sorted(retained):
        old_to_new[atom_index] = editable.AddAtom(
            _copy_atom(molecule.GetAtomWithIdx(atom_index))
        )
    for bond in molecule.GetBonds():
        begin = int(bond.GetBeginAtomIdx())
        end = int(bond.GetEndAtomIdx())
        if begin in retained and end in retained:
            editable.AddBond(
                old_to_new[begin], old_to_new[end], bond.GetBondType()
            )

    omitted_components = _omitted_components(molecule, retained)
    component_by_atom = {
        atom_index: component
        for component in omitted_components
        for atom_index in component
    }
    retained_component_by_atom = _connected_component_by_atom(
        molecule, retained
    )
    boundary_bonds: list[tuple[int, int, Any]] = []
    for bond in molecule.GetBonds():
        begin, end = int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())
        if begin in retained and end not in retained:
            boundary_bonds.append((begin, end, bond))
        elif end in retained and begin not in retained:
            boundary_bonds.append((end, begin, bond))
    boundary_retained_components: Dict[Tuple[int, ...], set[int]] = {}
    for retained_index, attachment_index, _ in boundary_bonds:
        omitted_component = component_by_atom[attachment_index]
        boundary_retained_components.setdefault(omitted_component, set()).add(
            retained_component_by_atom[retained_index]
        )
    for retained_index, attachment_index, bond in sorted(
        boundary_bonds,
        key=lambda value: (
            value[0],
            value[1],
        ),
    ):
        atom = molecule.GetAtomWithIdx(retained_index)
        should_cap = (
            retained_index in aromatic_policy_atom_indices
            and atom.GetIsAromatic()
            and atom.GetAtomicNum() == 6
            and len(
                boundary_retained_components[
                    component_by_atom[attachment_index]
                ]
            )
            == 1
        )
        omitted_atoms = component_by_atom[attachment_index]
        attachment_atom = molecule.GetAtomWithIdx(attachment_index)
        atom_maps = tuple(
            sorted(
                int(molecule.GetAtomWithIdx(value).GetAtomMapNum())
                for value in omitted_atoms
                if molecule.GetAtomWithIdx(value).GetAtomMapNum() > 0
            )
        )
        remote_class, continuity, correspondence_maps = _remote_metadata(
            core,
            side=side,
            component_index=component_index,
            atom_indices=omitted_atoms,
        )
        if not atom_maps:
            atom_maps = correspondence_maps
        substituent_profile = build_substituent_profile(
            molecule,
            fragment_atom_indices=omitted_atoms,
            attachment_atom_index=attachment_index,
            core_atom_index=retained_index,
            attachment_bond_order=str(bond.GetBondType()).upper(),
        )
        if remote_class is None:
            remote_class = substituent_profile.base_class
        dummy_atom_index = None
        if should_cap:
            removed_substituent_count += 1
            aromatic_carbons_to_cap.add(retained_index)
            boundary_action = "aromatic_hydrogen_cap"
        else:
            placeholder = Chem.Atom(0)
            placeholder.SetProp("atomLabel", "R")
            placeholder.SetProp("_displayLabel", "R")
            dummy_atom_index = editable.AddAtom(placeholder)
            editable.AddBond(
                old_to_new[retained_index],
                dummy_atom_index,
                bond.GetBondType(),
            )
            r_group_count += 1
            boundary_action = "r_placeholder"
        center_map = int(atom.GetAtomMapNum()) or None
        attachment_map = int(attachment_atom.GetAtomMapNum()) or None
        substituents.append(
            _SubstituentDraft(
                side=side,
                source_component_index=component_index,
                center_atom_index=retained_index,
                center_atom_map_number=center_map,
                center_element=str(atom.GetSymbol()),
                center_hybridization=str(atom.GetHybridization()),
                attachment_atom_index=attachment_index,
                attachment_atom_map_number=attachment_map,
                attachment_element=str(attachment_atom.GetSymbol()),
                attachment_bond_order=str(bond.GetBondType()).upper(),
                atom_indices=omitted_atoms,
                atom_map_numbers=atom_maps,
                fragment_smiles=_fragment_smiles(molecule, omitted_atoms),
                remote_class=remote_class,
                substituent_profile=substituent_profile,
                continuity=continuity,
                boundary_action=boundary_action,
                aromatic_relations=_aromatic_relations(
                    molecule,
                    attachment_ring_atom_index=retained_index,
                    aromatic_systems=aromatic_systems,
                    aromatic_seed_indices=aromatic_seed_indices,
                ),
                dummy_atom_index=dummy_atom_index,
            )
        )

    connector_ports: Dict[Tuple[int, ...], list[_SubstituentDraft]] = {}
    for substituent in substituents:
        if substituent.boundary_action == "r_placeholder":
            connector_ports.setdefault(substituent.atom_indices, []).append(
                substituent
            )
    connectors: list[_ConnectorDraft] = []
    for omitted_atoms, ports in sorted(connector_ports.items()):
        ports_by_retained_component: Dict[int, list[_SubstituentDraft]] = {}
        for port in ports:
            ports_by_retained_component.setdefault(
                retained_component_by_atom[port.center_atom_index], []
            ).append(port)
        if len(ports_by_retained_component) < 2:
            continue
        allowed_atoms = set(omitted_atoms)
        if len(ports_by_retained_component) == 2:
            first_component, second_component = sorted(
                ports_by_retained_component
            )
            candidates = []
            for first in ports_by_retained_component[first_component]:
                for second in ports_by_retained_component[second_component]:
                    path = _shortest_path_within(
                        molecule,
                        allowed_atom_indices=allowed_atoms,
                        start=first.attachment_atom_index,
                        end=second.attachment_atom_index,
                    )
                    candidates.append(
                        (
                            len(path),
                            first.center_atom_index,
                            first.attachment_atom_index,
                            second.center_atom_index,
                            second.attachment_atom_index,
                            first,
                            second,
                            path,
                        )
                    )
            best = min(candidates, key=lambda value: value[:5])
            ordered_ports = (best[5], best[6])
            shortest_path = best[7]
            shortest_path_bond_count: Optional[int] = len(shortest_path) - 1
        else:
            selected_ports = []
            for retained_component in sorted(ports_by_retained_component):
                candidates = []
                for port in ports_by_retained_component[retained_component]:
                    other_ports = (
                        other
                        for other_component, values in (
                            ports_by_retained_component.items()
                        )
                        if other_component != retained_component
                        for other in values
                    )
                    nearest = min(
                        len(
                            _shortest_path_within(
                                molecule,
                                allowed_atom_indices=allowed_atoms,
                                start=port.attachment_atom_index,
                                end=other.attachment_atom_index,
                            )
                        )
                        for other in other_ports
                    )
                    candidates.append(
                        (
                            nearest,
                            port.center_atom_index,
                            port.attachment_atom_index,
                            port,
                        )
                    )
                selected_ports.append(
                    min(candidates, key=lambda value: value[:3])[3]
                )
            ordered_ports = tuple(selected_ports)
            shortest_path = ()
            shortest_path_bond_count = None
        connectors.append(
            _ConnectorDraft(
                side=side,
                source_component_index=component_index,
                ports=ordered_ports,
                atom_indices=omitted_atoms,
                atom_map_numbers=ordered_ports[0].atom_map_numbers,
                fragment_smiles=ordered_ports[0].fragment_smiles,
                shortest_path_atom_indices=shortest_path,
                shortest_path_bond_count=shortest_path_bond_count,
            )
        )

    for atom_index in aromatic_carbons_to_cap:
        atom = editable.GetAtomWithIdx(old_to_new[atom_index])
        if atom.GetFormalCharge() == 0 and atom.GetNumRadicalElectrons() == 0:
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(False)

    display_molecule = editable.GetMol()
    try:
        Chem.SanitizeMol(display_molecule)
    except Exception as exc:
        raise ValueError(
            "display-minimized component could not be sanitized"
        ) from exc
    return _ComponentDraft(
        side=side,
        source_component_index=component_index,
        molecule=display_molecule,
        retained_atom_indices=tuple(sorted(retained)),
        active_atom_indices=tuple(sorted(active_atom_indices)),
        retained_aromatic_system_count=len(aromatic_systems),
        removed_substituent_count=removed_substituent_count,
        r_group_count=r_group_count,
        substituents=substituents,
        connectors=connectors,
    )


def _side_component_drafts(
    analysis: ReactionRenderContext,
    *,
    side: Literal["reactant", "product"],
) -> list[_ComponentDraft]:
    components = analysis.reactants if side == "reactant" else analysis.products
    active, aromatic_seeds = _active_coordinates(analysis, side=side)
    core = analysis.reaction_core
    explicit_remote: Dict[int, set[int]] = {}
    if core is not None:
        for subgraph in core.remote_subgraphs:
            continuity = str(subgraph.continuity)
            preserve_small_unresolved_heteroatom = (
                continuity == "unresolved"
                and str(subgraph.remote_class) == "heteroatom"
                and int(subgraph.fragment_heavy_atom_count) <= 1
            )
            if str(subgraph.side) == side and (
                continuity in {"departing", "appearing"}
                or preserve_small_unresolved_heteroatom
            ):
                explicit_remote.setdefault(
                    int(subgraph.component_index), set()
                ).update(int(value) for value in subgraph.atom_indices)
    output: list[_ComponentDraft] = []
    for component_index, active_indices in active.items():
        if component_index < 0 or component_index >= len(components):
            raise ValueError("reaction core references an unknown component")
        molecule = parse_smiles(components[component_index].input_smiles)
        if molecule is None:
            raise ValueError("reaction display component could not be reparsed")
        output.append(
            _display_component_draft(
                molecule,
                side=side,
                component_index=component_index,
                active_atom_indices=active_indices,
                aromatic_seed_indices=aromatic_seeds.get(
                    component_index, set()
                ),
                explicit_remote_atom_indices=explicit_remote.get(
                    component_index, set()
                ),
                core=core,
            )
        )
    return output


_SUPERSCRIPT_DIGITS = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def _placeholder_label(index: int, count: int) -> str:
    if count == 1:
        return "R"
    return f"R{str(index).translate(_SUPERSCRIPT_DIGITS)}"


def _placeholder_base_identity(draft: _SubstituentDraft) -> tuple[Any, ...]:
    return (
        "fragment",
        draft.fragment_smiles,
        draft.remote_class,
        draft.attachment_element,
        draft.center_element,
        draft.attachment_bond_order,
    )


def _assign_placeholder_labels(components: list[_ComponentDraft]) -> None:
    component_by_substituent = {
        id(substituent): component
        for component in components
        for substituent in component.substituents
    }
    placeholders = [
        substituent
        for component in components
        for substituent in component.substituents
        if substituent.boundary_action == "r_placeholder"
    ]
    grouped: Dict[tuple[str, tuple[Any, ...]], list[_SubstituentDraft]] = {}
    for draft in placeholders:
        grouped.setdefault(
            (draft.side, _placeholder_base_identity(draft)), []
        ).append(draft)
    instance_keys: Dict[int, tuple[Any, ...]] = {}
    for (side, base), values in grouped.items():
        del side
        ordered = sorted(
            values,
            key=lambda value: (
                value.source_component_index,
                value.center_atom_map_number or 0,
                value.center_atom_index,
                value.attachment_atom_map_number or 0,
                value.attachment_atom_index,
            ),
        )
        for ordinal, draft in enumerate(ordered, start=1):
            instance_keys[id(draft)] = (base, ordinal)
    ordered_keys = sorted(set(instance_keys.values()), key=repr)
    index_by_key = {
        identity: index for index, identity in enumerate(ordered_keys, start=1)
    }
    for draft in placeholders:
        index = index_by_key[instance_keys[id(draft)]]
        draft.placeholder_index = index
        draft.display_label = _placeholder_label(index, len(ordered_keys))
        if draft.dummy_atom_index is None:
            raise ValueError("R-group placeholder is missing its dummy atom")
        component = component_by_substituent[id(draft)]
        atom = component.molecule.GetAtomWithIdx(draft.dummy_atom_index)
        atom.SetAtomMapNum(index)
        atom.SetProp("atomLabel", draft.display_label)
        atom.SetProp("_displayLabel", draft.display_label)


def _connector_base_identity(draft: _ConnectorDraft) -> tuple[Any, ...]:
    placeholder_indices = tuple(
        sorted(
            int(port.placeholder_index)
            for port in draft.ports
            if port.placeholder_index is not None
        )
    )
    if len(placeholder_indices) != len(draft.ports):
        raise ValueError("hidden connector port is missing an R-group label")
    if draft.atom_map_numbers:
        return (
            "mapped_connector",
            draft.atom_map_numbers,
            placeholder_indices,
        )
    return (
        "connector",
        draft.fragment_smiles,
        placeholder_indices,
    )


def _assign_connector_labels(components: list[_ComponentDraft]) -> None:
    connectors = [
        connector
        for component in components
        for connector in component.connectors
    ]
    for draft in connectors:
        original_first = draft.ports[0]
        draft.ports = tuple(
            sorted(
                draft.ports,
                key=lambda port: int(port.placeholder_index or 0),
            )
        )
        if (
            len(draft.ports) == 2
            and draft.ports[0] is not original_first
            and draft.shortest_path_atom_indices
        ):
            draft.shortest_path_atom_indices = tuple(
                reversed(draft.shortest_path_atom_indices)
            )
    grouped: Dict[tuple[str, tuple[Any, ...]], list[_ConnectorDraft]] = {}
    for draft in connectors:
        grouped.setdefault(
            (draft.side, _connector_base_identity(draft)), []
        ).append(draft)
    instance_keys: Dict[int, tuple[Any, ...]] = {}
    for (side, base), values in grouped.items():
        del side
        ordered = sorted(
            values,
            key=lambda value: (
                value.source_component_index,
                value.atom_indices,
            ),
        )
        for ordinal, draft in enumerate(ordered, start=1):
            instance_keys[id(draft)] = (base, ordinal)
    ordered_keys = sorted(set(instance_keys.values()), key=repr)
    index_by_key = {
        identity: index for index, identity in enumerate(ordered_keys, start=1)
    }
    for draft in connectors:
        index = index_by_key[instance_keys[id(draft)]]
        draft.connector_index = index
        draft.display_label = (
            f"S{str(index).translate(_SUPERSCRIPT_DIGITS)}"
        )


def _public_substituent(
    draft: _SubstituentDraft,
) -> ReactionDisplaySubstituent:
    return ReactionDisplaySubstituent(
        substituent_id=_substituent_id(draft),
        side=draft.side,
        source_component_index=draft.source_component_index,
        center_atom_index=draft.center_atom_index,
        center_atom_map_number=draft.center_atom_map_number,
        center_element=draft.center_element,
        center_hybridization=draft.center_hybridization,
        attachment_atom_index=draft.attachment_atom_index,
        attachment_atom_map_number=draft.attachment_atom_map_number,
        attachment_element=draft.attachment_element,
        attachment_bond_order=draft.attachment_bond_order,
        atom_indices=draft.atom_indices,
        atom_map_numbers=draft.atom_map_numbers,
        fragment_smiles=draft.fragment_smiles,
        remote_class=draft.remote_class,
        substituent_profile=draft.substituent_profile,
        continuity=draft.continuity,
        boundary_action=draft.boundary_action,
        placeholder_index=draft.placeholder_index,
        display_label=draft.display_label,
        aromatic_relations=draft.aromatic_relations,
    )


def _public_connector(draft: _ConnectorDraft) -> ReactionDisplayConnector:
    if draft.connector_index is None or draft.display_label is None:
        raise ValueError("hidden connector has not been labeled")
    placeholder_indices = tuple(
        int(port.placeholder_index)
        for port in draft.ports
        if port.placeholder_index is not None
    )
    port_display_labels = tuple(
        str(port.display_label)
        for port in draft.ports
        if port.display_label is not None
    )
    if (
        len(placeholder_indices) != len(draft.ports)
        or len(port_display_labels) != len(draft.ports)
    ):
        raise ValueError("hidden connector has an unlabeled display port")
    return ReactionDisplayConnector(
        connector_id=f"S{draft.connector_index}",
        display_label=draft.display_label,
        side=draft.side,
        source_component_index=draft.source_component_index,
        port_substituent_ids=tuple(
            _substituent_id(port) for port in draft.ports
        ),
        placeholder_indices=placeholder_indices,
        port_display_labels=port_display_labels,
        attachment_atom_indices=tuple(
            port.attachment_atom_index for port in draft.ports
        ),
        atom_indices=draft.atom_indices,
        atom_map_numbers=draft.atom_map_numbers,
        fragment_smiles=draft.fragment_smiles,
        hidden_atom_count=len(draft.atom_indices),
        shortest_path_atom_indices=draft.shortest_path_atom_indices,
        shortest_path_bond_count=draft.shortest_path_bond_count,
    )


def _finalize_component(draft: _ComponentDraft) -> ReactionDisplayComponent:
    render_smiles = str(
        Chem.MolToSmiles(
            draft.molecule,
            canonical=True,
            isomericSmiles=True,
        )
    )
    generic = Chem.Mol(draft.molecule)
    for atom in generic.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomMapNum(0)
    display_smiles = str(
        Chem.MolToSmiles(generic, canonical=True, isomericSmiles=True)
    )
    substituents = tuple(
        sorted(
            (_public_substituent(value) for value in draft.substituents),
            key=lambda value: value.substituent_id,
        )
    )
    connectors = tuple(
        sorted(
            (_public_connector(value) for value in draft.connectors),
            key=lambda value: value.connector_id,
        )
    )
    return ReactionDisplayComponent(
        side=draft.side,
        source_component_index=draft.source_component_index,
        display_smiles=display_smiles,
        render_smiles=render_smiles,
        retained_atom_indices=draft.retained_atom_indices,
        active_atom_indices=draft.active_atom_indices,
        retained_aromatic_system_count=draft.retained_aromatic_system_count,
        removed_substituent_count=draft.removed_substituent_count,
        r_group_count=draft.r_group_count,
        substituents=substituents,
        connectors=connectors,
    )


def _finalize_components(
    drafts: list[_ComponentDraft],
) -> Tuple[ReactionDisplayComponent, ...]:
    return tuple(
        sorted(
            (_finalize_component(value) for value in drafts),
            key=lambda value: (
                value.display_smiles,
                value.source_component_index,
            ),
        )
    )


def _intramolecular_annotation(
    topology: Any,
    definition: Mapping[str, Any],
) -> str | None:
    if topology is None or str(topology.reaction_scope) != "intramolecular":
        return None
    ring_sizes = _formed_ring_sizes(topology)
    if len(ring_sizes) == 1:
        return str(definition["intramolecular_note_template"]).format(
            ring_size=ring_sizes[0]
        )
    if ring_sizes:
        sizes = ", ".join(str(value) for value in ring_sizes)
        return f"Intramolecular; forms rings of size {sizes}"
    return str(definition["intramolecular_note_without_size"])


def _formed_ring_sizes(topology: Any) -> Tuple[int, ...]:
    """Prefer explicit formed-cycle evidence over distance approximations."""
    if topology is None:
        return ()
    explicit = tuple(
        int(change.ring_size)
        for change in tuple(getattr(topology, "ring_changes", ()) or ())
        if str(change.change_type) == "formed"
    )
    if explicit:
        return tuple(sorted(set(explicit)))
    return tuple(
        sorted(
            set(
                int(value)
                for value in tuple(
                    getattr(topology, "formed_ring_sizes", ()) or ()
                )
            )
        )
    )


def build_reaction_display_projection(
    analysis: ReactionRenderContext,
) -> ReactionDisplayProjection:
    """Build a display-only minimum reaction SMILES from a reaction core."""
    if not isinstance(analysis, ReactionRenderContext):
        raise TypeError(
            "build_reaction_display_projection requires a ReactionRenderContext"
        )
    _require_rdkit()
    core = analysis.reaction_core
    if core is None:
        raise ValueError(
            "Reaction display minimization requires a ReactionCoreProjection."
        )
    definition = load_reaction_display_projection_definition()
    reactant_drafts = _side_component_drafts(analysis, side="reactant")
    product_drafts = _side_component_drafts(analysis, side="product")
    _assign_placeholder_labels(reactant_drafts + product_drafts)
    _assign_connector_labels(reactant_drafts + product_drafts)
    reactants = _finalize_components(reactant_drafts)
    products = _finalize_components(product_drafts)
    if not reactants or not products:
        raise ValueError("reaction display projection has no drawable reaction")
    reactant_smiles = ".".join(
        component.display_smiles for component in reactants
    )
    product_smiles = ".".join(
        component.display_smiles for component in products
    )
    render_reactant_smiles = ".".join(
        component.render_smiles for component in reactants
    )
    render_product_smiles = ".".join(
        component.render_smiles for component in products
    )
    substituents = tuple(
        sorted(
            (
                substituent
                for component in reactants + products
                for substituent in component.substituents
            ),
            key=lambda value: (
                value.side,
                value.source_component_index,
                value.substituent_id,
            ),
        )
    )
    connectors = tuple(
        sorted(
            (
                connector
                for component in reactants + products
                for connector in component.connectors
            ),
            key=lambda value: (
                value.side,
                value.source_component_index,
                value.connector_id,
            ),
        )
    )
    topology = analysis.reaction_topology
    annotation = _intramolecular_annotation(topology, definition)
    warnings = list(core.warnings)
    warnings.append("DISPLAY_PROJECTION_NOT_REACTION_IDENTITY")
    if connectors:
        warnings.append("HIDDEN_CONNECTOR_ABSTRACTED_IN_DISPLAY")
    if annotation is not None:
        warnings.append("INTRAMOLECULAR_TETHER_ABSTRACTED_IN_DISPLAY")
    return ReactionDisplayProjection(
        minimum_reaction_smiles=f"{reactant_smiles}>>{product_smiles}",
        render_reaction_smiles=(
            f"{render_reactant_smiles}>>{render_product_smiles}"
        ),
        reactants=reactants,
        products=products,
        substituents=substituents,
        connectors=connectors,
        reaction_scope=(
            str(topology.reaction_scope) if topology is not None else "unresolved"
        ),
        formed_ring_sizes=_formed_ring_sizes(topology),
        annotation=annotation,
        evidence_status=str(core.evidence_status),
        confidence=float(core.confidence),
        warnings=tuple(sorted(set(warnings))),
        definition_id=str(definition["definition_id"]),
        schema_version=str(definition["schema_version"]),
    )


__all__ = [
    "ReactionDisplayAromaticRelation",
    "ReactionDisplayComponent",
    "ReactionDisplayConnector",
    "ReactionDisplayProjection",
    "ReactionDisplaySubstituent",
    "build_reaction_display_projection",
    "load_reaction_display_projection_definition",
]
