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

from .chemistry.rdkit_utils import parse_smiles, rdkit_available
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
    continuity: str
    boundary_action: Literal["r_placeholder", "aromatic_hydrogen_cap"]
    placeholder_index: Optional[int]
    display_label: Optional[str]
    aromatic_relations: Tuple[ReactionDisplayAromaticRelation, ...]


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


@dataclass(frozen=True)
class ReactionDisplayProjection:
    """Deterministic minimum reaction SMILES intended only for rendering."""

    minimum_reaction_smiles: str
    render_reaction_smiles: str
    reactants: Tuple[ReactionDisplayComponent, ...]
    products: Tuple[ReactionDisplayComponent, ...]
    substituents: Tuple[ReactionDisplaySubstituent, ...]
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
        "definition_id": "reaction_display_projection.v1.2",
        "aromatic_system_policy": "retain_aromatic_bond_component",
        "multiple_bond_policy": "retain_contiguous_multiple_bond_unit",
        "active_heteroatom_shell_policy": (
            "retain_direct_noncarbon_neighbors"
        ),
        "saturated_carbon_policy": "retain_active_atom_only",
        "aromatic_carbon_boundary_policy": "remove_and_hydrogen_cap",
        "aromatic_heteroatom_boundary_policy": "retain_R_attachment",
        "nonaromatic_boundary_policy": "retain_R_attachment",
        "placeholder_smiles": "*",
        "placeholder_label_template": "R{index}",
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
    continuity: str
    boundary_action: Literal["r_placeholder", "aromatic_hydrogen_cap"]
    aromatic_relations: Tuple[ReactionDisplayAromaticRelation, ...]
    dummy_atom_index: Optional[int]
    placeholder_index: Optional[int] = None
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

    pending = [
        value
        for value in sorted(active_atom_indices)
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


def _fragment_smiles(molecule: Any, atom_indices: Tuple[int, ...]) -> str:
    copy = Chem.Mol(molecule)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return str(
        Chem.MolFragmentToSmiles(
            copy,
            atomsToUse=list(atom_indices),
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
) -> tuple[Optional[str], str]:
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
        return None, "unresolved"
    matches.sort(
        key=lambda value: (
            not value[1],
            -value[0],
            str(value[2].subgraph_id),
        )
    )
    _, exact, best = matches[0]
    remote_class = str(best.remote_class) if exact else None
    return remote_class, str(best.continuity)


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
    boundary_bonds: list[tuple[int, int, Any]] = []
    for bond in molecule.GetBonds():
        begin, end = int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())
        if begin in retained and end not in retained:
            boundary_bonds.append((begin, end, bond))
        elif end in retained and begin not in retained:
            boundary_bonds.append((end, begin, bond))
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
        remote_class, continuity = _remote_metadata(
            core,
            side=side,
            component_index=component_index,
            atom_indices=omitted_atoms,
        )
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
    if draft.atom_map_numbers:
        return ("mapped_fragment", draft.atom_map_numbers)
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
        continuity=draft.continuity,
        boundary_action=draft.boundary_action,
        placeholder_index=draft.placeholder_index,
        display_label=draft.display_label,
        aromatic_relations=draft.aromatic_relations,
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
    topology = analysis.reaction_topology
    annotation = _intramolecular_annotation(topology, definition)
    warnings = list(core.warnings)
    warnings.append("DISPLAY_PROJECTION_NOT_REACTION_IDENTITY")
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
    "ReactionDisplayProjection",
    "ReactionDisplaySubstituent",
    "build_reaction_display_projection",
    "load_reaction_display_projection_definition",
]
