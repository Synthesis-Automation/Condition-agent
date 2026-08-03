"""Build a display-only minimum reaction graph from observed reaction edits.

The projection is deliberately separate from :class:`ReactionCoreProjection`.
It is a chemist-facing drawing aid and must not be used as reaction identity or
as an input to recommendation, admission, or retrieval.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Literal, Mapping, Tuple

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
class ReactionDisplayComponent:
    """One source component after display-only graph abstraction."""

    side: Literal["reactant", "product"]
    source_component_index: int
    display_smiles: str
    retained_atom_indices: Tuple[int, ...]
    active_atom_indices: Tuple[int, ...]
    retained_aromatic_system_count: int
    removed_substituent_count: int
    r_group_count: int


@dataclass(frozen=True)
class ReactionDisplayProjection:
    """Deterministic minimum reaction SMILES intended only for rendering."""

    minimum_reaction_smiles: str
    reactants: Tuple[ReactionDisplayComponent, ...]
    products: Tuple[ReactionDisplayComponent, ...]
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
        "definition_id": "reaction_display_projection.v1.0",
        "aromatic_system_policy": "retain_aromatic_bond_component",
        "aromatic_carbon_boundary_policy": "remove_and_hydrogen_cap",
        "aromatic_heteroatom_boundary_policy": "retain_R_attachment",
        "nonaromatic_boundary_policy": "retain_R_attachment",
        "placeholder_smiles": "*",
        "placeholder_label": "R",
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


def _display_component(
    molecule: Any,
    *,
    side: Literal["reactant", "product"],
    component_index: int,
    active_atom_indices: set[int],
    aromatic_seed_indices: set[int],
    explicit_remote_atom_indices: set[int],
) -> ReactionDisplayComponent:
    retained = set(active_atom_indices).union(explicit_remote_atom_indices)
    aromatic_systems: set[tuple[int, ...]] = set()
    aromatic_policy_atom_indices: set[int] = set()
    for atom_index in sorted(aromatic_seed_indices):
        system = _aromatic_system(molecule, atom_index)
        if system:
            retained.update(system)
            aromatic_policy_atom_indices.update(system)
            aromatic_systems.add(tuple(sorted(system)))

    editable = Chem.RWMol()
    old_to_new: Dict[int, int] = {}
    removed_substituent_count = 0
    r_group_count = 0
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

    boundary_bonds = []
    for bond in molecule.GetBonds():
        begin = int(bond.GetBeginAtomIdx())
        end = int(bond.GetEndAtomIdx())
        if (begin in retained) == (end in retained):
            continue
        retained_index = begin if begin in retained else end
        boundary_bonds.append((retained_index, bond))
    for retained_index, bond in sorted(
        boundary_bonds,
        key=lambda value: (
            value[0],
            min(value[1].GetBeginAtomIdx(), value[1].GetEndAtomIdx()),
            max(value[1].GetBeginAtomIdx(), value[1].GetEndAtomIdx()),
        ),
    ):
        atom = molecule.GetAtomWithIdx(retained_index)
        if (
            retained_index in aromatic_policy_atom_indices
            and atom.GetIsAromatic()
            and atom.GetAtomicNum() == 6
        ):
            removed_substituent_count += 1
            aromatic_carbons_to_cap.add(retained_index)
            continue
        placeholder = Chem.Atom(0)
        placeholder.SetProp("atomLabel", "R")
        placeholder.SetProp("_displayLabel", "R")
        placeholder_index = editable.AddAtom(placeholder)
        editable.AddBond(
            old_to_new[retained_index],
            placeholder_index,
            bond.GetBondType(),
        )
        r_group_count += 1

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
    display_smiles = str(
        Chem.MolToSmiles(
            display_molecule,
            canonical=True,
            isomericSmiles=True,
        )
    )
    return ReactionDisplayComponent(
        side=side,
        source_component_index=component_index,
        display_smiles=display_smiles,
        retained_atom_indices=tuple(sorted(retained)),
        active_atom_indices=tuple(sorted(active_atom_indices)),
        retained_aromatic_system_count=len(aromatic_systems),
        removed_substituent_count=removed_substituent_count,
        r_group_count=r_group_count,
    )


def _side_components(
    analysis: ReactionRenderContext,
    *,
    side: Literal["reactant", "product"],
) -> Tuple[ReactionDisplayComponent, ...]:
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
    output = []
    for component_index, active_indices in active.items():
        if component_index < 0 or component_index >= len(components):
            raise ValueError("reaction core references an unknown component")
        molecule = parse_smiles(components[component_index].input_smiles)
        if molecule is None:
            raise ValueError("reaction display component could not be reparsed")
        output.append(
            _display_component(
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
            )
        )
    return tuple(
        sorted(
            output,
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
    reactants = _side_components(analysis, side="reactant")
    products = _side_components(analysis, side="product")
    if not reactants or not products:
        raise ValueError("reaction display projection has no drawable reaction")
    reactant_smiles = ".".join(
        component.display_smiles for component in reactants
    )
    product_smiles = ".".join(
        component.display_smiles for component in products
    )
    topology = analysis.reaction_topology
    annotation = _intramolecular_annotation(topology, definition)
    warnings = list(core.warnings)
    warnings.append("DISPLAY_PROJECTION_NOT_REACTION_IDENTITY")
    if annotation is not None:
        warnings.append("INTRAMOLECULAR_TETHER_ABSTRACTED_IN_DISPLAY")
    return ReactionDisplayProjection(
        minimum_reaction_smiles=f"{reactant_smiles}>>{product_smiles}",
        reactants=reactants,
        products=products,
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
    "ReactionDisplayComponent",
    "ReactionDisplayProjection",
    "build_reaction_display_projection",
    "load_reaction_display_projection_definition",
]
