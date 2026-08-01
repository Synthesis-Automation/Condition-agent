"""Render a mapped reaction core as a compact, scaffold-abstracted scheme."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple

from reactive_taxonomy.chemistry.rdkit_utils import parse_smiles, rdkit_available

from .rendering import apply_render_preset

try:  # pragma: no cover - public helpers exercise availability
    from rdkit import Chem
    from rdkit.Chem import Draw, rdChemReactions
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError:  # pragma: no cover
    Chem = None  # type: ignore[assignment]
    Draw = None  # type: ignore[assignment]
    rdChemReactions = None  # type: ignore[assignment]
    rdMolDraw2D = None  # type: ignore[assignment]


_DEFINITION_PATH = (
    Path(__file__).with_name("definitions")
    / "reaction_core_graphic.v1.json"
)
_SUPPORTED_FORMATS = frozenset({"png", "svg"})
_BOND_TYPES = {
    "SINGLE": "SINGLE",
    "DOUBLE": "DOUBLE",
    "TRIPLE": "TRIPLE",
    "AROMATIC": "AROMATIC",
}


@dataclass(frozen=True)
class ReactionCoreGraphicPlaceholder:
    """One exact remote subgraph represented by a short drawing label."""

    label: str
    remote_class: str
    fragment_smiles: str
    attachment_port_count: int


@dataclass(frozen=True)
class ReactionCoreGraphic:
    """Rendered minimized scheme plus the information hidden by placeholders."""

    image_bytes: bytes
    image_format: str
    core_id: str
    evidence_status: str
    confidence: float
    placeholders: Tuple[ReactionCoreGraphicPlaceholder, ...]
    warnings: Tuple[str, ...]
    definition_id: str
    schema_version: str


@dataclass(frozen=True)
class _MultisiteScaffoldCollapse:
    """One retained ring/scaffold shared by multiple active center atoms."""

    identity: tuple[Any, ...]
    side: str
    component_index: int
    core_atom_indices: Tuple[int, ...]
    subgraph_ids: Tuple[str, ...]
    remote_class: str
    fragment_smiles: str
    attachment_port_count: int


@lru_cache(maxsize=1)
def load_reaction_core_graphic_definition() -> Dict[str, Any]:
    """Load and validate the versioned placeholder-rendering definition."""
    with _DEFINITION_PATH.open("r", encoding="utf-8") as handle:
        definition = dict(json.load(handle))
    if str(definition.get("schema_version") or "") != "1.4":
        raise ValueError("unsupported reaction-core graphic schema")
    if str(definition.get("definition_id") or "") != (
        "reaction_core_graphic.v1.4"
    ):
        raise ValueError("unexpected reaction-core graphic definition ID")
    labels = definition.get("remote_class_labels")
    if not isinstance(labels, Mapping) or not labels:
        raise ValueError("reaction-core graphic requires remote-class labels")
    continuities = tuple(
        str(value)
        for value in definition.get("continuities_replaced_by_labels") or ()
    )
    if continuities != ("retained",):
        raise ValueError(
            "reaction-core graphic may abstract only retained subgraphs"
        )
    if definition.get("collapse_retained_multisite_scaffolds") is not True:
        raise ValueError(
            "reaction-core graphic must collapse retained multi-site scaffolds"
        )
    if definition.get("preserve_intramolecular_tethers") is not True:
        raise ValueError(
            "reaction-core graphic must preserve intramolecular tethers"
        )
    if definition.get("render_nonretained_subgraphs_explicitly") is not True:
        raise ValueError(
            "reaction-core graphic must render non-retained subgraphs"
        )
    explicit = definition.get("explicit_retained_subgraphs")
    if not isinstance(explicit, Mapping):
        raise ValueError(
            "reaction-core graphic requires an explicit-fragment policy"
        )
    max_heavy_atoms = explicit.get("max_heavy_atom_count")
    if not isinstance(max_heavy_atoms, int) or max_heavy_atoms < 0:
        raise ValueError(
            "explicit-fragment maximum heavy-atom count must be non-negative"
        )
    explicit_classes = explicit.get("remote_classes")
    if (
        not isinstance(explicit_classes, list)
        or not explicit_classes
        or not all(isinstance(value, str) and value for value in explicit_classes)
    ):
        raise ValueError(
            "explicit-fragment remote classes must be a non-empty string list"
        )
    template = str(definition.get("indexed_label_template") or "")
    if "{label}" not in template or "{index}" not in template:
        raise ValueError("reaction-core graphic requires an indexed label template")
    return definition


def _require_rdkit() -> None:
    if (
        not rdkit_available()
        or Chem is None
        or Draw is None
        or rdChemReactions is None
        or rdMolDraw2D is None
    ):
        raise RuntimeError(
            "RDKit is required for reaction-core graphic rendering."
        )


def _remote_identity(subgraph: Any) -> tuple[Any, ...]:
    """Pair the same retained graph across sides without using row order."""
    mapped_atoms = tuple(sorted(int(value) for value in subgraph.atom_map_numbers))
    core_maps = tuple(
        sorted(
            int(port.core_atom_map_number)
            for port in subgraph.attachment_ports
            if port.core_atom_map_number is not None
        )
    )
    if mapped_atoms:
        attachment_identity: object = ("mapped_atoms", mapped_atoms)
    else:
        attachment_identity = (
            "core_maps",
            core_maps,
            str(subgraph.fragment_smiles),
        )
    return (
        str(subgraph.remote_class),
        attachment_identity,
    )


def _render_remote_explicitly(subgraph: Any) -> bool:
    """Return whether a retained fragment must remain chemically visible."""
    definition = load_reaction_core_graphic_definition()
    policy = definition["explicit_retained_subgraphs"]
    return (
        subgraph.continuity == "retained"
        and int(subgraph.fragment_heavy_atom_count)
        <= int(policy["max_heavy_atom_count"])
        and str(subgraph.remote_class) in set(policy["remote_classes"])
    )


def _topology_protected_subgraph_ids(
    analysis: Any,
    *,
    side: str | None = None,
) -> set[str]:
    """Return retained tethers that cannot be safely contracted to dummies.

    A connected remote graph spanning multiple active atoms carries the path
    closed by an intramolecular formed bond. Replacing each attachment port by
    an unrelated dummy atom destroys that path and therefore the observed ring
    topology. Keep such graphs explicit whenever the analysis establishes a
    ring-forming intramolecular event.
    """
    topology = getattr(analysis, "reaction_topology", None)
    if topology is None or str(topology.reaction_scope) not in {
        "intramolecular",
        "mixed",
    }:
        return set()
    ring_sizes = tuple(getattr(topology, "formed_ring_sizes", ()) or ())
    ring_count_delta = getattr(topology, "ring_count_delta", None)
    if not ring_sizes and not (
        ring_count_delta is not None and int(ring_count_delta) > 0
    ):
        return set()
    core = getattr(analysis, "reaction_core", None)
    if core is None:
        return set()
    protected = set()
    for subgraph in core.remote_subgraphs:
        if (
            subgraph.continuity != "retained"
            or (side is not None and str(subgraph.side) != side)
        ):
            continue
        boundary = {
            (int(port.core_component_index), int(port.core_atom_index))
            for port in subgraph.attachment_ports
        }
        if len(boundary) >= 2:
            protected.add(str(subgraph.subgraph_id))
    return protected


def _multisite_scaffold_collapses(
    analysis: Any,
) -> Tuple[_MultisiteScaffoldCollapse, ...]:
    """Recognize one retained ring remainder split by several active sites."""
    from rdkit import Chem

    core = analysis.reaction_core
    topology_protected = _topology_protected_subgraph_ids(analysis)
    values = []
    for side, components in (
        ("reactant", analysis.reactants),
        ("product", analysis.products),
    ):
        component_by_index = {
            component.component_index: component for component in components
        }
        grouped: Dict[tuple[int, Tuple[int, ...]], list[Any]] = {}
        for subgraph in core.remote_subgraphs:
            if (
                subgraph.side == side
                and subgraph.continuity == "retained"
                and not _render_remote_explicitly(subgraph)
                and str(subgraph.subgraph_id) not in topology_protected
            ):
                core_atom_indices = tuple(
                    sorted(
                        {
                            int(port.core_atom_index)
                            for port in subgraph.attachment_ports
                            if port.core_component_index
                            == subgraph.component_index
                        }
                    )
                )
                if len(core_atom_indices) >= 2:
                    grouped.setdefault(
                        (subgraph.component_index, core_atom_indices),
                        [],
                    ).append(subgraph)
        for (component_index, core_atom_indices), subgraphs in grouped.items():
            if len(subgraphs) < 2:
                continue
            component = component_by_index.get(component_index)
            molecule = (
                parse_smiles(component.input_smiles)
                if component is not None
                else None
            )
            if molecule is None or not all(
                molecule.GetAtomWithIdx(atom_index).IsInRing()
                for atom_index in core_atom_indices
            ):
                continue
            scaffold_atom_indices = tuple(
                sorted(
                    set(core_atom_indices).union(
                        atom_index
                        for subgraph in subgraphs
                        for atom_index in subgraph.atom_indices
                    )
                )
            )
            mapped_atoms = tuple(
                sorted(
                    {
                        int(map_number)
                        for subgraph in subgraphs
                        for map_number in (
                            *subgraph.atom_map_numbers,
                            *(
                                port.core_atom_map_number
                                for port in subgraph.attachment_ports
                                if port.core_atom_index in core_atom_indices
                            ),
                        )
                        if map_number is not None and int(map_number) > 0
                    }
                )
            )
            if not mapped_atoms:
                continue
            remote_classes = {str(value.remote_class) for value in subgraphs}
            remote_class = (
                next(iter(remote_classes))
                if len(remote_classes) == 1
                else "generic_R"
            )
            copied = Chem.Mol(molecule)
            for atom in copied.GetAtoms():
                atom.SetAtomMapNum(0)
            fragment_smiles = str(
                Chem.MolFragmentToSmiles(
                    copied,
                    atomsToUse=list(scaffold_atom_indices),
                    canonical=True,
                    isomericSmiles=True,
                )
            )
            active_coordinates = _active_coordinates(core, side)
            external_neighbors = {
                int(neighbor.GetIdx())
                for atom_index in core_atom_indices
                for neighbor in molecule.GetAtomWithIdx(atom_index).GetNeighbors()
                if (
                    component_index,
                    int(neighbor.GetIdx()),
                )
                in active_coordinates
                and int(neighbor.GetIdx()) not in scaffold_atom_indices
            }
            identity = (
                remote_class,
                ("multisite_mapped_atoms", mapped_atoms),
            )
            values.append(
                _MultisiteScaffoldCollapse(
                    identity=identity,
                    side=side,
                    component_index=component_index,
                    core_atom_indices=core_atom_indices,
                    subgraph_ids=tuple(
                        sorted(str(value.subgraph_id) for value in subgraphs)
                    ),
                    remote_class=remote_class,
                    fragment_smiles=fragment_smiles,
                    attachment_port_count=len(external_neighbors),
                )
            )
    return tuple(
        sorted(
            values,
            key=lambda value: (
                value.side,
                value.component_index,
                value.identity,
            ),
        )
    )


def _placeholder_assignments(
    analysis: Any,
) -> tuple[
    Dict[tuple[Any, ...], str],
    Tuple[ReactionCoreGraphicPlaceholder, ...],
    Tuple[_MultisiteScaffoldCollapse, ...],
]:
    core = analysis.reaction_core
    definition = load_reaction_core_graphic_definition()
    labels = definition["remote_class_labels"]
    scaffold_collapses = _multisite_scaffold_collapses(analysis)
    topology_protected = _topology_protected_subgraph_ids(analysis)
    collapsed_subgraph_ids = {
        subgraph_id
        for collapse in scaffold_collapses
        for subgraph_id in collapse.subgraph_ids
    }
    retained = tuple(
        subgraph
        for subgraph in core.remote_subgraphs
        if (
            subgraph.continuity == "retained"
            and not _render_remote_explicitly(subgraph)
            and str(subgraph.subgraph_id) not in topology_protected
            and str(subgraph.subgraph_id) not in collapsed_subgraph_ids
        )
    )
    representative_by_identity: Dict[tuple[Any, ...], tuple[Any, ...]] = {}
    for subgraph in retained:
        representative_by_identity.setdefault(
            _remote_identity(subgraph),
            (
                str(subgraph.remote_class),
                str(subgraph.fragment_smiles),
                len(subgraph.attachment_ports),
            ),
        )
    for collapse in scaffold_collapses:
        representative_by_identity.setdefault(
            collapse.identity,
            (
                collapse.remote_class,
                collapse.fragment_smiles,
                collapse.attachment_port_count,
            ),
        )
    identities_by_base: Dict[str, list[tuple[Any, ...]]] = {}
    for identity, record in representative_by_identity.items():
        base = str(labels.get(record[0]) or "R")
        identities_by_base.setdefault(base, []).append(identity)
    assignments: Dict[tuple[Any, ...], str] = {}
    template = str(definition["indexed_label_template"])
    for base, identities in sorted(identities_by_base.items()):
        ordered = sorted(identities, key=repr)
        for index, identity in enumerate(ordered, start=1):
            assignments[identity] = (
                template.format(label=base, index=index)
                if len(ordered) > 1
                else base
            )
    placeholders = tuple(
        ReactionCoreGraphicPlaceholder(
            label=assignments[identity],
            remote_class=str(record[0]),
            fragment_smiles=str(record[1]),
            attachment_port_count=int(record[2]),
        )
        for identity, record in sorted(
            representative_by_identity.items(),
            key=lambda item: assignments[item[0]],
        )
    )
    return assignments, placeholders, scaffold_collapses


def _active_coordinates(core: Any, side: str) -> Dict[tuple[int, int], Any]:
    values: Dict[tuple[int, int], Any] = {}
    for transition in core.atom_transitions:
        state = (
            transition.before_state
            if side == "reactant"
            else transition.after_state
        )
        if state is not None:
            values[(state.component_index, state.atom_index)] = state
    return values


def _bond_type(token: str) -> Any:
    name = _BOND_TYPES.get(str(token).upper(), "SINGLE")
    return getattr(Chem.BondType, name)


def _embedded_core_placeholders(
    core: Any,
    *,
    side: str,
    component_index: int,
    active_atom_indices: set[int],
    assignments: Mapping[tuple[Any, ...], str],
    topology_protected: set[str],
) -> tuple[Dict[int, str], set[str]]:
    """Find ring/scaffold fragments that can absorb one embedded core atom.

    A retained fragment with multiple ports to the same active atom is the
    remainder of a ring or scaffold around that atom. Rendering a dummy for
    every port is safe but visually misleading. When exactly one such fragment
    surrounds an active atom, render the core atom itself as the fragment
    placeholder and preserve its bonds to the other active atoms.
    """
    candidates: Dict[int, list[Any]] = {}
    for subgraph in core.remote_subgraphs:
        if (
            subgraph.side != side
            or subgraph.continuity != "retained"
            or _render_remote_explicitly(subgraph)
            or str(subgraph.subgraph_id) in topology_protected
            or subgraph.component_index != component_index
        ):
            continue
        ports = tuple(
            port
            for port in subgraph.attachment_ports
            if (
                port.core_component_index == component_index
                and port.core_atom_index in active_atom_indices
            )
        )
        core_indices = {port.core_atom_index for port in ports}
        if len(ports) > 1 and len(core_indices) == 1:
            core_index = next(iter(core_indices))
            candidates.setdefault(core_index, []).append(subgraph)
    labels: Dict[int, str] = {}
    collapsed_subgraph_ids: set[str] = set()
    for core_index, subgraphs in candidates.items():
        if len(subgraphs) != 1:
            continue
        subgraph = subgraphs[0]
        labels[core_index] = assignments[_remote_identity(subgraph)]
        collapsed_subgraph_ids.add(str(subgraph.subgraph_id))
    return labels, collapsed_subgraph_ids


def _build_side_molecules(
    analysis: Any,
    *,
    side: str,
    assignments: Mapping[tuple[Any, ...], str],
    scaffold_collapses: Tuple[_MultisiteScaffoldCollapse, ...],
) -> Tuple[Any, ...]:
    core = analysis.reaction_core
    components = analysis.reactants if side == "reactant" else analysis.products
    active = _active_coordinates(core, side)
    active_by_component: Dict[int, set[int]] = {}
    for component_index, atom_index in active:
        active_by_component.setdefault(component_index, set()).add(atom_index)
    output = []
    topology_protected = _topology_protected_subgraph_ids(
        analysis,
        side=side,
    )
    for component_index, atom_indices in sorted(active_by_component.items()):
        source = parse_smiles(components[component_index].input_smiles)
        if source is None:
            raise ValueError("reaction-core component could not be reparsed")
        editable = Chem.RWMol()
        old_to_new: Dict[int, int] = {}
        component_scaffolds = tuple(
            collapse
            for collapse in scaffold_collapses
            if (
                collapse.side == side
                and collapse.component_index == component_index
            )
        )
        scaffold_by_atom = {
            atom_index: collapse
            for collapse in component_scaffolds
            for atom_index in collapse.core_atom_indices
        }
        scaffold_node_indices: Dict[tuple[Any, ...], int] = {}
        embedded_labels, collapsed_subgraph_ids = _embedded_core_placeholders(
            core,
            side=side,
            component_index=component_index,
            active_atom_indices=atom_indices,
            assignments=assignments,
            topology_protected=topology_protected,
        )
        collapsed_subgraph_ids.update(
            subgraph_id
            for collapse in component_scaffolds
            for subgraph_id in collapse.subgraph_ids
        )
        for atom_index in sorted(atom_indices):
            scaffold = scaffold_by_atom.get(atom_index)
            if scaffold is not None:
                if scaffold.identity not in scaffold_node_indices:
                    atom = Chem.Atom(0)
                    label = assignments[scaffold.identity]
                    atom.SetProp("atomLabel", label)
                    atom.SetProp("_displayLabel", label)
                    scaffold_node_indices[scaffold.identity] = editable.AddAtom(atom)
                old_to_new[atom_index] = scaffold_node_indices[scaffold.identity]
                continue
            if atom_index in embedded_labels:
                atom = Chem.Atom(0)
                label = embedded_labels[atom_index]
                atom.SetProp("atomLabel", label)
                atom.SetProp("_displayLabel", label)
            else:
                atom = Chem.Atom(source.GetAtomWithIdx(atom_index))
                atom.SetAtomMapNum(0)
            old_to_new[atom_index] = editable.AddAtom(atom)
        for bond in source.GetBonds():
            begin = int(bond.GetBeginAtomIdx())
            end = int(bond.GetEndAtomIdx())
            if begin in atom_indices and end in atom_indices:
                if old_to_new[begin] == old_to_new[end]:
                    continue
                if editable.GetBondBetweenAtoms(
                    old_to_new[begin], old_to_new[end]
                ) is not None:
                    continue
                editable.AddBond(
                    old_to_new[begin],
                    old_to_new[end],
                    bond.GetBondType(),
                )
        for subgraph in core.remote_subgraphs:
            if subgraph.side != side:
                continue
            if str(subgraph.subgraph_id) in collapsed_subgraph_ids:
                continue
            matching_ports = tuple(
                port
                for port in subgraph.attachment_ports
                if port.core_component_index == component_index
                and port.core_atom_index in atom_indices
            )
            if not matching_ports:
                continue
            if (
                subgraph.continuity != "retained"
                or _render_remote_explicitly(subgraph)
                or str(subgraph.subgraph_id) in topology_protected
            ):
                remote_old_to_new: Dict[int, int] = {}
                for atom_index in sorted(subgraph.atom_indices):
                    atom = Chem.Atom(source.GetAtomWithIdx(atom_index))
                    atom.SetAtomMapNum(0)
                    remote_old_to_new[atom_index] = editable.AddAtom(atom)
                remote_indices = set(subgraph.atom_indices)
                for bond in source.GetBonds():
                    begin = int(bond.GetBeginAtomIdx())
                    end = int(bond.GetEndAtomIdx())
                    if begin in remote_indices and end in remote_indices:
                        editable.AddBond(
                            remote_old_to_new[begin],
                            remote_old_to_new[end],
                            bond.GetBondType(),
                        )
                for port in matching_ports:
                    editable.AddBond(
                        old_to_new[port.core_atom_index],
                        remote_old_to_new[port.attachment_atom_index],
                        _bond_type(port.bond_order),
                    )
                continue
            label = assignments[_remote_identity(subgraph)]
            for port in matching_ports:
                placeholder = Chem.Atom(0)
                placeholder.SetProp("atomLabel", label)
                placeholder.SetProp("_displayLabel", label)
                placeholder_index = editable.AddAtom(placeholder)
                editable.AddBond(
                    old_to_new[port.core_atom_index],
                    placeholder_index,
                    _bond_type(port.bond_order),
                )
        molecule = editable.GetMol()
        try:
            Chem.SanitizeMol(molecule)
        except Exception:
            molecule.UpdatePropertyCache(strict=False)
        output.append(molecule)
    return tuple(output)


def _style_draw_options(
    options: Any,
    *,
    render_preset: str,
    molecules: tuple[Any, ...],
) -> None:
    apply_render_preset(
        options,
        render_preset,
        molecules=molecules,
        context="reaction_core",
    )
    if hasattr(options, "addAtomIndices"):
        options.addAtomIndices = False


def build_reaction_core_graphic(
    analysis: Any,
    *,
    size: tuple[int, int] = (960, 260),
    image_format: str = "svg",
    render_preset: str = "current",
) -> ReactionCoreGraphic:
    """Build a compact graphic from active atoms and retained remote groups."""
    _require_rdkit()
    core = getattr(analysis, "reaction_core", None)
    if core is None:
        raise ValueError(
            "Reaction minimization requires a ReactionCoreProjection."
        )
    width, height = size
    if width <= 0 or height <= 0:
        raise ValueError("reaction-core graphic size must be positive")
    normalized_format = str(image_format).casefold()
    if normalized_format not in _SUPPORTED_FORMATS:
        raise ValueError(
            f"image_format must be one of {sorted(_SUPPORTED_FORMATS)}"
        )
    assignments, placeholders, scaffold_collapses = _placeholder_assignments(
        analysis
    )
    reaction = rdChemReactions.ChemicalReaction()
    reactants = _build_side_molecules(
        analysis,
        side="reactant",
        assignments=assignments,
        scaffold_collapses=scaffold_collapses,
    )
    products = _build_side_molecules(
        analysis,
        side="product",
        assignments=assignments,
        scaffold_collapses=scaffold_collapses,
    )
    if not reactants or not products:
        raise ValueError("reaction-core graphic has no drawable reaction")
    for molecule in reactants:
        reaction.AddReactantTemplate(molecule)
    for molecule in products:
        reaction.AddProductTemplate(molecule)
    options = rdMolDraw2D.MolDrawOptions()
    _style_draw_options(
        options,
        render_preset=render_preset,
        molecules=tuple(reactants) + tuple(products),
    )
    panel_count = max(len(reactants) + len(products) + 1, 3)
    panel_size = (
        max(width // panel_count, 120),
        max(height, 140),
    )
    if normalized_format == "svg":
        drawing = Draw.ReactionToImage(
            reaction,
            subImgSize=panel_size,
            useSVG=True,
            drawOptions=options,
        )
    else:
        drawing = Draw.ReactionToImage(
            reaction,
            subImgSize=panel_size,
            returnPNG=True,
            drawOptions=options,
        )
    image_bytes = (
        drawing.encode("utf-8")
        if isinstance(drawing, str)
        else bytes(drawing)
    )
    definition = load_reaction_core_graphic_definition()
    return ReactionCoreGraphic(
        image_bytes=image_bytes,
        image_format=normalized_format,
        core_id=str(core.core_id),
        evidence_status=str(core.evidence_status),
        confidence=float(core.confidence),
        placeholders=placeholders,
        warnings=tuple(core.warnings),
        definition_id=str(definition["definition_id"]),
        schema_version=str(definition["schema_version"]),
    )


def render_reaction_core_image_bytes(
    analysis: Any,
    *,
    size: tuple[int, int] = (960, 260),
    image_format: str = "svg",
    render_preset: str = "current",
) -> bytes:
    """Return just the SVG or PNG bytes for a minimized reaction scheme."""
    return build_reaction_core_graphic(
        analysis,
        size=size,
        image_format=image_format,
        render_preset=render_preset,
    ).image_bytes


__all__ = [
    "ReactionCoreGraphic",
    "ReactionCoreGraphicPlaceholder",
    "build_reaction_core_graphic",
    "load_reaction_core_graphic_definition",
    "render_reaction_core_image_bytes",
]
