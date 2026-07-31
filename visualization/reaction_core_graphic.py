"""Render a mapped reaction core as a compact, scaffold-abstracted scheme."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple

from reactive_taxonomy.chemistry.rdkit_utils import parse_smiles, rdkit_available

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
    functional_group_ids: Tuple[str, ...]


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


@lru_cache(maxsize=1)
def load_reaction_core_graphic_definition() -> Dict[str, Any]:
    """Load and validate the versioned placeholder-rendering definition."""
    with _DEFINITION_PATH.open("r", encoding="utf-8") as handle:
        definition = dict(json.load(handle))
    if str(definition.get("schema_version") or "") != "1.0":
        raise ValueError("unsupported reaction-core graphic schema")
    if str(definition.get("definition_id") or "") != (
        "reaction_core_graphic.v1"
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


def _placeholder_assignments(
    core: Any,
) -> tuple[
    Dict[tuple[Any, ...], str],
    Tuple[ReactionCoreGraphicPlaceholder, ...],
]:
    definition = load_reaction_core_graphic_definition()
    labels = definition["remote_class_labels"]
    retained = tuple(
        subgraph
        for subgraph in core.remote_subgraphs
        if subgraph.continuity == "retained"
    )
    representative_by_identity: Dict[tuple[Any, ...], Any] = {}
    for subgraph in retained:
        representative_by_identity.setdefault(
            _remote_identity(subgraph),
            subgraph,
        )
    identities_by_base: Dict[str, list[tuple[Any, ...]]] = {}
    for identity, subgraph in representative_by_identity.items():
        base = str(labels.get(subgraph.remote_class) or "R")
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
            remote_class=str(subgraph.remote_class),
            fragment_smiles=str(subgraph.fragment_smiles),
            functional_group_ids=tuple(subgraph.functional_group_ids),
        )
        for identity, subgraph in sorted(
            representative_by_identity.items(),
            key=lambda item: assignments[item[0]],
        )
    )
    return assignments, placeholders


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


def _build_side_molecules(
    analysis: Any,
    *,
    side: str,
    assignments: Mapping[tuple[Any, ...], str],
) -> Tuple[Any, ...]:
    core = analysis.reaction_core
    components = analysis.reactants if side == "reactant" else analysis.products
    active = _active_coordinates(core, side)
    active_by_component: Dict[int, set[int]] = {}
    for component_index, atom_index in active:
        active_by_component.setdefault(component_index, set()).add(atom_index)
    output = []
    for component_index, atom_indices in sorted(active_by_component.items()):
        source = parse_smiles(components[component_index].input_smiles)
        if source is None:
            raise ValueError("reaction-core component could not be reparsed")
        editable = Chem.RWMol()
        old_to_new: Dict[int, int] = {}
        for atom_index in sorted(atom_indices):
            atom = Chem.Atom(source.GetAtomWithIdx(atom_index))
            atom.SetAtomMapNum(0)
            old_to_new[atom_index] = editable.AddAtom(atom)
        for bond in source.GetBonds():
            begin = int(bond.GetBeginAtomIdx())
            end = int(bond.GetEndAtomIdx())
            if begin in atom_indices and end in atom_indices:
                editable.AddBond(
                    old_to_new[begin],
                    old_to_new[end],
                    bond.GetBondType(),
                )
        for subgraph in core.remote_subgraphs:
            if subgraph.side != side or subgraph.continuity != "retained":
                continue
            matching_ports = tuple(
                port
                for port in subgraph.attachment_ports
                if port.core_component_index == component_index
                and port.core_atom_index in atom_indices
            )
            if not matching_ports:
                continue
            placeholder = Chem.Atom(0)
            label = assignments[_remote_identity(subgraph)]
            placeholder.SetProp("atomLabel", label)
            placeholder.SetProp("_displayLabel", label)
            placeholder_index = editable.AddAtom(placeholder)
            for port in matching_ports:
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


def _style_draw_options(options: Any) -> None:
    options.padding = 0.04
    options.bondLineWidth = 2.2
    if hasattr(options, "minFontSize"):
        options.minFontSize = 18
    if hasattr(options, "maxFontSize"):
        options.maxFontSize = 36
    if hasattr(options, "addAtomIndices"):
        options.addAtomIndices = False
    if hasattr(options, "useDefaultAtomPalette"):
        options.useDefaultAtomPalette()


def build_reaction_core_graphic(
    analysis: Any,
    *,
    size: tuple[int, int] = (960, 260),
    image_format: str = "svg",
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
    assignments, placeholders = _placeholder_assignments(core)
    reaction = rdChemReactions.ChemicalReaction()
    reactants = _build_side_molecules(
        analysis,
        side="reactant",
        assignments=assignments,
    )
    products = _build_side_molecules(
        analysis,
        side="product",
        assignments=assignments,
    )
    if not reactants or not products:
        raise ValueError("reaction-core graphic has no drawable reaction")
    for molecule in reactants:
        reaction.AddReactantTemplate(molecule)
    for molecule in products:
        reaction.AddProductTemplate(molecule)
    options = rdMolDraw2D.MolDrawOptions()
    _style_draw_options(options)
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
) -> bytes:
    """Return just the SVG or PNG bytes for a minimized reaction scheme."""
    return build_reaction_core_graphic(
        analysis,
        size=size,
        image_format=image_format,
    ).image_bytes


__all__ = [
    "ReactionCoreGraphic",
    "ReactionCoreGraphicPlaceholder",
    "build_reaction_core_graphic",
    "load_reaction_core_graphic_definition",
    "render_reaction_core_image_bytes",
]
