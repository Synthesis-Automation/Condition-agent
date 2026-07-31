"""Render molecular graphs and reactions with consistent RDKit styling.

The helpers are deliberately independent of the legacy :mod:`chemtools`
package. File-based functions remain useful for exports, while the byte-based
functions let application layers display drawings without temporary files.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple, Union

from reactive_taxonomy.chemistry.rdkit_utils import (
    parse_smiles,
    rdkit_available,
)

try:  # pragma: no cover - availability is exercised through public helpers
    from rdkit import Chem
    from rdkit.Chem import Draw, rdChemReactions, rdchem
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError:  # pragma: no cover - produces a clear runtime error
    Chem = None  # type: ignore[assignment]
    Draw = None  # type: ignore[assignment]
    rdChemReactions = None  # type: ignore[assignment]
    rdMolDraw2D = None  # type: ignore[assignment]
    rdchem = None  # type: ignore[assignment]

ImageSize = Tuple[int, int]
RenderedImage = Union[bytes, str]

DEFAULT_MOLECULE_SIZE: ImageSize = (480, 300)
DEFAULT_REACTION_SIZE: ImageSize = (960, 320)
SUPPORTED_FORMATS = frozenset({"png", "svg"})
RENDER_STYLE_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "render_styles.v1.json"
)

_HETEROATOM_PALETTE = {
    5: (0.85, 0.35, 0.10),  # B: brick orange
    7: (0.00, 0.00, 0.90),  # N: vivid blue
    8: (0.85, 0.05, 0.05),  # O: vivid red
    9: (0.15, 0.75, 0.15),  # F: green
    14: (0.50, 0.35, 0.05),  # Si: brown
    15: (1.00, 0.45, 0.00),  # P: orange
    16: (0.75, 0.60, 0.00),  # S: gold
    17: (0.00, 0.65, 0.00),  # Cl: medium green
    35: (0.55, 0.10, 0.10),  # Br: dark red
    53: (0.45, 0.00, 0.55),  # I: purple
}


@dataclass(frozen=True)
class RenderStyle:
    """Validated style options shared by molecule and reaction rendering."""

    size: ImageSize
    image_format: str = "png"
    kekulize: bool = True
    render_preset: str = "current"

    def normalized_format(self) -> str:
        """Return a supported lowercase output format."""
        image_format = str(self.image_format or "png").casefold()
        if image_format not in SUPPORTED_FORMATS:
            raise ValueError(
                "image_format must be one of "
                f"{sorted(SUPPORTED_FORMATS)}, got {image_format!r}"
            )
        return image_format

    def validated_size(self) -> ImageSize:
        """Return a positive integer image size."""
        if len(self.size) != 2:
            raise ValueError("size must contain width and height.")
        width, height = self.size
        if (
            isinstance(width, bool)
            or isinstance(height, bool)
            or not isinstance(width, int)
            or not isinstance(height, int)
            or width <= 0
            or height <= 0
        ):
            raise ValueError("size width and height must be positive integers.")
        return width, height

    def validated_preset(self) -> str:
        """Return a known versioned rendering preset identifier."""
        preset = str(self.render_preset or "current")
        definitions = load_render_style_definitions()
        if preset not in definitions["presets"]:
            raise ValueError(
                "render_preset must be one of "
                f"{sorted(definitions['presets'])}, got {preset!r}"
            )
        return preset


@lru_cache(maxsize=1)
def load_render_style_definitions() -> Dict[str, Any]:
    """Load and validate the versioned molecule-drawing presets."""
    with RENDER_STYLE_DEFINITION_PATH.open("r", encoding="utf-8") as handle:
        definition = dict(json.load(handle))
    if definition.get("schema_version") != "1.0":
        raise ValueError("unsupported render-style schema")
    if definition.get("definition_id") != "render_styles.v1":
        raise ValueError("unexpected render-style definition ID")
    presets = definition.get("presets")
    if not isinstance(presets, Mapping) or not presets:
        raise ValueError("render-style definition requires presets")
    if definition.get("default_preset") not in presets:
        raise ValueError("render-style default preset is not defined")
    for preset_id, preset in presets.items():
        if not isinstance(preset_id, str) or not preset_id:
            raise ValueError("render-style preset IDs must be non-empty strings")
        if not isinstance(preset, Mapping):
            raise ValueError(
                f"render-style preset {preset_id!r} must be an object"
            )
        if preset.get("mode") not in {"custom", "acs_1996"}:
            raise ValueError(
                f"render-style preset {preset_id!r} has invalid mode"
            )
        target_length = preset.get("target_bond_length_pixels")
        if target_length is not None and (
            not isinstance(target_length, (int, float)) or target_length <= 0
        ):
            raise ValueError(
                f"render-style preset {preset_id!r} has invalid bond length"
            )
        bond_length_scale = preset.get("bond_length_scale", 1.0)
        if (
            isinstance(bond_length_scale, bool)
            or not isinstance(bond_length_scale, (int, float))
            or not 0 < bond_length_scale <= 1
        ):
            raise ValueError(
                f"render-style preset {preset_id!r} has invalid bond-length scale"
            )
        for context in ("standard", "reaction_core"):
            if not isinstance(preset.get(context), Mapping):
                raise ValueError(
                    f"render-style preset {preset_id!r} requires {context} options"
                )
    return definition


def available_render_presets() -> Tuple[Tuple[str, str], ...]:
    """Return deterministic preset identifiers and display labels."""
    presets = load_render_style_definitions()["presets"]
    return tuple(
        (str(preset_id), str(preset.get("label") or preset_id))
        for preset_id, preset in presets.items()
    )


def _require_rdkit() -> None:
    if (
        not rdkit_available()
        or Chem is None
        or Draw is None
        or rdChemReactions is None
        or rdMolDraw2D is None
    ):
        raise RuntimeError("RDKit is required for molecule and reaction rendering.")


def _promote_hetero_hydrogens(molecule: Any) -> Any:
    """Materialize heteroatom hydrogen counts as readable atom labels."""
    editable = Chem.RWMol(molecule)
    changed = False
    for atom in editable.GetAtoms():
        if atom.GetAtomicNum() in (0, 1, 6):
            continue
        hydrogen_count = atom.GetTotalNumHs()
        if hydrogen_count:
            atom.SetNumExplicitHs(hydrogen_count)
            atom.SetNoImplicit(True)
            changed = True
    if changed:
        try:
            Chem.SanitizeMol(
                editable,
                sanitizeOps=rdchem.SanitizeFlags.SANITIZE_ADJUSTHS,
            )
        except Exception:
            Chem.SanitizeMol(editable)
    return editable.GetMol()


def _prepare_molecule(smiles: str, *, kekulize: bool) -> Any:
    molecule = parse_smiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES string: {smiles!r}")
    molecule = _promote_hetero_hydrogens(Chem.Mol(molecule))
    return rdMolDraw2D.PrepareMolForDrawing(
        molecule,
        kekulize=kekulize,
    )


def _split_reaction_sections(reaction_smiles: str) -> Tuple[str, str, str]:
    text = str(reaction_smiles or "").strip()
    if not text:
        raise ValueError("Reaction SMILES cannot be empty.")
    if ">>" in text:
        if text.count(">>") != 1:
            raise ValueError(f"Invalid reaction SMILES string: {reaction_smiles!r}")
        reactants, products = text.split(">>")
        agents = ""
    else:
        sections = text.split(">")
        if len(sections) != 3:
            raise ValueError(f"Invalid reaction SMILES string: {reaction_smiles!r}")
        reactants, agents, products = sections
    if not reactants.strip() or not products.strip():
        raise ValueError(
            "Reaction SMILES must contain at least one reactant and product."
        )
    return reactants, agents, products


def _split_side(section: str) -> Tuple[str, ...]:
    return tuple(part.strip() for part in section.split(".") if part.strip())


def _prepare_reaction(reaction_smiles: str, *, kekulize: bool) -> Any:
    reactants, agents, products = _split_reaction_sections(reaction_smiles)
    reaction = rdChemReactions.ChemicalReaction()
    for smiles in _split_side(reactants):
        reaction.AddReactantTemplate(
            _prepare_molecule(smiles, kekulize=kekulize)
        )
    for smiles in _split_side(agents):
        reaction.AddAgentTemplate(_prepare_molecule(smiles, kekulize=kekulize))
    for smiles in _split_side(products):
        reaction.AddProductTemplate(_prepare_molecule(smiles, kekulize=kekulize))
    return reaction


def _mean_bond_length(molecules: Tuple[Any, ...]) -> float:
    values = []
    for molecule in molecules:
        if molecule is None or molecule.GetNumBonds() <= 0:
            continue
        try:
            values.append(float(rdMolDraw2D.MeanBondLength(molecule)))
        except (RuntimeError, ValueError):
            continue
    return sum(values) / len(values) if values else 1.5


def apply_render_preset(
    options: Any,
    render_preset: str,
    *,
    molecules: Tuple[Any, ...],
    context: str = "standard",
) -> None:
    """Apply one validated definition-backed RDKit drawing preset."""
    definitions = load_render_style_definitions()
    if render_preset not in definitions["presets"]:
        raise ValueError(
            "render_preset must be one of "
            f"{sorted(definitions['presets'])}, got {render_preset!r}"
        )
    preset = definitions["presets"][render_preset]
    if context not in {"standard", "reaction_core"}:
        raise ValueError(f"unknown rendering context: {context!r}")
    context_options = preset[context]
    mean_bond_length = _mean_bond_length(molecules)
    if preset["mode"] == "acs_1996":
        rdMolDraw2D.SetACS1996Mode(options, mean_bond_length)
    if "padding" in context_options:
        padding = float(context_options["padding"])
        bond_length_scale = float(preset.get("bond_length_scale", 1.0))
        # RDKit's automatic fit has no relative bond-length option. Scaling the
        # drawable span preserves automatic fitting while shortening every bond
        # by the requested proportion: span = 1 - (2 * padding).
        options.padding = (
            1.0 - bond_length_scale * (1.0 - 2.0 * padding)
        ) / 2.0
    if "bond_line_width" in context_options:
        options.bondLineWidth = float(context_options["bond_line_width"])
    if "min_font_size" in context_options and hasattr(options, "minFontSize"):
        options.minFontSize = int(context_options["min_font_size"])
    if "max_font_size" in context_options and hasattr(options, "maxFontSize"):
        options.maxFontSize = int(context_options["max_font_size"])
    target_bond_length = preset.get("target_bond_length_pixels")
    if target_bond_length is not None and hasattr(options, "fixedBondLength"):
        options.fixedBondLength = float(target_bond_length) / mean_bond_length
    if hasattr(options, "explicitMethyl"):
        options.explicitMethyl = False
    if hasattr(options, "addStereoAnnotation"):
        options.addStereoAnnotation = True
    if preset["mode"] != "acs_1996" and hasattr(
        options, "useDefaultAtomPalette"
    ):
        options.useDefaultAtomPalette()
    if preset["mode"] != "acs_1996" and hasattr(options, "updateAtomPalette"):
        options.updateAtomPalette(_HETEROATOM_PALETTE)


def _make_molecule_drawer(style: RenderStyle, molecule: Any) -> Any:
    width, height = style.validated_size()
    if style.normalized_format() == "png":
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    else:
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    apply_render_preset(
        drawer.drawOptions(),
        style.validated_preset(),
        molecules=(molecule,),
    )
    try:
        drawer.SetBackgroundColour((1, 1, 1))
    except AttributeError:
        pass
    return drawer


def _drawing_bytes(drawing: RenderedImage) -> bytes:
    if isinstance(drawing, str):
        return drawing.encode("utf-8")
    return bytes(drawing)


def render_molecule_image_bytes(
    smiles: str,
    *,
    size: ImageSize = DEFAULT_MOLECULE_SIZE,
    image_format: str = "png",
    kekulize: bool = True,
    legend: str | None = None,
    render_preset: str = "current",
) -> bytes:
    """Render a molecule from SMILES and return PNG or UTF-8 SVG bytes."""
    _require_rdkit()
    style = RenderStyle(
        size=size,
        image_format=image_format,
        kekulize=kekulize,
        render_preset=render_preset,
    )
    molecule = _prepare_molecule(smiles, kekulize=style.kekulize)
    drawer = _make_molecule_drawer(style, molecule)
    drawer.DrawMolecule(molecule, legend=legend or "")
    drawer.FinishDrawing()
    return _drawing_bytes(drawer.GetDrawingText())


def render_reaction_image_bytes(
    reaction_smiles: str,
    *,
    size: ImageSize = DEFAULT_REACTION_SIZE,
    image_format: str = "png",
    kekulize: bool = True,
    render_preset: str = "current",
) -> bytes:
    """Render reaction SMILES and return PNG or UTF-8 SVG bytes."""
    _require_rdkit()
    style = RenderStyle(
        size=size,
        image_format=image_format,
        kekulize=kekulize,
        render_preset=render_preset,
    )
    width, height = style.validated_size()
    reaction = _prepare_reaction(
        reaction_smiles,
        kekulize=style.kekulize,
    )
    options = rdMolDraw2D.MolDrawOptions()
    molecules = (
        tuple(reaction.GetReactants())
        + tuple(reaction.GetAgents())
        + tuple(reaction.GetProducts())
    )
    apply_render_preset(
        options,
        style.validated_preset(),
        molecules=molecules,
    )
    panel_count = max(
        reaction.GetNumReactantTemplates()
        + reaction.GetNumAgentTemplates()
        + reaction.GetNumProductTemplates(),
        3,
    )
    panel_size = (max(width // panel_count, 100), max(height, 100))
    if style.normalized_format() == "svg":
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
            drawOptions=options,
            returnPNG=True,
        )
    return _drawing_bytes(drawing)


def _write_image(output_path: Union[str, Path], drawing: bytes) -> Path:
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(drawing)
    return destination


def render_molecule_image(
    smiles: str,
    output_path: Union[str, Path],
    *,
    size: ImageSize = DEFAULT_MOLECULE_SIZE,
    image_format: str = "png",
    kekulize: bool = True,
    legend: str | None = None,
    render_preset: str = "current",
) -> Path:
    """Render a molecule from SMILES to a PNG or SVG file."""
    drawing = render_molecule_image_bytes(
        smiles,
        size=size,
        image_format=image_format,
        kekulize=kekulize,
        legend=legend,
        render_preset=render_preset,
    )
    return _write_image(output_path, drawing)


def render_reaction_image(
    reaction_smiles: str,
    output_path: Union[str, Path],
    *,
    size: ImageSize = DEFAULT_REACTION_SIZE,
    image_format: str = "png",
    kekulize: bool = True,
    render_preset: str = "current",
) -> Path:
    """Render reaction SMILES to a PNG or SVG file."""
    drawing = render_reaction_image_bytes(
        reaction_smiles,
        size=size,
        image_format=image_format,
        kekulize=kekulize,
        render_preset=render_preset,
    )
    return _write_image(output_path, drawing)
