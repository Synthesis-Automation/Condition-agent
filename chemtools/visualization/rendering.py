"""High-level rendering helpers for molecules and reactions.

These utilities provide ChemTools-specific styling for RDKit renderings so we
can embed reaction or substrate depictions inside both the CLI and GUI tools.
They intentionally keep hydrogens implicit except when bound to heteroatoms,
matching the style in project documentation and the provided figure.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Union

from chemtools.core.rdkit import parse_smiles, rdkit_available

try:  # pragma: no cover - import guarded in tests
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions, rdchem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D
except Exception:  # pragma: no cover - fallback when RDKit missing
    Chem = None  # type: ignore
    rdChemReactions = None  # type: ignore
    rdMolDraw2D = None  # type: ignore
    rdchem = None  # type: ignore
    Draw = None  # type: ignore

ImageSize = Tuple[int, int]

DEFAULT_MOLECULE_SIZE: ImageSize = (480, 300)
DEFAULT_REACTION_SIZE: ImageSize = (960, 320)
SUPPORTED_FORMATS = {"png", "svg"}


@dataclass(frozen=True)
class RenderStyle:
    """Style options shared by molecule and reaction rendering."""

    size: ImageSize
    image_format: str = "png"
    kekulize: bool = True

    def normalized_format(self) -> str:
        fmt = (self.image_format or "png").lower()
        if fmt not in SUPPORTED_FORMATS:
            raise ValueError(
                f"image_format must be one of {sorted(SUPPORTED_FORMATS)}, got {fmt!r}"
            )
        return fmt


def _require_rdkit() -> None:
    if not rdkit_available() or Chem is None or rdMolDraw2D is None:
        raise RuntimeError(
            "RDKit is required for rendering. Install rdkit or set "
            "CHEMTOOLS_DISABLE_RDKIT=1 to skip image generation."
        )


def _ensure_output_path(path: Union[str, Path]) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def _promote_hetero_hydrogens(mol: "Chem.Mol") -> "Chem.Mol":
    """Materialize hydrogens attached to hetero atoms for drawing clarity."""
    rw = Chem.RWMol(mol)
    changed = False
    for atom in rw.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num in (0, 1, 6):
            continue
        total_h = atom.GetTotalNumHs()
        if total_h:
            atom.SetNumExplicitHs(total_h)
            atom.SetNoImplicit(True)
            changed = True
    if changed:
        try:
            Chem.SanitizeMol(rw, sanitizeOps=rdchem.SanitizeFlags.SANITIZE_ADJUSTHS)
        except Exception:
            Chem.SanitizeMol(rw)
    return rw.GetMol()


def _prepare_mol_for_drawing(
    smiles: str,
    *,
    kekulize: bool,
) -> "Chem.Mol":
    mol = parse_smiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles!r}")
    mol = Chem.Mol(mol)
    mol = _promote_hetero_hydrogens(mol)
    rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize)
    return mol


def _split_reaction_sections(reaction_smiles: str) -> Tuple[str, str, str]:
    txt = (reaction_smiles or "").strip()
    if not txt:
        raise ValueError("Reaction SMILES cannot be empty.")
    parts = txt.split(">>")
    if len(parts) == 3:
        reactants = parts[0]
        reagent = parts[1]
        products = parts[2]
        combined = f"{reactants}.{reagent}" if reagent else reactants
        combined = combined.strip(".")
        return combined, "", products
    if len(parts) == 2:
        return parts[0], "", parts[1]
    sections = txt.split(">")
    if len(sections) == 3:
        return tuple(sections)  # type: ignore[return-value]
    raise ValueError(f"Invalid reaction SMILES string: {reaction_smiles!r}")


def _split_side(section: str) -> Tuple[str, ...]:
    section = section.strip()
    if not section:
        return tuple()
    return tuple(part for part in section.split(".") if part)


def _prepare_reaction_for_drawing(
    reaction_smiles: str,
    *,
    kekulize: bool,
) -> "rdChemReactions.ChemicalReaction":
    if rdChemReactions is None:
        raise RuntimeError("RDKit reaction utilities are unavailable.")
    reactants, agents, products = _split_reaction_sections(reaction_smiles)
    rxn = rdChemReactions.ChemicalReaction()
    for smi in _split_side(reactants):
        rxn.AddReactantTemplate(_prepare_mol_for_drawing(smi, kekulize=kekulize))
    for smi in _split_side(products):
        rxn.AddProductTemplate(_prepare_mol_for_drawing(smi, kekulize=kekulize))
    for smi in _split_side(agents):
        rxn.AddAgentTemplate(_prepare_mol_for_drawing(smi, kekulize=kekulize))
    if rxn.GetNumReactantTemplates() == 0 or rxn.GetNumProductTemplates() == 0:
        raise ValueError("Reaction SMILES must contain at least one reactant and product.")
    return rxn


def _apply_heteroatom_palette(opts: "rdMolDraw2D.MolDrawOptions") -> None:
    """Apply CPK-like colours for heteroatoms so they stand out from carbon."""
    # Restore RDKit's built-in colour scheme first (undoes any prior BW call)
    if hasattr(opts, "useDefaultAtomPalette"):
        opts.useDefaultAtomPalette()

    # Fine-tune individual element colours: (R, G, B) in [0, 1]
    PALETTE: dict[int, tuple[float, float, float]] = {
        7:  (0.00, 0.00, 0.90),   # N  – vivid blue
        8:  (0.85, 0.05, 0.05),   # O  – vivid red
        9:  (0.15, 0.75, 0.15),   # F  – green
        15: (1.00, 0.45, 0.00),   # P  – orange
        16: (0.75, 0.60, 0.00),   # S  – gold / dark yellow
        17: (0.00, 0.65, 0.00),   # Cl – medium green
        35: (0.55, 0.10, 0.10),   # Br – dark red / brown
        53: (0.45, 0.00, 0.55),   # I  – purple
        14: (0.50, 0.35, 0.05),   # Si – brown
        5:  (0.85, 0.35, 0.10),   # B  – brick orange
    }
    if hasattr(opts, "atomColourPalette"):
        opts.atomColourPalette.update(PALETTE)


def _style_draw_options(opts: "rdMolDraw2D.MolDrawOptions") -> None:
    opts.padding = 0.08
    opts.bondLineWidth = 2.0
    if hasattr(opts, "explicitMethyl"):
        opts.explicitMethyl = False
    if hasattr(opts, "addStereoAnnotation"):
        opts.addStereoAnnotation = True
    if hasattr(opts, "minFontSize"):
        opts.minFontSize = 14
    if hasattr(opts, "maxFontSize"):
        opts.maxFontSize = 36
    _apply_heteroatom_palette(opts)


def _configure_draw_options(drawer: "rdMolDraw2D.MolDraw2D") -> None:
    opts = drawer.drawOptions()
    _style_draw_options(opts)
    try:
        drawer.SetBackgroundColour((1, 1, 1))
    except AttributeError:
        pass


def _make_drawer(style: RenderStyle) -> "rdMolDraw2D.MolDraw2D":
    width, height = style.size
    fmt = style.normalized_format()
    if fmt == "png":
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    else:
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    _configure_draw_options(drawer)
    return drawer


def _write_image(drawer: "rdMolDraw2D.MolDraw2D", output_path: Path) -> Path:
    drawer.FinishDrawing()
    data = drawer.GetDrawingText()
    if isinstance(data, str):
        output_path.write_text(data, encoding="utf-8")
    else:
        output_path.write_bytes(data)
    return output_path


def render_molecule_image(
    smiles: str,
    output_path: Union[str, Path],
    *,
    size: ImageSize = DEFAULT_MOLECULE_SIZE,
    image_format: str = "png",
    kekulize: bool = True,
    legend: str | None = None,
) -> Path:
    """
    Render a molecule from SMILES to an image file.

    Args:
        smiles: Compound SMILES string.
        output_path: Where to write the image.
        size: (width, height) in pixels.
        image_format: "png" or "svg".
        kekulize: Whether to kekulize aromatic systems before drawing.
        legend: Optional caption rendered below the molecule.

    Returns:
        Path to the written image.
    """
    _require_rdkit()
    style = RenderStyle(size=size, image_format=image_format, kekulize=kekulize)
    mol = _prepare_mol_for_drawing(smiles, kekulize=style.kekulize)
    drawer = _make_drawer(style)
    drawer.DrawMolecule(mol, legend=legend or "")
    return _write_image(drawer, _ensure_output_path(output_path))


def render_reaction_image(
    reaction_smiles: str,
    output_path: Union[str, Path],
    *,
    size: ImageSize = DEFAULT_REACTION_SIZE,
    image_format: str = "png",
    kekulize: bool = True,
) -> Path:
    """
    Render a reaction SMILES to an image with reagents, arrow, and products.

    Args:
        reaction_smiles: Reaction SMILES (reactants>agents>products).
        output_path: Where to write the image.
        size: (width, height) in pixels.
        image_format: "png" or "svg".
        kekulize: Whether to kekulize aromatic systems.

    Returns:
        Path to the written image.
    """
    _require_rdkit()
    if Draw is None:
        raise RuntimeError("RDKit drawing backend unavailable.")
    style = RenderStyle(size=size, image_format=image_format, kekulize=kekulize)
    rxn = _prepare_reaction_for_drawing(reaction_smiles, kekulize=style.kekulize)
    opts = rdMolDraw2D.MolDrawOptions()
    _style_draw_options(opts)

    fmt = style.normalized_format()
    panel_count = (
        max(rxn.GetNumReactantTemplates(), 1)
        + max(rxn.GetNumAgentTemplates(), 0)
        + max(rxn.GetNumProductTemplates(), 1)
    )
    panel_count = max(panel_count, 3)
    panel_width = max(int(style.size[0] / panel_count), 200)
    panel_height = max(style.size[1], 200)
    output = _ensure_output_path(output_path)

    if fmt == "svg":
        svg = Draw.ReactionToImage(
            rxn,
            subImgSize=(panel_width, panel_height),
            useSVG=True,
            drawOptions=opts,
        )
        output.write_text(svg, encoding="utf-8")
    else:
        png = Draw.ReactionToImage(
            rxn,
            subImgSize=(panel_width, panel_height),
            drawOptions=opts,
            returnPNG=True,
        )
        output.write_bytes(png)
    return output
