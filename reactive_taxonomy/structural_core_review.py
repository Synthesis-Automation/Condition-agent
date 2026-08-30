"""Self-contained review rendering for structural-core observations."""

from __future__ import annotations

import html
from pathlib import Path
from typing import Sequence

from .chemistry.rdkit_utils import parse_smiles
from .structural_cores import StructuralCoreAnalysis, StructuralCoreObservation


def _draw_analysis_svg(analysis: StructuralCoreAnalysis) -> str:
    from rdkit.Chem.Draw import rdMolDraw2D  # type: ignore

    molecule = parse_smiles(analysis.canonical_smiles or "")
    if molecule is None:
        return '<p class="invalid">Structure could not be rendered.</p>'
    focus_atoms = {
        index
        for observation in analysis.observations
        for index in observation.focus_atom_indices
    }
    focus_bond_pairs = {
        (int(parts[1]), int(parts[2]))
        for observation in analysis.observations
        for key in observation.focus_bond_keys
        if len(parts := key.split(":")) >= 3
    }
    focus_bonds = [
        int(bond.GetIdx())
        for left, right in sorted(focus_bond_pairs)
        if (bond := molecule.GetBondBetweenAtoms(left, right)) is not None
    ]
    highlighted = sorted(focus_atoms)
    colors = {index: (0.90, 0.20, 0.18) for index in highlighted}
    if not highlighted and analysis.observations:
        highlighted = list(analysis.observations[0].atom_indices)
        colors = {index: (0.96, 0.55, 0.16) for index in highlighted}
    drawer = rdMolDraw2D.MolDraw2DSVG(520, 300)
    drawer.drawOptions().addAtomIndices = True
    drawer.DrawMolecule(
        molecule,
        highlightAtoms=highlighted,
        highlightBonds=focus_bonds,
        highlightAtomColors=colors,
        highlightBondColors={index: (0.90, 0.20, 0.18) for index in focus_bonds},
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText().replace(
        "<?xml version='1.0' encoding='iso-8859-1'?>",
        "",
    )


def _observation_card(
    observation: StructuralCoreObservation,
) -> str:
    descriptors = ", ".join(
        f"{html.escape(name)}={value:g}"
        for name, value in observation.descriptor_values
    )
    attachments = ", ".join(observation.attachment_bond_keys) or "none"
    warnings = ", ".join(observation.warnings) or "none"
    evidence = ", ".join(observation.evidence_codes) or "none"
    topology = ", ".join(observation.topology_tokens) or "none"
    focus_bonds = ", ".join(observation.focus_bond_keys) or "none"
    return (
        '<article class="observation">'
        f"<h3>Observation {observation.rank}: "
        f"{html.escape(observation.kind)}</h3>"
        f"<p><code>{html.escape(observation.core_observation_id)}</code></p>"
        f"<p><b>Target atom indices:</b> "
        f"{html.escape(str(observation.atom_indices))}</p>"
        f"<p><b>Focus atom indices:</b> "
        f"{html.escape(str(observation.focus_atom_indices))}</p>"
        f"<p><b>Focus bonds:</b> {html.escape(focus_bonds)}</p>"
        f"<p><b>Attachment bonds:</b> {html.escape(attachments)}</p>"
        f"<p><b>Topology:</b> {html.escape(topology)}</p>"
        f"<p><b>Descriptors:</b> {descriptors}</p>"
        f"<p><b>Evidence:</b> {html.escape(evidence)}</p>"
        f"<p><b>Observation warnings:</b> {html.escape(warnings)}</p>"
        '<details><summary>Structural match keys</summary>'
        f"<p><b>Exact:</b> <code>{html.escape(observation.structural_exact_key)}</code></p>"
        f"<p><b>Typed:</b> <code>{html.escape(observation.structural_typed_key)}</code></p>"
        f"<p><b>Shape:</b> <code>{html.escape(observation.structural_shape_key)}</code></p>"
        "</details></article>"
    )


def render_structural_core_review_html(
    analyses: Sequence[StructuralCoreAnalysis],
    *,
    title: str = "Structural-core observation review",
) -> str:
    """Render analyses as a self-contained, answer-free HTML review packet."""

    target_cards: list[str] = []
    for sequence_index, analysis in enumerate(analyses, start=1):
        identity = analysis.molecule_id or "invalid target"
        warnings = ", ".join(analysis.warnings) or "none"
        if analysis.valid:
            observations = "".join(
                _observation_card(observation)
                for observation in analysis.observations
            )
            if not observations:
                observations = '<p class="empty">No bounded core was observed.</p>'
        else:
            observations = (
                '<p class="invalid">Invalid target: '
                + html.escape(analysis.error or "unknown error")
                + "</p>"
            )
        target_cards.append(
            '<section class="target">'
            f"<h2>Target {sequence_index}</h2>"
            f"<p><b>Canonical structure:</b> "
            f"<code>{html.escape(analysis.canonical_smiles or analysis.input_smiles)}</code></p>"
            f"<p><b>Molecule identity:</b> <code>{html.escape(identity)}</code></p>"
            f"<p><b>Analysis warnings:</b> {html.escape(warnings)}</p>"
            + (_draw_analysis_svg(analysis) if analysis.valid else "")
            + observations
            + "</section>"
        )
    escaped_title = html.escape(title)
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escaped_title}</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #17202a; }}
.notice {{ max-width: 75rem; color: #4d5656; }}
.target {{ border: 1px solid #aeb6bf; border-radius: 10px; padding: 1rem;
  margin: 1.25rem 0; break-inside: avoid; }}
.observation {{ border-top: 1px solid #d5d8dc; padding-top: 0.75rem;
  margin-top: 1rem; }}
svg {{ display: block; max-width: 520px; width: 100%; height: auto; }}
code {{ background: #f4f6f7; padding: 0.15rem 0.3rem; overflow-wrap: anywhere; }}
.invalid {{ color: #922b21; font-weight: 600; }}
.empty {{ color: #626567; font-style: italic; }}
</style></head><body>
<h1>{escaped_title}</h1>
<p class="notice">Red atoms and bonds are deterministic observation foci such as
balanced graph bridges, direct ring linkers, or stereocentres. Orange is used only
when an observation has no narrower focus. Atom numbers refer only to the normalized
target shown here; they are not reaction atom maps. This packet contains structural
observations and deliberately omits source answers, commercial availability, route
counts, and preferred synthetic strategies.</p>
{''.join(target_cards)}
</body></html>
"""


def write_structural_core_review_html(
    analyses: Sequence[StructuralCoreAnalysis],
    output_path: str | Path,
    *,
    title: str = "Structural-core observation review",
) -> Path:
    """Write a self-contained structural-core review packet and return its path."""

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        render_structural_core_review_html(analyses, title=title),
        encoding="utf-8",
    )
    return destination


__all__ = [
    "render_structural_core_review_html",
    "write_structural_core_review_html",
]
