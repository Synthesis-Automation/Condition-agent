"""Small JSON and HTML review output for synthetic partition landscapes."""

from __future__ import annotations

from html import escape
import json
from pathlib import Path
from typing import Iterable

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from .synthetic_partition import (
    SyntheticPartition,
    SyntheticPartitionLandscape,
    analyze_partition_target,
)


PARTITION_REVIEW_SCHEMA_VERSION = "1.0"
_COLORS = (
    (0.31, 0.59, 0.89),
    (0.91, 0.45, 0.35),
    (0.38, 0.72, 0.48),
    (0.73, 0.51, 0.86),
    (0.94, 0.70, 0.29),
    (0.30, 0.75, 0.75),
    (0.78, 0.64, 0.43),
    (0.65, 0.65, 0.65),
)


def render_partition_svg(partition: SyntheticPartition) -> str:
    """Render one target with role-neutral module coloring."""

    molecule = Chem.MolFromSmiles(partition.target_smiles)
    if molecule is None:
        return "<p>Target depiction unavailable.</p>"
    _, _, target_atoms, _ = analyze_partition_target(partition.target_smiles)
    index_by_target_map = {
        reference.atom_map: reference.atom_index for reference in target_atoms
    }
    highlight_atoms = []
    highlight_colors = {}
    for module_index, module in enumerate(partition.modules):
        color = _COLORS[module_index % len(_COLORS)]
        for target_map in module.target_atom_maps:
            atom_index = index_by_target_map.get(target_map)
            if atom_index is None or atom_index >= molecule.GetNumAtoms():
                continue
            highlight_atoms.append(atom_index)
            highlight_colors[atom_index] = color
    drawer = rdMolDraw2D.MolDraw2DSVG(560, 300)
    options = drawer.drawOptions()
    options.addAtomIndices = True
    drawer.DrawMolecule(
        molecule,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_colors,
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _items(values: Iterable[str]) -> str:
    items = tuple(values)
    if not items:
        return "<span class='muted'>none</span>"
    return "<ul>" + "".join(f"<li>{escape(item)}</li>" for item in items) + "</ul>"


def render_partition_landscape_html(
    landscape: SyntheticPartitionLandscape,
) -> str:
    """Render one static, chemist-readable partition review document."""

    cards = []
    for rank, partition in enumerate(landscape.partitions, start=1):
        module_rows = "".join(
            "<tr>"
            f"<td>M{index}</td>"
            f"<td>{len(module.target_atom_maps)}</td>"
            f"<td>{escape(', '.join(map(str, module.target_atom_maps)))}</td>"
            f"<td>{escape(', '.join(map(str, module.attachment_atom_maps)) or '—')}</td>"
            "</tr>"
            for index, module in enumerate(partition.modules, start=1)
        )
        interface_rows = "".join(
            "<tr>"
            f"<td>{escape(interface.interface_kind)}</td>"
            f"<td>{escape(interface.evidence_level)}</td>"
            f"<td>{escape(', '.join(interface.candidate_operator_ids) or '—')}</td>"
            f"<td>{escape('; '.join(f'{a}-{b}:{kind}' for a, b, kind in interface.target_bonds) or 'unary')}</td>"
            "</tr>"
            for interface in partition.interfaces
        )
        cards.append(
            "<section class='card'>"
            f"<h2>#{rank} · k={partition.k} · {escape(partition.evidence_level)}</h2>"
            f"<p><code>{escape(partition.partition_id)}</code></p>"
            f"<p>Source: <strong>{escape(partition.source_kind)}</strong>; "
            f"heuristic score: {partition.heuristic_score:.4f}; "
            f"realization: {escape(partition.realization_status)}</p>"
            f"<div class='depiction'>{render_partition_svg(partition)}</div>"
            "<h3>Role-neutral modules</h3>"
            "<table><thead><tr><th>Display</th><th>Atoms</th>"
            "<th>Target maps</th><th>Attachment maps</th></tr></thead>"
            f"<tbody>{module_rows}</tbody></table>"
            "<h3>Strategic interfaces</h3>"
            "<table><thead><tr><th>Kind</th><th>Evidence</th>"
            "<th>Operators</th><th>Target bonds</th></tr></thead>"
            f"<tbody>{interface_rows or '<tr><td colspan=4>none</td></tr>'}</tbody>"
            "</table>"
            "<h3>Warnings</h3>"
            f"{_items(partition.warnings)}"
            "</section>"
        )
    payload = json.dumps(
        landscape.to_dict(),
        sort_keys=True,
        separators=(",", ":"),
    ).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Synthetic partition landscape</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #20242a; background: #f5f6f8; }}
header, .card {{ max-width: 980px; margin: 0 auto 1.2rem; background: white; padding: 1.2rem; border-radius: 10px; box-shadow: 0 1px 5px #0002; }}
h1, h2, h3 {{ margin-top: 0.4rem; }}
table {{ border-collapse: collapse; width: 100%; margin-bottom: 1rem; }}
th, td {{ border: 1px solid #d8dce2; padding: 0.45rem; text-align: left; vertical-align: top; }}
th {{ background: #eef1f5; }}
code {{ overflow-wrap: anywhere; }}
.depiction {{ overflow-x: auto; text-align: center; }}
.muted {{ color: #6c737f; }}
.warning {{ color: #8b4a00; }}
</style>
</head>
<body>
<header>
<h1>Synthetic partition landscape</h1>
<p><code>{escape(landscape.target_smiles)}</code></p>
<p>{len(landscape.partitions)} retained views across k={escape(', '.join(map(str, landscape.searched_k_values)))}.</p>
<p>Abstained: <strong>{str(landscape.abstained).lower()}</strong></p>
<div class="warning">{_items(landscape.warnings + landscape.abstention_reasons)}</div>
</header>
{''.join(cards)}
<script type="application/json" id="partition-landscape-data">{payload}</script>
</body>
</html>
"""


def write_partition_landscape_review(
    landscape: SyntheticPartitionLandscape,
    json_path: str | Path,
    html_path: str | Path,
) -> None:
    """Write deterministic JSON and a static HTML review."""

    output_json = Path(json_path)
    output_html = Path(html_path)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(landscape.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    output_html.write_text(
        render_partition_landscape_html(landscape),
        encoding="utf-8",
    )


__all__ = [
    "PARTITION_REVIEW_SCHEMA_VERSION",
    "render_partition_svg",
    "render_partition_landscape_html",
    "write_partition_landscape_review",
]
