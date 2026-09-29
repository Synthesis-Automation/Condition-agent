"""Highlight provenance remains aligned with the canonical target atom IDs."""

import xml.etree.ElementTree as ET

import pytest

from reactive_taxonomy.search_fragments import suggest_search_fragments
from visualization.rendering import render_molecule_image_bytes


@pytest.mark.parametrize("smiles", ["CC(=O)c1ccc2c(c1)COc1ccccc1-2", "CC[C@H](O)CCCCCCN"])
def test_svg_highlights_exact_selected_input_atoms(smiles):
    result = suggest_search_fragments(smiles)
    chosen = result.candidates[0]
    svg = render_molecule_image_bytes(result.target_smiles, image_format="svg",
                                      highlight_atom_indices=chosen.target_atom_ids)
    drawing = ET.fromstring(svg)
    highlighted = {int(node.attrib["class"].removeprefix("atom-"))
                   for node in drawing.iter("{http://www.w3.org/2000/svg}ellipse")}
    assert highlighted == set(chosen.target_atom_ids)


def test_invalid_highlight_and_unsupported_canvas_do_not_silently_draw_wrong_atoms():
    with pytest.raises(ValueError, match="index"):
        render_molecule_image_bytes("CO", highlight_atom_indices=(10,))
    with pytest.raises(ValueError, match="fixed canvas"):
        render_molecule_image_bytes("CO", image_format="svg", expand_canvas=True, highlight_atom_indices=(0,))
