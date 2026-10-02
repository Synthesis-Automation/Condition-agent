"""Compact annotated SVGs preserve explicit chemistry and complete provenance."""

import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from visualization import (
    SchemeAnnotation, SchemeMolecule, load_annotated_scheme_style,
    render_annotated_scheme_svg,
)


def test_scheme_has_compact_visible_labels_and_complete_attributed_description() -> None:
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("Bromoethane", "CCBr"), SchemeMolecule("Ammonia", "N")),
        (SchemeMolecule("Ethylamine", "CCN"),), title="Illustrative substitution", basis="proposed",
        conditions=(SchemeAnnotation("NaOH, EtOH", "reported"),),
        yield_info=SchemeAnnotation("81% isolated yield after workup", "reported"),
    ))
    ns = "{http://www.w3.org/2000/svg}"
    visible = " ".join(node.text or "" for node in root.findall(ns + "text"))
    description = root.find(ns + "desc").text
    assert "Bromoethane" not in visible and "Ethylamine" not in visible
    assert "Bromoethane: CCBr" in description and "Ethylamine: CCN" in description
    assert not root.findall(ns + "text[@data-role='molecule-name']")
    assert "NaOH" in visible and "EtOH" in visible and "81%" in visible
    assert all(text not in visible for text in ("Reported", "Proposed", "Yield", "transformation", "agent-authored"))
    assert "reported: NaOH, EtOH" in description
    assert "Yield (reported): 81% isolated yield after workup" in description
    assert "Proposed transformation" in description
    assert "Not a feasibility assessment" in description
    assert len(root.findall(ns + "g" + "[@data-role='molecule']")) == 3
    assert len(root.findall(".//" + ns + "path")) > 10
    assert root.get("data-definition") == "annotated_scheme.v1"
    assert root.get("data-schema-version") == "1.2"
    assert not root.findall(".//" + ns + "unknown")


def test_missing_conditions_are_not_invented_and_bad_structures_fail() -> None:
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("Input", "CCO"),), (SchemeMolecule("Output", "CC=O"),),
        title="Hypothesis", basis="proposed",
    ))
    ns = "{http://www.w3.org/2000/svg}"
    description = root.find(ns + "desc").text
    assert "Conditions not supplied" in description and "Yield not supplied" in description
    assert not root.findall(ns + "text[@data-role='conditions']")
    assert not root.findall(ns + "text[@data-role='yield']")
    with pytest.raises(ValueError, match="Invalid SMILES"):
        render_annotated_scheme_svg(
            (SchemeMolecule("Bad", "C1CC"),), (SchemeMolecule("Output", "CC=O"),),
            title="Invalid", basis="unknown",
        )
    with pytest.raises(ValueError, match="explicit reactants"):
        render_annotated_scheme_svg((), (), title="Empty", basis="unknown")


def test_long_annotations_are_bounded_visually_and_complete_in_description() -> None:
    long = "NaOH, EtOH, DMF, DMSO, THF, K2CO3, Cs2CO3, acetonitrile, triethylamine"
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("A long molecule name " * 8, "C[C@H](O)Cl"),),
        (SchemeMolecule("Product", "CC(=O)Cl"),), title="Long labels", basis="proposed",
        conditions=(SchemeAnnotation(long, "proposed"),),
    ))
    ns = "{http://www.w3.org/2000/svg}"
    assert long in root.find(ns + "desc").text
    assert "C[C@H](O)Cl" in root.find(ns + "desc").text
    visible = " ".join(node.text or "" for node in root.findall(ns + "text"))
    assert "…" in visible
    assert "More conditions" not in visible
    assert float(root.get("height")) < 430
    style = load_annotated_scheme_style()
    assert len(root.findall(ns + "text[@data-role='conditions']")) <= style["max_annotation_lines"]
    assert not root.findall(ns + "text[@data-role='molecule-name']")
    width, height = float(root.get("width")), float(root.get("height"))
    for node in root.findall(ns + "text"):
        assert 0 <= float(node.get("x")) <= width
        size = float(node.get("font-size"))
        assert size <= float(node.get("y")) <= height - size / 4
        # Wrapped labels respect the versioned width estimate and canvas height.
        if node.get("data-role") == "conditions":
            assert len(node.text) * size * style["character_width_em"] <= style["arrow_width"]


def test_unknown_annotations_stay_in_metadata_without_visible_placeholders() -> None:
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("Input", "CCO"),), (SchemeMolecule("Output", "CC=O"),),
        title="Hypothesis", basis="proposed",
        conditions=(SchemeAnnotation("Conditions not supplied", "unknown"),),
        yield_info=SchemeAnnotation("Not established", "unknown"),
    ))
    ns = "{http://www.w3.org/2000/svg}"
    assert not root.findall(ns + "text[@data-role='conditions']")
    assert not root.findall(ns + "text[@data-role='yield']")
    assert "unknown: Conditions not supplied" in root.find(ns + "desc").text
    assert "Yield (unknown): Not established" in root.find(ns + "desc").text


def test_compact_layout_shortens_arrow_without_scaling_molecule_panels() -> None:
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("Input", "CCO"),), (SchemeMolecule("Output", "CC=O"),),
        title="Compact", basis="reported",
    ))
    style = load_annotated_scheme_style()
    assert style["arrow_width"] == 240
    assert style["molecule_preset"] == "web_consistent"
    assert style["molecule_width"] == 100 and style["molecule_height"] == 100
    assert float(root.get("width")) == 2 * 12 + 2 * 100 + 2 * 20 + 240
    assert float(root.get("height")) == 124


def test_hydrolysis_scheme_omits_names_substrate_ids_and_operating_details() -> None:
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("B: ethyl 2-acetamido-1,3-oxazole-4-carboxylate", "CCOC(=O)c1coc(NC(C)=O)n1", "B"),
         SchemeMolecule("Water", "O", "water")),
        (SchemeMolecule("C: 2-acetamido-1,3-oxazole-4-carboxylic acid", "O=C(O)c1coc(NC(C)=O)n1"),),
        title="Ester hydrolysis", basis="proposed",
        conditions=(SchemeAnnotation("B, lithium hydroxide (2 equiv), substrate concentration 0.1 M, 25 °C, 2 h", "proposed"),),
    ))
    ns = "{http://www.w3.org/2000/svg}"
    labels = root.findall(ns + "text[@data-role='conditions']")
    assert " ".join(node.text for node in labels) == "lithium hydroxide"
    assert not root.findall(ns + "text[@data-role='molecule-name']")
    assert "substrate concentration 0.1 M" in root.find(ns + "desc").text
    assert "ethyl 2-acetamido-1,3-oxazole-4-carboxylate" in root.find(ns + "desc").text
    assert float(root.get("width")) < 900
    assert float(root.get("height")) < 200


def test_common_ligand_formula_remains_intact_on_compact_arrow() -> None:
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("Input", "CCO"),), (SchemeMolecule("Output", "CC=O"),),
        title="Supplied catalyst", basis="reported",
        conditions=(SchemeAnnotation("Pd[P(t-Bu)3]2, K2CO3, THF", "reported"),),
    ))
    labels = root.findall("{http://www.w3.org/2000/svg}text[@data-role='conditions']")
    assert any("Pd[P(t-Bu)3]2" in (node.text or "") for node in labels)


@pytest.mark.parametrize("key,value", [
    ("font_size", 0), ("gap", True), ("schema_version", "2.0"),
    ("bottom_padding", 0), ("max_annotation_lines", False), ("line_height", 1), ("arrow_width", 1),
    ("character_width_em", True), ("character_width_em", 2),
])
def test_scheme_definitions_are_validated(monkeypatch, key, value) -> None:
    definition = load_annotated_scheme_style()
    definition[key] = value
    monkeypatch.setattr(Path, "read_text", lambda *a, **k: json.dumps(definition))
    with pytest.raises(ValueError):
        load_annotated_scheme_style()
