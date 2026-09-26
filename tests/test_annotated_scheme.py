"""Annotated SVG display preserves explicit chemistry and epistemic labels."""

import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from visualization import (
    SchemeAnnotation, SchemeMolecule, load_annotated_scheme_style,
    render_annotated_scheme_svg,
)


def test_scheme_has_named_structures_conditions_and_separate_yield_status() -> None:
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("Bromoethane", "CCBr"), SchemeMolecule("Ammonia", "N")),
        (SchemeMolecule("Ethylamine", "CCN"),), title="Illustrative substitution", basis="proposed",
        conditions=(SchemeAnnotation("Solvent <unknown> & temperature unspecified", "unknown"),),
        yield_info=SchemeAnnotation("Not established", "unknown"),
    ))
    text = " ".join(root.itertext())
    assert "Bromoethane" in text and "Ethylamine" in text
    assert "Unknown: Solvent <unknown>" in text
    assert "Yield (unknown):" in text
    assert "Proposed transformation" in text
    assert "Not a feasibility assessment" in text
    ns = "{http://www.w3.org/2000/svg}"
    assert len(root.findall(ns + "g" + "[@data-role='molecule']")) == 3
    assert len(root.findall(".//" + ns + "path")) > 10
    assert root.get("data-definition") == "annotated_scheme.v1"
    assert not root.findall(".//" + ns + "unknown")


def test_missing_conditions_are_not_invented_and_bad_structures_fail() -> None:
    text = render_annotated_scheme_svg(
        (SchemeMolecule("Input", "CCO"),), (SchemeMolecule("Output", "CC=O"),),
        title="Hypothesis", basis="proposed",
    ).decode()
    assert "Conditions not supplied" in text and "Yield not supplied" in text
    with pytest.raises(ValueError, match="Invalid SMILES"):
        render_annotated_scheme_svg(
            (SchemeMolecule("Bad", "C1CC"),), (SchemeMolecule("Output", "CC=O"),),
            title="Invalid", basis="unknown",
        )
    with pytest.raises(ValueError, match="explicit reactants"):
        render_annotated_scheme_svg((), (), title="Empty", basis="unknown")


def test_long_annotations_are_bounded_visually_and_complete_in_description() -> None:
    long = "Long unverified condition annotation. " * 90
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("A long molecule name " * 8, "C[C@H](O)Cl"),),
        (SchemeMolecule("Product", "CC(=O)Cl"),), title="Long labels", basis="proposed",
        conditions=(SchemeAnnotation(long, "proposed"),),
    ))
    ns = "{http://www.w3.org/2000/svg}"
    assert long in root.find(ns + "desc").text
    assert "C[C@H](O)Cl" in root.find(ns + "desc").text
    assert "More conditions in step details" in " ".join(root.itertext())
    assert float(root.get("height")) < 900
    width, height = float(root.get("width")), float(root.get("height"))
    for node in root.findall(ns + "text"):
        assert 0 <= float(node.get("x")) <= width
        assert 0 <= float(node.get("y")) <= height


@pytest.mark.parametrize("key,value", [("font_size", 0), ("gap", True), ("schema_version", "2.0")])
def test_scheme_definitions_are_validated(monkeypatch, key, value) -> None:
    definition = load_annotated_scheme_style()
    definition[key] = value
    monkeypatch.setattr(Path, "read_text", lambda *a, **k: json.dumps(definition))
    with pytest.raises(ValueError):
        load_annotated_scheme_style()
