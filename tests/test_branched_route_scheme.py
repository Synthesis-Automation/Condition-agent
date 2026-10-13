"""Branched chemical drawings preserve supplied structures and dependency edges."""

from dataclasses import asdict, replace
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import pytest

from visualization import (
    BranchedRouteStep, SchemeAnnotation, SchemeMolecule,
    load_branched_route_style, render_branched_route_scheme_svg,
)


NS = "{http://www.w3.org/2000/svg}"


def route() -> tuple[BranchedRouteStep, ...]:
    """Two independent preparations converge after one extra branch step."""
    alcohol = SchemeMolecule("Alcohol", "CCO", "alcohol")
    aldehyde = SchemeMolecule("Aldehyde", "CC=O", "aldehyde")
    nitro = SchemeMolecule("Nitro", "C[N+](=O)[O-]", "nitro")
    amine = SchemeMolecule("Amine", "CN", "amine")
    ketone = SchemeMolecule("Ketone", "CC(C)=O", "ketone")
    product = SchemeMolecule("Product", "CC(C)NC", "product")
    return (
        BranchedRouteStep("amine", (nitro,), (amine,), title="Prepare amine"),
        BranchedRouteStep("aldehyde", (alcohol,), (aldehyde,), title="Prepare aldehyde"),
        BranchedRouteStep("ketone", (aldehyde,), (ketone,), ("aldehyde",)),
        BranchedRouteStep("join", (amine, ketone), (product,), ("amine", "ketone")),
    )


def test_convergence_retains_all_reactions_and_only_declared_dependencies() -> None:
    steps = route()
    original = [asdict(s) for s in steps]
    svg = render_branched_route_scheme_svg(steps)
    assert svg == render_branched_route_scheme_svg(steps)
    root = ET.fromstring(svg)
    assert root.get("data-layout") == "layered-dependencies"
    panels = root.findall(f"{NS}g[@data-role='route-step']")
    assert [p.get("data-depth") for p in panels] == ["0", "0", "1", "2"]
    dependencies = root.findall(f"{NS}g[@data-role='dependency']")
    assert {(e.get("data-parent"), e.get("data-child")) for e in dependencies} == {
        ("aldehyde", "ketone"), ("amine", "join"), ("ketone", "join"),
    }
    assert len(root.findall(f".//{NS}path[@data-role='reaction-arrow']")) == 4
    # One vector drawing per supplied component, including repeated intermediates.
    assert len(root.findall(f".//{NS}g[@data-role='molecule']")) == 9
    assert not root.findall(f".//{NS}image")
    assert json.loads(root.find(NS + "metadata").text)["steps"] == json.loads(json.dumps(original))
    assert [asdict(s) for s in steps] == original
    for panel in panels:
        x, y = map(float, re.fullmatch(
            r"translate\(([^,]+),([^\)]+)\)", panel.get("transform"),
        ).groups())
        assert 0 <= x <= x + float(panel.get("data-width")) <= float(root.get("width"))
        assert 0 <= y <= y + float(panel.get("data-height")) <= float(root.get("height"))
    # The dependency spanning a layer goes through the left exterior gutter.
    long_edge = next(e for e in dependencies if e.get("data-parent") == "amine")
    assert "H 24 " in long_edge.find(NS + "path").get("d")


def test_unordered_input_and_fanout_keep_topological_layout() -> None:
    steps = route()
    fork = replace(steps[2], id="alternate")
    root = ET.fromstring(render_branched_route_scheme_svg((steps[3], fork, *steps[:3])))
    depths = {p.get("data-step-id"): p.get("data-depth")
              for p in root.findall(f"{NS}g[@data-role='route-step']")}
    assert depths == {"join": "2", "alternate": "1", "amine": "0", "aldehyde": "0", "ketone": "1"}
    assert len(root.findall(f"{NS}g[@data-parent='aldehyde']")) == 2


@pytest.mark.parametrize("case", ["duplicate", "unknown", "unshared", "cycle", "conflict", "invalid"])
def test_invalid_or_conflicting_graph_is_not_drawn_as_a_route(case) -> None:
    steps = list(route())
    if case == "duplicate":
        steps.append(steps[0])
    elif case == "unknown":
        steps[3] = replace(steps[3], after_step_ids=("missing",))
    elif case == "unshared":
        steps[1] = replace(steps[1], after_step_ids=("amine",))
    elif case == "cycle":
        steps[0] = replace(steps[0], reactants=steps[3].products, after_step_ids=("join",))
    elif case == "conflict":
        steps[3] = replace(steps[3], reactants=(replace(steps[0].products[0], smiles="CO"),))
    else:
        steps[0] = replace(steps[0], reactants=(SchemeMolecule("Invalid", "C1"),))
    with pytest.raises(ValueError):
        render_branched_route_scheme_svg(steps)


def test_labels_are_escaped_and_reagents_only_keeps_full_attribution() -> None:
    steps = list(route())
    steps[0] = replace(
        steps[0], title="<script>long title & label</script> " * 10,
        basis="proposed", conditions=(SchemeAnnotation("Procedure only", "reported"),),
        reagents=(SchemeAnnotation("Pd/C", "proposed"),),
        yield_info=SchemeAnnotation("70%", "reported"),
    )
    root = ET.fromstring(render_branched_route_scheme_svg(steps, reagents_only=True))
    assert not root.findall(f".//{NS}script")
    text = " ".join(n.text or "" for n in root.iter(NS + "text"))
    assert "Pd/C" in text and "Procedure only" not in text and "70%" not in text
    metadata = json.loads(root.find(NS + "metadata").text)
    assert metadata["steps"][0]["conditions"][0]["basis"] == "reported"
    assert metadata["steps"][0]["yield_info"]["text"] == "70%"
    assert metadata["steps"][0]["title"] == steps[0].title


@pytest.mark.parametrize("key,value", [
    ("schema_version", "2"), ("definition_id", "other"), ("padding", 0),
    ("row_gap", True), ("lane_gap", 1), ("line_height", 1), ("unexpected", 1),
])
def test_definition_validation(monkeypatch, key, value) -> None:
    definition = asdict(load_branched_route_style())
    definition[key] = value
    monkeypatch.setattr(Path, "read_text", lambda *args, **kwargs: json.dumps(definition))
    with pytest.raises(ValueError):
        load_branched_route_style()
