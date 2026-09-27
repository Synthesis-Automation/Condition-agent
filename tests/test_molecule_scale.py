"""Molecule size must reflect graph extent, not independent box-filling zoom."""

import base64
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import pytest
from fastapi.testclient import TestClient

from app.web_api.main import create_app
from app.web_api.runtime import LocalRecommendationRuntime
from app.web_api.scientific_presentation import present_message
from core_retrosynthesis.html_report import molecule_svg, reaction_svg
from visualization import (
    SchemeMolecule, render_annotated_scheme_svg, render_molecule_image_bytes,
    render_reaction_image_bytes,
)


APIXABAN = "COc1ccc(-n2nc(C(N)=O)c3c2C(=O)N(c2ccc(N4CCCCC4=O)cc2)CC3)cc1"
EMTRICITABINE = "Nc1nc(=O)n([C@@H]2CS[C@H](CO)O2)cc1F"
NS = "{http://www.w3.org/2000/svg}"


def longest_bond(root: ET.Element) -> float:
    """Measure complete carbon bond segments, excluding atom glyphs."""
    lengths = []
    for node in root.iter(NS + "path"):
        if not node.get("class", "").startswith("bond-"):
            continue
        match = re.match(r"M ([\d.-]+),([\d.-]+) L ([\d.-]+),([\d.-]+)", node.get("d", ""))
        if match:
            values = tuple(map(float, match.groups()))
            lengths.append(math.dist(values[:2], values[2:]))
    return max(lengths)


def test_expanding_canvas_keeps_bonds_and_atom_glyphs_at_the_same_scale() -> None:
    for smiles in ("CCl", APIXABAN, EMTRICITABINE):
        roots = [ET.fromstring(render_molecule_image_bytes(
            smiles, size=size, image_format="svg", render_preset="web_consistent",
            expand_canvas=True,
        )) for size in ((10, 10), (800, 600))]
        paths = [[node.attrib for node in root.iter(NS + "path")] for root in roots]
        assert paths[0] == paths[1], "canvas size must not rescale bonds or labels"
        assert float(roots[0].get("width").removesuffix("px")) > 10


def test_workspace_and_example_drugs_use_the_same_bond_scale() -> None:
    view = present_message(f"`{APIXABAN}` and `{EMTRICITABINE}`", "a" * 32)
    assert len(view["structures"]) == 2
    for item in view["structures"]:
        root = ET.fromstring(base64.b64decode(item["image_url"].split(",", 1)[1]))
        assert longest_bond(root) == pytest.approx(30, abs=1.1)
    gallery = (Path(__file__).resolve().parents[1] / "examples" /
               "top_20_selling_small_molecule_drugs_smiles.html").read_text("utf-8")
    for name in ("Apixaban", "Emtricitabine"):
        svg = re.search(rf'aria-label="2D molecular structure of {name}">\s*(<svg.*?</svg>)', gallery, re.S)
        assert svg is not None
        assert longest_bond(ET.fromstring(svg.group(1))) == pytest.approx(30, abs=1.1)


def test_annotated_scheme_grows_molecule_panels_without_shrinking_compounds() -> None:
    root = ET.fromstring(render_annotated_scheme_svg(
        (SchemeMolecule("Apixaban", APIXABAN),),
        (SchemeMolecule("Emtricitabine", EMTRICITABINE),),
        title="Rendering scale fixture, not a chemical transformation", basis="unknown",
    ))
    groups = root.findall(NS + "g[@data-role='molecule']")
    assert len(groups) == 2
    for group in groups:
        assert longest_bond(group) == pytest.approx(30, abs=1.1)
        rect = group.find(NS + "rect")
        x, y = map(float, re.search(r"translate\(([^,]+),([^\)]+)\)", group.get("transform")).groups())
        assert x >= 0 and y >= 0
        assert x + float(rect.get("width")) <= float(root.get("width"))
        assert y + float(rect.get("height")) <= float(root.get("height"))


def test_reaction_canvas_grows_to_fit_tall_molecules_and_agents() -> None:
    root = ET.fromstring(render_reaction_image_bytes(
        f"{APIXABAN}>{APIXABAN}>{APIXABAN}", size=(400, 100),
        image_format="svg", render_preset="web_consistent",
    ))
    height = float(root.get("height").removesuffix("px"))
    for group in root.findall(NS + "g"):
        x, y = map(float, re.search(r"translate\(([^ ]+) ([^\)]+)\)", group.get("transform")).groups())
        rect = group.find(NS + "rect")
        assert x >= 0 and y >= 0
        assert y + float(rect.get("height").removesuffix("px")) <= height


@pytest.mark.parametrize("options", [{"image_format": "png", "render_preset": "web_consistent"},
                                     {"image_format": "svg", "render_preset": "current"}])
def test_expanding_canvas_requires_a_vector_drawing_and_defined_scale(options) -> None:
    with pytest.raises(ValueError, match="expand_canvas requires"):
        render_molecule_image_bytes("CCO", expand_canvas=True, **options)


def test_html_report_helpers_keep_intrinsic_scale_despite_report_css() -> None:
    for smiles in (APIXABAN, EMTRICITABINE):
        for document in (molecule_svg(smiles, width=160, height=100),
                         reaction_svg(f"{smiles}>>{smiles}")):
            container = ET.fromstring(document)
            assert container.get("tabindex") == "0"
            assert "overflow:auto" in container.get("style")
            drawing = container.find(NS + "svg")
            assert longest_bond(drawing) == pytest.approx(30, abs=1.1)
            assert "max-width:none;max-height:none" in drawing.get("style")
            assert f"width:{float(drawing.get('width').removesuffix('px')):g}px" in drawing.get("style")
    assert "drawing failed" in molecule_svg("not-a-smiles")


def test_web_molecule_endpoint_expands_small_requested_canvases() -> None:
    web = TestClient(create_app(runtime=LocalRecommendationRuntime()))
    for smiles in (APIXABAN, EMTRICITABINE):
        response = web.post("/api/v1/render/molecule", json={
            "molecule_smiles": smiles, "width": 160, "height": 100,
        })
        assert response.status_code == 200
        drawing = ET.fromstring(response.content)
        assert longest_bond(drawing) == pytest.approx(30, abs=1.1)
        assert float(drawing.get("width").removesuffix("px")) > 160


def test_scientific_workspace_api_keeps_scale_in_saved_cards_and_schemes() -> None:
    attribution = {"basis": "unknown", "source_ids": [], "limitations": []}
    answer = {
        "schema_version": "scientific_answer.v2", "answer_markdown": "Rendering fixture only",
        "evidence_refs": [], "uncertainties": [], "needs_user_input": False,
        "sources": [], "target_molecule_ids": [], "routes": [], "claims": [],
        "molecules": [{"id": name, "name": name, "smiles": smiles, **attribution}
                      for name, smiles in (("apixaban", APIXABAN), ("emtricitabine", EMTRICITABINE))],
        "steps": [{"id": "layout", "title": "Layout fixture, not a chemical transformation",
                   "reactant_ids": ["apixaban", "emtricitabine"], "product_ids": ["apixaban"],
                   "after_step_ids": [], "conditions": [], "yield_info": None, **attribution}],
    }

    class SavedService:
        def get(self, identity):
            return {"id": identity, "turns": [{"question": f"`{EMTRICITABINE}`",
                                               "status": "completed", "answer": answer}]}

    web = TestClient(create_app(runtime=object(), scientific_service=SavedService(), recommendation_only=False),
                     base_url="http://127.0.0.1")
    response = web.get("/api/v1/scientific/conversations/" + "a" * 32)
    assert response.status_code == 200
    turn = response.json()["turns"][0]
    view = turn["structured_presentation"]
    for item in [*turn["question_presentation"]["structures"], *view["molecules"]]:
        root = ET.fromstring(base64.b64decode(item["image_url"].split(",", 1)[1]))
        assert longest_bond(root) == pytest.approx(30, abs=1.1)
    step = view["steps"][0]
    root = ET.fromstring(base64.b64decode(step["image_url"].split(",", 1)[1]))
    assert float(root.get("width")) == step["scheme_width"]
    for group in root.findall(NS + "g[@data-role='molecule']"):
        assert longest_bond(group) == pytest.approx(30, abs=1.1)
