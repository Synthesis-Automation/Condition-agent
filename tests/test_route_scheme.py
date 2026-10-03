"""Fixed route geometry, component retention and CLI regressions."""

import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import pytest

from visualization import RouteSchemeStep, SchemeAnnotation, render_route_scheme_svg
from visualization.route_cli import main
from visualization.route_scheme import load_route_scheme_style


NS = "{http://www.w3.org/2000/svg}"
ROUTE = (
    "CCO>>CC=O",
    "CC=O.N>>CCN",
    "CCN.CC(=O)Cl>>CCNC(C)=O",
    "CCNC(C)=O>>CCN(C)C(C)=O",
)


@pytest.mark.parametrize("count", [1, 2, 3, 4])
def test_grid_order_each_intermediate_once_and_no_captions(count) -> None:
    svg = render_route_scheme_svg(ROUTE[:count])
    assert svg == render_route_scheme_svg(ROUTE[:count])
    root = ET.fromstring(svg)
    blocks = root.findall(NS + "g")
    assert len(blocks) == 2 * count + 1
    for i, block in enumerate(blocks):
        assert int(block.get("data-row")) == i // 3
        assert int(block.get("data-column")) == i % 3
        assert block.get("data-role") == ("arrow-block" if i % 2 else "route-molecule")
    assert not root.findall(".//" + NS + "g[@data-role='molecule-name']")
    assert len(root.findall(".//" + NS + "path[@data-role='reaction-arrow']")) == count
    assert not root.findall(".//" + NS + "image")
    assert "scale(" not in svg.decode()
    assert root.get("data-definition") == "route_scheme.v1"


def test_partners_agents_and_other_products_are_retained() -> None:
    root = ET.fromstring(
        render_route_scheme_svg(
            (
                RouteSchemeStep("N.CCBr>O>CCN.Br", product_index=0),
                "CCN.CC(=O)Cl>>CCNC(C)=O",
            ),
            main_reactant_index=1,
        )
    )
    expected = {
        "partners": ["N", "CC(=O)Cl"],
        "agents": ["O"],
        "other-products": ["Br"],
    }
    for role, smiles in expected.items():
        groups = root.findall(
            f".//{NS}g[@data-role='{role}']/{NS}g[@data-role='molecule']"
        )
        assert [node.get("data-smiles") for node in groups] == smiles
    start = root.find(NS + "g[@data-role='route-molecule']")
    assert start.find(NS + "g").get("data-smiles") == "CCBr"


def test_default_start_preserves_all_reactants_without_guessing_main_material() -> None:
    root = ET.fromstring(render_route_scheme_svg(("CCBr.N>>CCN",)))
    start = root.find(NS + "g[@data-role='route-molecule']")
    assert len(start.findall(NS + "g[@data-role='molecule']")) == 2
    assert not root.findall(".//" + NS + "g[@data-role='partners']")


def test_large_structures_long_annotations_and_attribution_fit_without_truncation() -> (
    None
):
    large = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO"
    note = "Long supplied annotation & <unsafe> " * 15
    svg = render_route_scheme_svg(
        (
            RouteSchemeStep(
                f"CCO>>{large}",
                above=(SchemeAnnotation(note, "proposed"),),
                below=(SchemeAnnotation("THF, 25 °C, 16 h", "reported"),),
                yield_info=SchemeAnnotation("70% (illustrative)", "unknown"),
            ),
        ),
        title="<script>alert(1)</script>",
    )
    root = ET.fromstring(svg)
    assert not root.findall(".//" + NS + "script")
    metadata = json.loads(root.find(NS + "metadata").text)
    assert metadata["steps"][0]["above"][0] == {"text": note, "basis": "proposed"}
    assert metadata["steps"][0]["below"][0]["basis"] == "reported"
    texts = root.findall(".//" + NS + "g[@data-role='conditions-above']/" + NS + "text")
    assert texts
    assert "".join("".join(node.text.split()) for node in texts) == "".join(
        note.split()
    )
    assert float(root.get("data-column-width")) > 700
    for block in root.findall(NS + "g"):
        x, y = map(
            float,
            re.fullmatch(
                r"translate\(([^,]+),([^\)]+)\)", block.get("transform")
            ).groups(),
        )
        assert x >= 0 and y - float(block.get("data-above")) >= 0
        assert x + float(block.get("data-width")) <= float(root.get("width"))
        assert y + float(block.get("data-below")) <= float(root.get("height"))
    # Small and large structures retain equal carbon bond lengths.
    bond_lengths = []
    for molecule in root.findall(".//" + NS + "g[@data-role='molecule']"):
        lengths = []
        for node in molecule.iter(NS + "path"):
            if not node.get("class", "").startswith("bond-"):
                continue
            match = re.match(
                r"M ([\d.-]+),([\d.-]+) L ([\d.-]+),([\d.-]+)", node.get("d", "")
            )
            if match:
                coords = tuple(map(float, match.groups()))
                lengths.append(math.dist(coords[:2], coords[2:]))
        bond_lengths.append(max(lengths))
    assert bond_lengths[0] == pytest.approx(bond_lengths[1], abs=0.15)


@pytest.mark.parametrize(
    "key,value",
    [
        ("columns", 4),
        ("columns", True),
        ("schema_version", "2.0"),
        ("row_gap", 0),
        ("line_height", 1),
        ("minimum_column_width", 1),
        ("minimum_arrow_width", 1),
        ("minimum_arrow_width", True),
        ("molecule_preset", "current"),
        ("molecule_preset", []),
        ("character_width_em", False),
        ("arrow_length_scale", False),
        ("arrow_length_scale", 0),
        ("arrow_length_scale", 1.1),
        ("unexpected", 1),
    ],
)
def test_definition_validation(monkeypatch, key, value) -> None:
    from dataclasses import asdict

    definition = asdict(load_route_scheme_style())
    definition[key] = value
    monkeypatch.setattr(Path, "read_text", lambda *a, **kw: json.dumps(definition))
    with pytest.raises(ValueError):
        load_route_scheme_style()


def test_cli_round_trip_and_invalid_input_does_not_overwrite(tmp_path) -> None:
    source, output = tmp_path / "route.json", tmp_path / "scheme.svg"
    source.write_text(
        json.dumps(
            {
                "steps": [
                    {
                        "reaction_smiles": ROUTE[0],
                        "above": [{"text": "Test", "basis": "proposed"}],
                    },
                    ROUTE[1],
                ]
            }
        ),
        encoding="utf-8",
    )
    assert main([str(source), str(output)]) == 0
    original = output.read_bytes()
    assert ET.fromstring(original).get("data-layout") == "three-column-row-major"
    source.write_text('["C>>N", "O>>CO"]', encoding="utf-8")
    with pytest.raises(SystemExit):
        main([str(source), str(output)])
    assert output.read_bytes() == original


def test_reagents_only_preserves_procedures_in_metadata() -> None:
    step = RouteSchemeStep(
        "CCBr.N>O>CCN",
        above=(SchemeAnnotation("Add slowly; heat for 16 h", "reported"),),
        below=(SchemeAnnotation("Extract and purify", "reported"),),
        yield_info=SchemeAnnotation("70%", "reported"),
        reagents=(SchemeAnnotation("K2CO3", "proposed"),),
    )
    root = ET.fromstring(render_route_scheme_svg((step,), reagents_only=True))
    visible = " ".join(node.text or "" for node in root.iter(NS + "text"))
    assert visible == "+ K2CO3"
    assert root.findall(f".//{NS}g[@data-role='agents']")
    metadata = json.loads(root.find(NS + "metadata").text)
    assert metadata["steps"][0]["above"][0]["text"] == step.above[0].text
    assert metadata["steps"][0]["below"][0]["text"] == step.below[0].text
    assert metadata["steps"][0]["yield_info"]["text"] == "70%"
    assert metadata["steps"][0]["reagents"][0]["basis"] == "proposed"


def test_arrow_spacing_does_not_inherit_large_molecule_width() -> None:
    large = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO"
    root = ET.fromstring(render_route_scheme_svg((f"CCO>>{large}",)))
    start, arrow, product = root.findall(NS + "g")
    style = load_route_scheme_style()
    assert float(arrow.get("data-width")) == style.minimum_arrow_width
    assert float(product.get("data-width")) > 700

    def left(node):
        return float(node.get("transform").split("(")[1].split(",")[0])

    assert left(product) - left(arrow) == pytest.approx(
        style.minimum_arrow_width + style.column_gap
    )
    expected = (
        2 * style.padding + style.minimum_column_width
        + style.minimum_arrow_width + float(product.get("data-width"))
        + 2 * style.column_gap
    )
    assert float(root.get("width")) == expected


def test_cli_accepts_explicit_reagent_labels(tmp_path) -> None:
    source, output = tmp_path / "route.json", tmp_path / "scheme.svg"
    source.write_text(json.dumps({
        "reagents_only": True,
        "steps": [{"reaction_smiles": "CCO>>CC=O", "reagents": ["PCC"],
                   "above": ["Procedure retained in metadata"]}],
    }), encoding="utf-8")
    assert main([str(source), str(output)]) == 0
    root = ET.fromstring(output.read_bytes())
    assert [node.text for node in root.iter(NS + "text")] == ["PCC"]
