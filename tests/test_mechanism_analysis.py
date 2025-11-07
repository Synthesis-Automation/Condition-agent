from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemtools.mechanism import classify_mechanism_simple
from chemtools.mechanism.renderer import build_mechanism_narrative
from chem_assistant.chemtools_wrapper import analyze_mechanism_tool


GOLDEN_PATH = Path(__file__).parent / "data" / "mechanism_golden_set.json"


@pytest.fixture(scope="module")
def golden_reactions() -> list[dict]:
    with GOLDEN_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def test_classify_mechanism_simple_matches_golden(golden_reactions: list[dict]) -> None:
    for entry in golden_reactions:
        result = classify_mechanism_simple(entry["expected_family"])
        assert result["mechanism_type"] == entry["expected_mechanism"]
        assert result["steps"], f"No steps generated for {entry['name']}"
        assert result["confidence"] > 0.3


def test_classify_mechanism_simple_with_bond_evidence() -> None:
    bond_changes = {
        "broken_bonds": [[1, 7]],
        "formed_bonds": [[1, 8]],
        "changed_atoms": {1, 7, 8},
        "leaving_groups": ["Br"],
    }

    result = classify_mechanism_simple(
        "suzuki_miyaura",
        bond_changes=bond_changes,
        detail_level="high",
    )

    assert result["steps"][-1]["title"] == "Bond change summary"
    sources = {ref["source"] for ref in result["evidence_refs"]}
    assert "bond_analysis" in sources
    assert "reaction_detection" in sources
    assert not result["warnings"]


def test_analyze_mechanism_tool_returns_structured_payload(golden_reactions: list[dict]) -> None:
    payload = {
        "reaction_smiles": golden_reactions[0]["reaction_smiles"],
        "detail_level": "summary",
    }

    first = analyze_mechanism_tool.invoke(payload)
    assert first["success"] is True
    assert first["mechanism_type"] == golden_reactions[0]["expected_mechanism"]
    assert first["steps"], "Mechanism steps missing"
    assert first["telemetry"]["cache_hit"] is False
    assert first["narrative"]
    assert first["highlights"]
    assert first["normalized_reaction"]

    second = analyze_mechanism_tool.invoke(payload)
    assert second["success"] is True
    assert second["telemetry"]["cache_hit"] is True
    assert second["mechanism_type"] == first["mechanism_type"]
    assert second["narrative"] == first["narrative"]


def test_renderer_adds_context_and_warnings() -> None:
    mechanism = classify_mechanism_simple("snar")
    detection = {"family": "snar", "status": "conflict"}
    bond_changes = {"broken_bonds": [], "formed_bonds": []}
    context = {"solvent": "DMF", "temperature_c": 90}

    rendered = build_mechanism_narrative(
        mechanism,
        detection=detection,
        bond_changes=bond_changes,
        context=context,
    )

    assert "DMF" in rendered["narrative"]
    assert rendered["highlights"], "Expected highlights for renderer"
    assert any("disagreement" in warning.lower() for warning in rendered["warnings"])
