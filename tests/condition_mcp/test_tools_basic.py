from __future__ import annotations

import json
import pathlib

import pytest

from condition_mcp.tools import detect_family, featurize_substrates, normalize_reaction


@pytest.fixture(scope="module")
def sample_reactants() -> list[str]:
    # Typical Buchwald/Ullmann-style coupling example
    return ["Clc1ccc(Cl)cc1", "Nc1ccccc1"]


def test_normalize_reaction_returns_canonical_payload():
    reaction = "Clc1ccc(Cl)cc1.Nc1ccccc1>>Clc1ccc(Nc2ccccc2)cc1"
    result = normalize_reaction({"smiles_rxn": reaction})

    assert result["schema_version"] == "0.1.0"
    assert result["normalized"].count(">") == 2
    assert len(result["reactants"]) == 2
    assert "Cl" in result["reactants"][0]
    assert result["products"][0]
    assert result["errors"] == []


def test_detect_family_identifies_cn_coupling(sample_reactants: list[str]):
    outcome = detect_family({"reactants": sample_reactants})

    assert outcome["schema_version"] == "0.1.0"
    assert outcome["family"] in {"Ullmann_CN", "Buchwald_CN"}
    assert 0.0 <= outcome["confidence"] <= 1.0
    assert isinstance(outcome["hits"], dict)


def test_featurize_substrates_produces_descriptors(sample_reactants: list[str]):
    payload = featurize_substrates({"reactants": sample_reactants})

    assert payload["schema_version"] == "0.1.0"
    assert payload["electrophile"]
    assert "descriptors" in payload
    assert "bin" in payload["descriptors"]


def test_condition_set_schema_is_well_formed():
    schema_path = pathlib.Path(__file__).parents[2] / "condition_mcp" / "resources" / "schemas" / "condition_set.json"
    data = json.loads(schema_path.read_text(encoding="utf-8"))

    assert data["title"] == "ConditionSet"
    assert "core" in data["properties"]
    assert data["properties"]["components"]["items"]["required"] == ["role", "identifier"]
