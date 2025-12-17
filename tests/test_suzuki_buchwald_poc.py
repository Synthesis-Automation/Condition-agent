from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemtools.taxonomy.new.feature_engine import (
    FeatureValidationError,
    classify_reaction_smiles,
    compute_atomic_features,
    compute_derived_features,
    load_feature_definitions,
    validate_features,
)
from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available


pytestmark = pytest.mark.skipif(not rdkit_available(), reason="RDKit not available or disabled")


POC_DIR = Path(__file__).resolve().parents[1] / "chemtools" / "taxonomy" / "new"


@pytest.fixture(scope="module")
def poc_defs() -> dict[str, object]:
    atomic_defs, derived_defs = load_feature_definitions(
        atomic_path=POC_DIR / "calculable_features.atomic.suzuki_buchwald.poc.json",
        derived_path=POC_DIR / "calculable_features.derived.suzuki_buchwald.poc.json",
    )
    reactant_types = json.loads((POC_DIR / "reactant_types.suzuki_buchwald.poc.json").read_text(encoding="utf-8"))
    reaction_types = json.loads((POC_DIR / "reaction_types.suzuki_buchwald.poc.json").read_text(encoding="utf-8"))
    templates = json.loads((POC_DIR / "smarts_templates.suzuki_buchwald.poc.json").read_text(encoding="utf-8"))

    validate_features(
        atomic_defs,
        derived_defs,
        reactant_types=reactant_types,
        reaction_types=reaction_types,
    )

    return {
        "atomic_defs": atomic_defs,
        "derived_defs": derived_defs,
        "reactant_types": reactant_types,
        "reaction_types": reaction_types,
        "templates": templates,
    }


def test_atomic_feature_matches(poc_defs: dict[str, object]) -> None:
    atomic_defs = poc_defs["atomic_defs"]
    assert isinstance(atomic_defs, list)

    cases = [
        ("Brc1ccccc1", "aromatic_bromide_present", True),
        ("Clc1ccccc1", "aromatic_chloride_present", True),
        ("Ic1ccccc1", "aromatic_iodide_present", True),
        ("C=CBr", "vinyl_bromide_present", True),
        ("OB(O)c1ccccc1", "aromatic_boronic_acid_present", True),
        ("Nc1ccccc1", "amine_with_nh_present", True),
        ("CCN(CC)CC", "tertiary_amine_present", True),
        ("CCN(CC)CC", "amine_with_nh_present", False),
    ]

    for smiles, token, expected in cases:
        mol = parse_smiles(smiles)
        assert mol is not None, f"Failed to parse SMILES: {smiles}"
        features = compute_atomic_features(mol, atomic_defs)
        assert features.get(token, False) is expected, f"{smiles} -> {token}"


def test_derived_logic_sanity() -> None:
    derived_defs = [
        {"token": "x_present", "type": "bool", "derive": {"any_of": ["a_present", "b_present"]}},
        {"token": "y_present", "type": "bool", "derive": {"all_of": ["x_present"], "none_of": ["c_present"]}},
    ]
    base = {"a_present": False, "b_present": True, "c_present": False}
    out = compute_derived_features(base, derived_defs)
    assert out["x_present"] is True
    assert out["y_present"] is True


def test_derived_cycle_detection() -> None:
    derived_defs = [
        {"token": "a_present", "type": "bool", "derive": {"any_of": ["b_present"]}},
        {"token": "b_present", "type": "bool", "derive": {"any_of": ["a_present"]}},
    ]
    with pytest.raises(FeatureValidationError):
        validate_features([], derived_defs)


def test_reaction_classification_suzuki_vs_buchwald(poc_defs: dict[str, object]) -> None:
    atomic_defs = poc_defs["atomic_defs"]
    derived_defs = poc_defs["derived_defs"]
    reaction_types = poc_defs["reaction_types"]
    assert isinstance(atomic_defs, list)
    assert isinstance(derived_defs, list)
    assert isinstance(reaction_types, list)

    suzuki = classify_reaction_smiles(
        ["Brc1ccccc1", "OB(O)c1ccccc1"],
        atomic_defs=atomic_defs,
        derived_defs=derived_defs,
        reaction_types=reaction_types,
    )
    suzuki_ids = sorted([m.id for m in (suzuki["reaction_type_matches"] or [])])
    assert suzuki_ids == ["SuzukiMiyaura"]

    buchwald = classify_reaction_smiles(
        ["Brc1ccccc1", "Nc1ccccc1"],
        atomic_defs=atomic_defs,
        derived_defs=derived_defs,
        reaction_types=reaction_types,
    )
    buchwald_ids = sorted([m.id for m in (buchwald["reaction_type_matches"] or [])])
    assert buchwald_ids == ["BuchwaldHartwigAmination"]

    none = classify_reaction_smiles(
        ["OB(O)c1ccccc1", "Nc1ccccc1"],
        atomic_defs=atomic_defs,
        derived_defs=derived_defs,
        reaction_types=reaction_types,
    )
    none_ids = sorted([m.id for m in (none["reaction_type_matches"] or [])])
    assert none_ids == []

    tertiary_amine = classify_reaction_smiles(
        ["Brc1ccccc1", "CCN(CC)CC"],
        atomic_defs=atomic_defs,
        derived_defs=derived_defs,
        reaction_types=reaction_types,
    )
    tertiary_ids = sorted([m.id for m in (tertiary_amine["reaction_type_matches"] or [])])
    assert tertiary_ids == []


def test_template_provenance_matches_atomic_file(poc_defs: dict[str, object]) -> None:
    atomic_defs = poc_defs["atomic_defs"]
    templates = poc_defs["templates"]
    assert isinstance(atomic_defs, list)
    assert isinstance(templates, dict)

    atomic_by_token = {str(e.get("token")): e for e in atomic_defs if isinstance(e, dict)}
    generated = templates.get("generated_atomic_features", [])
    assert isinstance(generated, list)

    for entry in generated:
        assert isinstance(entry, dict)
        token = str(entry.get("token", ""))
        smarts = str(entry.get("smarts", ""))
        assert token in atomic_by_token
        detect = atomic_by_token[token].get("detect") or {}
        smarts_any = detect.get("smarts_any") or []
        assert smarts in smarts_any

