"""External proposals preserve chemistry uncertainty, rank and graph identity."""

import json

import pytest

from core_retrosynthesis.supplier_predictions import normalize_supplier_predictions


def row(reactants, product="COC", **kwargs):
    return {"product": product, "reactants": reactants, **kwargs}


def test_normalization_deduplicates_order_and_mapping_but_preserves_multiplicity():
    result = normalize_supplier_predictions("COC", [
        row(["[CH3:1]Br", "CO"], score=2, probability=0.8),
        row(["OC", "BrC"]), row(["CO", "CBr", "CO"]),
    ])
    assert len(result["candidates"]) == 2
    first, second = result["candidates"]
    assert [first["supplier_rank"], second["supplier_rank"]] == [1, 3]
    assert first["model_probability"] == 0.8
    assert first["evidence_kind"] == "predicted"
    assert first["validation_status"] == "structures_parsed_reaction_unverified"
    assert len(second["precursor_smiles"]) == 3
    assert result["duplicates"][0]["candidate_id"] == first["candidate_id"]
    reordered = normalize_supplier_predictions("COC", [row(["OC", "BrC"])])
    assert reordered["candidates"][0]["candidate_id"] == first["candidate_id"]


def test_invalid_and_conflicting_results_remain_auditable():
    result = normalize_supplier_predictions("COC", [
        row(["bad"]), row(["COC", "O"]), row(["C", "O"], product="CO"),
        row(["C"], probability=float("nan")), row([]),
        row(["C"], individual_ranks={"model": -1}),
    ])
    assert result["candidates"] == []
    assert len(result["rejected"]) == 6
    assert all("raw" in item and "reason" in item for item in result["rejected"])
    json.dumps(result, allow_nan=False)


def test_unsupported_chemistry_remains_unverified_and_stereo_distinct():
    result = normalize_supplier_predictions("CCO", [
        row(["N"], product="CCO"),
        row(["F[C@H](Cl)Br"], product="CCO"),
        row(["F[C@@H](Cl)Br"], product="CCO"),
    ])
    assert len(result["candidates"]) == 3
    assert len({item["candidate_id"] for item in result["candidates"]}) == 3
    assert all(item["validation_status"].endswith("unverified") for item in result["candidates"])


@pytest.mark.parametrize("target", ["", "not smiles", "[CH5]"])
def test_invalid_target_is_rejected(target):
    with pytest.raises(ValueError):
        normalize_supplier_predictions(target, [])
