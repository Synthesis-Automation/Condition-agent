"""Forward library, product prediction, and competition regressions."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis import build_generic_library
from forward_synthesis import (
    ForwardOperatorLibrary,
    assess_proposed_step,
    build_forward_library,
    load_forward_library,
    predict_products,
    save_forward_library,
)


S_ALKYLATION = "Cc1[nH]cnc1CCl.NCCS>>NCCSCc1c(C)[nH]cn1"
N_ALKYLATION = "Cc1[nH]cnc1CCl.NCCS>>SCCNCc1c(C)[nH]cn1"


def _library(*reactions: str) -> ForwardOperatorLibrary:
    generic = build_generic_library(
        (
            {
                "reaction_id": f"reaction-{index}",
                "reference_id": f"reference-{index}",
                "reaction_smiles": reaction,
            }
            for index, reaction in enumerate(reactions)
        ),
        levels=("L1", "L2"),
        admission_mode="data_driven",
    )
    return build_forward_library(generic)


def test_forward_library_requires_independent_source_round_trip() -> None:
    library = _library("CCBr.N>>CCN")

    assert library.source_template_count == 2
    assert library.admitted_operator_count == 2
    assert library.rejection_counts == {}
    assert all(operator.forward_smarts for operator in library.operators)


def test_forward_library_round_trips_through_json_gzip(tmp_path) -> None:
    library = _library("CCBr.N>>CCN")
    path = tmp_path / "forward-library.json.gz"

    save_forward_library(library, path)
    restored = load_forward_library(path)

    assert restored == library


def test_blind_prediction_recovers_product_without_target_input() -> None:
    library = _library("CCBr.N>>CCN")

    result = predict_products("N.CCBr", library)

    assert result.valid
    assert result.status == "predicted"
    assert tuple(candidate.product_smiles for candidate in result.candidates) == (
        "CCN",
    )
    assert result.candidates[0].operator_edit_agreement
    assert result.candidates[0].reverse_round_trip
    assert len(result.candidates[0].alternative_pathway_ids) == 1
    assert result.diagnostics.valid_pathway_count == 2
    assert "CONDITIONS_NOT_SUPPLIED_PRODUCTS_ARE_POSSIBILITIES" in result.warnings


def test_competing_endpoint_products_are_retained_and_route_is_competitive() -> None:
    library = _library(S_ALKYLATION, N_ALKYLATION)
    starting_materials = "Cc1[nH]cnc1CCl.NCCS"

    prediction = predict_products(starting_materials, library)
    assessment = assess_proposed_step(
        starting_materials,
        "NCCSCc1c(C)[nH]cn1",
        library,
    )

    assert prediction.diagnostics.unique_product_count == 2
    assert len(prediction.candidates) == 2
    assert {group.competition_level for group in prediction.competition_groups} == {
        "operator",
        "site",
        "product",
    }
    assert assessment.intended_product_rank == 2
    assert assessment.disposition == "competitive"
    assert assessment.best_competitor_product is not None


def test_targeted_replay_does_not_rescue_blind_operator_mismatch() -> None:
    library = _library("CCBr.N>>CCN")
    corrupted = replace(
        library.operators[0],
        edit_tokens=("formed:C-O:NONE>SINGLE",),
    )
    bad_library = replace(
        library,
        operators=(corrupted,),
        admitted_operator_count=1,
        precursor_index=replace(
            library.precursor_index,
            component_count_to_operator_ids={2: (corrupted.forward_operator_id,)},
            operator_required_atomic_numbers={
                corrupted.forward_operator_id: library.precursor_index.operator_required_atomic_numbers[
                    library.operators[0].forward_operator_id
                ]
            },
        ),
    )

    result = predict_products("CCBr.N", bad_library)

    assert not result.candidates
    assert result.diagnostics.operator_edit_mismatch_count == 1


def test_invalid_query_returns_typed_failure() -> None:
    result = predict_products("not-smiles", _library("CCBr.N>>CCN"))

    assert not result.valid
    assert result.error == "INVALID_STARTING_MATERIALS"


def test_recipe_compatibility_is_audited_not_presented_as_yield() -> None:
    result = predict_products(
        "CCBr.N",
        _library("CCBr.N>>CCN"),
        recipe={"temperature_c": 25.0, "atmosphere": "nitrogen"},
    )

    assert result.conditions_supplied
    assert result.candidates[0].recipe_evidence.evaluated
    assert result.candidates[0].recipe_evidence.compatible
    assert result.candidates[0].recipe_evidence.score is not None
