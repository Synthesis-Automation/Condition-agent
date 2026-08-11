"""Dataset-driven condition-aware selectivity POC regressions."""

from __future__ import annotations

import pytest

from core_retrosynthesis_poc import (
    ConditionalEditChoiceModel,
    build_reaction_choice_set,
    build_reaction_choice_set_from_record,
    condition_tokens_from_mapping,
    condition_tokens_from_recipe,
    detect_functional_group_competition,
)


S_ALKYLATION = (
    "Cc1[nH]cnc1CCl.NCCS>>NCCSCc1c(C)[nH]cn1"
)
N_ALKYLATION = (
    "Cc1[nH]cnc1CCl.NCCS>>SCCNCc1c(C)[nH]cn1"
)
ACIDIC_RECIPE = {
    "medium": "acidic",
    "salt_state": "hydrochloride",
    "temperature_c": 100.0,
}
NEUTRAL_RECIPE = {
    "medium": "neutral",
    "salt_state": "free_base",
    "temperature_c": 25.0,
}


def test_condition_tokens_are_stable_and_order_invariant() -> None:
    left = condition_tokens_from_mapping(
        {"solvents": ["water", "ethanol"], "temperature_c": 25.0}
    )
    right = condition_tokens_from_mapping(
        {"temperature_c": 25.0, "solvents": ["ethanol", "water"]}
    )

    assert left == right
    assert "condition.temperature_c=25" in left


def test_resolved_recipe_projection_excludes_provenance_and_ids() -> None:
    tokens = condition_tokens_from_recipe(
        {
            "recipe_id": "must-not-leak",
            "temperature_c": 80.0,
            "bases": [
                {
                    "substance_id": "cas:584-08-7",
                    "primary_role": "base",
                    "source_field": "must-not-leak",
                    "provenance": {"document": "must-not-leak"},
                }
            ],
        }
    )

    assert "condition.bases.identity=cas:584-08-7" in tokens
    assert "condition.bases.role=base" in tokens
    assert "condition.temperature_c=80" in tokens
    assert all("must-not-leak" not in token for token in tokens)


def test_canonical_converted_record_adapter_uses_resolved_recipe() -> None:
    choice_set = build_reaction_choice_set_from_record(
        {
            "reaction_smiles": S_ALKYLATION,
            "reference_id": "reference-1",
            "resolved_recipe": ACIDIC_RECIPE,
        },
        label_strength=0.8,
    )

    assert choice_set.reference_id == "reference-1"
    assert choice_set.label_strength == 0.8
    assert "condition.temperature_c=100" in choice_set.condition_tokens


def test_cysteamine_choice_set_enumerates_s_and_n_outcomes() -> None:
    choice_set = build_reaction_choice_set(
        S_ALKYLATION,
        ACIDIC_RECIPE,
        reference_id="patent-example",
        label_strength=1.0,
    )

    assert choice_set.selected_candidate.element == "S"
    assert {candidate.element for candidate in choice_set.candidates} == {"N", "S"}
    assert len({candidate.product_smiles for candidate in choice_set.candidates}) == 2
    assert all(candidate.product_smiles for candidate in choice_set.candidates)
    assert choice_set.schema_version == "1.0"


def test_cysteamine_competition_warning_is_review_only_and_explained() -> None:
    warning = detect_functional_group_competition(S_ALKYLATION)

    assert warning is not None
    assert warning.code == "POSSIBLE_FUNCTIONAL_GROUP_COMPETITION"
    assert warning.selected_outcome.element == "S"
    assert {outcome.element for outcome in warning.competing_outcomes} == {"N"}
    assert warning.ranking_impact == "none"
    assert warning.conditions_evaluated is False
    assert "conditions and relative outcome probabilities" in warning.message


def test_single_endpoint_reaction_has_no_competition_warning() -> None:
    assert detect_functional_group_competition("CCBr.N>>CCN") is None


def test_same_enumerator_handles_n_o_competition() -> None:
    choice_set = build_reaction_choice_set(
        "CCBr.NCCO>>CCOCCN",
        {"solvent": "water"},
        reference_id="aminoethanol-example",
    )

    assert choice_set.selected_candidate.element == "O"
    assert {candidate.element for candidate in choice_set.candidates} == {"N", "O"}


def _training_rows() -> tuple:
    rows = []
    for index in range(4):
        rows.append(
            build_reaction_choice_set(
                S_ALKYLATION,
                ACIDIC_RECIPE,
                reference_id=f"acidic-{index}",
                label_strength=1.0,
            )
        )
        rows.append(
            build_reaction_choice_set(
                N_ALKYLATION,
                NEUTRAL_RECIPE,
                reference_id=f"neutral-{index}",
                label_strength=1.0,
            )
        )
    return tuple(rows)


def test_listwise_model_learns_condition_dependent_site_preference() -> None:
    rows = _training_rows()
    model = ConditionalEditChoiceModel(feature_dimension=1024)

    report = model.fit(rows, epochs=200, learning_rate=0.15)
    acidic_query = build_reaction_choice_set(S_ALKYLATION, ACIDIC_RECIPE)
    neutral_query = build_reaction_choice_set(N_ALKYLATION, NEUTRAL_RECIPE)
    acidic_assessment = model.assess(acidic_query)
    neutral_assessment = model.assess(neutral_query)

    assert report.final_loss < report.initial_loss
    assert acidic_assessment.desired_probability > 0.95
    assert acidic_assessment.ranked_outcomes[0].element == "S"
    assert acidic_assessment.exact_condition_reference_support == 4
    assert neutral_assessment.desired_probability > 0.95
    assert neutral_assessment.ranked_outcomes[0].element == "N"
    assert neutral_assessment.exact_condition_reference_support == 4


def test_model_exposes_condition_mismatch_instead_of_a_specific_rule() -> None:
    model = ConditionalEditChoiceModel(feature_dimension=1024)
    model.fit(_training_rows(), epochs=200, learning_rate=0.15)
    unsupported_s_query = build_reaction_choice_set(
        S_ALKYLATION,
        NEUTRAL_RECIPE,
    )

    assessment = model.assess(unsupported_s_query)

    assert assessment.desired_probability < 0.05
    assert assessment.probability_margin < 0.0
    assert assessment.ranked_outcomes[0].element == "N"
    assert assessment.exact_condition_reference_support == 0
    assert assessment.evidence_status == "model_only"


def test_training_and_serialization_are_deterministic() -> None:
    rows = _training_rows()
    first = ConditionalEditChoiceModel(feature_dimension=256)
    second = ConditionalEditChoiceModel(feature_dimension=256)

    first_report = first.fit(rows, epochs=80, learning_rate=0.1)
    second_report = second.fit(reversed(rows), epochs=80, learning_rate=0.1)
    query = build_reaction_choice_set(S_ALKYLATION, ACIDIC_RECIPE)
    restored = ConditionalEditChoiceModel.from_dict(first.to_dict())

    assert first_report == second_report
    assert first.weights == second.weights
    assert restored.predict_probabilities(query) == first.predict_probabilities(query)


def test_unsupported_multi_event_reaction_is_explicit() -> None:
    with pytest.raises(ValueError, match="one formed and one broken"):
        build_reaction_choice_set(
            "CCBr.NCCO>>CCOCCN.CC",
            {"solvent": "water"},
        )
