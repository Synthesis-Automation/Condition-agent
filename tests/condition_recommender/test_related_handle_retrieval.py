"""Graph and evidence regressions for controlled Br/I analogue retrieval."""

from dataclasses import asdict

import pytest

from reactive_taxonomy import featurize_reaction
from condition_registry.constraints import ConditionConstraintSet, normalize_condition_constraint
from condition_recommender import ChemistRankingPreferences, GenericConditionRecommender
from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.sqlite_indexing import load_sqlite_generic_index, save_sqlite_generic_index
from condition_recommender.related_handles import load_related_handle_rules
from condition_recommender.ranking_preferences import available_ranking_profiles


QUERY = "Ic1ccccc1.CN>>CNc1ccccc1"
BROMIDE = "Brc1ccccc1.CN>>CNc1ccccc1"


def record(number: int, smiles: str, *, reference: str = "", temperature: int = 25) -> dict:
    """Build a current, structurally analysed precedent with independent provenance."""
    analysis = featurize_reaction(smiles)
    assert analysis.valid and analysis.reaction_signature and analysis.reaction_core
    return {
        "schema_version": "10.3",
        "converter_definition_version": "generic_conversion.v10.3",
        "admission_tier": "verified", "index_eligibility": "eligible",
        "precedent_tier": "trusted", "core_eligibility": "trusted_core",
        "core_eligibility_definition_version": "core_eligibility.v1@1.0",
        "chemistry_status": "verified", "condition_status": "resolved_complete",
        "condition_stage_status": "single_stage", "outcome_status": "usable",
        "reaction_id": f"reaction-{number}", "observation_id": f"obs-{number}",
        "reaction_smiles": smiles, "yield_pct": 50 + number,
        "source_dataset": "test", "reference_id": reference or f"REF1:{number}",
        "reference_condition_series_id": f"series-{number}",
        "reaction_signature": asdict(analysis.reaction_signature),
        "reaction_core": asdict(analysis.reaction_core),
        "fallback_descriptor": asdict(analysis.fallback_descriptor),
        "resolved_recipe_id": f"recipe-{number}",
        "resolved_recipe_core_id": f"core-{number}",
        "resolved_recipe": {
            "recipe_id": f"recipe-{number}", "recipe_core_id": f"core-{number}",
            "temperature_c": temperature,
        },
        "condition_resolution": {"has_uncertainty": False},
    }


@pytest.mark.parametrize("storage", ["memory", "sqlite"])
def test_automatic_finds_bromide_without_rebuilding_index(storage, tmp_path):
    index = build_generic_index([record(1, BROMIDE)])
    if storage == "sqlite":
        path = tmp_path / "index.sqlite"
        save_sqlite_generic_index(index, path)
        index = load_sqlite_generic_index(path)
    recommender = GenericConditionRecommender(index=index)
    strict = recommender.recommend(QUERY, search_scope="same_handle")
    assert not strict.recommendations
    result = recommender.recommend(QUERY)
    assert result.valid and result.query_reaction_smiles == QUERY
    assert result.search_scope == "automatic"
    assert "related_handle_retrieval.v1@1.0" in result.retrieval_definition_version
    item = result.recommendations[0]
    assert item.match_level == 3 and item.match_label == "Related handle"
    assert item.precedent_reaction_ids == ("reaction-1",)
    assert item.match_details == ("Query Ar-I; precedent Ar-Br",)
    assert any("not a validated transfer" in value for value in item.cautions)
    assert any("below the target" in value for value in item.explanation)


def test_support_uses_independent_references_even_after_recipe_target_reached():
    index = build_generic_index([
        record(1, QUERY, reference="REF1:repeated"),
        record(2, QUERY, reference="REF1:repeated"), record(3, BROMIDE),
    ])
    result = GenericConditionRecommender(index=index).recommend(QUERY, top_k=2)
    exact = result.retrieval_trace[0]
    assert exact.independent_compatible_candidate_count == 1
    assert exact.status == "selected_limited_support"
    assert any(trace.level == "related_handle_I_to_Br" for trace in result.retrieval_trace)
    assert all(item.match_level == 1 for item in result.recommendations)


def test_automatic_stops_at_sufficient_support_broad_checks_analogues():
    index = build_generic_index([record(1, QUERY), record(2, QUERY), record(3, BROMIDE)])
    recommender = GenericConditionRecommender(index=index)
    automatic = recommender.recommend(QUERY, top_k=3)
    broad = recommender.recommend(QUERY, top_k=3, search_scope="broad")
    assert not any(trace.level.startswith("related_handle_") for trace in automatic.retrieval_trace)
    assert [item.match_level for item in broad.recommendations] == [1, 1, 3]
    assert "requested" in " ".join(broad.recommendations[-1].explanation)


def test_original_condition_constraints_filter_analogues_before_support():
    index = build_generic_index([record(1, BROMIDE, temperature=150)])
    constraint = normalize_condition_constraint(
        "maximum_temperature_c", "80", provenance="explicit_user",
    ).constraint
    assert constraint is not None
    result = GenericConditionRecommender(index=index).recommend(
        QUERY, condition_constraints=ConditionConstraintSet((constraint,)),
    )
    assert not result.recommendations


@pytest.mark.parametrize("smiles", [
    "Clc1ccccc1.CN>>CNc1ccccc1",
    "Brc1ccccc1.CNC>>CN(C)c1ccccc1",
    "Brc1ccccc1.CO>>COc1ccccc1",
    "CCBr.CN>>CCNC",
])
def test_unapproved_handles_and_different_partner_chemistry_are_excluded(smiles):
    result = GenericConditionRecommender(index=build_generic_index([record(1, smiles)])).recommend(QUERY)
    assert not result.recommendations


def test_partner_order_invariance_and_reverse_direction():
    index = build_generic_index([record(1, QUERY)])
    result = GenericConditionRecommender(index=index).recommend("CN.Brc1ccccc1>>CNc1ccccc1")
    assert result.valid
    assert result.recommendations[0].match_details == ("Query Ar-Br; precedent Ar-I",)


def test_rank_profile_cannot_enable_excluded_analogue():
    recommender = GenericConditionRecommender(index=build_generic_index([record(1, BROMIDE)]))
    for profile in available_ranking_profiles():
        result = recommender.recommend(
            QUERY, search_scope="same_handle",
            ranking_preferences=ChemistRankingPreferences(profile_id=profile["profile_id"]),
        )
        assert not result.recommendations
    assert load_related_handle_rules()["calibration_status"] == "chemist_prior_pending_validation"


def test_unknown_search_scope_rejected():
    with pytest.raises(ValueError, match="search scope"):
        GenericConditionRecommender(index=build_generic_index([])).recommend(QUERY, search_scope="anything")


def test_supported_analogue_tier_does_not_expand_to_fill_display_limit():
    result = GenericConditionRecommender(index=build_generic_index([
        record(1, BROMIDE), record(2, BROMIDE),
    ])).recommend(QUERY, top_k=5)
    assert len(result.recommendations) == 2
    assert all(item.match_level == 3 for item in result.recommendations)
    assert result.retrieval_trace[-1].level == "related_handle_I_to_Br"
