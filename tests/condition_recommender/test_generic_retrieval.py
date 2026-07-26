from pathlib import Path

import pytest

from reactive_taxonomy import (
    REACTION_SIGNATURE_SCHEMA_VERSION,
    reaction_signature_definition_versions,
)

from condition_recommender import (
    GenericConditionRecommender,
    recommend_generic_conditions,
)
import condition_recommender.generic_api as generic_api_module
from condition_recommender.conversion.engine import convert_datasets
from condition_recommender.generic_api import recommend_indexed_signature
from condition_recommender.generic_indexing import (
    build_generic_index,
    save_generic_index,
)
from condition_recommender.generic_retrieval import (
    generic_signature_similarity,
    load_generic_retrieval_rules,
    retrieve_compatible_generic_pool,
    retrieve_compatible_generic_pool_with_trace,
    retrieve_generic_pool,
    retrieve_generic_pool_with_trace,
)
from condition_recommender.recipe_ranking import (
    load_generic_ranking_rules,
    validate_generic_ranking_rules,
)
from condition_recommender.similarity import (
    assess_signature_similarity,
    load_generic_similarity_rules,
    validate_generic_similarity_rules,
)


def _signature(
    token: str,
    *,
    exact: str = "exact-a",
    handles: str = "handles-a",
    transformation: str = "transformation-a",
    bond: str = "bond-a",
    environment: str = "environment-a",
    family: str | None = "family-a",
    family_confidence: float = 1.0,
    reaction_scope: str = "intermolecular",
    ring_size: int | None = None,
) -> dict:
    return {
        "schema_version": REACTION_SIGNATURE_SCHEMA_VERSION,
        "definition_versions": reaction_signature_definition_versions(),
        "signature_id": token,
        "exact_signature_key": exact,
        "handle_signature_key": handles,
        "transformation_signature_key": transformation,
        "bond_edit_signature_key": bond,
        "environment_signature_key": environment,
        "formed_bond_types": ["C-N:SINGLE"],
        "broken_bond_types": ["B-C:SINGLE"],
        "order_changes": [],
        "event_count": 1,
        "event_scope": "single_event",
        "events": [
            {
                "formed_bond_types": ["C-N:SINGLE"],
                "broken_bond_types": ["B-C:SINGLE"],
                "order_changes": [],
                "transformation_class": "generic_c_n_coupling",
                "topology": {"reaction_scope": reaction_scope},
            }
        ],
        "partners": [
            {
                "handle_tokens": ["B(OH)2"],
                "anchor_contexts": ["Ar"],
                "steric": {"class": "open"},
                "electronic": {"class": "neutral"},
                "nearby_groups": [],
            },
            {
                "handle_tokens": ["N-H"],
                "anchor_contexts": ["Ar"],
                "steric": {"class": "open"},
                "electronic": {"class": "neutral"},
                "nearby_groups": [],
            },
        ],
        "spectator_groups": [],
        "transformation_class": "generic_c_n_coupling",
        "named_family": family,
        "family_confidence": family_confidence,
        "topology": {
            "reaction_scope": reaction_scope,
            "formed_bond_scopes": [reaction_scope]
            if reaction_scope in {"intramolecular", "intermolecular"}
            else [],
            "formed_ring_sizes": [ring_size] if ring_size is not None else [],
            "ring_count_delta": 1 if ring_size is not None else 0,
        },
    }


def _record(index: int, signature: dict, *, tier: str = "verified") -> dict:
    recipe_id = f"RCR1:{index % 2}"
    recipe_core_id = f"RCORE1:{index % 2}"
    return {
        "schema_version": "1.9",
        "converter_definition_version": "generic_conversion.v1.6",
        "admission_tier": tier,
        "index_eligibility": "eligible" if tier == "verified" else "review_only",
        "chemistry_status": "verified",
        "condition_status": "resolved_complete",
        "outcome_status": "usable",
        "reaction_id": f"reaction-{index}",
        "observation_id": f"observation-{index}",
        "reaction_smiles": "C.N>>CN",
        "yield_pct": 60.0 + index,
        "source_dataset": f"dataset-{index % 2}",
        "reference_id": f"REF1:{index}",
        "reaction_signature": signature,
        "reference_condition_series_id": f"RCS1:{index}",
        "resolved_recipe_id": recipe_id,
        "resolved_recipe_core_id": recipe_core_id,
        "resolved_recipe": {
            "recipe_id": recipe_id,
            "recipe_core_id": recipe_core_id,
        },
        "condition_resolution": {"has_uncertainty": False},
    }


def test_generic_index_admits_only_usable_verified_records() -> None:
    tier_only = _record(4, _signature("four"))
    tier_only.pop("index_eligibility")
    index = build_generic_index(
        [
            _record(1, _signature("one")),
            _record(2, _signature("two"), tier="review"),
            {**_record(3, _signature("three")), "resolved_recipe": None},
            tier_only,
        ]
    )
    assert len(index.rows) == 1
    assert index.exact["exact-a"] == (0,)


def test_retrieval_uses_exact_signature_when_supported() -> None:
    query = _signature("query")
    records = [_record(index, _signature(str(index))) for index in range(3)]
    records.append(
        _record(9, _signature("other", exact="exact-b", handles="handles-b"))
    )
    level, pool = retrieve_generic_pool(query, build_generic_index(records))
    assert level == "exact_signature"
    assert len(pool) == 3


def test_retrieval_falls_back_exact_to_handles_then_family() -> None:
    query = _signature("query")
    records = [
        _record(
            index,
            _signature(str(index), exact=f"other-{index}"),
        )
        for index in range(3)
    ]
    level, _ = retrieve_generic_pool(query, build_generic_index(records))
    assert level == "handle_signature"

    family_records = [
        _record(
            index,
            _signature(
                str(index),
                exact=f"other-{index}",
                handles=f"other-{index}",
                transformation=f"other-{index}",
            ),
        )
        for index in range(3)
    ]
    level, pool = retrieve_generic_pool(query, build_generic_index(family_records))
    assert level == "named_family"
    assert len(pool) == 3


def test_low_confidence_family_uses_transformation_fallback() -> None:
    query = _signature("query", family_confidence=0.5)
    records = [
        _record(
            index,
            _signature(
                str(index),
                exact=f"exact-{index}",
                handles=f"handles-{index}",
                family="different-family",
            ),
        )
        for index in range(3)
    ]
    level, pool = retrieve_generic_pool(query, build_generic_index(records))
    assert level == "transformation_signature"
    assert len(pool) == 3


def test_family_fallback_never_crosses_bond_edit_gate() -> None:
    query = _signature("query")
    records = [
        _record(
            index,
            _signature(
                str(index),
                exact=f"exact-{index}",
                handles=f"handles-{index}",
                transformation=f"transformation-{index}",
                bond="incompatible-bond",
            ),
        )
        for index in range(3)
    ]
    level, pool = retrieve_generic_pool(query, build_generic_index(records))
    assert level == "no_compatible_bond_edit"
    assert pool == ()


def test_compatibility_exclusion_continues_to_relaxed_tier() -> None:
    query = _signature("query")
    query["spectator_groups"] = [
        {"group_id": "aldehyde", "tags": ["oxidation_sensitive"]}
    ]
    unsafe = _record(0, _signature("unsafe"))
    unsafe["resolved_recipe"] = {
        "recipe_id": unsafe["resolved_recipe_id"],
        "oxidants": [
            {"roles": [{"family_id": "peroxides_hydroperoxides", "role_id": "oxidant"}]}
        ],
    }
    safe_records = [
        _record(
            index,
            _signature(str(index), exact=f"safe-exact-{index}"),
        )
        for index in range(1, 4)
    ]
    level, accepted, raw_count, excluded_count = retrieve_compatible_generic_pool(
        query,
        build_generic_index([unsafe, *safe_records]),
        minimum_pool_size=2,
    )
    assert level == "handle_signature"
    assert len(accepted) == 3
    assert raw_count == 4
    assert excluded_count == 1

    traced = retrieve_compatible_generic_pool_with_trace(
        query,
        build_generic_index([unsafe, *safe_records]),
        minimum_pool_size=2,
    )
    assert traced.level == "handle_signature"
    assert traced.independent_candidate_count == 4
    assert traced.independent_compatible_candidate_count == 3
    assert traced.trace[-1].status == "selected"
    assert traced.trace[-1].excluded_candidate_count == 1


def test_support_threshold_counts_one_reference_once() -> None:
    query = _signature("query")
    records = [_record(index, _signature(str(index))) for index in range(3)]
    for record in records:
        record["reference_id"] = "REF1:one-paper"

    level, rows, trace = retrieve_generic_pool_with_trace(
        query,
        build_generic_index(records),
        minimum_pool_size=3,
    )

    assert level == "exact_signature_limited_support"
    assert len(rows) == 3
    exact = next(item for item in trace if item.level == "exact_signature")
    assert exact.candidate_count == 3
    assert exact.independent_candidate_count == 1
    assert exact.status == "selected_limited_support"


def test_environment_neighbors_narrow_the_bond_edit_pool() -> None:
    query = _signature(
        "query",
        exact="query-exact",
        handles="query-handles",
        transformation="query-transformation",
        family=None,
    )
    records = [
        _record(
            index,
            _signature(
                str(index),
                exact=f"exact-{index}",
                handles=f"handles-{index}",
                transformation=f"transformation-{index}",
                family=None,
            ),
        )
        for index in range(4)
    ]
    records[-1]["reaction_signature"]["partners"][0]["steric"]["class"] = "crowded"
    records[-1]["reaction_signature"]["partners"][0]["electronic"]["class"] = (
        "electron_poor"
    )
    records[-1]["reaction_signature"]["partners"][1]["steric"]["class"] = "crowded"
    records[-1]["reaction_signature"]["partners"][1]["electronic"]["class"] = (
        "electron_poor"
    )

    level, pool, trace = retrieve_generic_pool_with_trace(
        query,
        build_generic_index(records),
        minimum_pool_size=3,
    )

    assert level == "environment_neighbors"
    assert len(pool) == 3
    environment = next(
        item for item in trace if item.level == "environment_neighbors"
    )
    assert environment.status == "selected"
    assert environment.candidate_count == 3
    assert all(item.level != "bond_edit_signature" for item in trace)


def test_bond_edit_fallback_remains_reachable_without_environment_overlap() -> None:
    query = _signature(
        "query",
        exact="query-exact",
        handles="query-handles",
        transformation="query-transformation",
        family=None,
    )
    records = [
        _record(
            index,
            _signature(
                str(index),
                exact=f"exact-{index}",
                handles=f"handles-{index}",
                transformation=f"transformation-{index}",
                family=None,
            ),
        )
        for index in range(3)
    ]
    for record in records:
        for partner in record["reaction_signature"]["partners"]:
            partner["steric"]["class"] = "crowded"
            partner["electronic"]["class"] = "electron_poor"

    level, pool, trace = retrieve_generic_pool_with_trace(
        query,
        build_generic_index(records),
        minimum_pool_size=3,
    )

    assert level == "bond_edit_signature"
    assert len(pool) == 3
    assert next(
        item for item in trace if item.level == "environment_neighbors"
    ).status == "empty"
    assert trace[-1].level == "bond_edit_signature"
    assert trace[-1].status == "selected"


def test_similarity_does_not_treat_missing_features_as_matches() -> None:
    score, components = generic_signature_similarity(
        {
            "formed_bond_types": ["C-N:SINGLE"],
            "partners": [],
        },
        {
            "formed_bond_types": ["C-N:SINGLE"],
            "partners": [],
        },
    )
    assert components["edit_topology"] == 1.0
    assert components["handles"] == 0.0
    assert components["environment"] == 0.0
    assert 0.0 < score < 1.0

    assessment = assess_signature_similarity(
        {"formed_bond_types": ["C-N:SINGLE"], "partners": []},
        {"formed_bond_types": ["C-N:SINGLE"], "partners": []},
    )
    assert assessment.score == score
    assert round(sum(assessment.contributions.values()), 6) == score
    assert assessment.definition_version == "1.0"


def test_similarity_prefers_matching_reaction_topology() -> None:
    query = _signature("query", reaction_scope="intramolecular", ring_size=5)
    matched = _signature("matched", reaction_scope="intramolecular", ring_size=5)
    different_ring = _signature(
        "different-ring", reaction_scope="intramolecular", ring_size=7
    )
    intermolecular = _signature("intermolecular")

    matched_score, matched_components = generic_signature_similarity(query, matched)
    ring_score, ring_components = generic_signature_similarity(query, different_ring)
    intermolecular_score, intermolecular_components = generic_signature_similarity(
        query, intermolecular
    )

    assert matched_components["reaction_topology"] == 1.0
    assert 0.0 < ring_components["reaction_topology"] < 1.0
    assert intermolecular_components["reaction_topology"] == 0.0
    assert matched_score > ring_score > intermolecular_score


def test_similarity_preserves_reaction_event_multiplicity() -> None:
    query = _signature("query")
    query["events"] = [*query["events"], *query["events"]]
    query["event_count"] = 2
    query["event_scope"] = "multi_event"
    exact = {**_signature("exact"), "events": list(query["events"])}
    partial = _signature("partial")

    exact_score, exact_components = generic_signature_similarity(query, exact)
    partial_score, partial_components = generic_signature_similarity(query, partial)

    assert exact_components["reaction_events"] == 1.0
    assert partial_components["reaction_events"] == 0.5
    assert exact_score > partial_score


def test_generic_fallback_discloses_reaction_scope_mismatch() -> None:
    query = _signature(
        "query",
        exact="query-exact",
        handles="query-handles",
        transformation="query-transformation",
        family=None,
        reaction_scope="intramolecular",
        ring_size=5,
    )
    precedents = [
        _record(
            index,
            _signature(
                str(index),
                exact=f"other-exact-{index}",
                handles=f"other-handles-{index}",
                transformation=f"other-transformation-{index}",
                family=None,
                reaction_scope="intermolecular",
            ),
        )
        for index in range(3)
    ]

    result = recommend_indexed_signature(
        query,
        build_generic_index(precedents),
        minimum_pool_size=1,
    )

    assert result.valid
    assert result.retrieval_level == "environment_neighbors"
    assert result.retrieval_definition_version == "1.5"
    assert "REACTION_TOPOLOGY_FALLBACK_USED" in result.warnings
    assert any(
        caution.startswith("Reaction-scope mismatch:")
        for caution in result.recommendations[0].cautions
    )


def test_aggregation_counts_references_as_independent_evidence() -> None:
    query = _signature("query")
    records = [_record(index, _signature(str(index))) for index in range(4)]
    for record in records:
        record["resolved_recipe_id"] = "RCR1:shared"
        record["resolved_recipe_core_id"] = "RCORE1:shared"
        record["resolved_recipe"] = {
            "recipe_id": "RCR1:shared",
            "recipe_core_id": "RCORE1:shared",
        }
    for record in records[:3]:
        record["reference_id"] = "REF1:shared"
        record["reference_condition_series_id"] = "RCS1:shared"

    result = recommend_indexed_signature(
        query,
        build_generic_index(records),
        minimum_pool_size=1,
    )

    recommendation = result.recommendations[0]
    assert recommendation.observation_support == 4
    assert recommendation.reference_support == 2
    assert recommendation.precedent_reference_ids == (
        "REF1:3",
        "REF1:shared",
    )
    assert recommendation.condition_series_support == 2
    assert any(
        "one independent evidence unit" in caution
        for caution in recommendation.cautions
    )


def test_recommendation_can_report_unknown_expected_yield() -> None:
    records = [_record(index, _signature(str(index))) for index in range(3)]
    for record in records:
        record["yield_pct"] = None
        record["outcome_status"] = "missing"

    result = recommend_indexed_signature(
        _signature("query"),
        build_generic_index(records),
        minimum_pool_size=1,
    )

    assert result.valid
    recommendation = result.recommendations[0]
    assert recommendation.expected_yield_pct is None
    assert recommendation.score_trace.ranking_components["yield"] is None
    assert recommendation.score_trace.applied_ranking_weights["yield"] == 0.0
    assert round(
        sum(recommendation.score_trace.ranking_contributions.values()),
        6,
    ) == recommendation.score
    assert recommendation.score_trace.observed_outcome_count == 0
    assert any("No usable yield evidence" in item for item in recommendation.cautions)


def test_condition_uncertainty_is_an_explicit_ranking_component() -> None:
    certain = _record(0, _signature("certain"))
    uncertain = _record(1, _signature("uncertain"))
    certain["yield_pct"] = uncertain["yield_pct"] = 75.0
    uncertain["condition_resolution"] = {"has_uncertainty": True}

    result = recommend_indexed_signature(
        _signature("query"),
        build_generic_index([certain, uncertain]),
        minimum_pool_size=1,
    )

    assert result.valid
    assert len(result.recommendations) == 2
    assert result.recommendations[0].recipe_core_id == "RCORE1:0"
    assert (
        result.recommendations[0].score_trace.ranking_components[
            "condition_certainty"
        ]
        == 1.0
    )
    assert (
        result.recommendations[1].score_trace.ranking_components[
            "condition_certainty"
        ]
        == 0.0
    )


def test_generic_retrieval_weights_are_normalized() -> None:
    rules = load_generic_retrieval_rules()
    similarity = load_generic_similarity_rules()
    ranking = load_generic_ranking_rules()
    assert round(sum(similarity["weights"].values()), 10) == 1.0
    assert round(sum(ranking["weights"].values()), 10) == 1.0
    assert rules["retrieval_ladder"][-2:] == [
        "environment_neighbors",
        "bond_edit_signature",
    ]


def test_similarity_and_ranking_definitions_reject_stale_schemas() -> None:
    similarity = {**load_generic_similarity_rules(), "schema_version": "stale"}
    ranking = {**load_generic_ranking_rules(), "schema_version": "stale"}

    with pytest.raises(ValueError, match="similarity definition schema"):
        validate_generic_similarity_rules(similarity)
    with pytest.raises(ValueError, match="ranking definition schema"):
        validate_generic_ranking_rules(ranking)


def test_preloaded_recommender_loads_the_index_once(
    tmp_path: Path,
    monkeypatch,
) -> None:
    path = tmp_path / "index.json"
    save_generic_index(
        build_generic_index([_record(1, _signature("one"))]),
        path,
    )
    original_loader = generic_api_module.load_generic_index
    calls = []

    def counting_loader(source):
        calls.append(str(source))
        return original_loader(source)

    monkeypatch.setattr(generic_api_module, "load_generic_index", counting_loader)
    recommender = GenericConditionRecommender.from_path(path)
    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"

    recommender.recommend(reaction, minimum_pool_size=1)
    recommender.recommend(reaction, minimum_pool_size=1)

    assert calls == [str(path)]


def test_real_pilot_returns_resolved_recipe(tmp_path: Path) -> None:
    output = tmp_path / "generic_conversion_chan_lam_pilot"
    convert_datasets(
        "data-processor/reaction_dataset/ChanLam_Narylation.csv",
        output,
        max_rows=100,
    )
    path = output / "records.jsonl"
    result = recommend_generic_conditions(
        "c1ccc2[nH]cnc2c1.COc1ccc(B(O)O)cc1>>COc1ccc(-n2cnc3ccccc32)cc1",
        records_path=path,
        top_k=3,
        minimum_pool_size=1,
    )
    assert result.valid
    assert result.retrieval_level == "exact_signature"
    assert result.candidate_count == 1
    assert result.compatible_candidate_count == 1
    assert result.excluded_candidate_count == 0
    assert result.schema_version == "1.5"
    assert result.retrieval_trace[-1].status == "selected"
    assert result.recommendations
    assert result.recommendations[0].recipe_id.startswith("RCR1:")
    assert result.recommendations[0].resolved_recipe["recipe_id"].startswith("RCR1:")
    assert 0.0 <= result.recommendations[0].compatibility_score <= 1.0
