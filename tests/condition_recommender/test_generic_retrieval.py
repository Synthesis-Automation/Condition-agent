from dataclasses import asdict, replace
import json
from pathlib import Path

import pytest

from reactive_taxonomy import (
    REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
    REACTION_CORE_PROJECTION_SCHEMA_VERSION,
    REACTION_SIGNATURE_SCHEMA_VERSION,
    featurize_reaction,
    reaction_signature_definition_versions,
)

from condition_recommender import (
    ChemistRankingPreferences,
    GenericConditionRecommender,
    recommend_generic_conditions,
)
import condition_recommender.generic_api as generic_api_module
import condition_recommender.generic_retrieval as generic_retrieval_module
from condition_recommender.conversion.engine import convert_datasets
from condition_recommender.generic_api import recommend_indexed_signature
from condition_recommender.generic_indexing import (
    build_generic_index,
    load_generic_index,
)
from condition_recommender.sqlite_indexing import save_sqlite_generic_index
from condition_recommender.core_retrieval import (
    retrieve_core_pool_with_trace,
)
from condition_recommender.hypothesis_retrieval import (
    load_edit_hypothesis_retrieval_rules,
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
from condition_recommender.ranking_preferences import RANKING_COMPONENTS
from condition_recommender.similarity import (
    assess_signature_similarity,
    load_generic_similarity_rules,
    validate_generic_similarity_rules,
)


def _aromatic_profile(
    *,
    access: str = "open",
    electronic: str = "balanced",
) -> dict:
    return {
        "schema_version": "1.0",
        "context_kind": "aromatic",
        "context": {
            "context_kind": "aromatic",
            "ring_family": "benzene",
            "ring_sizes": [6],
            "ortho_occupancy_count": 0,
            "ortho_capacity": 2,
            "ortho_burden_class": "none",
            "heteroatoms": [],
        },
        "steric": {
            "accessibility_class": access,
            "accessibility_score": 0.0 if access == "open" else 0.9,
        },
        "electronic": {
            "activation_axis": "electronic_demand",
            "activation_class": electronic,
            "activation_score": 0.0 if electronic == "balanced" else 0.8,
        },
        "reactive_center": {},
        "modifiers": [],
    }


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
                "role": "transfer_partner",
                "handle_tokens": ["B(OH)2"],
                "anchor_contexts": ["Ar"],
                "reactivity_profile": _aromatic_profile(),
                "nearby_groups": [],
            },
            {
                "role": "nucleophile",
                "handle_tokens": ["N-H"],
                "anchor_contexts": ["Ar"],
                "reactivity_profile": _aromatic_profile(),
                "nearby_groups": [],
            },
        ],
        "spectator_groups": [],
        "transformation_class": "generic_c_n_coupling",
        "named_family": family,
        "family_confidence": family_confidence,
        "topology": {
            "schema_version": "2.0",
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
        "schema_version": "10.1",
        "converter_definition_version": "generic_conversion.v10.1",
        "admission_tier": tier,
        "index_eligibility": "eligible" if tier == "verified" else "review_only",
        "precedent_tier": "trusted" if tier == "verified" else None,
        "core_eligibility": "trusted_core",
        "core_eligibility_definition_version": "core_eligibility.v1@1.0",
        "chemistry_status": "verified",
        "condition_status": "resolved_complete",
        "condition_stage_status": "single_stage",
        "outcome_status": "usable",
        "reaction_id": f"reaction-{index}",
        "observation_id": f"observation-{index}",
        "reaction_smiles": "C.N>>CN",
        "reaction_label": {
            "concise": f"Precedent reaction {index}",
            "detailed": f"Precedent reaction {index}",
            "status": "observed_edits",
        },
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


def _core(
    token: str,
    *,
    shape: str = "RSH2:shared",
    center: str = "RCS2:shared",
    evidence_status: str = "verified",
) -> dict:
    return {
        "schema_version": REACTION_CORE_PROJECTION_SCHEMA_VERSION,
        "algorithm_version": REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
        "core_id": f"RCP2:{token}",
        "exact_core_key": f"RCX2:{token}",
        "typed_core_key": f"RCT2:{token}",
        "shape_core_key": shape,
        "center_transition_key": center,
        "mapping_equivalence_key": f"RME1:{token}",
        "event_count": 1,
        "evidence_status": evidence_status,
        "warnings": [],
    }


def _fischer_signature(token: str) -> dict:
    signature = _signature(
        token,
        exact=f"exact-{token}",
        handles=f"handles-{token}",
        transformation=f"transformation-{token}",
        bond=f"bond-{token}",
        environment=f"environment-{token}",
        family=None,
        family_confidence=0.0,
        reaction_scope="intermolecular",
        ring_size=5,
    )
    signature.update(
        {
            "formed_bond_types": ["C-C:AROMATIC", "C-N:AROMATIC"],
            "broken_bond_types": ["C-O:DOUBLE", "N-N:SINGLE"],
            "order_changes": [
                "C-C:SINGLE>AROMATIC",
                "C-N:SINGLE>AROMATIC",
            ],
            "transformation_class": "generic_multi_event_graph_transformation",
        }
    )
    return signature


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


def test_additive_protocol_output_accepts_v10_index_rows() -> None:
    record = _record(1, _signature("compatible-v10"))
    record["schema_version"] = "10.0"
    record["converter_definition_version"] = "generic_conversion.v10.0"

    index = build_generic_index([record])

    assert len(index.rows) == 1
    assert index.record_schema_versions == ("10.0",)


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


def test_related_edit_graph_can_cross_exact_bond_key_gate() -> None:
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
    assert level == "edit_graph_neighbors"
    assert len(pool) == 3


def test_retrieval_uses_context_core_after_exact_edit_tiers() -> None:
    query = _signature(
        "query",
        exact="query-exact",
        handles="query-handles",
        transformation="query-transformation",
        bond="query-bond",
        environment="query-environment",
        family=None,
    )
    records = []
    for index in range(2):
        record = _record(
            index,
            _signature(
                str(index),
                exact=f"exact-{index}",
                handles=f"handles-{index}",
                transformation=f"transformation-{index}",
                bond=f"bond-{index}",
                environment=f"environment-{index}",
                family=None,
            ),
        )
        record["reaction_core"] = _core(str(index))
        records.append(record)

    level, pool = retrieve_generic_pool(
        query,
        build_generic_index(records),
        reaction_core=_core("query", evidence_status="external"),
    )

    assert level == "reaction_core_context"
    assert len(pool) == 2


def test_exact_retrieval_does_not_compute_broader_fallbacks(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    signature = _signature("query")
    index = build_generic_index(
        [_record(1, signature), _record(2, signature)]
    )

    def unexpected_fallback(*args: object, **kwargs: object) -> set[int]:
        raise AssertionError("broader fallback was computed before it was needed")

    monkeypatch.setattr(
        generic_retrieval_module,
        "_environment_neighbor_positions",
        unexpected_fallback,
    )
    monkeypatch.setattr(
        generic_retrieval_module,
        "_core_level_positions",
        unexpected_fallback,
    )
    monkeypatch.setattr(
        generic_retrieval_module,
        "_edit_graph_neighbor_positions",
        unexpected_fallback,
    )

    result = retrieve_compatible_generic_pool_with_trace(signature, index)

    assert result.level == "exact_signature"
    assert len(result.pool) == 2
    assert tuple(trace.level for trace in result.trace) == ("exact_signature",)


def test_approximate_tiers_bulk_load_candidate_rows() -> None:
    records = [_record(1, _signature("one")), _record(2, _signature("two"))]
    records[0]["reaction_core"] = _core("one")
    records[1]["reaction_core"] = _core("two")
    index = build_generic_index(records)

    class BulkOnlyRows:
        def __init__(self, rows: object) -> None:
            self.rows = rows
            self.select_calls = 0

        def __len__(self) -> int:
            return len(self.rows)  # type: ignore[arg-type]

        def __getitem__(self, position: int) -> object:
            raise AssertionError(f"candidate row {position} was loaded individually")

        def select(self, positions: object) -> tuple[object, ...]:
            self.select_calls += 1
            return tuple(
                self.rows[position]  # type: ignore[index]
                for position in positions  # type: ignore[union-attr]
            )

        def edit_graph_candidate_positions(self, prototype: object) -> tuple[int, ...]:
            return tuple(range(len(self)))

    bulk_rows = BulkOnlyRows(index.rows)
    bulk_index = replace(index, rows=bulk_rows)  # type: ignore[arg-type]

    edit_positions = generic_retrieval_module._edit_graph_neighbor_positions(
        _signature("query", bond="different-bond"),
        bulk_index,
        exclude=set(),
    )
    assert edit_positions == {0, 1}
    assert bulk_rows.select_calls == 1

    core_positions = generic_retrieval_module._core_level_positions(
        _core("one"),
        bulk_index,
        level="reaction_core_exact",
        query_eligible=True,
    )
    assert core_positions == {0}
    assert bulk_rows.select_calls == 2


def test_core_ladder_prefers_exact_then_local_before_context() -> None:
    exact_records = []
    for index in range(2):
        record = _record(index, _signature(str(index)))
        record["reaction_core"] = _core("shared-exact")
        exact_records.append(record)
    exact = retrieve_core_pool_with_trace(
        _core("shared-exact", evidence_status="external"),
        _signature("compatibility"),
        build_generic_index(exact_records),
    )

    local_records = []
    for index in range(2):
        record = _record(index, _signature(str(index)))
        core = _core(f"candidate-{index}")
        core["typed_core_key"] = "RCT2:shared-local"
        record["reaction_core"] = core
        local_records.append(record)
    local_query = _core("query", evidence_status="external")
    local_query["typed_core_key"] = "RCT2:shared-local"
    local = retrieve_core_pool_with_trace(
        local_query,
        _signature("compatibility"),
        build_generic_index(local_records),
    )

    assert exact.level == "reaction_core_exact"
    assert [value.status for value in exact.trace] == ["selected"]
    assert local.level == "reaction_core_local"
    assert [value.status for value in local.trace] == ["empty", "selected"]


def test_blocked_core_quality_cannot_retrieve() -> None:
    record = _record(1, _signature("one"))
    record["reaction_core"] = _core("shared")
    query = _core("shared", evidence_status="external")
    query["quality"] = {
        "status": "blocked",
        "blocking_reasons": ["formed_bond_state_inconsistent"],
    }

    result = retrieve_core_pool_with_trace(
        query,
        _signature("compatibility"),
        build_generic_index([record]),
    )

    assert result.level == "reaction_core_quality_blocked"
    assert not result.pool


def test_mapping_variants_do_not_inflate_core_support() -> None:
    records = []
    for index in range(2):
        record = _record(index, _signature(str(index)))
        record["reference_id"] = ""
        record["canonical_reaction_id"] = f"CRX1:map-variant-{index}"
        core = _core("shared")
        core["mapping_equivalence_key"] = "RME1:same-chemistry"
        record["reaction_core"] = core
        records.append(record)
    query = _core("shared", evidence_status="external")
    query["mapping_equivalence_key"] = "RME1:same-chemistry"

    result = retrieve_core_pool_with_trace(
        query,
        _signature("compatibility"),
        build_generic_index(records),
        minimum_pool_size=2,
    )

    assert result.level == "reaction_core_exact_limited_support"
    assert result.candidate_count == 2
    assert result.independent_candidate_count == 1
    assert result.independent_compatible_candidate_count == 1


def test_core_center_key_alone_cannot_retrieve_precedents() -> None:
    records = []
    for index in range(2):
        record = _record(index, _signature(str(index)))
        record["reaction_core"] = _core(
            str(index),
            shape=f"RSH2:candidate-{index}",
            center="RCS2:coarse-shared",
        )
        records.append(record)
    index = build_generic_index(records)

    result = retrieve_core_pool_with_trace(
        _core(
            "query",
            shape="RSH2:query-only",
            center="RCS2:coarse-shared",
            evidence_status="external",
        ),
        _signature("compatibility"),
        index,
    )

    assert result.level == "no_reaction_core_precedent"
    assert not result.pool
    assert index.core_centers["RCS2:coarse-shared"] == (0, 1)


def test_reaction_core_with_unresolved_remote_continuity_uses_review_tier() -> None:
    record = _record(1, _signature("one"))
    record["reaction_core"] = _core("one")
    query_core = _core("query", evidence_status="external")
    query_core["warnings"] = [
        "REACTION_CORE_REMOTE_CONTINUITY_UNRESOLVED"
    ]

    result = retrieve_core_pool_with_trace(
        query_core,
        _signature("compatibility"),
        build_generic_index([record]),
    )

    assert result.level == "reaction_core_context_limited_support"
    assert result.pool


def test_unsigned_mapped_core_query_returns_review_qualified_result() -> None:
    reaction = (
        "[Br:19][c:5]1[cH:4][cH:1][c:8]([Br:20])[cH:7][cH:6]1."
        "O[B:21](O)[c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1"
        ">>[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][cH:7]"
        "[c:8](-[c:9]3[cH:10][cH:11][cH:12][cH:13][cH:14]3)"
        "[cH:15][cH:16]2)[cH:17][cH:18]1"
    )
    analysis = featurize_reaction(reaction)
    assert analysis.reaction_signature is None
    assert analysis.reaction_core is not None
    assert analysis.partial_product_transformation is None
    core = asdict(analysis.reaction_core)
    fallback = (
        asdict(analysis.fallback_descriptor)
        if analysis.fallback_descriptor is not None
        else {}
    )
    records = []
    for index in range(2):
        record = _record(index, _signature(str(index)))
        record["reaction_core"] = core
        record["fallback_descriptor"] = fallback
        records.append(record)

    result = GenericConditionRecommender(
        build_generic_index(records)
    ).recommend(reaction)

    assert result.valid
    assert result.recommendation_mode == "reaction_core_review"
    assert result.retrieval_strategy == "reaction_core_ladder"
    assert result.retrieval_level == "reaction_core_exact"
    assert result.query_signature_id is None
    assert result.query_reaction_core_id == core["core_id"]
    assert "QUERY_TRANSFORMATION_NOT_VERIFIED" in result.warnings
    assert "REACTION_CORE_RETRIEVAL_USED" in result.warnings


def test_compatibility_exclusion_continues_to_relaxed_tier() -> None:
    query = _signature("query")
    query["spectator_groups"] = [
        {"group_id": "aldehyde", "tags": ["oxidation_sensitive"]}
    ]
    unsafe = _record(0, _signature("unsafe"))
    unsafe["resolved_recipe"] = {
        "recipe_id": unsafe["resolved_recipe_id"],
        "oxidants": [
            {
                "identity_status": "resolved",
                "role_status": "assigned",
                "primary_role": "oxidant",
                "roles": [{"role_id": "oxidant"}],
            }
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
    for partner in records[-1]["reaction_signature"]["partners"]:
        partner["reactivity_profile"] = _aromatic_profile(
            access="severe",
            electronic="electron_poor",
        )

    level, pool, trace = retrieve_generic_pool_with_trace(
        query,
        build_generic_index(records),
        minimum_pool_size=3,
    )

    assert level == "environment_neighbors"
    assert len(pool) == 3
    environment = next(item for item in trace if item.level == "environment_neighbors")
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
            partner["reactivity_profile"] = _aromatic_profile(
                access="severe",
                electronic="electron_poor",
            )

    level, pool, trace = retrieve_generic_pool_with_trace(
        query,
        build_generic_index(records),
        minimum_pool_size=3,
    )

    assert level == "bond_edit_signature"
    assert len(pool) == 3
    assert (
        next(item for item in trace if item.level == "environment_neighbors").status
        == "empty"
    )
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
    assert result.retrieval_definition_version == "1.8"
    assert "REACTION_TOPOLOGY_FALLBACK_USED" in result.warnings
    assert any(
        caution.startswith("Reaction-scope mismatch:")
        for caution in result.recommendations[0].cautions
    )


def test_recommendation_result_preserves_structured_query_context() -> None:
    query = _signature("query")
    query["spectator_groups"] = [
        {
            "group_id": "nitrile",
            "chemist_label": "R–C≡N",
            "component_index": 0,
            "graph_distance": 3,
        }
    ]

    precedent_signature = _signature("precedent")
    precedent_signature["spectator_groups"] = [
        {
            "group_id": "ether",
            "chemist_label": "R–O–R",
            "component_index": 1,
            "graph_distance": 2,
        }
    ]
    result = recommend_indexed_signature(
        query,
        build_generic_index([_record(1, precedent_signature)]),
        reaction_label={
            "concise": "Ar1–Br + Ar2–B(OH)2 → Ar1–Ar2",
            "detailed": "Ar1–Br + Ar2–B(OH)2 → Ar1–Ar2",
            "status": "observed_edits",
        },
        minimum_pool_size=1,
    )

    assert result.valid
    assert result.reaction_label["concise"] == "Ar1–Br + Ar2–B(OH)2 → Ar1–Ar2"
    assert result.reaction_label["status"] == "observed_edits"
    assert result.spectator_groups[0]["group_id"] == "nitrile"
    assert (
        result.reaction_partners[0]["reactivity_profile"]["steric"][
            "accessibility_class"
        ]
        == "open"
    )
    assert (
        result.reaction_partners[0]["reactivity_profile"]["electronic"][
            "activation_class"
        ]
        == "balanced"
    )
    hit_context = result.recommendations[0].precedent_reaction_contexts[0]
    assert hit_context["reaction_label"]["concise"] == "Precedent reaction 1"
    assert hit_context["reaction_label"]["status"] == "observed_edits"
    assert hit_context["spectator_groups"][0]["group_id"] == "ether"
    assert (
        hit_context["reaction_partners"][0]["reactivity_profile"]["steric"][
            "accessibility_class"
        ]
        == "open"
    )


def test_inferred_correspondence_precedent_discloses_review_caution() -> None:
    query = _signature("query")
    inferred = _signature("inferred")
    inferred["evidence_quality"] = "global_atom_correspondence"
    record = _record(1, inferred)
    record["chemistry_status"] = "review"

    result = recommend_indexed_signature(
        query,
        build_generic_index([record]),
        minimum_pool_size=1,
    )

    assert result.valid
    assert any(
        "inferred atom correspondence" in caution
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
    assert (
        round(
            sum(recommendation.score_trace.ranking_contributions.values()),
            6,
        )
        == recommendation.score
    )
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
        result.recommendations[0].score_trace.ranking_components["condition_certainty"]
        == 1.0
    )
    assert (
        result.recommendations[1].score_trace.ranking_components["condition_certainty"]
        == 0.0
    )


def test_chemist_reactant_category_priority_reranks_with_audit_trace() -> None:
    query_reaction = "Brc1ccccc1.CCNCC>>CCN(CC)c1ccccc1"
    aryl_amine_reaction = (
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    )
    query_analysis = featurize_reaction(query_reaction)
    aryl_analysis = featurize_reaction(aryl_amine_reaction)
    assert query_analysis.reaction_signature is not None
    assert query_analysis.reaction_core is not None
    assert aryl_analysis.reaction_signature is not None
    assert aryl_analysis.reaction_core is not None

    def precedent(
        index: int,
        analysis,
        *,
        recipe_core_id: str,
        yield_pct: float,
    ) -> dict:
        record = _record(index, asdict(analysis.reaction_signature))
        recipe_id = recipe_core_id.replace("RCORE", "RCR")
        record.update(
            {
                "reaction_smiles": (
                    query_reaction
                    if analysis is query_analysis
                    else aryl_amine_reaction
                ),
                "reaction_core": asdict(analysis.reaction_core),
                "yield_pct": yield_pct,
                "resolved_recipe_id": recipe_id,
                "resolved_recipe_core_id": recipe_core_id,
                "resolved_recipe": {
                    "recipe_id": recipe_id,
                    "recipe_core_id": recipe_core_id,
                },
            }
        )
        return record

    records = [
        precedent(
            0,
            query_analysis,
            recipe_core_id="RCORE:aliphatic",
            yield_pct=10.0,
        ),
        *[
            precedent(
                index,
                aryl_analysis,
                recipe_core_id="RCORE:aryl",
                yield_pct=99.0,
            )
            for index in range(1, 4)
        ],
    ]
    result = recommend_indexed_signature(
        asdict(query_analysis.reaction_signature),
        build_generic_index(records),
        reaction_core=asdict(query_analysis.reaction_core),
        minimum_pool_size=2,
        ranking_preferences=ChemistRankingPreferences(
            profile_id="reactant_category"
        ),
    )

    assert result.valid
    assert result.retrieval_level == "progressive_reaction_facets"
    assert result.ranking_preferences["profile_id"] == "reactant_category"
    recommendation = result.recommendations[0]
    assert recommendation.recipe_core_id == "RCORE:aliphatic"
    assert recommendation.retrieval_level == "reaction_facet_exact"
    assert result.recommendations[1].retrieval_level == "transformation_signature"
    assert (
        recommendation.score_trace.ranking_components["partner_category"]
        == 1.0
    )
    evidence = recommendation.factor_evidence["partner_category"]
    assert "secondary aliphatic amine" in evidence["query_categories"]
    assert recommendation.score_trace.default_ranking_contributions


def test_unavailable_exclusive_custom_factor_has_zero_contribution() -> None:
    weights = {name: 0.0 for name in RANKING_COMPONENTS}
    weights["functional_group_tolerance"] = 1.0
    result = recommend_indexed_signature(
        _signature("query"),
        build_generic_index([_record(0, _signature("precedent"))]),
        minimum_pool_size=1,
        ranking_preferences=ChemistRankingPreferences(
            profile_id="tolerance-only",
            weights=weights,
            customized=True,
        ),
    )

    assert result.valid
    recommendation = result.recommendations[0]
    assert recommendation.score == 0.0
    assert not any(recommendation.score_trace.applied_ranking_weights.values())
    assert any(
        "No unchanged query functional groups" in caution
        for caution in recommendation.cautions
    )


def test_unassigned_multistage_ingredient_set_is_ranked_with_caution() -> None:
    record = _record(0, _signature("unassigned-stage"))
    record["condition_stage_status"] = "unassigned_multistage"

    result = recommend_indexed_signature(
        _signature("query"),
        build_generic_index([record]),
        minimum_pool_size=1,
    )

    assert result.valid
    recommendation = result.recommendations[0]
    assert recommendation.score_trace.ranking_components["condition_certainty"] == 0.0
    assert any(
        "assignment to ordered reaction stages is unavailable" in caution
        for caution in recommendation.cautions
    )


def test_generic_retrieval_weights_are_normalized() -> None:
    rules = load_generic_retrieval_rules()
    similarity = load_generic_similarity_rules()
    ranking = load_generic_ranking_rules()
    assert round(sum(similarity["weights"].values()), 10) == 1.0
    assert round(sum(ranking["weights"].values()), 10) == 1.0
    assert rules["retrieval_ladder"][-6:] == [
        "environment_neighbors",
        "bond_edit_signature",
        "reaction_core_exact",
        "reaction_core_local",
        "reaction_core_context",
        "edit_graph_neighbors",
    ]
    hypothesis_rules = load_edit_hypothesis_retrieval_rules()
    assert round(
        sum(hypothesis_rules["similarity_weights"].values()),
        10,
    ) == 1.0
    assert (
        hypothesis_rules["consensus_policy"]
        == "intersection_across_all_hypotheses"
    )


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
    path = tmp_path / "index.sqlite"
    save_sqlite_generic_index(
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


def test_review_mode_uses_current_conversion_report_when_artifact_report_is_stale(
    tmp_path: Path,
) -> None:
    path = tmp_path / "generic_index.sqlite"
    save_sqlite_generic_index(
        build_generic_index([_record(1, _signature("one"))]),
        path,
    )
    (tmp_path / "recommendation_artifacts_report.json").write_text(
        json.dumps(
            {
                "trusted_precedent_count": 20_000,
                "unrestricted_precedent_count": 20_000,
                "review_core_precedent_count": 0,
                "review_index_reuses_trusted": True,
                "artifacts": {"fast_index": {"path": str(path)}},
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "conversion_report.json").write_text(
        json.dumps(
            {
                "output_row_count": 1,
                "failed_shard_count": 0,
                "precedent_tier_counts": {"trusted": 1},
                "integrity": {"valid": True, "verified_row_count": 1},
            }
        ),
        encoding="utf-8",
    )

    recommender = GenericConditionRecommender.from_path(path, include_review=True)

    assert len(recommender.index.rows) == 1
    assert recommender.review_index_reuses_trusted


def test_reaction_core_index_round_trip_preserves_lookup_maps(
    tmp_path: Path,
) -> None:
    record = _record(1, _signature("one"))
    record["reaction_core"] = _core("one")
    path = tmp_path / "index.sqlite"

    save_sqlite_generic_index(build_generic_index([record]), path)
    restored = load_generic_index(path)

    assert restored.rows[0].reaction_core["core_id"] == "RCP2:one"
    assert restored.core_shapes["RSH2:shared"] == (0,)
    assert restored.reaction_core_schema_version == (
        REACTION_CORE_PROJECTION_SCHEMA_VERSION
    )
    assert restored.reaction_core_algorithm_version == (
        REACTION_CORE_PROJECTION_ALGORITHM_VERSION
    )


def test_ambiguous_query_retrieves_only_edit_hypothesis_consensus() -> None:
    index = build_generic_index(
        [
            _record(201, _fischer_signature("fischer-a")),
            _record(202, _fischer_signature("fischer-b")),
            _record(203, _signature("unrelated")),
        ]
    )
    recommender = GenericConditionRecommender(index)

    result = recommender.recommend(
        "O=C1CCCCC1.Cl.NNc1ccc(F)cc1"
        ">>Fc1ccc2[nH]c3c(c2c1)CCCC3"
    )

    assert result.valid
    assert result.recommendation_mode == "ambiguous_edit_hypotheses"
    assert result.retrieval_strategy == "edit_hypothesis_consensus"
    assert result.retrieval_level == "edit_hypothesis_consensus"
    assert len(result.query_edit_hypothesis_ids) == 2
    assert result.compatible_candidate_count == 2
    assert result.independent_compatible_candidate_count == 2
    assert "AMBIGUOUS_EDIT_HYPOTHESIS_RETRIEVAL_USED" in result.warnings
    assert result.recommendations
    assert all(
        any(
            caution.startswith("Query atom correspondence is ambiguous")
            for caution in recommendation.cautions
        )
        for recommendation in result.recommendations
    )


def test_ambiguous_query_abstains_without_independent_consensus_support() -> None:
    index = build_generic_index(
        [_record(201, _fischer_signature("fischer-a"))]
    )
    recommender = GenericConditionRecommender(index)

    result = recommender.recommend(
        "O=C1CCCCC1.Cl.NNc1ccc(F)cc1"
        ">>Fc1ccc2[nH]c3c(c2c1)CCCC3"
    )

    assert not result.valid
    assert result.recommendation_mode == "ambiguous_edit_hypotheses"
    assert (
        result.retrieval_level
        == "insufficient_edit_hypothesis_consensus_support"
    )
    assert result.recommendations == ()
    assert result.error == "NO_ROBUST_EDIT_HYPOTHESIS_PRECEDENT"


def test_real_pilot_returns_resolved_recipe(tmp_path: Path) -> None:
    query_reaction = (
        "c1ccc2[nH]cnc2c1.COc1ccc(B(O)O)cc1>>"
        "COc1ccc(-n2cnc3ccccc32)cc1"
    )
    output = tmp_path / "generic_conversion_chan_lam_pilot"
    convert_datasets(
        "raw_dataset/literature_reaction_dataset/ChanLam_Narylation.csv",
        output,
        max_rows=100,
    )
    path = output / "records.jsonl"
    result = recommend_generic_conditions(
        query_reaction,
        records_path=path,
        top_k=3,
        minimum_pool_size=1,
    )
    assert result.valid
    assert result.retrieval_level == "reaction_facet_exact"
    assert result.candidate_count >= 1
    assert result.compatible_candidate_count >= 1
    assert result.compatible_candidate_count <= result.candidate_count
    assert result.excluded_candidate_count == 0
    assert result.schema_version == "3.4"
    assert result.retrieval_trace[-1].status == "selected_target_reached"
    assert result.recommendations
    assert result.recommendations[0].recipe_id.startswith("RCR2:")
    assert result.recommendations[0].resolved_recipe["recipe_id"].startswith("RCR2:")
    assert result.recommendations[0].precedent_reaction_smiles
    protocol = result.recommendations[0].synthesis_protocol
    assert protocol["reaction_smiles"] == query_reaction
    assert sum(
        material["category"] == "reaction_input"
        for material in protocol["materials"]
    ) == 2
    assert any(
        material["category"] == "condition" and material["cas"]
        for material in protocol["materials"]
    )
    assert 0.0 <= result.recommendations[0].compatibility_score <= 1.0
