"""Dataset binding, dual-channel parity and condition qualification tests."""

from dataclasses import asdict
import json
import sqlite3

import pytest

from reactive_taxonomy import featurize_reaction
from condition_registry import ConditionConstraintSet
from condition_registry.constraints import normalize_condition_constraint
from condition_recommender import GenericConditionRecommender
from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.sqlite_indexing import (
    save_sqlite_generic_index,
    load_sqlite_generic_index,
)
from condition_recommender.shared_core_index import (
    build_shared_core_index,
    load_shared_core_index,
)

QUERY = "Ic1ccccc1.N#C[Cu]>>N#Cc1ccccc1"
BROMIDE = "Brc1ccccc1.N#C[Cu]>>N#Cc1ccccc1"


def record(number, reaction, *, reference=None, temperature=25):
    analysis = featurize_reaction(reaction)
    assert analysis.valid and analysis.reaction_signature and analysis.reaction_core
    return {
        "schema_version": "10.3",
        "converter_definition_version": "generic_conversion.v10.3",
        "admission_tier": "verified",
        "index_eligibility": "eligible",
        "precedent_tier": "trusted",
        "core_eligibility": "trusted_core",
        "core_eligibility_definition_version": "core_eligibility.v1@1.0",
        "chemistry_status": "verified",
        "condition_status": "resolved_complete",
        "condition_stage_status": "single_stage",
        "outcome_status": "usable",
        "reaction_id": f"reaction-{number}",
        "observation_id": f"obs-{number}",
        "reaction_smiles": reaction,
        "yield_pct": 70,
        "source_dataset": "controlled-test-fixture",
        "reference_id": reference or f"REF1:{number}",
        "reaction_signature": asdict(analysis.reaction_signature),
        "reaction_core": asdict(analysis.reaction_core),
        "fallback_descriptor": asdict(analysis.fallback_descriptor),
        "resolved_recipe_id": f"recipe-{number}",
        "resolved_recipe_core_id": f"core-{number}",
        "resolved_recipe": {
            "recipe_id": f"recipe-{number}",
            "recipe_core_id": f"core-{number}",
            "temperature_c": temperature,
        },
        "condition_resolution": {"has_uncertainty": False},
    }


def engine(tmp_path, records, *, sqlite=False):
    index = build_generic_index(records)
    if sqlite:
        source = tmp_path / "generic_index.sqlite"
        save_sqlite_generic_index(index, source)
        index = load_sqlite_generic_index(source)
    path = tmp_path / "generic_index.shared_core.sqlite"
    build_shared_core_index(index, path)
    return GenericConditionRecommender(
        index, shared_core_index=load_shared_core_index(path, index)
    )


@pytest.mark.parametrize("sqlite", [False, True])
def test_generalized_core_retrieves_other_leaving_groups_and_sources(tmp_path, sqlite):
    recommender = engine(
        tmp_path,
        [record(1, BROMIDE), record(2, "CS(=O)(=O)Oc1ccccc1.C#N>>N#Cc1ccccc1")],
        sqlite=sqlite,
    )
    result = recommender.recommend(QUERY, search_scope="broad")
    assert result.valid and len(result.recommendations) == 2
    assert result.query_reaction_smiles == QUERY
    assert all(
        item.match_label == "L1: shared transformation"
        for item in result.recommendations
    )
    assert all(
        item.evidence_relation == "analogue_evidence" for item in result.recommendations
    )
    assert all(item.source_input_requirements for item in result.recommendations)
    assert not recommender.recommend(QUERY, search_scope="same_handle").recommendations


def test_product_channel_and_retro_seeds_have_identical_qualification(tmp_path):
    recommender = engine(tmp_path, [record(1, BROMIDE)])
    normal = recommender.recommend(QUERY)
    seeded = recommender.recommend(
        QUERY, preferred_reaction_ids=("reaction-1", "missing-id")
    )
    assert normal.valid and seeded.valid
    first = next(
        x for x in normal.shared_core_trace if x.get("reaction_id") == "reaction-1"
    )
    second = next(
        x for x in seeded.shared_core_trace if x.get("reaction_id") == "reaction-1"
    )
    assert first["comparison"] == second["comparison"]
    assert first["condition_compatibility"] == second["condition_compatibility"]
    assert set(first["channels"]) == {"direct", "product_side", "reactant_side"}
    assert set(second["channels"]) == {
        "direct",
        "product_side",
        "reactant_side",
        "retro_precedent",
    }
    assert seeded.independent_compatible_candidate_count == 1
    assert seeded.candidate_count == 1
    assert any(
        x.get("reason") == "RETRO_PRECEDENT_NOT_IN_CONDITION_INDEX"
        for x in seeded.shared_core_trace
    )


def test_same_product_wrong_transformation_is_diagnosed_and_excluded(tmp_path):
    # Amide dehydration reaches the same nitrile through a different edit graph.
    recommender = engine(tmp_path, [record(1, "NC(=O)c1ccccc1>>N#Cc1ccccc1")])
    result = recommender.recommend(QUERY, preferred_reaction_ids=("reaction-1",))
    assert not result.recommendations
    assert result.shared_core_trace[0]["comparison"]["reasons"] == (
        "DIFFERENT_TRANSFORMATION_SAME_PRODUCT",
    )


def test_constraints_filter_before_support_and_seed_cannot_bypass(tmp_path):
    recommender = engine(tmp_path, [record(1, BROMIDE, temperature=150)])
    constraint = normalize_condition_constraint(
        "maximum_temperature_c", "80", provenance="explicit_user"
    ).constraint
    result = recommender.recommend(
        QUERY,
        preferred_reaction_ids=("reaction-1",),
        condition_constraints=ConditionConstraintSet((constraint,)),
    )
    assert not result.recommendations and result.excluded_candidate_count == 1
    assert not result.shared_core_trace[0]["condition_compatibility"]["compatible"]


def test_duplicate_references_and_retrieval_channels_count_once(tmp_path):
    recommender = engine(
        tmp_path,
        [
            record(1, BROMIDE, reference="REF1:same"),
            record(2, BROMIDE, reference="REF1:same"),
        ],
    )
    result = recommender.recommend(
        QUERY, preferred_reaction_ids=("reaction-1", "reaction-2")
    )
    assert (
        result.candidate_count == 2
        and result.independent_compatible_candidate_count == 1
    )


def test_whole_reaction_precedes_analogue_and_stops_on_independent_support(tmp_path):
    recommender = engine(
        tmp_path, [record(1, QUERY), record(2, QUERY), record(3, BROMIDE)]
    )
    result = recommender.recommend(QUERY, top_k=5)
    assert len(result.recommendations) == 2
    assert all(item.match_label == "Whole reaction" for item in result.recommendations)
    assert all(
        item.candidate_channels == ("direct",) for item in result.recommendations
    )
    broad = recommender.recommend(QUERY, top_k=5, search_scope="broad")
    assert [item.match_level for item in broad.recommendations] == [1, 1, 3]


def test_same_recipe_different_required_sources_remains_separate(tmp_path):
    first, second = record(1, BROMIDE), record(2, "Brc1ccccc1.C#N>>N#Cc1ccccc1")
    for key in ("resolved_recipe_id", "resolved_recipe_core_id", "resolved_recipe"):
        second[key] = first[key]
    recommender = engine(tmp_path, [first, second])
    result = recommender.recommend(QUERY, search_scope="broad")
    assert len(result.recommendations) == 2
    assert all(item.support == 1 for item in result.recommendations)


def test_same_recipe_source_context_aggregates_remote_scaffold_variants(tmp_path):
    first = record(1, "Brc1ncccc1.N#C[Cu]>>N#Cc1ncccc1")
    second = record(2, "Brc1nccnc1.N#C[Cu]>>N#Cc1nccnc1")
    for key in ("resolved_recipe_id", "resolved_recipe_core_id", "resolved_recipe"):
        second[key] = first[key]
    result = engine(tmp_path, [first, second]).recommend(QUERY, search_scope="broad")
    assert len(result.recommendations) == 1
    item = result.recommendations[0]
    assert item.support == 2
    assert set(item.precedent_reaction_ids) == {"reaction-1", "reaction-2"}
    assert len(item.source_input_requirements) == 2
    assert item.evidence_relation == "analogue_evidence"


def test_general_edit_retrieval_and_product_seeds_share_qualification(tmp_path):
    recommender = engine(tmp_path, [record(1, "CCC=O>>CCCO"), record(2, "CCCO>>CCC=O")])
    result = recommender.recommend("CC=O>>CCO", preferred_reaction_ids=("reaction-2",))
    assert result.valid
    assert result.recommendations[0].precedent_reaction_ids == ("reaction-1",)
    assert any(
        entry.get("reaction_id") == "reaction-2" and not entry["comparison"]["eligible"]
        for entry in result.shared_core_trace
    )


def test_artifact_binding_and_corruption_are_explicit_failures(tmp_path):
    recommender = engine(tmp_path, [record(1, BROMIDE)])
    wrong = build_generic_index([record(1, QUERY)])
    with pytest.raises(ValueError, match="ARTIFACT_MISMATCH"):
        recommender.shared_core_index.validate_binding(wrong)
    with sqlite3.connect(recommender.shared_core_index.path) as connection:
        connection.execute("UPDATE projection SET payload='{}'")
    with pytest.raises(ValueError, match="OBSERVATION_MISMATCH"):
        recommender.recommend(QUERY)


def test_unknown_projection_version_is_rejected(tmp_path):
    recommender = engine(tmp_path, [record(1, BROMIDE)])
    path = recommender.shared_core_index.path
    with sqlite3.connect(path) as connection:
        payload = json.loads(
            connection.execute("SELECT payload FROM metadata").fetchone()[0]
        )
        payload["definition_hash"] = "old"
        connection.execute("UPDATE metadata SET payload=?", (json.dumps(payload),))
    with pytest.raises(ValueError, match="ARTIFACT_MISMATCH"):
        load_shared_core_index(path, recommender.index)


def test_cannot_overwrite_source_index(tmp_path):
    recommender = engine(tmp_path, [record(1, BROMIDE)], sqlite=True)
    with pytest.raises(ValueError, match="overwrite"):
        build_shared_core_index(recommender.index, tmp_path / "generic_index.sqlite")


def test_atomic_failure_preserves_existing_artifact(tmp_path, monkeypatch):
    import condition_recommender.shared_core_index as module

    recommender = engine(tmp_path, [record(1, BROMIDE)])
    path = recommender.shared_core_index.path
    before = path.read_bytes()

    def fail(*args, **kwargs):
        raise RuntimeError("simulated build interruption")

    monkeypatch.setattr(module, "build_shared_reaction_core", fail)
    with pytest.raises(RuntimeError, match="interruption"):
        build_shared_core_index(recommender.index, path)
    assert path.read_bytes() == before
    assert not list(tmp_path.glob("*.tmp"))


def test_cancelled_build_preserves_existing_artifact(tmp_path):
    from condition_recommender.shared_core_index import SharedCoreBuildCancelled

    recommender = engine(tmp_path, [record(1, BROMIDE)])
    path = recommender.shared_core_index.path
    before = path.read_bytes()
    with pytest.raises(SharedCoreBuildCancelled):
        build_shared_core_index(recommender.index, path, cancel_check=lambda: True)
    assert path.read_bytes() == before
    assert not list(tmp_path.glob("*.tmp"))


def test_parallel_build_matches_serial_payloads_and_lookups(tmp_path):
    recommender = engine(tmp_path, [record(1, BROMIDE), record(2, QUERY)], sqlite=True)
    serial = recommender.shared_core_index.path
    parallel = tmp_path / "parallel.sqlite"
    build_shared_core_index(recommender.index, parallel, workers=2, batch_size=1)
    from contextlib import closing

    with (
        closing(sqlite3.connect(serial)) as left,
        closing(sqlite3.connect(parallel)) as right,
    ):
        for table, order in (
            ("projection", "position"),
            ("lookup", "kind,key,position"),
        ):
            assert (
                left.execute(f"SELECT * FROM {table} ORDER BY {order}").fetchall()
                == right.execute(f"SELECT * FROM {table} ORDER BY {order}").fetchall()
            )


def test_resume_reuses_committed_prefix_and_rejects_wrong_source(tmp_path, monkeypatch):
    from condition_recommender.shared_core_index import SharedCoreBuildCancelled
    import condition_recommender.shared_core_index as module

    recommender = engine(tmp_path, [record(1, BROMIDE), record(2, QUERY)], sqlite=True)
    target = tmp_path / "resumed.sqlite"
    progress = []
    with pytest.raises(SharedCoreBuildCancelled):
        build_shared_core_index(
            recommender.index,
            target,
            batch_size=1,
            resume=True,
            progress_callback=progress.append,
            cancel_check=lambda: bool(progress and progress[-1] >= 1),
        )
    assert not target.exists()
    checkpoint = target.with_name(target.name + ".building")
    with pytest.raises(ValueError, match="manifest"):
        load_shared_core_index(checkpoint, recommender.index)
    wrong = build_generic_index([record(3, BROMIDE)])
    with pytest.raises(ValueError, match="CHECKPOINT_MISMATCH"):
        build_shared_core_index(wrong, target, resume=True)
    projected = []
    original = module._project_rows

    def tracked(start, rows):
        projected.append(start)
        return original(start, rows)

    monkeypatch.setattr(module, "_project_rows", tracked)
    build_shared_core_index(recommender.index, target, batch_size=1, resume=True)
    assert projected == [1]
    assert not checkpoint.exists()
    assert load_shared_core_index(target, recommender.index).row_count == 2


def test_candidate_budget_is_bounded_and_reported(tmp_path, monkeypatch):
    import condition_recommender.shared_core_retrieval as module

    recommender = engine(tmp_path, [record(1, BROMIDE), record(2, BROMIDE)])
    original = module.load_shared_retrieval_rules()
    monkeypatch.setattr(
        module,
        "load_shared_retrieval_rules",
        lambda: {**original, "candidate_limit": 1},
    )
    result = recommender.recommend(QUERY, search_scope="broad")
    assert result.candidate_count == 1
    assert "SHARED_CORE_CANDIDATE_BUDGET_REACHED" in result.warnings


def test_reactant_channel_alone_recovers_qualified_precedents(tmp_path, monkeypatch):
    from condition_recommender.shared_core_index import SharedCoreIndex

    recommender = engine(tmp_path, [record(1, BROMIDE)])
    original = SharedCoreIndex.lookup

    def only_reactants(self, kind, key, limit):
        return (
            original(self, kind, key, limit)
            if kind.startswith("reactant_")
            else ((), False)
        )

    monkeypatch.setattr(SharedCoreIndex, "lookup", only_reactants)
    result = recommender.recommend(QUERY)
    assert result.valid and result.candidate_count == 1
    assert result.recommendations[0].candidate_channels == ("reactant_side",)
    assert result.recommendations[0].evidence_relation == "analogue_evidence"


def test_same_reactants_wrong_product_cannot_supply_conditions(tmp_path):
    query = "[CH3:1][CH2:2][Br:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]"
    elimination = "[CH3:1][CH2:2][Br:3].[NH3:4]>>[CH2:1]=[CH2:2]"
    result = engine(tmp_path, [record(1, elimination)]).recommend(query)
    assert not result.recommendations
    entry = next(
        x for x in result.shared_core_trace if x.get("reaction_id") == "reaction-1"
    )
    assert "reactant_side" in entry["channels"]
    assert not entry["comparison"]["eligible"]


def test_reactant_seed_cannot_bypass_condition_constraints(tmp_path, monkeypatch):
    from condition_recommender.shared_core_index import SharedCoreIndex

    recommender = engine(tmp_path, [record(1, BROMIDE, temperature=150)])
    original = SharedCoreIndex.lookup
    monkeypatch.setattr(
        SharedCoreIndex,
        "lookup",
        lambda self, kind, key, limit: original(self, kind, key, limit)
        if kind.startswith("reactant_")
        else ((), False),
    )
    constraint = normalize_condition_constraint(
        "maximum_temperature_c", "80", provenance="explicit_user"
    ).constraint
    result = recommender.recommend(
        QUERY, condition_constraints=ConditionConstraintSet((constraint,))
    )
    assert not result.recommendations and result.excluded_candidate_count == 1
    assert not result.shared_core_trace[0]["condition_compatibility"]["compatible"]


def test_auxiliary_channels_share_reserved_budget(tmp_path, monkeypatch):
    import condition_recommender.shared_core_retrieval as module
    from condition_recommender.shared_core_index import SharedCoreIndex

    recommender = engine(tmp_path, [record(i, QUERY) for i in range(6)])
    rules = module.load_shared_retrieval_rules()
    monkeypatch.setattr(
        module,
        "load_shared_retrieval_rules",
        lambda: {
            **rules,
            "candidate_limit": 4,
            "direct_candidate_limit": 2,
        },
    )

    def lookup(self, kind, key, limit):
        # Existing direct hits must not consume the product channel's turns.
        values = {
            "whole_reaction": (0, 1),
            "product_identity": (0, 1, 2, 4),
            "reactant_identity": (3, 5),
        }
        return values.get(kind, ()), False

    monkeypatch.setattr(SharedCoreIndex, "lookup", lookup)
    result = recommender.recommend(QUERY, search_scope="broad", top_k=6)
    entries = {
        x["reaction_id"]: x for x in result.shared_core_trace if "reaction_id" in x
    }
    assert set(entries) == {"reaction-0", "reaction-1", "reaction-2", "reaction-3"}
    assert entries["reaction-2"]["channels"] == ["product_side"]
    assert entries["reaction-3"]["channels"] == ["reactant_side"]
    assert "SHARED_CORE_CANDIDATE_BUDGET_REACHED" in result.warnings


def test_artifact_without_reactant_projection_contract_requires_rebuild(tmp_path):
    recommender = engine(tmp_path, [record(1, BROMIDE)])
    with sqlite3.connect(recommender.shared_core_index.path) as connection:
        metadata = json.loads(
            connection.execute("SELECT payload FROM metadata").fetchone()[0]
        )
        metadata.pop("reactant_projection_hash")
        connection.execute("UPDATE metadata SET payload=?", (json.dumps(metadata),))
    with pytest.raises(ValueError, match="artifact manifest"):
        load_shared_core_index(recommender.shared_core_index.path, recommender.index)


def test_web_runtime_defaults_to_shared_core_and_response_contract(
    tmp_path, monkeypatch
):
    from fastapi.testclient import TestClient
    from app.web_api.main import create_app
    from app.web_api.runtime import LocalRecommendationRuntime

    monkeypatch.delenv("CONDITION_SHARED_CORE_EXPERIMENTAL", raising=False)
    engine(tmp_path, [record(1, BROMIDE)], sqlite=True)
    runtime = LocalRecommendationRuntime(index_path=tmp_path / "generic_index.sqlite")
    assert runtime.capabilities()["recommendation_engine"] == "shared_reaction_core.v2"
    with TestClient(create_app(runtime=runtime, recommendation_only=False)) as client:
        response = client.post(
            "/api/v1/recommendations",
            json={"reaction_smiles": QUERY, "library_mode": "full"},
        )
    assert response.status_code == 200
    result = response.json()["data"]
    assert result["recommendation_mode"] == "experimental_shared_core"
    assert result["recommendations"][0]["match_namespace"] == "shared_reaction_core.v2"
    assert result["shared_core_trace"]


def test_default_loader_requires_bound_companion_and_baseline_is_explicit(tmp_path):
    engine(tmp_path, [record(1, BROMIDE)], sqlite=True)
    source = tmp_path / "generic_index.sqlite"
    default = GenericConditionRecommender.from_path(source)
    assert (
        default.recommend(QUERY).recommendations[0].match_namespace
        == "shared_reaction_core.v2"
    )
    baseline = GenericConditionRecommender.from_path(source, use_shared_core=False)
    assert baseline.shared_core_index is None
    with pytest.raises(ValueError, match="requires use_shared_core"):
        GenericConditionRecommender.from_path(
            source,
            use_shared_core=False,
            shared_core_path=source.with_suffix(".shared_core.sqlite"),
        )
    source.with_suffix(".shared_core.sqlite").unlink()
    with pytest.raises(FileNotFoundError, match="Shared-core artifact is unavailable"):
        GenericConditionRecommender.from_path(source)
    assert (
        GenericConditionRecommender.from_path(
            source, use_shared_core=False
        ).shared_core_index
        is None
    )


def test_web_runtime_baseline_override_and_constructor_precedence(
    tmp_path, monkeypatch
):
    from app.web_api.runtime import LocalRecommendationRuntime

    engine(tmp_path, [record(1, BROMIDE)], sqlite=True)
    source = tmp_path / "generic_index.sqlite"
    monkeypatch.setenv("CONDITION_SHARED_CORE_EXPERIMENTAL", "0")
    runtime = LocalRecommendationRuntime(index_path=source)
    assert runtime.capabilities()["recommendation_engine"] == "baseline"
    assert (
        runtime._get_recommender(
            library_mode="full", use_rxnmapper=False, include_review=False
        ).shared_core_index
        is None
    )
    explicit = LocalRecommendationRuntime(index_path=source, shared_core_enabled=True)
    assert (
        explicit._get_recommender(
            library_mode="full", use_rxnmapper=False, include_review=False
        ).shared_core_index
        is not None
    )


def test_default_review_loader_uses_review_bound_projection(tmp_path, monkeypatch):
    from app.web_api.runtime import LocalRecommendationRuntime

    monkeypatch.delenv("CONDITION_SHARED_CORE_EXPERIMENTAL", raising=False)
    engine(tmp_path, [record(1, BROMIDE)], sqlite=True)
    source = tmp_path / "generic_index.sqlite"
    review_source = tmp_path / "generic_review_index.sqlite"
    review_index = build_generic_index(
        [record(1, BROMIDE), record(2, QUERY)], include_review=True
    )
    save_sqlite_generic_index(review_index, review_source)
    review_index = load_sqlite_generic_index(review_source)
    companion = review_source.with_suffix(".shared_core.sqlite")
    build_shared_core_index(review_index, companion)
    runtime = LocalRecommendationRuntime(index_path=source)
    recommender = runtime._get_recommender(
        library_mode="full", use_rxnmapper=False, include_review=True
    )
    assert len(recommender.index.rows) == 2
    assert recommender.shared_core_index.path == companion.resolve()
    assert (
        GenericConditionRecommender.from_path(
            source, include_review=True
        ).shared_core_index.path
        == companion.resolve()
    )


def test_default_loader_rejects_stale_projection(tmp_path):
    engine(tmp_path, [record(1, BROMIDE)], sqlite=True)
    source = tmp_path / "generic_index.sqlite"
    save_sqlite_generic_index(build_generic_index([record(2, QUERY)]), source)
    with pytest.raises(ValueError, match="SHARED_CORE_ARTIFACT_MISMATCH"):
        GenericConditionRecommender.from_path(source)


def test_shared_review_builds_projection_artifact_from_training_rows_only(tmp_path):
    from condition_recommender.chemist_review import generate_chemist_review_packet

    engine(
        tmp_path,
        [record(i, BROMIDE if i % 2 else QUERY) for i in range(1, 9)],
        sqlite=True,
    )
    report = generate_chemist_review_packet(
        tmp_path / "generic_index.sqlite",
        tmp_path / "review",
        experimental_shared_core=True,
        max_cases=3,
    )
    assert (
        report["shared_core_manifest"]["row_count"]
        == report["split"]["train_row_count"]
    )
    assert report["split"]["leakage_group_count"] == 0
    assert report["shared_core_artifact_sha256"]
