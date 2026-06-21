from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest


def test_plan_sources_analysis_only_no_dataset_load():
    from chemtools.recommend.models import RecommendationRequest
    from chemtools.recommend.planner import plan_sources

    plan = plan_sources(
        RecommendationRequest(reaction_smiles="CCO>>CC=O", analysis_only=True)
    )
    assert plan.needs_hte_data is False
    assert plan.sources_to_run == ()
    assert "no dataset-backed recommendation" in plan.notes[0]


def test_plan_sources_precedent_only_single_pass_coerces_literature():
    from chemtools.recommend.models import OutputView, RecommendationRequest, SourceGroup
    from chemtools.recommend.planner import plan_sources

    plan = plan_sources(
        RecommendationRequest(
            reaction_smiles="A.B>>P",
            source_group=SourceGroup.RULES,
            output_view=OutputView.PRECEDENT_ONLY,
        )
    )
    assert plan.needs_hte_data is True
    assert plan.single_run_source_group == "literature"
    assert plan.needs_precedent_data is True


def test_analyze_recommendation_query_uses_featurizers(monkeypatch):
    from chemtools.recommend.models import RecommendationRequest
    from chemtools.recommend.query_analysis import analyze_recommendation_query

    fake_smiles = types.SimpleNamespace(
        normalize_reaction=lambda txt: {
            "normalized": "A.B>>P",
            "reactants": [{"smiles_norm": "A"}, {"smiles_norm": "B"}],
            "agents": [],
            "products": [{"smiles_norm": "P"}],
        }
    )
    monkeypatch.setitem(sys.modules, "chemtools.core.smiles", fake_smiles)
    monkeypatch.setattr(
        "chemtools.recommend.query_analysis.pick_electrophile_nucleophile",
        lambda reactants: (reactants[0], reactants[1]),
    )
    monkeypatch.setattr(
        "chemtools.recommend.query_analysis._featurize_query_reaction",
        lambda _reaction_smiles: {
            "reaction": {
                "reaction_key": "CRK-demo",
                "aggregates": {
                    "reacted_motifs": ["Ar-Br", "B(OH)2"],
                    "formed_motifs": ["biaryl"],
                    "spectator_motifs": ["OMe"],
                    "spectator_groups_ranked": ["ether"],
                },
                "reaction_type": {
                    "reaction_type": "Suzuki",
                    "confidence": 0.91,
                },
            }
        },
    )
    monkeypatch.setattr(
        "chemtools.recommend.query_analysis._resolve_reaction_type",
        lambda label: (
            "suzuki_miyaura" if str(label).lower() == "suzuki" else str(label),
            "suzuki_miyaura" if str(label).lower() == "suzuki" else None,
            "Suzuki-Miyaura" if str(label).lower() == "suzuki" else None,
            "cross_coupling" if str(label).lower() == "suzuki" else None,
        ),
    )

    qa = analyze_recommendation_query(
        RecommendationRequest(reaction_smiles="A.B>>P", reaction_type_filter="Suzuki")
    )
    assert qa.reaction_smiles_normalized == "A.B>>P"
    assert qa.reactant_a_smiles == "A"
    assert qa.reactant_b_smiles == "B"
    assert qa.product_smiles == "P"
    assert qa.reaction_key == "CRK-demo"
    assert qa.detected_reaction_type_id == "suzuki_miyaura"
    assert qa.detected_reaction_type_category == "cross_coupling"
    assert qa.requested_reaction_type_filter_canonical == "suzuki_miyaura"


def test_data_manager_caches_recommender_instances(monkeypatch, tmp_path: Path):
    from chemtools.recommend.data_manager import RecommendationDataManager
    import chemtools.recommend.recommender as recommender_mod

    calls = {"n": 0}

    class FakeRecommender:
        def __init__(self, db_path: str):
            calls["n"] += 1
            self.db_path = db_path

    monkeypatch.setattr(recommender_mod, "HTERecommender", FakeRecommender)
    dm = RecommendationDataManager(base_db_path=str(tmp_path))
    r1, i1 = dm.get_recommender(source_group="rules")
    r2, i2 = dm.get_recommender(source_group="rules")

    assert r1 is r2
    assert calls["n"] == 1
    assert i1.cache_hit is False
    assert i2.cache_hit is True


def test_api_recommend_analysis_only_skips_data_loading(monkeypatch):
    from chemtools.recommend.api import recommend
    from chemtools.recommend.models import RecommendationRequest

    monkeypatch.setattr(
        "chemtools.recommend.api.analyze_recommendation_query",
        lambda req: types.SimpleNamespace(reaction_smiles_normalized="A>>B"),
    )

    class _DM:
        def get_recommender(self, **kwargs):  # pragma: no cover - should never run
            raise AssertionError("get_recommender should not be called")

    out = recommend(
        RecommendationRequest(reaction_smiles="A>>B", analysis_only=True),
        data_manager=_DM(),
    )
    assert out.plan.needs_hte_data is False
    assert out.recommendation is None
    assert out.loaded_resources["skipped"] == "no_dataset_load_required"


def test_api_recommend_per_source_keeps_results_by_source_without_merge(monkeypatch):
    from chemtools.recommend.api import recommend
    from chemtools.recommend.models import OutputView, RecommendationRequest, RunStrategy
    from chemtools.recommend.data_manager import LoadedResourceInfo

    monkeypatch.setattr(
        "chemtools.recommend.api.analyze_recommendation_query",
        lambda req: types.SimpleNamespace(reactant_a_smiles="A", reactant_b_smiles="B", product_smiles="P"),
    )

    class FakeRec:
        def __init__(self, tag: str):
            self.tag = tag
            self.recommendations = [f"{tag}-combined"]
            self.recommendations_by_source = {tag: [f"{tag}-1", f"{tag}-2"]}

    class FakeRecommender:
        def recommend(self, **kwargs):
            src = kwargs.get("source_group") or "all"
            if src == "literature":
                r = FakeRec("literature")
                r.recommendations_by_source["precedent"] = ["prec-1"]
                return r
            return FakeRec(src)

    class FakeDM:
        def __init__(self):
            self.calls = []

        def get_recommender(self, **kwargs):
            src = kwargs.get("source_group") or "all"
            self.calls.append(src)
            return FakeRecommender(), LoadedResourceInfo(
                cache_key=src,
                db_path=f"/tmp/{src}",
                source_group=str(src),
                cache_hit=False,
            )

    dm = FakeDM()
    out = recommend(
        RecommendationRequest(
            reaction_smiles="A.B>>P",
            run_strategy=RunStrategy.PER_SOURCE,
            output_view=OutputView.BY_SOURCE,
        ),
        data_manager=dm,
    )

    assert out.plan.run_strategy.value == "per_source"
    assert out.recommendation is not None
    assert out.recommendation.recommendations == []
    source_map = out.recommendation.recommendations_by_source
    assert "literature" in source_map
    assert "experiments" in source_map
    assert "rules" in source_map
    assert "precedent" in source_map
    assert dm.calls[:3] == ["literature", "experiments", "rules"]


def test_api_recommend_precedent_only_trims_output(monkeypatch):
    from chemtools.recommend.api import recommend
    from chemtools.recommend.data_manager import LoadedResourceInfo
    from chemtools.recommend.models import OutputView, RecommendationRequest

    monkeypatch.setattr(
        "chemtools.recommend.api.analyze_recommendation_query",
        lambda req: types.SimpleNamespace(reactant_a_smiles="A", reactant_b_smiles="B", product_smiles="P"),
    )

    class FakeRec:
        def __init__(self):
            self.recommendations = ["full-1"]
            self.recommendations_by_source = {"literature": ["lit-1"], "precedent": ["prec-1", "prec-2"]}

    class FakeRecommender:
        def recommend(self, **kwargs):  # noqa: ARG002
            return FakeRec()

    class FakeDM:
        def get_recommender(self, **kwargs):
            return FakeRecommender(), LoadedResourceInfo(
                cache_key="k", db_path="/tmp/db", source_group="literature", cache_hit=True
            )

    out = recommend(
        RecommendationRequest(reaction_smiles="A.B>>P", output_view=OutputView.PRECEDENT_ONLY),
        data_manager=FakeDM(),
    )
    assert out.plan.single_run_source_group == "literature"
    assert out.recommendation.recommendations == ["prec-1", "prec-2"]
    assert list(out.recommendation.recommendations_by_source.keys()) == ["precedent"]


def test_plan_sources_motif_strategy_auto_loads_experiments_and_rules():
    from chemtools.recommend.models import RecommendationRequest
    from chemtools.recommend.planner import plan_sources

    plan = plan_sources(RecommendationRequest(reaction_smiles="A.B>>P", strategy="motif"))
    assert plan.recommendation_strategy == "motif"
    assert plan.run_strategy.value == "per_source"
    assert plan.sources_to_run == ("experiments", "rules")
    assert plan.single_run_source_group is None


def test_plan_sources_similarity_strategy_pins_literature():
    from chemtools.recommend.models import RecommendationRequest
    from chemtools.recommend.planner import plan_sources

    plan = plan_sources(
        RecommendationRequest(reaction_smiles="A.B>>P", strategy="similarity", source_group="rules")
    )
    assert plan.recommendation_strategy == "similarity"
    assert plan.single_run_source_group == "literature"
    assert plan.needs_precedent_data is True


def test_api_recommend_similarity_strategy_filters_to_similarity(monkeypatch):
    from chemtools.recommend.api import recommend
    from chemtools.recommend.models import RecommendationRequest

    monkeypatch.setattr(
        "chemtools.recommend.api.analyze_recommendation_query",
        lambda req: types.SimpleNamespace(reactant_a_smiles="A", reactant_b_smiles="B", product_smiles="P"),
    )

    class FakeRec:
        def __init__(self):
            self.recommendations = ["combined-lit"]
            self.recommendations_by_source = {
                "literature": ["lit-1"],
                "precedent": ["prec-1", "prec-2"],
            }

    class FakeDM:
        def __init__(self):
            self.calls = []

        def get_recommender(self, **kwargs):
            src = kwargs.get("source_group") or "all"
            self.calls.append(src)
            raise AssertionError("similarity strategy should use fast path and skip data-manager loads")

    monkeypatch.setattr(
        "chemtools.recommend.api._run_similarity_fast_pass",
        lambda req, analysis: (FakeRec(), {"similarity_fast_path": True}),
    )

    dm = FakeDM()
    out = recommend(
        RecommendationRequest(reaction_smiles="A.B>>P", strategy="similarity"),
        data_manager=dm,
    )
    assert dm.calls == []
    assert list(out.recommendation.recommendations_by_source.keys()) == ["similarity"]
    assert out.recommendation.recommendations == ["prec-1", "prec-2"]


def test_api_recommend_motif_strategy_auto_loads_and_filters_sources(monkeypatch):
    from chemtools.recommend.api import recommend
    from chemtools.recommend.data_manager import LoadedResourceInfo
    from chemtools.recommend.models import RecommendationRequest

    monkeypatch.setattr(
        "chemtools.recommend.api.analyze_recommendation_query",
        lambda req: types.SimpleNamespace(reactant_a_smiles="A", reactant_b_smiles="B", product_smiles="P"),
    )

    class FakeRec:
        def __init__(self, tag: str):
            self.recommendations = [f"{tag}-combined"]
            self.recommendations_by_source = {tag: [f"{tag}-1", f"{tag}-2"]}

    class FakeRecommender:
        def recommend(self, **kwargs):
            src = kwargs.get("source_group") or "all"
            if src == "experiments":
                return FakeRec("experiments")
            if src == "rules":
                rec = FakeRec("rules")
                rec.recommendations_by_source["precedent"] = ["should-drop"]
                return rec
            return FakeRec(str(src))

    class FakeDM:
        def __init__(self):
            self.calls = []

        def get_recommender(self, **kwargs):
            src = kwargs.get("source_group") or "all"
            self.calls.append(src)
            return FakeRecommender(), LoadedResourceInfo(
                cache_key=str(src), db_path=f"/tmp/{src}", source_group=str(src), cache_hit=False
            )

    dm = FakeDM()
    out = recommend(
        RecommendationRequest(reaction_smiles="A.B>>P", strategy="motif", top_k=3),
        data_manager=dm,
    )
    assert dm.calls == ["experiments", "rules"]
    assert set(out.recommendation.recommendations_by_source.keys()) == {"experiments", "rules"}
    assert out.recommendation.recommendations[0] == "experiments-1"
    assert "should-drop" not in out.recommendation.recommendations
