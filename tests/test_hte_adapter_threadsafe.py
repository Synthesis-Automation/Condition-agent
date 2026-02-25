from __future__ import annotations

import sys
import threading
import time
import types


def test_get_default_recommender_is_threadsafe(monkeypatch):
    from chemtools.recommend import hte_adapter

    calls = {"n": 0}

    class FakeRecommender:
        def __init__(self, path: str):
            time.sleep(0.05)  # widen the race window
            calls["n"] += 1
            self.path = path

    fake_module = types.SimpleNamespace(HTERecommender=FakeRecommender)
    monkeypatch.setitem(sys.modules, "chemtools.recommend.recommender", fake_module)
    monkeypatch.setattr(hte_adapter, "_DEFAULT_RECOMMENDER", None)
    monkeypatch.setattr(hte_adapter, "_DEFAULT_DB_PATH", None)

    barrier = threading.Barrier(2)
    results = []
    errors = []

    def worker():
        try:
            barrier.wait(timeout=2.0)
            results.append(hte_adapter._get_default_recommender("data/HTE_db"))
        except Exception as exc:  # pragma: no cover - diagnostic path
            errors.append(exc)

    t1 = threading.Thread(target=worker)
    t2 = threading.Thread(target=worker)
    t1.start()
    t2.start()
    t1.join(timeout=3.0)
    t2.join(timeout=3.0)

    assert not errors
    assert len(results) == 2
    assert calls["n"] == 1
    assert results[0] is results[1]


def test_recommend_from_reaction_emits_timing_breakdown(monkeypatch):
    from chemtools.recommend import hte_adapter

    class _FakeResult:
        recommendations = []
        recommendations_by_source = {}
        predicted_reaction_type = "suzuki_miyaura"
        reaction_type_confidence = 0.8
        reactant_a_type = None
        reactant_b_type = None
        reactant_a_category = None
        reactant_b_category = None
        total_matching_experiments = 0
        database_coverage = 0.0
        is_fallback_match = False
        matched_motifs = ()
        reacted_motifs = ()
        formed_motifs = ()
        spectator_motifs = ()

    class _FakeRecommender:
        def recommend(self, **kwargs):  # noqa: ARG002
            return _FakeResult()

    monkeypatch.setattr(hte_adapter, "_extract_reaction_parts", lambda r: (["BrC1=CC=CC=C1", "B(O)C1=CC=CC=C1"], ["prod"], "rxn_norm"))
    monkeypatch.setattr(hte_adapter, "_select_reactants", lambda reactants, a, b: (reactants[0], reactants[1]))
    monkeypatch.setattr(hte_adapter, "_get_default_recommender", lambda p=None: _FakeRecommender())
    monkeypatch.setattr(hte_adapter, "_build_recommendation_entries", lambda result, rt: [])  # noqa: ARG005

    # Success path uses 10 perf_counter calls.
    ticks = iter([0.0, 0.0, 0.01, 0.01, 0.03, 0.03, 0.08, 0.08, 0.09, 0.10])
    monkeypatch.setattr(hte_adapter.time, "perf_counter", lambda: next(ticks))

    payload = hte_adapter.recommend_from_reaction("A>>B")
    timing = payload["meta"]["timing_ms"]

    assert payload["meta"]["status"] == "success"
    assert timing["input_parse_ms"] == 10.0
    assert timing["recommender_get_ms"] == 20.0
    assert timing["recommend_compute_ms"] == 50.0
    assert timing["postprocess_ms"] == 10.0
    assert timing["total_ms"] == 100.0
    assert payload["extras"]["hte"]["timing_ms"] == timing
