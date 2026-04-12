from __future__ import annotations

from types import SimpleNamespace

from app import HTE_recommender_gui as gui


def _fake_result(
    *,
    recommendations_by_source=None,
    recommendations=None,
    predicted_reaction_type=None,
    reaction_type_confidence=0.0,
):
    return SimpleNamespace(
        recommendations_by_source=recommendations_by_source or {},
        recommendations=recommendations or [],
        predicted_reaction_type=predicted_reaction_type,
        reaction_type_confidence=reaction_type_confidence,
    )


def test_full_mode_reuses_precedent_for_similarity_when_kmn_missing(monkeypatch) -> None:
    worker = gui.RecommendationWorker(
        db_path="data/HTE_db",
        reaction_smiles="A.B>>P",
        top_k=10,
        min_exp=1,
        reaction_filter="",
        catalyst_filter="",
        strategy="Full (all 4 modes)",
        source_override="",
        use_aryl_weighting=False,
        prefer_mixfp_similarity=False,
    )

    calls = []

    def fake_run_single(recommender, reactant_a, reactant_b, product, source_group):
        calls.append(source_group)
        if source_group is None:
            return _fake_result(
                recommendations_by_source={"precedent": ["prec-1"]},
                recommendations=["base-1"],
                predicted_reaction_type="Suzuki_miyaura",
                reaction_type_confidence=0.92,
            )
        if source_group == "literature":
            return _fake_result(
                recommendations_by_source={"literature": ["lit-1"], "precedent": ["prec-1", "prec-2"]},
                recommendations=["lit-1"],
            )
        if source_group == "motif":
            return _fake_result(recommendations_by_source={"motif": ["motif-1"]}, recommendations=["motif-1"])
        if source_group == "rules":
            return _fake_result(recommendations_by_source={"rules": ["rules-1"]}, recommendations=["rules-1"])
        raise AssertionError(f"unexpected source_group: {source_group}")

    def fake_run_similarity_only(*args, **kwargs):
        raise AssertionError("dedicated similarity pass should be skipped when KMN is missing")

    monkeypatch.setattr(worker, "_run_single", fake_run_single)
    monkeypatch.setattr(worker, "_run_similarity_only", fake_run_similarity_only)
    monkeypatch.setattr(gui, "_kmn_index_built", lambda: False)

    result = worker._run_all_recommendations(object(), "A", "B", "P")

    assert calls == [None, "literature", "motif", "rules"]
    assert result.recommendations_by_source["similarity"] == ["prec-1", "prec-2"]
    assert result.recommendations_by_source["precedent"] == ["prec-1", "prec-2"]


def test_similarity_hint_from_result_requires_confident_prediction() -> None:
    low_conf = _fake_result(predicted_reaction_type="Suzuki_miyaura", reaction_type_confidence=0.2)
    high_conf = _fake_result(predicted_reaction_type="Suzuki_miyaura", reaction_type_confidence=0.8)

    assert gui._similarity_hint_from_result(low_conf) is None
    assert gui._similarity_hint_from_result(high_conf) == "Suzuki_miyaura"
