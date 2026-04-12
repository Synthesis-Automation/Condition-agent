from __future__ import annotations

from chemtools.recommend.api import _run_similarity_fast_pass
from chemtools.recommend.models import QueryAnalysis, RecommendationRequest


def test_similarity_fast_pass_uses_detected_family_when_confident(monkeypatch) -> None:
    captured = {}

    def fake_knn(
        reactant_a_smiles,
        reactant_b_smiles,
        product_smiles,
        reaction_type_filter,
        top_k,
        **kwargs,
    ):
        captured["reaction_type_filter"] = reaction_type_filter
        captured["top_k"] = top_k
        return [{"reaction_id": "prec-1", "similarity": 0.9}]

    monkeypatch.setattr("chemtools.recommend.recommender._run_precedent_knn", fake_knn)

    req = RecommendationRequest(reaction_smiles="A.B>>P", strategy="similarity", top_k=7)
    analysis = QueryAnalysis(
        reaction_smiles_input="A.B>>P",
        reactant_a_smiles="A",
        reactant_b_smiles="B",
        product_smiles="P",
        detected_reaction_type="Suzuki",
        detected_reaction_type_id="suzuki_miyaura",
        reaction_type_confidence=0.91,
    )

    rec_obj, loaded = _run_similarity_fast_pass(req, analysis)

    assert captured["reaction_type_filter"] == "suzuki_miyaura"
    assert captured["top_k"] == 7
    assert loaded["similarity_family_filter"] == "suzuki_miyaura"
    assert loaded["similarity_family_filter_source"] == "detected"
    assert len(rec_obj.recommendations) == 1


def test_similarity_fast_pass_keeps_cross_family_when_detection_is_weak(monkeypatch) -> None:
    captured = {}

    def fake_knn(
        reactant_a_smiles,
        reactant_b_smiles,
        product_smiles,
        reaction_type_filter,
        top_k,
        **kwargs,
    ):
        captured["reaction_type_filter"] = reaction_type_filter
        return []

    monkeypatch.setattr("chemtools.recommend.recommender._run_precedent_knn", fake_knn)

    req = RecommendationRequest(reaction_smiles="A.B>>P", strategy="similarity")
    analysis = QueryAnalysis(
        reaction_smiles_input="A.B>>P",
        reactant_a_smiles="A",
        reactant_b_smiles="B",
        product_smiles="P",
        detected_reaction_type="Suzuki",
        detected_reaction_type_id="suzuki_miyaura",
        reaction_type_confidence=0.2,
    )

    _rec_obj, loaded = _run_similarity_fast_pass(req, analysis)

    assert captured["reaction_type_filter"] is None
    assert loaded["similarity_family_filter"] == "cross_family"
    assert loaded["similarity_family_filter_source"] == "cross_family"
