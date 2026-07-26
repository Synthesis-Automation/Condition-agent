from pathlib import Path

from app import generic_recommender_gui as gui
from condition_recommender.models import (
    GenericConditionRecommendation,
    GenericRecommendationResult,
    RecommendationScoreTrace,
)


def _recommendation() -> GenericConditionRecommendation:
    trace = RecommendationScoreTrace(
        similarity_components={"bond_edits": 1.0},
        similarity_contributions={"bond_edits": 0.3},
        ranking_components={"similarity": 0.8},
        ranking_contributions={"similarity": 0.4},
        applied_ranking_weights={"similarity": 0.5},
        independent_evidence_count=2,
        observed_outcome_count=2,
        pool_yield_prior_pct=75.0,
        definition_versions={"generic_ranking.v1": "1.0"},
    )
    return GenericConditionRecommendation(
        rank=1,
        recipe_id="RCR1:variant",
        recipe_core_id="RCORE1:core",
        recipe_variant_ids=("RCR1:variant",),
        resolved_recipe={
            "recipe_id": "RCR1:variant",
            "recipe_core_id": "RCORE1:core",
            "catalysts": [{"canonical_name": "Palladium catalyst"}],
            "bases": [{"canonical_name": "Potassium carbonate"}],
            "solvents": [{"canonical_name": "1,4-Dioxane"}],
            "temperature_c": 80.0,
            "time_h": 12.0,
        },
        score=0.82,
        similarity_score=0.91,
        compatibility_score=1.0,
        expected_yield_pct=78.5,
        support=3,
        observation_support=4,
        reference_support=2,
        condition_series_support=2,
        dataset_support=2,
        retrieval_level="exact_signature",
        precedent_reaction_ids=("reaction-1",),
        precedent_reference_ids=("REF1:one",),
        explanation=("Exact bond-edit and handle match",),
        score_trace=trace,
    )


def _result() -> GenericRecommendationResult:
    return GenericRecommendationResult(
        query_reaction_smiles="BrC.B(O)O>>CC",
        valid=True,
        query_signature_id="RS1:query",
        named_family="suzuki_miyaura",
        transformation_class="c_c_transfer_coupling",
        retrieval_level="exact_signature",
        candidate_count=4,
        independent_candidate_count=3,
        compatible_candidate_count=3,
        independent_compatible_candidate_count=2,
        excluded_candidate_count=1,
        recommendations=(_recommendation(),),
    )


def test_window_uses_literature_recommendation_data_by_default(
    qtbot,
) -> None:
    window = gui.GenericRecommenderWindow()
    qtbot.addWidget(window)

    assert Path(window.data_path_edit.text()) == (
        gui.default_recommendation_data_path()
    )
    assert window.top_k_spin.value() == 5
    assert window.reaction_edit.metaObject().className() == "QLineEdit"
    assert window.results_table.columnCount() == 9
    assert not window.export_button.isEnabled()


def test_worker_reuses_recommender_contract(monkeypatch) -> None:
    expected = _result()
    calls = []

    class FakeRecommender:
        def recommend(self, reaction, *, top_k, minimum_pool_size):
            calls.append((reaction, top_k, minimum_pool_size))
            return expected

    monkeypatch.setattr(
        gui,
        "_get_cached_recommender",
        lambda path: FakeRecommender(),
    )
    worker = gui.GenericRecommendationWorker(
        "index.json.gz",
        "A.B>>P",
        top_k=3,
        minimum_pool_size=2,
    )
    progress = []
    finished = []
    worker.progress.connect(progress.append)
    worker.finished.connect(
        lambda success, result, error: finished.append(
            (success, result, error)
        )
    )

    worker.run()

    assert calls == [("A.B>>P", 3, 2)]
    assert len(progress) == 2
    assert finished == [(True, expected, "")]


def test_window_renders_recipe_summary_and_details(qtbot) -> None:
    window = gui.GenericRecommenderWindow()
    qtbot.addWidget(window)

    window._render_result(_result())

    assert "Suzuki Miyaura" in window.summary_box.toPlainText()
    assert window.results_table.rowCount() == 1
    assert "Palladium catalyst" in window.results_table.item(0, 7).text()
    assert "Potassium carbonate" in window.details_box.toPlainText()
    assert "Exact bond-edit and handle match" in (
        window.details_box.toPlainText()
    )
    assert window.status_label.text() == "Done — 1 recipe(s)"


def test_recipe_summary_keeps_roles_and_operating_conditions() -> None:
    summary = gui.format_recipe_summary(_recommendation().resolved_recipe)

    assert "Catalyst: Palladium catalyst" in summary
    assert "Base: Potassium carbonate" in summary
    assert "Solvent: 1,4-Dioxane" in summary
    assert "80 °C" in summary
    assert "12 h" in summary
