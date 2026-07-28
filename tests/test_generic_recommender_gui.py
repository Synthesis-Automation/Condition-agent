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
        precedent_reaction_ids=("reaction-1", "reaction-2"),
        precedent_reaction_smiles=(
            "BrC.B(O)O>>CC",
            "ClC.B(O)O>>CC",
        ),
        precedent_reference_ids=("REF1:one",),
        explanation=("Exact bond-edit and handle match",),
        score_trace=trace,
    )


def _result() -> GenericRecommendationResult:
    return GenericRecommendationResult(
        query_reaction_smiles="BrC.B(O)O>>CC",
        valid=True,
        query_signature_id="RS2:query",
        reaction_label="Ar1–Br + Ar2–B(OH)2 → Ar1–Ar2",
        reaction_label_status="exact_product",
        named_family="suzuki_miyaura",
        transformation_class="c_c_transfer_coupling",
        spectator_groups=(
            {
                "group_id": "nitrile",
                "chemist_label": "R–C≡N",
                "component_index": 0,
                "graph_distance": 3,
            },
        ),
        reaction_partners=(
            {
                "role": "electrophile",
                "chemist_label": "Ar–Br",
                "anchor_contexts": ("Ar",),
                "steric": {
                    "class": "ortho_hindered",
                    "ortho_substituent_count": 2,
                    "local_heavy_atoms_r2": 6,
                },
                "electronic": {
                    "class": "electron_poor",
                    "qualitative_sum": 0.75,
                },
            },
        ),
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
    assert window.summary_box.height() == 168
    assert window.data_row_layout.indexOf(window.data_label) == 0
    assert window.data_row_layout.indexOf(window.data_path_edit) == 1
    assert window.query_summary_layout.count() == 2
    assert window.reaction_image_label.objectName() == "queryReactionGraph"
    assert window.selected_details_layout.count() == 2
    assert (
        window.selected_reaction_image_label.objectName()
        == "selectedPrecedentReactionGraph"
    )
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

    summary = window.summary_box.toPlainText()
    assert "Ar1–Br + Ar2–B(OH)2 → Ar1–Ar2" in summary
    assert "Label evidence: Exact Product" in summary
    assert "Suzuki Miyaura" in summary
    assert "R–C≡N [nitrile] (reactant 1, d=3)" in summary
    assert "steric ortho hindered, 2 ortho substituent(s)" in summary
    assert "electronic electron poor (local score +0.75)" in summary
    reaction_pixmap = window.reaction_image_label.pixmap()
    assert reaction_pixmap is not None
    assert not reaction_pixmap.isNull()
    precedent_pixmap = window.selected_reaction_image_label.pixmap()
    assert precedent_pixmap is not None
    assert not precedent_pixmap.isNull()
    assert (
        window.selected_reaction_image_label.toolTip()
        == "BrC.B(O)O>>CC"
    )
    assert window.results_table.rowCount() == 1
    assert "Palladium catalyst" in window.results_table.item(0, 7).text()
    assert "Potassium carbonate" in window.details_box.toPlainText()
    assert "Exact bond-edit and handle match" in (
        window.details_box.toPlainText()
    )
    assert window.status_label.text() == "Done — 1 recipe(s)"


def test_window_requests_vector_reaction_drawing(qtbot, monkeypatch) -> None:
    calls = []

    def render(reaction_smiles, *, size, image_format):
        calls.append((reaction_smiles, size, image_format))
        return (
            b'<svg xmlns="http://www.w3.org/2000/svg" '
            b'width="560" height="142"></svg>'
        )

    monkeypatch.setattr(gui, "render_reaction_image_bytes", render)
    window = gui.GenericRecommenderWindow()
    qtbot.addWidget(window)

    window._render_reaction_graph("CCO>>CC=O")

    assert calls == [
        ("CCO>>CC=O", gui.QUERY_REACTION_IMAGE_SIZE, "svg")
    ]
    assert window.reaction_image_label._source_svg


def test_reaction_label_trims_empty_svg_canvas_space(qtbot) -> None:
    label = gui.ReactionImageLabel()
    qtbot.addWidget(label)
    image = gui.QtGui.QImage(
        120,
        80,
        gui.QtGui.QImage.Format.Format_RGB32,
    )
    image.fill(gui.QtCore.Qt.GlobalColor.white)
    painter = gui.QtGui.QPainter(image)
    painter.fillRect(
        gui.QtCore.QRect(30, 20, 60, 40),
        gui.QtCore.Qt.GlobalColor.black,
    )
    painter.end()

    trimmed = label._trim_white_space(image, margin=3)

    assert trimmed.width() < image.width()
    assert trimmed.height() < image.height()
    assert trimmed.width() >= 60
    assert trimmed.height() >= 40


def test_main_window_starts_maximized() -> None:
    events = []

    class FakeWindow:
        def showMaximized(self):
            events.append(("show", "maximized"))

    gui._show_main_window(FakeWindow())

    assert ("show", "maximized") in events


def test_recipe_summary_keeps_roles_and_operating_conditions() -> None:
    summary = gui.format_recipe_summary(_recommendation().resolved_recipe)

    assert "Catalyst: Palladium catalyst" in summary
    assert "Base: Potassium carbonate" in summary
    assert "Solvent: 1,4-Dioxane" in summary
    assert "80 °C" in summary
    assert "12 h" in summary
