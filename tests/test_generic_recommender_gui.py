from dataclasses import replace
from pathlib import Path

from app import reaction_recommender_gui as gui
from condition_recommender.models import (
    GenericConditionRecommendation,
    GenericRecommendationResult,
    RecommendationScoreTrace,
)
from condition_recommender.discovery_models import (
    DiscoveryScoreTrace,
    ReactionDiscoveryHit,
    ReactionDiscoveryResult,
)


def _aromatic_profile(
    *,
    access: str,
    electronic: str,
    ortho_count: int,
    burden: str,
    steric_score: float,
    electronic_score: float,
) -> dict:
    return {
        "schema_version": "1.0",
        "context_kind": "aromatic",
        "context": {
            "context_kind": "aromatic",
            "ring_family": "benzene",
            "ring_sizes": (6,),
            "ortho_occupancy_count": ortho_count,
            "ortho_capacity": 2,
            "ortho_burden_class": burden,
            "heteroatoms": (),
        },
        "steric": {
            "accessibility_class": access,
            "accessibility_score": steric_score,
        },
        "electronic": {
            "activation_axis": "electronic_demand",
            "activation_class": electronic,
            "activation_score": electronic_score,
        },
        "reactive_center": {},
        "modifiers": (),
    }


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
        precedent_reaction_contexts=(
            {
                "reaction_id": "reaction-1",
                "reaction_smiles": "BrC.B(O)O>>CC",
                "reaction_label": {
                    "text": "Precedent Ar–Br coupling",
                    "status": "structural_equation",
                    "basis": "reaction_sites",
                },
                "spectator_groups": (
                    {
                        "group_id": "ether",
                        "chemist_label": "R–O–R",
                        "component_index": 1,
                        "graph_distance": 2,
                    },
                ),
                "reaction_partners": (
                    {
                        "role": "transfer_partner",
                        "chemist_label": "Ar–B(OH)2",
                        "anchor_contexts": ("Ar",),
                        "reactivity_profile": _aromatic_profile(
                            access="open",
                            electronic="electron_rich",
                            ortho_count=0,
                            burden="none",
                            steric_score=0.0,
                            electronic_score=-0.5,
                        ),
                    },
                ),
            },
        ),
        precedent_reference_ids=("REF1:one",),
        explanation=("Exact bond-edit and handle match",),
        score_trace=trace,
    )


def _result() -> GenericRecommendationResult:
    return GenericRecommendationResult(
        query_reaction_smiles="IC.B(O)O>>CC",
        valid=True,
        query_signature_id="RS3:query",
        reaction_label={
            "text": "Ar¹–Br + Ar²–B(OH)₂ → Ar¹–Ar²",
            "status": "structural_equation",
            "basis": "reaction_sites",
        },
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
                "reactivity_profile": _aromatic_profile(
                    access="hindered",
                    electronic="electron_poor",
                    ortho_count=2,
                    burden="high",
                    steric_score=0.7,
                    electronic_score=0.75,
                ),
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
    assert not window.unrestricted_fallback_check.isChecked()
    assert window.use_rxnmapper_check.isChecked()
    assert window.use_rxnmapper_check.objectName() == "useRxnMapper"
    assert window.unrestricted_fallback_check.objectName() == "unrestrictedFallback"
    assert window.reaction_edit.metaObject().className() == "QLineEdit"
    assert not window.completed_query_edit.isVisible()
    assert window.mode_combo.parentWidget() == window.reaction_row_label
    assert window.reaction_row_label.layout().indexOf(window.mode_combo) == 0
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
    assert window.results_table.columnCount() == 10
    assert not window.export_button.isEnabled()
    assert not any(
        label.text() == "Reaction Condition Recommender"
        for label in window.findChildren(gui.QtWidgets.QLabel)
    )


def _discovery_result() -> ReactionDiscoveryResult:
    trace = DiscoveryScoreTrace(
        components={
            "edit_similarity": 1.0,
            "reaction_center": 0.8,
            "local_environment": 0.7,
            "partner_category": 0.5,
            "spectator_groups": None,
            "reaction_topology": 1.0,
            "reactive_scaffold": 0.0,
        },
        contributions={"edit_similarity": 0.4},
        configured_weights={
            "edit_similarity": 0.35,
            "reaction_center": 0.20,
            "local_environment": 0.15,
            "partner_category": 0.10,
            "spectator_groups": 0.08,
            "reaction_topology": 0.07,
            "reactive_scaffold": 0.05,
        },
        effective_weights={"edit_similarity": 0.4},
        matches=("retrieval:bond_edit_signature",),
        mismatches=("reactive_scaffold=0.000",),
        definition_id="discovery_retrieval.v1",
        definition_version="1.0",
    )
    hit = ReactionDiscoveryHit(
        rank=1,
        reaction_id="discovery-1",
        observation_id="observation-1",
        canonical_reaction_id="canonical-1",
        reaction_smiles="BrC.B(O)O>>CC",
        reaction_label={"text": "Related coupling"},
        relation_class="direct_edit_analogue",
        relation_tiers=("bond_edit_signature",),
        discovery_score=0.81,
        yield_pct=12.0,
        outcome_status="usable",
        evidence_tier="trusted",
        chemistry_status="verified",
        source_dataset="literature",
        reference_id="REF1:discovery",
        resolved_recipe={
            "catalysts": [{"canonical_name": "Discovery catalyst"}],
            "temperature_c": 60.0,
        },
        recipe_id="RCR1:discovery",
        recipe_core_id="RCORE1:discovery",
        hypothesis_id=None,
        score_trace=trace,
        insights=("Observed conditions from literature",),
        cautions=("Low observed yield (12.0%) may reveal a failure boundary",),
    )
    return ReactionDiscoveryResult(
        query_reaction_smiles="IC.B(O)O>>CC",
        valid=True,
        reaction_label={"text": "Query coupling", "status": "structural_equation"},
        transformation_class="c_c_transfer_coupling",
        candidate_count=4,
        relation_counts={"direct_edit_analogue": 4},
        hits=(hit,),
        warnings=("DISCOVERY_CONDITIONS_ARE_OBSERVED_NOT_RECOMMENDED",),
    )


def test_worker_reuses_recommender_contract(monkeypatch) -> None:
    expected = _result()
    calls = []

    class FakeRecommender:
        def recommend(
            self,
            reaction,
            *,
            top_k,
            minimum_pool_size,
            unrestricted_fallback,
            ranking_preferences,
            completion_selections,
        ):
            calls.append(
                (
                    reaction,
                    top_k,
                    minimum_pool_size,
                    unrestricted_fallback,
                    ranking_preferences,
                    completion_selections,
                )
            )
            return expected

    monkeypatch.setattr(
        gui,
        "_get_cached_recommender",
        lambda path, *, use_rxnmapper, include_review: (
            calls.append(("index_options", use_rxnmapper, include_review))
            or FakeRecommender()
        ),
    )
    worker = gui.GenericRecommendationWorker(
        "index.sqlite",
        "A.B>>P",
        top_k=3,
        minimum_pool_size=2,
        unrestricted_fallback=True,
        use_rxnmapper=True,
    )
    progress = []
    finished = []
    worker.progress.connect(progress.append)
    worker.finished.connect(
        lambda success, result, error: finished.append((success, result, error))
    )

    worker.run()

    assert calls == [
        ("index_options", True, True),
        ("A.B>>P", 3, 2, True, None, ()),
    ]
    assert len(progress) == 2
    assert finished == [(True, expected, "")]


def test_discovery_worker_uses_separate_explorer_contract(monkeypatch) -> None:
    expected = _discovery_result()
    calls = []

    class FakeExplorer:
        def discover(self, reaction, **options):
            calls.append((reaction, options))
            return expected

    monkeypatch.setattr(
        gui,
        "_get_cached_explorer",
        lambda path, *, use_rxnmapper, include_review: (
            calls.append((path, use_rxnmapper, include_review)) or FakeExplorer()
        ),
    )
    worker = gui.ReactionDiscoveryWorker(
        "index.sqlite",
        "A.B>>P",
        top_k=8,
        view="failure_informed",
        include_low_yield=True,
        include_unreported_outcomes=False,
        use_rxnmapper=False,
    )
    finished = []
    worker.finished.connect(
        lambda success, result, error: finished.append((success, result, error))
    )

    worker.run()

    assert calls[0] == ("index.sqlite", False, False)
    assert calls[1] == (
        "A.B>>P",
        {
            "top_k": 8,
            "view": "failure_informed",
            "include_low_yield": True,
            "include_unreported_outcomes": False,
        },
    )
    assert finished == [(True, expected, "")]


def test_window_renders_recipe_summary_and_details(qtbot) -> None:
    window = gui.GenericRecommenderWindow()
    qtbot.addWidget(window)

    window._render_result(_result())

    summary = window.summary_box.toPlainText()
    assert "Ar¹–Br + Ar²–B(OH)₂ → Ar¹–Ar²" in summary
    assert "Label evidence: Structural Equation" in summary
    assert "Suzuki Miyaura" in summary
    assert "Mode: Verified Signature" in summary
    assert "R–C≡N [nitrile] (reactant 1, d=3)" in summary
    assert "benzene (6-membered)" in summary
    assert "ortho burden high (2/2)" in summary
    assert "electron demand electron poor" in summary
    reaction_pixmap = window.reaction_image_label.pixmap()
    assert reaction_pixmap is not None
    assert not reaction_pixmap.isNull()
    precedent_pixmap = window.selected_reaction_image_label.pixmap()
    assert precedent_pixmap is not None
    assert not precedent_pixmap.isNull()
    assert window.selected_reaction_image_label.toolTip() == "BrC.B(O)O>>CC"
    assert window.results_table.rowCount() == 1
    assert "Palladium catalyst" in window.results_table.item(0, 8).text()
    details = window.details_box.toPlainText()
    assert "Potassium carbonate" in details
    assert "Reaction label: Precedent Ar–Br coupling" in details
    assert "Ar¹–Br + Ar²–B(OH)₂ → Ar¹–Ar²" not in details
    assert "Selected hit reaction SMILES: BrC.B(O)O>>CC" in details
    assert "Selected hit reaction SMILES: IC.B(O)O>>CC" not in details
    assert "Spectator groups: R–O–R [ether] (reactant 2, d=2)" in details
    assert "R–C≡N [nitrile]" not in details
    assert "Reactivity profile:" in details
    assert "ortho burden none (0/2)" in details
    assert "electron demand electron rich" in details
    assert "ortho burden high" not in details
    assert "Recipe core:" not in details
    assert "Observed recipe variants:" not in details
    assert "RCORE1:core" not in details
    assert "RCR1:variant" not in details
    assert "Exact bond-edit and handle match" in (details)
    assert window.status_label.text() == "Done — 1 recipe(s)"


def test_window_switches_to_discovery_and_renders_evidence(qtbot) -> None:
    window = gui.GenericRecommenderWindow()
    qtbot.addWidget(window)

    window.mode_combo.setCurrentIndex(window.mode_combo.findData("discovery"))
    window._render_discovery_result(_discovery_result())

    assert window.run_button.text() == "Discover Related Reactions"
    assert window.results_heading.text() == "Related reaction precedents"
    assert window.results_table.horizontalHeaderItem(8).text() == "Observed conditions"
    assert window.results_table.rowCount() == 1
    assert "Discovery catalyst" in window.results_table.item(0, 8).text()
    assert "not condition recipes" in window.summary_box.toPlainText()
    details = window.details_box.toPlainText()
    assert "Observed conditions (precedent evidence, not a recommendation)" in details
    assert "configured weight 0.350" in details
    assert "failure boundary" in details
    assert window.status_label.text() == "Done — 1 analogue(s)"


def test_row_header_selection_moves_stale_current_cell(qtbot) -> None:
    window = gui.GenericRecommenderWindow()
    qtbot.addWidget(window)
    recommendations = tuple(
        replace(
            _recommendation(),
            rank=rank,
            recipe_id=f"RCR1:variant-{rank}",
            recipe_core_id=f"RCORE1:core-{rank}",
        )
        for rank in range(1, 4)
    )
    window._render_result(replace(_result(), recommendations=recommendations))
    table = window.results_table
    table.setCurrentCell(
        0,
        5,
        gui.QtCore.QItemSelectionModel.SelectionFlag.NoUpdate,
    )

    selection_model = table.selectionModel()
    selection_model.select(
        table.model().index(2, 0),
        gui.QtCore.QItemSelectionModel.SelectionFlag.ClearAndSelect
        | gui.QtCore.QItemSelectionModel.SelectionFlag.Rows,
    )

    assert table.currentRow() == 2
    assert table.currentColumn() == 0
    assert {index.row() for index in table.selectedIndexes()} == {2}


def test_query_summary_discloses_external_mapping_provenance() -> None:
    result = replace(
        _result(),
        external_mapping_status="external_mapping_internal_consensus",
        external_mapping_provider="rxnmapper",
        external_mapping_confidence=0.656,
        recommendation_mode="external_mapping_internal_consensus",
    )

    summary = gui.format_query_summary(result)

    assert "External Mapping Internal Consensus" in summary
    assert "provider rxnmapper" in summary
    assert "confidence 0.656" in summary


def test_completion_dialog_records_confirmed_and_edited_sources(qtbot) -> None:
    reaction = (
        "C=Cn1cc[n+](Cc2c(F)c(F)c(F)c(F)c2F)c1.[Br-]"
        ">>C=Cn1cc[n+](Cc2c(F)c(F)c(N=[N+]=[N-])c(F)c2F)c1.[Br-]"
    )
    proposal = gui.propose_reaction_completion(reaction)
    dialog = gui.CompletionSourceDialog(proposal)
    qtbot.addWidget(dialog)
    requirement = proposal.requirements[0]
    combo = dialog.source_combos[requirement.requirement_id]

    sodium_index = next(
        index
        for index, option in enumerate(requirement.options)
        if option.substance_id == "cas:26628-22-8"
    )
    combo.setCurrentIndex(sodium_index)
    confirmed = dialog.selections[0]
    assert confirmed.provenance == "user_confirmed"
    assert confirmed.substance_id == "cas:26628-22-8"

    combo.setEditText("26628-22-8")
    edited = dialog.selections[0]
    assert edited.provenance == "user_edited"
    assert edited.substance_id == "cas:26628-22-8"


def test_query_summary_discloses_completion_provenance() -> None:
    result = replace(
        _result(),
        completion_proposal={
            "requirements": (
                {"rooted_fragment_smiles": "*N=[N+]=[N-]"},
            )
        },
        completion_selections=(
            {
                "display_name": "NaN3",
                "provenance": "user_confirmed",
            },
        ),
        effective_query_reaction_smiles="A>[Na+].[N-]=[N+]=[N-]>P",
    )

    summary = gui.format_query_summary(result)

    assert "Missing source requirement: *N=[N+]=[N-]" in summary
    assert "Source completion: NaN3 [User Confirmed]" in summary
    assert "Completed query used: A>[Na+].[N-]=[N+]=[N-]>P" in summary


def test_stale_completion_result_reports_rebuild_status(qtbot) -> None:
    window = gui.GenericRecommenderWindow()
    qtbot.addWidget(window)
    result = replace(
        _result(),
        valid=False,
        effective_query_reaction_smiles="A>[Na+].[N-]=[N+]=[N-]>P",
        recommendations=(),
        error="RECOMMENDATION_INDEX_REBUILD_REQUIRED_FOR_COMPLETION",
    )

    window._render_result(result)

    assert window.status_label.text() == "Completed — index rebuild required"
    assert "confirmed source was applied" in window.summary_box.toPlainText()
    assert window.completed_query_edit.text() == "A>[Na+].[N-]=[N+]=[N-]>P"
    assert not window.completed_query_edit.isHidden()


def test_window_requests_vector_reaction_drawing(qtbot, monkeypatch) -> None:
    calls = []

    def render(reaction_smiles, *, size, image_format):
        calls.append((reaction_smiles, size, image_format))
        return (
            b'<svg xmlns="http://www.w3.org/2000/svg" width="560" height="142"></svg>'
        )

    monkeypatch.setattr(gui, "render_reaction_image_bytes", render)
    window = gui.GenericRecommenderWindow()
    qtbot.addWidget(window)

    window._render_reaction_graph("CCO>>CC=O")

    assert calls == [("CCO>>CC=O", gui.QUERY_REACTION_IMAGE_SIZE, "svg")]
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
