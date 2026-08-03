"""Offscreen checks for the auto-detecting Qt6 featurizer."""

from __future__ import annotations

import os

import pytest


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
QtWidgets = pytest.importorskip("PyQt6.QtWidgets")

from app.featurizer_gui import (  # noqa: E402
    REACTION_EXAMPLE,
    ReactiveTaxonomyWindow,
    detect_input_kind,
)
import app.featurizer_gui as gui  # noqa: E402


def test_detect_input_kind() -> None:
    assert detect_input_kind("Brc1ccccc1") == "molecule"
    assert detect_input_kind("CCBr.N>>CCN") == "reaction"
    assert detect_input_kind("CCBr>N>CCN") == "reaction"


def test_window_analyzes_reaction_and_molecule() -> None:
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReactiveTaxonomyWindow()
    try:
        assert "background: #f4f7f9" not in window.styleSheet()
        assert "background: #ffffff" not in window.styleSheet()
        assert "background-color: #0078d7" in window.styleSheet()
        assert window.use_rxnmapper_check.isChecked()
        assert window.use_rxnmapper_check.objectName() == "useRxnMapper"
        assert not window.force_core_mapping_check.isChecked()
        assert window.force_core_mapping_check.objectName() == "forceCoreMapping"
        assert window.render_style_combo.objectName() == "renderStylePreset"
        assert window.render_style_combo.currentData() == "current"

        window.input_edit.setText(REACTION_EXAMPLE)
        assert window.kind_label.text() == "Detected: reaction"
        window.analyze()
        reaction_output = window.output.toPlainText()
        priority_review = window.review_output.toPlainText()
        assert reaction_output.startswith("REACTION FEATURIZATION")
        assert (
            "RXNMapper: not_requested_resolved_internal_evidence"
            in reaction_output
        )
        assert "Reaction: Ar–Br + Ar–B(OH)₂ → Ar–Ar" in reaction_output
        assert "Evidence: global_atom_correspondence" in reaction_output
        assert priority_review.startswith(
            "Reaction label: Ar–Br + Ar–B(OH)₂ → Ar–Ar"
        )
        assert "Spectators: None detected" in priority_review
        assert "Electronic / steric analysis: Unavailable" in priority_review
        assert "reaction input · valid" in window.status_label.text()
        reaction_pixmap = window.structure_image_label.pixmap()
        assert reaction_pixmap is not None
        assert not reaction_pixmap.isNull()
        assert window.graph_heading.text() == "Reaction graph"
        assert window.structure_image_label.toolTip() == REACTION_EXAMPLE
        assert window.graph_tabs.currentIndex() == 1
        core_pixmap = window.core_image_label.pixmap()
        assert core_pixmap is not None
        assert not core_pixmap.isNull()

        window.input_edit.setText("Brc1ccccc1")
        assert window.kind_label.text() == "Detected: molecule"
        window.analyze()
        molecule_output = window.output.toPlainText()
        assert molecule_output.startswith("MOLECULE FEATURIZATION")
        assert "Reactive-site hypotheses:" in molecule_output
        assert window.review_output.toPlainText() == (
            "Reaction review applies to reaction SMILES."
        )
        assert "Ar–Br — leaving_group, available" in molecule_output
        assert "molecule input · valid" in window.status_label.text()
        molecule_pixmap = window.structure_image_label.pixmap()
        assert molecule_pixmap is not None
        assert not molecule_pixmap.isNull()
        assert window.graph_heading.text() == "Compound graph"
        assert window.structure_image_label.toolTip() == "Brc1ccccc1"
        assert not window.graph_tabs.isTabEnabled(1)
    finally:
        window.close()
        application.processEvents()


def test_window_reports_empty_input_without_crashing() -> None:
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReactiveTaxonomyWindow()
    try:
        window.analyze()
        assert "Enter a molecule or reaction SMILES." in window.output.toPlainText()
        assert window.status_label.text() == "Analysis failed"
        assert window.structure_image_label.text() == (
            "Structure graph unavailable."
        )
    finally:
        window.close()
        application.processEvents()


def test_window_explains_ambiguous_reaction_evidence() -> None:
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReactiveTaxonomyWindow()
    reaction = (
        "O=C1CCCCC1.Cl.NNc1ccc(F)cc1"
        ">>Fc1ccc2[nH]c3c(c2c1)CCCC3"
    )
    try:
        window.use_rxnmapper_check.setChecked(False)
        window.input_edit.setText(reaction)
        window.analyze()
        output = window.output.toPlainText()

        assert "RXNMapper: disabled" in output
        assert "Correspondence ambiguity: 2 distinct edit hypotheses" in output
        assert output.count("REH1:") == 2
        assert "4 correspondences; unverified" in output
        assert "Atoms not in the main product: Cl × 1, N × 1, O × 1" in output
        assert "Net bond inventory (unmapped, not verified edits):" in output
        assert "Retrieval: not eligible" in output
        assert "ambiguous edit hypotheses" in output
    finally:
        window.close()
        application.processEvents()


def test_window_displays_mapped_reaction_minimization() -> None:
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReactiveTaxonomyWindow()
    reaction = (
        "[CH3:1][OH:2].O[CH3:5]."
        "[CH:3](=[O:4])[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
        ">>[CH3:1][O:2][CH:3]([O:4][CH3:5])"
        "[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
    )
    try:
        window.use_rxnmapper_check.setChecked(False)
        window.input_edit.setText(reaction)
        window.analyze()
        output = window.output.toPlainText()

        assert "Reaction minimization:" in output
        assert "Reaction: R–CH=O + Alk–OH" in output
        assert "Core evidence: verified (validated_atom_mapping" in output
        assert "Core shape (retrieval): RSH2:" in output
        assert "Center transition (diagnostic only): RCS2:" in output
        assert "retained aryl [Fc1ccccc1] (1 port)" in output
        core_pixmap = window.core_image_label.pixmap()
        assert core_pixmap is not None
        assert not core_pixmap.isNull()
        assert window.graph_tabs.currentIndex() == 1
        assert "Ar = Fc1ccccc1" in window.core_graphic_note.text()
        assert "Alk = C" in window.core_graphic_note.text()
    finally:
        window.close()
        application.processEvents()


def test_drawing_style_switch_rerenders_without_reanalysis(monkeypatch) -> None:
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    calls = []
    original_renderer = gui.render_molecule_image_bytes

    def recording_renderer(*args, **kwargs):
        calls.append(kwargs.get("render_preset"))
        return original_renderer(*args, **kwargs)

    monkeypatch.setattr(gui, "render_molecule_image_bytes", recording_renderer)
    window = ReactiveTaxonomyWindow()
    try:
        window.input_edit.setText("CCO")
        window.analyze()
        assert calls == ["current"]
        assert window.structure_image_label._trim_svg_white_space
        assert window.structure_image_label._svg_max_upscale is None
        current_pixmap = window.structure_image_label.pixmap()
        current_width = current_pixmap.width() / current_pixmap.devicePixelRatio()

        compact_index = window.render_style_combo.findData("acs_1996_compact")
        window.render_style_combo.setCurrentIndex(compact_index)

        assert calls == ["current", "acs_1996_compact"]
        assert window.structure_image_label._trim_svg_white_space
        assert window.structure_image_label._svg_max_upscale == 6.0
        compact_pixmap = window.structure_image_label.pixmap()
        compact_width = compact_pixmap.width() / compact_pixmap.devicePixelRatio()
        assert compact_width > current_width * 0.15
        assert compact_width < current_width * 0.60
        assert window.output.toPlainText().startswith("MOLECULE FEATURIZATION")
    finally:
        window.close()
        application.processEvents()


def test_main_window_starts_maximized() -> None:
    events = []

    class FakeWindow:
        def showMaximized(self) -> None:
            events.append("maximized")

    gui._show_main_window(FakeWindow())

    assert events == ["maximized"]
