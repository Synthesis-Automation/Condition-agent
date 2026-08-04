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
    format_core_graph_analysis,
)
import app.featurizer_gui as gui  # noqa: E402
from reactive_taxonomy import featurize_reaction  # noqa: E402


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
        assert window.core_analysis_heading.text() == "Core graph analysis"
        assert (
            window.core_analysis_heading.objectName()
            == "coreGraphAnalysisHeading"
        )
        assert (
            window.core_analysis_output.objectName()
            == "coreGraphAnalysisOutput"
        )
        assert not hasattr(window, "review_output")

        window.input_edit.setText(REACTION_EXAMPLE)
        assert window.kind_label.text() == "Detected: reaction"
        window.analyze()
        reaction_output = window.output.toPlainText()
        core_analysis = window.core_analysis_output.toPlainText()
        assert reaction_output.startswith("REACTION FEATURIZATION")
        assert (
            "RXNMapper: not needed; internal correspondence resolved"
            in reaction_output
        )
        assert (
            "Observation evidence: whole-reaction atom correspondence"
            in reaction_output
        )
        assert "Generic signature: available" in reaction_output
        assert "Primary pattern: organoboron C–C coupling" in reaction_output
        assert "RS3:" not in reaction_output
        assert "Reaction minimization:" not in reaction_output
        assert "Product connection:" not in reaction_output
        assert core_analysis.startswith(
            "Reaction: Ar–Br + Ar–B(OH)₂ → Ar–Ar"
        )
        assert "Evidence:" in core_analysis
        assert "Bond changes:" in core_analysis
        assert "R-group attachment profiles:" in core_analysis
        assert "RSH2:" not in core_analysis
        assert "RSE1:" not in core_analysis
        assert "reaction input · valid" in window.status_label.text()
        reaction_pixmap = window.structure_image_label.pixmap()
        assert reaction_pixmap is not None
        assert not reaction_pixmap.isNull()
        assert window.graph_heading.text() == "Reaction graph"
        assert window.structure_image_label.toolTip() == REACTION_EXAMPLE
        graph_layout = window.full_structure_panel.parentWidget().layout()
        assert graph_layout.indexOf(window.full_structure_panel) < (
            graph_layout.indexOf(window.minimized_panel)
        )
        assert window.full_structure_panel.title() == "Full structure"
        assert window.minimized_panel.title() == "Minimized reaction"
        result_layout = (
            window.full_structure_panel.parentWidget().parentWidget().layout()
        )
        assert result_layout.stretch(0) == 9
        assert result_layout.stretch(1) == 11
        core_pixmap = window.core_image_label.pixmap()
        assert core_pixmap is not None
        assert not core_pixmap.isNull()

        window.input_edit.setText("Brc1ccccc1")
        assert window.kind_label.text() == "Detected: molecule"
        window.analyze()
        molecule_output = window.output.toPlainText()
        assert molecule_output.startswith("MOLECULE FEATURIZATION")
        assert "Reactive-site hypotheses:" in molecule_output
        assert window.core_analysis_output.toPlainText() == (
            "Core graph analysis applies to reaction SMILES."
        )
        assert "Ar–Br — leaving_group, available" in molecule_output
        assert "molecule input · valid" in window.status_label.text()
        molecule_pixmap = window.structure_image_label.pixmap()
        assert molecule_pixmap is not None
        assert not molecule_pixmap.isNull()
        assert window.graph_heading.text() == "Compound graph"
        assert window.structure_image_label.toolTip() == "Brc1ccccc1"
        assert window.core_image_label.text() == (
            "Reaction minimization applies only to reactions."
        )
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
        assert "Generic signature: unavailable" in output
        assert "Completeness:" in output
        assert "Warnings:" in output
        assert "REH1:" not in output
        assert "Net bond inventory" not in output
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
        core_analysis = window.core_analysis_output.toPlainText()

        assert "Reaction minimization:" not in output
        assert "Reaction: R–CH=O + Alk–OH" in core_analysis
        assert "Evidence: Verified from validated atom mapping" in core_analysis
        assert "RSH2:" not in core_analysis
        assert "Bond changes:" in core_analysis
        assert "R-group attachment profiles:" in core_analysis
        core_pixmap = window.core_image_label.pixmap()
        assert core_pixmap is not None
        assert not core_pixmap.isNull()
        assert "Display-only minimum:" in window.core_graphic_note.text()
        assert "R groups:" in window.core_graphic_note.text()
        assert "removed unchanged aromatic substituents:" in (
            window.core_graphic_note.text()
        )
    finally:
        window.close()
        application.processEvents()


def test_window_uses_r_group_display_minimization() -> None:
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReactiveTaxonomyWindow()
    try:
        window.use_rxnmapper_check.setChecked(False)
        window.input_edit.setText("C=CC1=CC=CC=C1>>CCc1ccccc1")
        window.analyze()

        core_pixmap = window.core_image_label.pixmap()
        assert core_pixmap is not None
        assert not core_pixmap.isNull()
        assert "Display-only minimum: *C=C>>*CC" in (
            window.core_graphic_note.text()
        )
        assert "R groups: 2" in window.core_graphic_note.text()
        assert "R-group attachment profiles:" in (
            window.core_analysis_output.toPlainText()
        )
    finally:
        window.close()
        application.processEvents()


def test_window_reports_multisite_hidden_connector() -> None:
    application = (
        QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    )
    window = ReactiveTaxonomyWindow()
    reaction = (
        "Cc1c([N+](=O)[O-])cnc2c1c("
        "C1=CCN(C(=O)OC(C)(C)C)C(C)(C)C1)cn2C"
        ">>Cc1c(N)cnc2c1c("
        "C1CCN(C(=O)OC(C)(C)C)C(C)(C)C1)cn2C"
    )
    try:
        window.use_rxnmapper_check.setChecked(False)
        window.input_edit.setText(reaction)
        window.analyze()

        core_pixmap = window.core_image_label.pixmap()
        assert core_pixmap is not None
        assert not core_pixmap.isNull()
        assert "hidden connectors: S¹ R³⋯R⁴ " in (
            window.core_graphic_note.text()
        )
        assert "(4 hidden-path atoms)" in window.core_graphic_note.text()
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


def test_core_graph_analysis_exposes_shared_r_ports_and_site_context() -> None:
    reaction = (
        "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O."
        "CC(C)c1nc(Cl)nc(Cl)c1Br>>"
        "CC(C)c1nc(Cl)nc(N2CCN(c3cc4c(cc3F)c(=O)c(C(=O)O)"
        "cn4C3CC3)CC2)c1Br"
    )

    text = format_core_graph_analysis(featurize_reaction(reaction))

    assert "R¹/R² — retained; primary carbon attachment; saturated ring" in text
    assert "R¹ and R² are connected through the same omitted scaffold." in text
    assert "N two bonds away" in text
    assert "break C–Cl; form C–N; remove N–H" in text
    assert "Unchanged functional groups on the R¹/R² scaffold:" in text
    assert "tertiary amine (2 bonds from each R attachment" in text
    assert "aryl halide (4 bonds from each R attachment" in text
    assert "carboxylic acid (9 bonds from each R attachment" in text
    assert "Active-site steric/electronic context:" in text
    assert "secondary N" in text
    assert "hindered access" in text
    assert "high lone-pair availability" in text
    assert "RSH2:" not in text
    assert "RSE1:" not in text
    assert "REACTION_CORE" not in text


def test_main_window_starts_maximized() -> None:
    events = []

    class FakeWindow:
        def showMaximized(self) -> None:
            events.append("maximized")

    gui._show_main_window(FakeWindow())

    assert events == ["maximized"]
