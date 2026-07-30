"""Offscreen checks for the auto-detecting Qt6 featurizer."""

from __future__ import annotations

import os

import pytest


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
QtWidgets = pytest.importorskip("PyQt6.QtWidgets")

from app.reactive_taxonomy_gui import (  # noqa: E402
    REACTION_EXAMPLE,
    ReactiveTaxonomyWindow,
    detect_input_kind,
)


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

        window.input_edit.setText(REACTION_EXAMPLE)
        assert window.kind_label.text() == "Detected: reaction"
        window.analyze()
        reaction_output = window.output.toPlainText()
        assert reaction_output.startswith("REACTION FEATURIZATION")
        assert "Reaction: Ar1–Br + Ar2–B(OH)2 → Ar1–Ar2" in reaction_output
        assert "Evidence: exact_product_reconstruction" in reaction_output
        assert "reaction input · valid" in window.status_label.text()

        window.input_edit.setText("Brc1ccccc1")
        assert window.kind_label.text() == "Detected: molecule"
        window.analyze()
        molecule_output = window.output.toPlainText()
        assert molecule_output.startswith("MOLECULE FEATURIZATION")
        assert "Reactive sites:" in molecule_output
        assert "Ar–Br — leaving_group, available" in molecule_output
        assert "molecule input · valid" in window.status_label.text()
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
        window.input_edit.setText(reaction)
        window.analyze()
        output = window.output.toPlainText()

        assert "Correspondence ambiguity: 2 distinct edit hypotheses" in output
        assert "Atoms not in the main product: Cl × 1, N × 1, O × 1" in output
        assert "Net bond inventory (unmapped, not verified edits):" in output
        assert "Retrieval: not eligible" in output
        assert "does not invalidate the product structure" in output
    finally:
        window.close()
        application.processEvents()
