"""Minimal offscreen smoke coverage for the Qt6 template wrapper."""

from __future__ import annotations

import os

import pytest


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
QtWidgets = pytest.importorskip("PyQt6.QtWidgets")

from app.reaction_template_registry_gui import (  # noqa: E402
    ReactionTemplateRegistryWindow,
)


def test_qt6_registry_window_exposes_authoring_and_query_controls() -> None:
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReactionTemplateRegistryWindow()
    try:
        assert window.findChild(QtWidgets.QLineEdit, "registryPath") is not None
        assert window.findChild(QtWidgets.QLineEdit, "templateId") is not None
        assert (
            window.findChild(
                QtWidgets.QLineEdit, "mappedReferenceReaction"
            )
            is not None
        )
        assert window.findChild(QtWidgets.QLineEdit, "queryReaction") is not None
        assert window.findChild(QtWidgets.QPushButton, "importTemplate") is not None
        assert (
            window.findChild(QtWidgets.QPushButton, "featurizeQuery")
            is not None
        )
        assert window.findChild(QtWidgets.QPushButton, "matchQuery") is not None
        assert "background-color: #ffffff" not in window.styleSheet()
        assert "color: #23313f" in window.status_label.styleSheet()
    finally:
        window.close()
        application.processEvents()


def test_qt6_registry_window_featurizes_reaction_test_input() -> None:
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReactionTemplateRegistryWindow()
    try:
        reaction = "CCO.COc1cccc(C=O)c1>>CCOC(OCC)c1cccc(OC)c1"
        window.query_edit.setText(reaction)
        window.featurize_query()

        details = window.details.toPlainText()
        assert details.startswith("REACTION FEATURIZATION")
        assert f"Input: {reaction}" in details
        assert "Status: valid" in details
        assert "Status: incomplete" in details
        assert "RS3 signature: unavailable" in details
        assert "TEMPLATE REGISTRY MATCHES (1)" in details
        assert "Family: acetalization" in details
        assert "Match: multiplicity-assisted exact reconstruction" in details
        assert "Confidence: 0.85" in details
        assert "Reaction label: Acetal formation" in details
        assert (
            "Structural label: R–CH=O + 2 × R–OH → acetal"
            in details
        )
        assert "carbonyl: R–CH=O" in details
        assert "Steric: moderate (0.28); electronic: slightly_poor" in details
        assert "alcohol × 2: R–OH" in details
        assert "Steric: open (0.10); electronic: medium" in details
        assert (
            "Reconstructed product: CCOC(OCC)c1cccc(OC)c1"
            in details
        )
        assert "Inferred repeated participant: yes" in details
        assert not details.lstrip().startswith("{")
        assert "Featurization: valid" in window.status_label.text()
        assert "completeness incomplete" in window.status_label.text()
        assert (
            "acetalization (multiplicity-assisted exact)"
            in window.status_label.text()
        )
    finally:
        window.close()
        application.processEvents()
