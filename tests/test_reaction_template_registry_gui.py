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
                QtWidgets.QPlainTextEdit, "mappedReferenceReaction"
            )
            is not None
        )
        assert window.findChild(QtWidgets.QPushButton, "importTemplate") is not None
        assert window.findChild(QtWidgets.QPushButton, "matchQuery") is not None
    finally:
        window.close()
        application.processEvents()
