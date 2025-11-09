"""Entry point for the ChemTools Qt assistant UI."""

from __future__ import annotations

import sys

from PyQt6.QtWidgets import QApplication

from .main_window import ChemAssistantWindow


def main() -> None:
    """Launch the Qt application."""
    app = QApplication(sys.argv)
    window = ChemAssistantWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
