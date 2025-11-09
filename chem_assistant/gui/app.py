"""Entry point for the ChemTools Qt assistant UI."""

from __future__ import annotations

import sys
from pathlib import Path

# Add project root to path for direct execution
if __name__ == "__main__":
    project_root = Path(__file__).parent.parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

from PyQt6.QtWidgets import QApplication

# Handle both relative and absolute imports
try:
    from .main_window import ChemAssistantWindow
except ImportError:
    from chem_assistant.gui.main_window import ChemAssistantWindow


def main() -> None:
    """Launch the Qt application."""
    app = QApplication(sys.argv)
    window = ChemAssistantWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
