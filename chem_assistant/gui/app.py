"""Entry point for the ChemTools Qt assistant UI."""

from __future__ import annotations

import sys
from pathlib import Path

# Add project root to path for direct execution
if __name__ == "__main__":
    project_root = Path(__file__).parent.parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    # Enable relative imports when executed as a script (bypass package __init__)
    if __package__ is None:
        __package__ = "chem_assistant.gui"

from PyQt6.QtWidgets import QApplication

# Handle both relative and absolute imports
try:
    from .main_window import ChemAssistantWindow
except ImportError:
    from chem_assistant.gui.main_window import ChemAssistantWindow


def main() -> None:
    """Launch the Qt application."""
    app = QApplication(sys.argv)

    startup_message = None
    try:
        from chemtools.taxonomy import load_registry

        registry = load_registry()
        term_count = len(list(registry.iter_chem_terms()))
        startup_message = (
            f"Taxonomy v{registry.manifest.taxonomy_version} "
            f"(schema {registry.manifest.schema_version}) | "
            f"chem terms: {term_count}"
        )
    except Exception as exc:
        startup_message = f"Taxonomy unavailable: {exc}"

    window = ChemAssistantWindow(startup_message=startup_message)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
