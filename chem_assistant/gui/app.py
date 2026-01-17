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

    startup_parts = []
    try:
        from chemtools.taxonomy.reaction_catalog import load_reaction_catalog
        from chemtools.reagent.reagent_v2 import ReagentTaxonomyV2

        # Load new taxonomy v2
        rxn_defs, _ = load_reaction_catalog()
        reagent_tax = ReagentTaxonomyV2.from_path()

        rxn_count = len(rxn_defs)
        reagent_count = len(list(reagent_tax.iter_families()))

        startup_parts.append(
            f"Taxonomy v2 | Reactions: {rxn_count} | Reagent Families: {reagent_count}"
        )
    except Exception as exc:
        startup_parts.append(f"Taxonomy v2 unavailable: {exc}")

    try:
        from chem_assistant.chemtools_wrapper import CHEMTOOLS_TOOLS

        tool_names = {tool.name for tool in CHEMTOOLS_TOOLS}
        hte_status = "HTE tools: ready" if "hte_recommend_conditions" in tool_names else "HTE tools: missing"
        startup_parts.append(hte_status)
    except Exception as exc:
        startup_parts.append(f"HTE tools unavailable: {exc}")

    startup_message = " | ".join(startup_parts).strip() or None
    window = ChemAssistantWindow(startup_message=startup_message)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
