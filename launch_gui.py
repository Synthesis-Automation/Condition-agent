"""
Launcher script for ChemAssistant GUI Application.

This script provides an easy way to launch the ChemTools assistant
graphical user interface.
"""

import sys
from pathlib import Path

# Ensure project root is in path
project_root = Path(__file__).parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from chem_assistant.gui.app import main

if __name__ == "__main__":
    main()
