import os
import sys

# Ensure project root is importable for local package imports
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Speed up tests by disabling RDKit-heavy paths; the code has graceful fallbacks
os.environ.setdefault("CHEMTOOLS_DISABLE_RDKIT", "1")

