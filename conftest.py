"""
Root conftest.py — ensures the project root is on sys.path so that both
`pytest` and `python -m pytest` find local packages (chemtools, chem_coworker,
app, llmtools, etc.) without requiring an editable install.
"""
import sys
from pathlib import Path

# Insert project root at the front of sys.path if it isn't already there.
_ROOT = str(Path(__file__).parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
