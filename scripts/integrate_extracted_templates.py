"""
Integrate high-quality extracted retro templates into the production
hte_templates.json and retron_patterns.json files.

This script is a thin wrapper — the implementation now lives in
``chemtools.retro.integrate_templates``.

Usage:
    python scripts/integrate_extracted_templates.py

Or use the module directly:
    python -m chemtools.retro.integrate_templates
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from chemtools.retro.integrate_templates import integrate

if __name__ == "__main__":
    n_templates, n_retrons = integrate()
