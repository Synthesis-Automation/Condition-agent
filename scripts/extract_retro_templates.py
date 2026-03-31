"""
Extract retrosynthetic templates from HTE reaction data.

This script is a thin wrapper — the implementation now lives in
``chemtools.retro.extract_templates``.

Usage:
    python scripts/extract_retro_templates.py                    # small test (3 files)
    python scripts/extract_retro_templates.py --all              # all literature CSVs
    python scripts/extract_retro_templates.py --file suzuki      # specific file pattern
    python scripts/extract_retro_templates.py --limit 200        # cap rows per file
    python scripts/extract_retro_templates.py --min-confidence 0.6  # confidence filter
    python scripts/extract_retro_templates.py --min-count 3      # frequency filter

Or use the module directly:
    python -m chemtools.retro.extract_templates --all --limit 0
"""
from __future__ import annotations

import sys
from pathlib import Path

# Ensure project root is importable
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from chemtools.retro.extract_templates import main

if __name__ == "__main__":
    main()
