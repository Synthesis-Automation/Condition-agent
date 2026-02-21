"""
Retron SMARTS Pattern Library.

A retron is a substructure in a target molecule that implies a specific
retrosynthetic transform — i.e., a bond whose formation defines the
"last step" in a synthetic route.

Each entry contains:
  name            : unique identifier
  product_smarts  : SMARTS that matches the bond/motif in the TARGET molecule
  reaction_name   : maps to the existing chemtools reaction taxonomy ID
  taxonomy_id     : canonical taxonomy cross-reference (reaction_registry.py)
  difficulty      : 0.0 (trivial) → 1.0 (heroic); guides ranking
  description     : human-readable retrosynthetic description
  notes           : chemistry notes and caveats
  precursor_hints : list of precursor type names (for LLM context)
  category        : reaction class grouping (for browsing / filtering)

Data is loaded from chemtools/retro/data/retron_patterns.json at import
time so the library is easy to extend without editing Python source.
"""
from __future__ import annotations

import json
import pathlib
from typing import Any, Dict, List

_DATA_FILE = pathlib.Path(__file__).parent.parent / "taxonomy" / "data" / "retron_patterns.json"


def _load() -> List[Dict[str, Any]]:
    with open(_DATA_FILE, encoding="utf-8") as fh:
        return json.load(fh)["retrons"]


RETRONS: List[Dict[str, Any]] = _load()


# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------

def get_retron_by_name(name: str) -> Dict[str, Any]:
    """Return a retron dict by its name, or empty dict if not found."""
    for r in RETRONS:
        if r["name"] == name:
            return r
    return {}


def get_retrons_by_reaction(reaction_name: str) -> List[Dict[str, Any]]:
    """Return all retrons that map to a given reaction taxonomy name."""
    return [r for r in RETRONS if r.get("reaction_name") == reaction_name]


def get_retrons_by_difficulty(max_difficulty: float = 1.0) -> List[Dict[str, Any]]:
    """Return all retrons with difficulty <= max_difficulty, sorted ascending."""
    filtered = [r for r in RETRONS if r.get("difficulty", 0.5) <= max_difficulty]
    return sorted(filtered, key=lambda r: r.get("difficulty", 0.5))
