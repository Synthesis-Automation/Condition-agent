"""
Retron SMARTS Pattern Library.

A retron is a substructure in a target molecule that implies a specific
retrosynthetic transform — i.e., a bond whose formation defines the
"last step" in a synthetic route.

Each entry contains:
  name            : unique identifier
  product_smarts  : SMARTS that matches the bond/motif in the TARGET molecule
  reaction_name   : legacy retro transform label (kept for compatibility)
  retro_transform_id : preferred explicit retro transform identifier (v2)
  taxonomy_family_id : preferred canonical taxonomy family link (v2)
  taxonomy_id     : legacy compatibility field (normalized to canonical family ID at load)
  difficulty      : 0.0 (trivial) → 1.0 (heroic); guides ranking
  description     : human-readable retrosynthetic description
  notes           : chemistry notes and caveats
  precursor_hints : list of precursor type names (for LLM context)
  category        : reaction class grouping (for browsing / filtering)

Data is loaded from the retro-owned file
``chemtools/retro/data/retron_patterns.json`` at import time so the library
is easy to extend without editing Python source.
"""
from __future__ import annotations

from typing import Any, Dict, List


def get_retron_taxonomy_id(retron: Dict[str, Any]) -> str:
    """Return canonical taxonomy family ID for a retron (v2-first, v1 fallback)."""
    if not isinstance(retron, dict):
        return ""
    return str(retron.get("taxonomy_family_id") or retron.get("taxonomy_id") or "")


def _load() -> List[Dict[str, Any]]:
    from chemtools.taxonomy import loader as taxonomy_loader

    payload = taxonomy_loader.load_retron_patterns()
    entries = payload.get("retrons") if isinstance(payload, dict) else None
    if not isinstance(entries, list):
        return []

    normalized: List[Dict[str, Any]] = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        item = dict(entry)
        canonical_tid = get_retron_taxonomy_id(item)
        if canonical_tid:
            # Keep both v2 and legacy link fields aligned during migration.
            item["taxonomy_family_id"] = canonical_tid
            item["taxonomy_id"] = canonical_tid
        normalized.append(item)
    return normalized


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
