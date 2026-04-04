"""
HTE-backed Retrosynthetic Template Library.

Each template is a retrosynthetic SMARTS in the form:
    [product_pattern] >> [precursor_1] . [precursor_2]

Applied via AllChem.RunReactants((target_mol,)) — the target is the "reactant"
and the output tuples are the precursor pairs.  Atom-mapped atoms (:1, :2, ...)
route the correct molecular context to each precursor fragment; the surrounding
un-mapped atoms follow their mapped neighbours automatically.

Each entry contains:
  name          : unique identifier
  hte_families  : list of HTE CSV family stem(s) this template covers
  taxonomy_family_id : preferred canonical taxonomy family link (v2)
  taxonomy_id   : legacy compatibility field (normalized to canonical family ID at load)
  retro_smarts  : retrosynthetic SMARTS (product >> precursor(s))
  description   : human-readable
  difficulty    : 0.0–1.0
  n_precursors  : 1 or 2 (output fragments from RunReactants)
  notes         : chemistry notes and caveats (optional)
  category      : reaction class grouping (for browsing / filtering)

Data is loaded from the retro-owned file
``chemtools/retro/data/hte_templates.json`` at import time so the library
is easy to extend without editing Python source.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional


def get_template_taxonomy_id(template: Dict[str, Any]) -> str:
    """Return canonical taxonomy family ID for a template (v2-first, v1 fallback)."""
    if not isinstance(template, dict):
        return ""
    return str(template.get("taxonomy_family_id") or template.get("taxonomy_id") or "")


def _load() -> List[Dict[str, Any]]:
    from chemtools.taxonomy import loader as taxonomy_loader

    payload = taxonomy_loader.load_hte_templates()
    entries = payload.get("templates") if isinstance(payload, dict) else None
    if not isinstance(entries, list):
        return []

    normalized: List[Dict[str, Any]] = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        item = dict(entry)
        canonical_tid = get_template_taxonomy_id(item)
        if canonical_tid:
            # Keep legacy field populated with canonical family ID for callers that
            # still read `taxonomy_id` only during migration.
            item["taxonomy_family_id"] = canonical_tid
            item["taxonomy_id"] = canonical_tid
        normalized.append(item)
    return normalized


HTE_TEMPLATES: List[Dict[str, Any]] = _load()


# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------

def get_template_by_name(name: str) -> Optional[Dict[str, Any]]:
    """Return a template dict by name, or None if not found."""
    for t in HTE_TEMPLATES:
        if t["name"] == name:
            return t
    return None


def get_templates_for_family(hte_family: str) -> List[Dict[str, Any]]:
    """Return all templates whose hte_families list contains the given family stem."""
    family_lower = hte_family.lower()
    return [
        t for t in HTE_TEMPLATES
        if any(f.lower() == family_lower for f in t.get("hte_families", []))
    ]


def get_all_template_names() -> List[str]:
    return [t["name"] for t in HTE_TEMPLATES]
