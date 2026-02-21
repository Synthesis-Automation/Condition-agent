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
  taxonomy_id   : canonical taxonomy cross-reference (reaction_registry.py)
  retro_smarts  : retrosynthetic SMARTS (product >> precursor(s))
  description   : human-readable
  difficulty    : 0.0–1.0
  n_precursors  : 1 or 2 (output fragments from RunReactants)
  notes         : chemistry notes and caveats (optional)
  category      : reaction class grouping (for browsing / filtering)

Data is loaded from chemtools/retro/data/hte_templates.json at import
time so the library is easy to extend without editing Python source.
"""
from __future__ import annotations

import json
import pathlib
from typing import Any, Dict, List, Optional

_DATA_FILE = pathlib.Path(__file__).parent.parent / "taxonomy" / "data" / "hte_templates.json"


def _load() -> List[Dict[str, Any]]:
    with open(_DATA_FILE, encoding="utf-8") as fh:
        return json.load(fh)["templates"]


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
