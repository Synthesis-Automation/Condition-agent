"""
Unified Reaction Registry
=========================

Central lookup table linking the three retrosynthesis knowledge systems:
  - Retron patterns (SMARTS-based matching, chemtools/retro/retron_patterns.py)
  - HTE templates  (RDKit RunReactants, chemtools/retro/hte_templates.py)
  - Taxonomy       (reaction classification, chemtools/taxonomy/data/)

Every retron and HTE template carries a ``taxonomy_id`` field that maps to a
canonical taxonomy ID.  This module aggregates those links at import time and
exposes bidirectional lookups:

  taxonomy_id → retrons, templates, hte_families, difficulty, description
  retron_name → taxonomy_id
  template_name → taxonomy_id
  retron_name → hte_families (for precedent search)

Usage
-----
    from chemtools.retro.reaction_registry import (
        get_by_taxonomy_id,
        get_taxonomy_id_for_retron,
        get_taxonomy_id_for_template,
        get_hte_families_for_retron,
        validate_registry,
    )

    entry = get_by_taxonomy_id("Suzuki_miyaura")
    entry.retron_names   # ["biaryl_suzuki"]
    entry.template_names # ["suzuki_miyaura"]
    entry.hte_families   # ["suzuki_miyaura"]
    entry.difficulty     # 0.15

    warnings = validate_registry()
    # [] if all cross-references are consistent
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class ReactionEntry:
    """Unified metadata for a single reaction type, aggregated from all three systems."""

    taxonomy_id: str
    """Canonical taxonomy identifier (matches ``id`` in reaction_types.v4.0.json)."""

    retron_names: List[str] = field(default_factory=list)
    """Names of retron patterns that map to this reaction type."""

    template_names: List[str] = field(default_factory=list)
    """Names of HTE templates that map to this reaction type."""

    hte_families: List[str] = field(default_factory=list)
    """HTE database family names aggregated from matching templates."""

    difficulty: float = 0.5
    """
    Unified difficulty score (0 = trivial, 1 = heroic).
    Set to the minimum of all contributing retron/template difficulty values.
    """

    description: str = ""
    """Human-readable description from the taxonomy entry (if available)."""


# ---------------------------------------------------------------------------
# Registry construction
# ---------------------------------------------------------------------------

def _build_registry() -> Dict[str, ReactionEntry]:
    """
    Build the unified reaction registry at import time.

    1. Load retrons and templates (each now carry a ``taxonomy_id`` field).
    2. Group by taxonomy_id.
    3. Resolve difficulty as min(all contributing values).
    4. Optionally enrich description from the taxonomy JSON.
    """
    from .retron_patterns import RETRONS
    from .hte_templates import HTE_TEMPLATES

    registry: Dict[str, ReactionEntry] = {}

    # -- Retrons --
    for retron in RETRONS:
        tid = retron.get("taxonomy_id")
        if not tid:
            continue
        if tid not in registry:
            registry[tid] = ReactionEntry(taxonomy_id=tid)
        entry = registry[tid]
        name = retron["name"]
        if name not in entry.retron_names:
            entry.retron_names.append(name)
        diff = retron.get("difficulty", 0.5)
        entry.difficulty = min(entry.difficulty, diff)
        if not entry.description and retron.get("description"):
            entry.description = retron["description"]

    # -- HTE Templates --
    for tpl in HTE_TEMPLATES:
        tid = tpl.get("taxonomy_id")
        if not tid:
            continue
        if tid not in registry:
            registry[tid] = ReactionEntry(taxonomy_id=tid)
        entry = registry[tid]
        name = tpl["name"]
        if name not in entry.template_names:
            entry.template_names.append(name)
        for fam in tpl.get("hte_families", []):
            if fam not in entry.hte_families:
                entry.hte_families.append(fam)
        diff = tpl.get("difficulty", 0.5)
        entry.difficulty = min(entry.difficulty, diff)
        if not entry.description and tpl.get("description"):
            entry.description = tpl["description"]

    # -- Enrich descriptions from taxonomy loader (best-effort) --
    try:
        from chemtools.taxonomy.loader import load_reaction_types
        tax_entries = load_reaction_types()
        # Build a quick id→description map
        desc_map: Dict[str, str] = {}
        for te in tax_entries:
            desc = getattr(te, "description", None) or ""
            if te.id:
                desc_map[te.id] = desc
        for tid, entry in registry.items():
            if not entry.description and desc_map.get(tid):
                entry.description = desc_map[tid]
    except Exception as exc:
        logger.debug("Could not enrich registry from taxonomy loader: %s", exc)

    return registry


# Registry singleton — built once at import time
REACTION_REGISTRY: Dict[str, ReactionEntry] = _build_registry()


# ---------------------------------------------------------------------------
# Reverse-lookup indexes (built from the registry)
# ---------------------------------------------------------------------------

# retron_name → taxonomy_id
_RETRON_TO_TAXONOMY: Dict[str, str] = {}
# template_name → taxonomy_id
_TEMPLATE_TO_TAXONOMY: Dict[str, str] = {}

for _entry in REACTION_REGISTRY.values():
    for _rn in _entry.retron_names:
        _RETRON_TO_TAXONOMY[_rn] = _entry.taxonomy_id
    for _tn in _entry.template_names:
        _TEMPLATE_TO_TAXONOMY[_tn] = _entry.taxonomy_id


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_by_taxonomy_id(taxonomy_id: str) -> Optional[ReactionEntry]:
    """Return the registry entry for *taxonomy_id*, or None if not found."""
    return REACTION_REGISTRY.get(taxonomy_id)


def get_taxonomy_id_for_retron(retron_name: str) -> Optional[str]:
    """Return the canonical taxonomy_id for a retron, or None if unmapped."""
    return _RETRON_TO_TAXONOMY.get(retron_name)


def get_taxonomy_id_for_template(template_name: str) -> Optional[str]:
    """Return the canonical taxonomy_id for an HTE template, or None if unmapped."""
    return _TEMPLATE_TO_TAXONOMY.get(template_name)


def get_hte_families_for_retron(retron_name: str) -> List[str]:
    """
    Return HTE database family names applicable to a retron.

    Resolves: retron_name → taxonomy_id → ReactionEntry.hte_families.

    Falls back to the retron's own ``reaction_name`` if no taxonomy entry is
    found (maintains backward compatibility).
    """
    tid = _RETRON_TO_TAXONOMY.get(retron_name)
    if tid:
        entry = REACTION_REGISTRY.get(tid)
        if entry and entry.hte_families:
            return entry.hte_families

    # Fallback: look up directly by reaction_name (legacy behaviour)
    try:
        from .retron_patterns import RETRONS
        for r in RETRONS:
            if r["name"] == retron_name:
                rxn_name = r.get("reaction_name", "")
                # Check if rxn_name maps to any template families
                for entry in REACTION_REGISTRY.values():
                    if entry.taxonomy_id.lower() == rxn_name.lower():
                        return entry.hte_families
                return [rxn_name] if rxn_name else []
    except Exception:
        pass

    return []


def get_hte_families_for_template(template_name: str) -> List[str]:
    """Return HTE database family names for an HTE template."""
    tid = _TEMPLATE_TO_TAXONOMY.get(template_name)
    if tid:
        entry = REACTION_REGISTRY.get(tid)
        if entry:
            return entry.hte_families
    # Fallback: direct lookup from HTE_TEMPLATES
    try:
        from .hte_templates import HTE_TEMPLATES
        for t in HTE_TEMPLATES:
            if t["name"] == template_name:
                return t.get("hte_families", [])
    except Exception:
        pass
    return []


def list_all_taxonomy_ids() -> List[str]:
    """Return all taxonomy IDs present in the registry."""
    return sorted(REACTION_REGISTRY.keys())


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_registry() -> List[str]:
    """
    Validate the registry for consistency issues.

    Checks:
    1. All retron taxonomy_id values exist in the registry.
    2. All HTE template taxonomy_id values exist in the registry.
    3. No retron or template has a missing taxonomy_id field.
    4. No registry entry is empty (no retrons AND no templates).

    Returns a list of warning strings (empty = all good).
    """
    warnings: List[str] = []

    try:
        from .retron_patterns import RETRONS
        for retron in RETRONS:
            tid = retron.get("taxonomy_id")
            name = retron.get("name", "?")
            if not tid:
                warnings.append(
                    f"Retron '{name}' is missing 'taxonomy_id' field"
                )
            elif tid not in REACTION_REGISTRY:
                warnings.append(
                    f"Retron '{name}' has taxonomy_id '{tid}' not found in registry"
                )
    except Exception as exc:
        warnings.append(f"Could not load retron_patterns: {exc}")

    try:
        from .hte_templates import HTE_TEMPLATES
        for tpl in HTE_TEMPLATES:
            tid = tpl.get("taxonomy_id")
            name = tpl.get("name", "?")
            if not tid:
                warnings.append(
                    f"HTE template '{name}' is missing 'taxonomy_id' field"
                )
            elif tid not in REACTION_REGISTRY:
                warnings.append(
                    f"HTE template '{name}' has taxonomy_id '{tid}' not found in registry"
                )
    except Exception as exc:
        warnings.append(f"Could not load hte_templates: {exc}")

    for tid, entry in REACTION_REGISTRY.items():
        if not entry.retron_names and not entry.template_names:
            warnings.append(
                f"Registry entry '{tid}' has no retrons and no templates"
            )

    if warnings:
        for w in warnings:
            logger.warning("reaction_registry: %s", w)
    else:
        logger.debug("reaction_registry: all %d entries validated OK", len(REACTION_REGISTRY))

    return warnings
