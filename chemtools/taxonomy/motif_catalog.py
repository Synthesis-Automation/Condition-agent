"""Canonical motif catalog accessors."""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, Iterable, List, Optional

from .compound_catalog import build_documented_compound_catalog
from .models import MotifDefinition


def _entry_id(entry: dict) -> str:
    explicit = str(entry.get("id") or entry.get("compound_id") or "").strip()
    if explicit:
        return explicit
    group_a = str(entry.get("A") or "").strip()
    group_b = str(entry.get("B") or "").strip()
    if group_a and group_b:
        return f"{group_a}{'' if group_b.startswith('-') else '-'}{group_b}"
    return ""


@lru_cache(maxsize=1)
def load_motif_catalog() -> Dict[str, MotifDefinition]:
    """Load documented motifs keyed by canonical motif id."""
    payload = build_documented_compound_catalog()
    entries = payload.get("compounds", []) if isinstance(payload, dict) else []
    out: Dict[str, MotifDefinition] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        motif_id = _entry_id(entry)
        if motif_id:
            out[motif_id] = MotifDefinition(motif_id=motif_id, payload=entry)
    return out


def get_motif(motif_id: str) -> Optional[MotifDefinition]:
    """Return a motif definition by canonical id."""
    return load_motif_catalog().get(str(motif_id or "").strip())


def list_motif_ids() -> List[str]:
    """Return all known motif ids."""
    return sorted(load_motif_catalog())


def known_motif_ids(motif_ids: Iterable[str]) -> List[str]:
    """Filter an iterable down to known motif ids."""
    catalog = load_motif_catalog()
    return [item for item in motif_ids if str(item).strip() in catalog]


__all__ = ["get_motif", "known_motif_ids", "list_motif_ids", "load_motif_catalog"]
