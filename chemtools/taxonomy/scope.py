"""Scope expansion services for taxonomy motifs."""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, Iterable, List, Set

from . import loader as taxonomy_loader


@lru_cache(maxsize=1)
def motif_scope_children() -> Dict[str, Set[str]]:
    """Return parent motif -> transitive child motifs."""
    payload = taxonomy_loader.load_motif_scope_index()
    raw_scope = payload.get("scope_map", {}) if isinstance(payload, dict) else {}
    if not isinstance(raw_scope, dict):
        return {}

    direct: Dict[str, Set[str]] = {}
    for key, value in raw_scope.items():
        parent = str(key).strip()
        if not parent:
            continue
        children = {
            str(item).strip()
            for item in (value or [])
            if str(item).strip() and str(item).strip() != parent
        }
        if children:
            direct[parent] = children

    transitive: Dict[str, Set[str]] = {}
    for parent in direct:
        seen: Set[str] = set()
        stack = list(direct.get(parent, set()))
        while stack:
            child = stack.pop()
            if child in seen:
                continue
            seen.add(child)
            stack.extend(direct.get(child, set()))
        if seen:
            transitive[parent] = seen
    return transitive


@lru_cache(maxsize=1)
def motif_scope_parents() -> Dict[str, Set[str]]:
    """Return child motif -> transitive parent motifs."""
    parents: Dict[str, Set[str]] = {}
    for parent, children in motif_scope_children().items():
        for child in children:
            parents.setdefault(child, set()).add(parent)
    return parents


def expand_motif_scope(motif_ids: Iterable[str], *, include_self: bool = True) -> List[str]:
    """Expand motifs to include known taxonomy scope descendants."""
    expanded: Set[str] = set()
    children = motif_scope_children()
    for motif_id in motif_ids:
        token = str(motif_id or "").strip()
        if not token:
            continue
        if include_self:
            expanded.add(token)
        expanded.update(children.get(token, set()))
    return sorted(expanded)


__all__ = ["expand_motif_scope", "motif_scope_children", "motif_scope_parents"]
