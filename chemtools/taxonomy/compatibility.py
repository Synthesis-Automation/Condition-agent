"""Taxonomy compatibility checks used by reaction and recommendation domains."""

from __future__ import annotations

from typing import Iterable, Set

from .scope import motif_scope_children, motif_scope_parents


def motif_tokens_compatible(left: str, right: str) -> bool:
    """Return whether two motif tokens are equal or in a parent/child scope relation."""
    left_text = str(left or "").strip()
    right_text = str(right or "").strip()
    if not left_text or not right_text:
        return False
    if left_text == right_text:
        return True

    children = motif_scope_children()
    if right_text in children.get(left_text, set()):
        return True
    if left_text in children.get(right_text, set()):
        return True

    parents = motif_scope_parents()
    if right_text in parents.get(left_text, set()):
        return True
    if left_text in parents.get(right_text, set()):
        return True
    return False


def motif_sets_compatible(required: Iterable[str], available: Iterable[str]) -> bool:
    """Return True when every required motif has a compatible available motif."""
    required_set: Set[str] = {str(item).strip() for item in required if str(item).strip()}
    available_set: Set[str] = {str(item).strip() for item in available if str(item).strip()}
    if not required_set:
        return True
    if not available_set:
        return False
    return all(
        any(motif_tokens_compatible(req, candidate) for candidate in available_set)
        for req in required_set
    )


__all__ = ["motif_sets_compatible", "motif_tokens_compatible"]
