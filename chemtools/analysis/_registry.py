from __future__ import annotations

"""
Internal helpers for accessing the shared taxonomy registry.

This module mirrors the lightweight caching that previously lived in
``chemtools.router`` so that the new analysis package can reuse the same
registry instance without introducing circular imports.
"""

from typing import Optional

from ..taxonomy.archive import load_registry
from ..taxonomy.archive.registry import TaxonomyRegistry

_REGISTRY_CACHE: Optional[TaxonomyRegistry] = None
_REGISTRY_LOAD_FAILED = False


def get_registry() -> Optional[TaxonomyRegistry]:
    """Return the cached taxonomy registry, loading it on first access."""
    global _REGISTRY_CACHE, _REGISTRY_LOAD_FAILED
    if _REGISTRY_CACHE is not None:
        return _REGISTRY_CACHE
    if _REGISTRY_LOAD_FAILED:
        return None
    try:
        _REGISTRY_CACHE = load_registry()
        return _REGISTRY_CACHE
    except Exception:
        _REGISTRY_LOAD_FAILED = True
        return None


def clear_registry_cache() -> None:
    """Clear the cached registry (mainly used in tests)."""
    global _REGISTRY_CACHE, _REGISTRY_LOAD_FAILED
    _REGISTRY_CACHE = None
    _REGISTRY_LOAD_FAILED = False


__all__ = ["get_registry", "clear_registry_cache"]
