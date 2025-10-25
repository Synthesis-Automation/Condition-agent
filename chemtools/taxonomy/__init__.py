"""
Unified taxonomy registry loader.

This package provides the public API for accessing the centralised taxonomy
definitions that unify reaction, reactant, and reagent metadata.

The registry lazily loads JSON assets from ``chemtools/taxonomy/data`` (or a
custom path when provided) and exposes typed lookup helpers. The conversion
script in ``scripts/taxonomy/convert_legacy_taxonomies.py`` can be used to
regenerate the data files from the legacy sources.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .registry import TaxonomyRegistry

__all__ = ["TaxonomyRegistry", "load_registry", "reset_registry"]

_CACHE: dict[Path, TaxonomyRegistry] = {}


def load_registry(root: Optional[Path] = None) -> TaxonomyRegistry:
    """
    Load (or retrieve from cache) the taxonomy registry.

    Args:
        root: Optional path to the directory containing taxonomy data files.
              Defaults to ``chemtools/taxonomy/data``.

    Returns:
        TaxonomyRegistry: Shared registry instance rooted at ``root``.

    Raises:
        FileNotFoundError: If the manifest or required JSON files are missing.
        ValueError: If the manifest fails validation.
    """
    target = root or Path(__file__).resolve().parent / "data"
    if target in _CACHE:
        return _CACHE[target]
    registry = TaxonomyRegistry.from_path(target)
    _CACHE[target] = registry
    return registry


def reset_registry(root: Optional[Path] = None) -> None:
    """
    Remove a cached registry instance.

    Primarily useful for tests that need to reload the taxonomy after mutating
    the underlying files in a temporary directory.
    """
    target = root or Path(__file__).resolve().parent / "data"
    _CACHE.pop(target, None)
