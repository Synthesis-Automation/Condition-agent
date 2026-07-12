"""Central SMARTS compilation cache owned by the standalone package."""

from __future__ import annotations

from functools import lru_cache
from typing import Any, Dict, Optional

from .rdkit_utils import _chem_module


@lru_cache(maxsize=1024)
def compile_smarts(smarts: str, *, validate: bool = False) -> Optional[Any]:
    """Compile and cache one SMARTS expression."""
    chem = _chem_module()
    if chem is None:
        if validate:
            raise ValueError("RDKit is unavailable; cannot compile SMARTS")
        return None
    try:
        pattern = chem.MolFromSmarts(str(smarts or ""))
    except Exception as exc:
        if validate:
            raise ValueError(f"Invalid SMARTS pattern: {smarts}") from exc
        return None
    if pattern is None and validate:
        raise ValueError(f"Invalid SMARTS pattern: {smarts}")
    return pattern


def compile_smarts_batch(
    patterns: Dict[str, str], *, validate: bool = False
) -> Dict[str, Optional[Any]]:
    """Compile a named collection through the shared cache."""
    return {
        name: compile_smarts(smarts, validate=validate)
        for name, smarts in patterns.items()
    }


def clear_smarts_cache() -> None:
    """Clear cached SMARTS, primarily for taxonomy validation tests."""
    compile_smarts.cache_clear()


__all__ = ["clear_smarts_cache", "compile_smarts", "compile_smarts_batch"]
