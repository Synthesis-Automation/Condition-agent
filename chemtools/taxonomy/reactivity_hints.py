"""Loader for reactivity hint metadata.

Provides likely bond-breaking hints for common functional groups.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional

_HINTS_FILE = "reactivity_hints.v1.json"


def _resolve_data_root(root: Optional[Path] = None) -> Path:
    if root is not None:
        return Path(root).resolve()
    return Path(__file__).resolve().parent / "data"


@lru_cache(maxsize=4)
def load_reactivity_hints(root: Optional[Path] = None) -> Dict[str, Any]:
    """Load reactivity hint metadata from taxonomy data.

    Returns a dict with keys:
    - group_hints: group_id -> hint metadata
    - reaction_overrides: reaction_id -> group_id -> override metadata
    """
    data_root = _resolve_data_root(root)
    path = data_root / _HINTS_FILE
    if not path.exists():
        return {"group_hints": {}, "reaction_overrides": {}}
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)
    except Exception:
        return {"group_hints": {}, "reaction_overrides": {}}
    if not isinstance(payload, dict):
        return {"group_hints": {}, "reaction_overrides": {}}
    payload.setdefault("group_hints", {})
    payload.setdefault("reaction_overrides", {})
    return payload


__all__ = ["load_reactivity_hints"]
