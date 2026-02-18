"""
Shared helper utilities for all ChemCoworker tool wrappers.

Follows the same pattern as reaction_agent/reasoning_tools.py and
chem_assistant/chemtools_wrapper.py for consistency.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional


def _to_jsonable(value: Any) -> Any:
    """Recursively convert to JSON-serializable types."""
    from dataclasses import asdict, is_dataclass
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, dict):
        return {k: _to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_to_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        return _to_jsonable(to_dict())
    to_list = getattr(value, "tolist", None)
    if callable(to_list):
        return _to_jsonable(to_list())
    return value


def _success(data: Any) -> Dict[str, Any]:
    """Wrap a successful result."""
    payload: Dict[str, Any] = {"success": True}
    if isinstance(data, dict):
        payload.update(_to_jsonable(data))
    else:
        payload["data"] = _to_jsonable(data)
    return payload


def _error(message: str, extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Wrap an error result."""
    payload: Dict[str, Any] = {"success": False, "error": message}
    if extra:
        payload.update(_to_jsonable(extra))
    return payload
