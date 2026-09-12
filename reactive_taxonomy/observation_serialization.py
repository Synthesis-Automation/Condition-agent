"""Typed restoration of schema-checked structural observations for reannotation."""

from __future__ import annotations

import types
from dataclasses import fields, is_dataclass
from functools import lru_cache
from typing import Any, Literal, Mapping, Union, get_args, get_origin, get_type_hints

from .reaction_models import ReactionObservation


@lru_cache(maxsize=None)
def _hints(cls: type) -> dict[str, Any]:
    return get_type_hints(cls)


def _restore(value: Any, target: Any) -> Any:
    """Traverse only types referenced by the explicit observation contract."""
    if target is Any:
        return value
    origin, args = get_origin(target), get_args(target)
    if origin in (Union, types.UnionType):
        if value is None and type(None) in args:
            return None
        for choice in args:
            if choice is type(None):
                continue
            try:
                return _restore(value, choice)
            except (ValueError, TypeError):
                continue
        raise ValueError("Observation value does not match its typed contract")
    if origin is Literal:
        if value not in args:
            raise ValueError(f"Invalid observation vocabulary value: {value!r}")
        return value
    if origin is tuple:
        if len(args) == 2 and args[1] is Ellipsis:
            return tuple(_restore(item, args[0]) for item in value)
        if len(value) != len(args):
            raise ValueError("Invalid observation tuple length")
        return tuple(_restore(item, kind) for item, kind in zip(value, args))
    if origin is dict:
        return {
            _restore(key, args[0]): _restore(item, args[1])
            for key, item in value.items()
        }
    if origin is list:
        return [_restore(item, args[0]) for item in value]
    if is_dataclass(target):
        if not isinstance(value, Mapping):
            raise TypeError("Observation object must be a mapping")
        hints = _hints(target)
        allowed = {item.name for item in fields(target) if item.init}
        if set(value) - allowed:
            raise ValueError(f"Unexpected observation fields for {target.__name__}")
        return target(
            **{key: _restore(item, hints[key]) for key, item in value.items()}
        )
    if target is float and isinstance(value, (int, float)):
        return float(value)
    if target in (str, int, bool, float) and not isinstance(value, target):
        raise TypeError("Incorrect primitive observation type")
    return value


def restore_reaction_observation(value: Mapping[str, Any]) -> ReactionObservation:
    """Restore structural facts; executable type names never come from JSON."""
    if value.get("schema_version") != "3.0":
        raise ValueError("Unsupported stored reaction observation schema")
    return _restore(value, ReactionObservation)
