"""Public registry API used by converters and recommenders."""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, Iterable, Optional, Tuple

from .models import ResolutionResult
from .resolver import ConditionRegistry


@lru_cache(maxsize=1)
def get_registry() -> ConditionRegistry:
    return ConditionRegistry()


def resolve_substance(*, cas: Optional[str] = None, name: Optional[str] = None) -> ResolutionResult:
    return get_registry().resolve(cas=cas, name=name)


def resolve_substance_id(substance_id: str) -> ResolutionResult:
    """Resolve a stable registry identity."""
    return get_registry().resolve_id(substance_id)


def resolve_condition_components(**source_fields: Iterable[str]) -> Dict[str, Tuple[ResolutionResult, ...]]:
    return {
        field: tuple(resolve_substance(cas=value) for value in values if str(value).strip())
        for field, values in source_fields.items()
    }


__all__ = [
    "get_registry",
    "resolve_condition_components",
    "resolve_substance",
    "resolve_substance_id",
]
