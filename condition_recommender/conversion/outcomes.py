"""Shared measurement contract for reported single-step yields."""

from __future__ import annotations

from typing import Any, Mapping


STEP_YIELD_TYPES = frozenset({
    "reported_yield_pct", "reported_product_yield_pct", "uv_area_yield_pct",
})


def has_step_yield_scope(outcome: Mapping[str, Any]) -> bool:
    """Recognize supported measurement types and honor explicit source scope."""
    metadata = outcome.get("metadata") or {}
    return (
        outcome.get("outcome_type") in STEP_YIELD_TYPES
        and metadata.get("primary_eligible", True) is True
        and metadata.get("scope") in {None, "", "step", "single_step", "product"}
    )
