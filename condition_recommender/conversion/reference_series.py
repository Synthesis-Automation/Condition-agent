"""Reference-local condition-series identity."""

from __future__ import annotations

import hashlib
import json
from typing import Optional


def _stage_class(stages: str) -> str:
    value = str(stages or "").strip()
    if not value:
        return "unknown"
    try:
        count = int(value)
    except ValueError:
        return "invalid"
    return f"count:{count}"


def reference_condition_series_id(
    *,
    reference_id: str,
    recipe_core_id: str,
    chemistry_key: str,
    stages: str,
    steps: str = "",
    temperature_c: Optional[float] = None,
    time_h: Optional[float] = None,
) -> str:
    """Group one publication's compatible recipe core and process regime."""
    if not reference_id or not recipe_core_id or not chemistry_key:
        return ""
    payload = json.dumps(
        {
            "reference_id": reference_id,
            "recipe_core_id": recipe_core_id,
            "chemistry_key": chemistry_key,
            "stage_class": _stage_class(stages),
            "step_class": _stage_class(steps),
            "temperature_c": temperature_c,
            "time_h": time_h,
            "schema_version": "1.2",
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return "RCS1:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


__all__ = ["reference_condition_series_id"]
