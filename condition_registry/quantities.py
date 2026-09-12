"""Deterministic quantity normalization and source-conflict preservation."""

from __future__ import annotations

import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

from .models import ResolvedConditionComponent


@lru_cache(maxsize=1)
def quantity_rules() -> dict[str, Any]:
    """Load the closed unit-conversion and duplicate policy definition."""
    path = Path(__file__).with_name("definitions") / "quantity_normalization.v1.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    if (
        value.get("schema_version") != "1.0"
        or value.get("duplicate_policy") != "preserve_conflicting_observations"
    ):
        raise ValueError("Unsupported condition quantity definition")
    for target, factor in value["units"].values():
        if not target or not math.isfinite(float(factor)) or float(factor) <= 0:
            raise ValueError("Invalid condition quantity conversion")
    return value


def quantity_definition_version() -> str:
    """Identify the exact rules used to resolve recipe quantities."""
    rules = quantity_rules()
    digest = hashlib.sha256(json.dumps(rules, sort_keys=True).encode()).hexdigest()[:16]
    return f"{rules['definition_version']}@sha256:{digest}"


def merge_quantities(
    components: Iterable[ResolvedConditionComponent],
) -> tuple[float | None, str | None, str, tuple[dict[str, Any], ...]]:
    """Resolve equivalent reports; sum only explicitly identified same-stage additions."""
    observations = []
    for component in components:
        if component.quantity_observations:
            observations.extend(
                dict(value) for value in component.quantity_observations
            )
            continue
        if component.amount is None:
            continue
        raw_unit = component.amount_unit
        unit = (
            str(raw_unit or "").strip().casefold().replace("µ", "u").replace("μ", "u")
        )
        normalized, factor = quantity_rules()["units"].get(unit, (unit or None, 1.0))
        observations.append(
            {
                "amount": round(float(component.amount) * factor, 12),
                "unit": normalized,
                "raw_amount": component.amount,
                "raw_unit": raw_unit,
                "source_field": component.source_field,
                "raw_identifier": component.raw_identifier,
                "provenance": dict(component.provenance),
            }
        )
    ordered = tuple(
        sorted(
            observations,
            key=lambda value: json.dumps(value, sort_keys=True, default=str),
        )
    )
    if not ordered:
        return None, None, "unreported", ()
    values = {(value["amount"], value["unit"]) for value in ordered}
    if any(
        not math.isfinite(value["amount"]) or value["amount"] < 0 for value in ordered
    ):
        return None, None, "conflicting", ordered
    provenances = [value.get("provenance") or {} for value in ordered]
    additions = [value.get("addition_id") for value in provenances]
    stages = [value.get("stage_index") for value in provenances]
    if (
        len(ordered) > 1
        and all(additions)
        and len(set(additions)) == len(additions)
        and all(stage is not None for stage in stages)
        and len(set(stages)) == 1
        and len({value["unit"] for value in ordered}) == 1
        and all(
            value.get("quantity_relationship")
            == quantity_rules()["separate_addition_provenance"]
            for value in provenances
        )
    ):
        return (
            round(sum(value["amount"] for value in ordered), 12),
            ordered[0]["unit"],
            "reported",
            ordered,
        )
    if len(values) == 1:
        amount, unit = next(iter(values))
        return amount, unit, "reported", ordered
    return None, None, "conflicting", ordered
