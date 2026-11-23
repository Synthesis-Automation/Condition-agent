"""Compare the explicitly stated components and operating conditions of two recipes."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, fields
from functools import lru_cache
import json
from math import isfinite
from pathlib import Path
from typing import Any

from .models import CONDITION_RECIPE_COMPONENT_BUCKETS as BUCKETS
from .models import ConditionComponentInput, ConditionProcessStage, ResolvedConditionRecipe
from .recipes import build_resolved_recipe_from_inputs


@lru_cache(maxsize=1)
def _rules() -> dict[str, Any]:
    path = Path(__file__).with_name("definitions") / "recipe_evidence_comparison.v1.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    if (data["definition_version"] != "recipe_evidence_comparison.v1@1.0"
            or not set(data["role_omissions"].values()) <= set(BUCKETS)
            or set(data["illumination_on"]) & set(data["illumination_off"])):
        raise ValueError("Invalid recipe evidence comparison definition")
    return data


_NUMERIC = ("temperature_c", "time_h", "concentration_m")
_OPERATING = (*_NUMERIC, "atmosphere", "illumination")


@dataclass(frozen=True)
class RecipeComparison:
    """Agreement is scoped to the conditions explicitly supplied in the query."""

    status: str
    agreements: tuple[str, ...]
    differences: tuple[str, ...]
    unknowns: tuple[str, ...]
    source_only: tuple[str, ...]
    definition_version: str


def _sequence(value: Any, key: str) -> Sequence:
    if not isinstance(value, (tuple, list)):
        raise ValueError(f"{key} must be an array")
    return value


def _mapping(value: Any, key: str) -> dict:
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be an object")
    return dict(value)


def _number(key: str, value: Any) -> None:
    if type(value) not in (int, float) or not isfinite(value):
        raise ValueError(f"{key} must be finite and numeric")
    if key in ("time_h", "concentration_m", "amount") and value < 0:
        raise ValueError(f"{key} cannot be negative")


def _validate_amount(component: Mapping[str, Any]) -> None:
    amount, unit = component.get("amount"), component.get("amount_unit")
    if amount is not None:
        _number("amount", amount)
        if amount == 0:
            raise ValueError("Zero component amount cannot establish presence; use declared_absences for an omitted role")
        if not isinstance(unit, str) or not unit.strip():
            raise ValueError("A component amount requires amount_unit")
    elif unit is not None:
        raise ValueError("A component amount_unit requires an amount")


def _normalized(key: str, value: Any) -> Any:
    if key in _NUMERIC:
        _number(key, value)
        return value
    if key == "illumination" and isinstance(value, bool):
        return "on" if value else "off"
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{key} requires a nonempty string")
    value = " ".join(value.casefold().split())
    rules = _rules()
    if value in rules["unreported_values"]:
        return None
    if key == "atmosphere":
        return rules["atmosphere_aliases"].get(value, value)
    if value in rules["illumination_off"]:
        return "off"
    if value in rules["illumination_on"]:
        return "on"
    return rules["illumination_aliases"].get(value, value)


def normalize_evidence_recipe(recipe: ResolvedConditionRecipe | Mapping[str, Any]) -> dict:
    """Accept an exported resolved recipe or resolve compact typed component inputs."""
    if isinstance(recipe, ResolvedConditionRecipe):
        recipe = recipe.to_dict()
    raw = _mapping(recipe, "recipe")
    absences = _sequence(raw.get("declared_absences", ()), "declared_absences")
    if any(not isinstance(x, str) or not x.strip() for x in absences):
        raise ValueError("declared_absences must contain nonempty strings")
    for key in _OPERATING:
        if raw.get(key) is not None:
            _normalized(key, raw[key])
    stage_rows = []
    for stage in _sequence(raw.get("stages", ()), "stages"):
        stage = _mapping(stage, "stage")
        if type(stage.get("stage_index")) is not int or stage["stage_index"] < 1:
            raise ValueError("stage_index must be a positive integer")
        if set(stage) - {"stage_index", "time_h", "temperature_c", "atmosphere", "illumination", "provenance"}:
            raise ValueError("Unrecognized stage field")
        provenance = _mapping(stage.get("provenance", {}), "stage provenance")
        if "illumination" in stage:
            provenance["illumination"] = stage.pop("illumination")
        for key in ("time_h", "temperature_c", "atmosphere", "illumination"):
            value = provenance.get(key) if key == "illumination" else stage.get(key)
            if value is not None:
                _normalized(key, value)
        stage_rows.append(ConditionProcessStage(**{**stage, "provenance": provenance}))
    if len({s.stage_index for s in stage_rows}) != len(stage_rows):
        raise ValueError("Duplicate recipe stage indices")
    if "components" in raw or (not set(raw).intersection(BUCKETS) and "recipe_id" not in raw):
        allowed = {"components", "stages", "declared_absences", *_OPERATING}
        if set(raw) - allowed:
            raise ValueError(f"Unrecognized compact recipe fields: {sorted(set(raw)-allowed)}")
        inputs = []
        for item in _sequence(raw.get("components", ()), "components"):
            item = _mapping(item, "component")
            if set(item) - {f.name for f in fields(ConditionComponentInput)}:
                raise ValueError("Unrecognized component field")
            if not isinstance(item.get("raw_identifier"), str) or not item["raw_identifier"].strip():
                raise ValueError("Each component requires raw_identifier")
            for key in ("identifier_type", "source_role_hint", "source_field"):
                if key in item and item[key] is not None and not isinstance(item[key], str):
                    raise ValueError(f"{key} must be text")
            if "provenance" in item:
                item["provenance"] = _mapping(item["provenance"], "component provenance")
            _validate_amount(item)
            item.setdefault("source_field", "recipe_query")
            inputs.append(ConditionComponentInput(**item))
        resolved = build_resolved_recipe_from_inputs(
            inputs, stages=stage_rows, declared_absences=absences,
            **{k:raw[k] for k in (*_NUMERIC, "atmosphere") if k in raw},
        ).to_dict()
        if raw.get("illumination") is not None:
            resolved["illumination"] = raw["illumination"]
    else:
        allowed = {f.name for f in fields(ResolvedConditionRecipe)} | {"illumination"}
        if set(raw) - allowed:
            raise ValueError(f"Unrecognized resolved recipe fields: {sorted(set(raw)-allowed)}")
        resolved = raw
        for bucket in BUCKETS:
            for component in _sequence(resolved.get(bucket, ()), bucket):
                component = _mapping(component, "resolved component")
                if not (component.get("substance_id") or component.get("raw_identifier")):
                    raise ValueError("A resolved component requires an identity")
                _validate_amount(component)
        if stage_rows:
            resolved["stages"] = [asdict(s) for s in stage_rows]
    facts = _facts(resolved)
    if not (facts.values or facts.uncertain):
        raise ValueError("Supply a component, declared absence, or operating condition")
    if facts.conflicts:
        raise ValueError(f"Contradictory recipe conditions: {sorted(facts.conflicts)}")
    return resolved


@dataclass
class _Facts:
    values: dict[str, Any]
    uncertain: set[str]
    conflicts: set[str]
    bucket_ids: dict[str, set[str]]


def _facts(recipe: Mapping[str, Any]) -> _Facts:
    values: dict[str, Any] = {}
    uncertain: set[str] = set()
    conflicts: set[str] = set()
    bucket_ids: dict[str, set[str]] = {k: set() for k in BUCKETS}
    for bucket in BUCKETS:
        for c in recipe.get(bucket) or ():
            identity = c.get("substance_id")
            name = c.get("canonical_name") or c.get("raw_identifier") or identity
            if c.get("amount") == 0 or (c.get("quantity_status") == "conflicting" and any(
                    v.get("amount") == 0 for v in c.get("quantity_observations") or ())):
                uncertain.add(f"{bucket}/{name} presence unconfirmed by zero dose")
                continue
            if not identity or c.get("identity_status") not in (None, "resolved"):
                uncertain.add(f"{bucket}/{name}")
                continue
            values[f"ingredient/{identity}"] = name
            if c.get("role_status") in ("ambiguous", "conflicting", "unassigned"):
                uncertain.add(f"role/{identity} ({name}, {c['role_status']})")
                if c.get("source_role_hint"):
                    values[f"reported_role_hint/{identity}"] = c["source_role_hint"]
            else:
                key = f"{bucket}/{identity}"
                values[key] = name
                bucket_ids[bucket].add(key)
            key = f"amount/{identity}"
            if c.get("quantity_status") == "conflicting":
                uncertain.add(key)
            elif c.get("amount") is not None:
                amount = c["amount"]
                if (type(amount) in (float, int) and isfinite(amount) and amount >= 0
                        and isinstance(c.get("amount_unit"), str) and c["amount_unit"].strip()):
                    values[key] = (c["amount"], c["amount_unit"])
                else:
                    uncertain.add(key)
    stages = sorted(recipe.get("stages") or (), key=lambda s: s["stage_index"])
    if stages:
        values["stage_count"] = len(stages)
    for key in _OPERATING:
        raw = recipe.get(key)
        scalar = _normalized(key, raw) if raw is not None else None
        if scalar is not None:
            values[key] = scalar
        elif raw is not None:
            uncertain.add(key)
        stage_values = []
        for s in stages:
            value = (s.get("provenance") or {}).get(key) if key == "illumination" else s.get(key)
            if value is None:
                continue
            v = _normalized(key, value)
            if v is None:
                uncertain.add(f"stage[{s['stage_index']}].{key}")
                continue
            stage_values.append(v)
            if len(stages) > 1:
                values[f"stage[{s['stage_index']}].{key}"] = v
        if len(stages) > 1 and key == "time_h":
            # A segment duration does not on its own establish the total elapsed time.
            continue
        if scalar is not None and any(v != scalar for v in stage_values):
            conflicts.add(key)
        if len(stage_values) == len(stages) and stage_values and len(set(stage_values)) == 1:
            if key not in conflicts:
                values[key] = stage_values[0]
    for absence in recipe.get("declared_absences") or ():
        absence = absence.strip().casefold()
        if absence in ("light", "illumination"):
            if any(value != "off" for key, value in values.items()
                   if key == "illumination" or key.endswith("].illumination")):
                conflicts.add("illumination")
            values["illumination"] = "off"
        else:
            bucket = _rules()["role_omissions"].get(absence, absence)
            values[f"absence/{bucket}"] = True
            if bucket_ids.get(bucket):
                conflicts.add(f"absence/{bucket}")
    return _Facts(values, uncertain, conflicts, bucket_ids)


def compare_evidence_recipes(query: Mapping[str, Any], source: Mapping[str, Any] | None) -> RecipeComparison:
    """Compare recorded conditions, preserving missing fields and explicit omissions."""
    q, s = _facts(query), _facts(source or {})
    same, different, unknown = [], [], sorted(q.uncertain)
    for key, value in sorted(q.values.items()):
        if key in q.conflicts | s.conflicts:
            unknown.append(f"{key} conflicting in record")
            continue
        if key in s.values:
            observed = s.values[key]
            if key.startswith(tuple(b + "/" for b in (*BUCKETS, "ingredient"))):
                same.append(f"{key} ({value})")  # Identified by registry ID, not display name.
            elif (key.endswith("illumination") and value == "on"
                  and observed in _rules()["illumination_sources"]):
                same.append(f"{key} = {observed}")
            elif key.endswith("illumination") and value == "on" and observed not in ("on", "off"):
                unknown.append(f"{key} source setting {observed!r} has unknown illumination status")
            elif key.endswith("illumination") and observed == "on" and value not in ("on", "off"):
                unknown.append(f"{key} source reports illumination without its type")
            elif value == observed:
                same.append(f"{key} = {observed}")
            elif key.startswith("amount/") and value[1] != observed[1]:
                unknown.append(f"{key} incomparable units {value[1]} and {observed[1]}")
            else:
                different.append(f"{key} requested {value}, reported {observed}")
            continue
        category, _, rest = key.partition("/")
        if category == "absence" and s.bucket_ids.get(rest):
            different.append(f"{key} requested, present in source")
        elif category in BUCKETS and f"absence/{category}" in s.values:
            different.append(f"{key} requested, {category} explicitly absent")
        elif category in BUCKETS and f"ingredient/{rest}" in s.values:
            unknown.append(f"{key} ingredient present, requested role unconfirmed in source")
        elif (category in BUCKETS and s.bucket_ids[category]
              and not any(v.startswith(category + "/") or v.startswith("other_components/") for v in s.uncertain)):
            different.append(f"{key} requested, source reports {', '.join(sorted(s.bucket_ids[category]))}")
        else:
            unknown.append(f"{key} unconfirmed in source")
    extra = tuple(f"{k} = {v}" for k, v in sorted(s.values.items()) if k not in q.values)
    extra += tuple(f"{k} unresolved in source" for k in sorted(s.uncertain) if k not in q.values)
    status = "different" if different else "matches_stated_fields"
    if not different and (unknown or not same):
        status = "partially_reported" if same else "unknown"
    return RecipeComparison(status, tuple(same), tuple(different), tuple(unknown), extra, _rules()["definition_version"])
