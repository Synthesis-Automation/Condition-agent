"""Load and validate package-owned condition-recipe templates."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Tuple

from .loader import load_taxonomy
from .resolver import ConditionRegistry
from .template_models import (
    ConditionRecipeTemplate,
    ConditionRecipeTemplateSet,
    RecipeTemplateOption,
    RecipeTemplateSlot,
)


RECIPE_TEMPLATES_PATH = (
    Path(__file__).with_name("definitions") / "recipe_templates.v1.json"
)
_ALLOWED_STATUSES = {"draft", "active", "retired"}
_ALLOWED_SELECTION_POLICIES = {"present_alternatives", "select_one"}
_ROOT_KEYS = {"schema_version", "definition_id", "templates"}
_TEMPLATE_KEYS = {
    "template_id",
    "status",
    "slots",
    "temperature_c",
    "time_h",
    "concentration_m",
    "pressure_bar",
    "atmosphere",
    "forbidden_substance_ids",
    "notes",
    "provenance",
}
_SLOT_KEYS = {
    "slot_id",
    "role_id",
    "required",
    "alternatives",
    "selection_policy",
    "notes",
}
_OPTION_KEYS = {"substance_id", "preference", "notes"}


def _tuple_strings(values: Iterable[Any]) -> Tuple[str, ...]:
    return tuple(str(value) for value in values)


def _unknown_keys(value: Mapping[str, Any], allowed: set[str]) -> Tuple[str, ...]:
    return tuple(sorted(str(key) for key in set(value) - allowed))


def _option(payload: Mapping[str, Any]) -> RecipeTemplateOption:
    return RecipeTemplateOption(
        substance_id=str(payload.get("substance_id") or ""),
        preference=int(payload.get("preference") or 0),
        notes=_tuple_strings(payload.get("notes") or ()),
    )


def _slot(payload: Mapping[str, Any]) -> RecipeTemplateSlot:
    return RecipeTemplateSlot(
        slot_id=str(payload.get("slot_id") or ""),
        role_id=str(payload.get("role_id") or ""),
        required=bool(payload.get("required", True)),
        alternatives=tuple(
            _option(value) for value in payload.get("alternatives") or ()
        ),
        selection_policy=str(
            payload.get("selection_policy") or "present_alternatives"
        ),  # type: ignore[arg-type]
        notes=_tuple_strings(payload.get("notes") or ()),
    )


def _template(payload: Mapping[str, Any], schema_version: str) -> ConditionRecipeTemplate:
    return ConditionRecipeTemplate(
        template_id=str(payload.get("template_id") or ""),
        status=str(payload.get("status") or ""),  # type: ignore[arg-type]
        slots=tuple(_slot(value) for value in payload.get("slots") or ()),
        temperature_c=float(payload["temperature_c"])
        if payload.get("temperature_c") is not None
        else None,
        time_h=float(payload["time_h"])
        if payload.get("time_h") is not None
        else None,
        concentration_m=float(payload["concentration_m"])
        if payload.get("concentration_m") is not None
        else None,
        pressure_bar=float(payload["pressure_bar"])
        if payload.get("pressure_bar") is not None
        else None,
        atmosphere=str(payload["atmosphere"])
        if payload.get("atmosphere")
        else None,
        forbidden_substance_ids=_tuple_strings(
            payload.get("forbidden_substance_ids") or ()
        ),
        notes=_tuple_strings(payload.get("notes") or ()),
        provenance=dict(payload.get("provenance") or {}),
        schema_version=schema_version,
    )


def validate_recipe_template_payload(payload: Mapping[str, Any]) -> Tuple[str, ...]:
    """Return deterministic errors for malformed or unresolved templates."""
    errors = []
    for key in _unknown_keys(payload, _ROOT_KEYS):
        errors.append(f"unknown_root_key:{key}")
    schema_version = str(payload.get("schema_version") or "")
    if schema_version != "1.0":
        errors.append(f"unsupported_schema_version:{schema_version}")
    if not str(payload.get("definition_id") or ""):
        errors.append("missing_definition_id")
    templates = payload.get("templates")
    if not isinstance(templates, list):
        return tuple(sorted((*errors, "templates_must_be_array")))

    registry = ConditionRegistry()
    taxonomy = load_taxonomy()
    roles = {str(value.get("id") or "") for value in taxonomy.get("roles") or ()}
    template_ids = []
    for template_index, template in enumerate(templates):
        prefix = f"templates[{template_index}]"
        if not isinstance(template, Mapping):
            errors.append(f"{prefix}:must_be_object")
            continue
        for key in _unknown_keys(template, _TEMPLATE_KEYS):
            errors.append(f"{prefix}:unknown_key:{key}")
        template_id = str(template.get("template_id") or "")
        template_ids.append(template_id)
        if not template_id:
            errors.append(f"{prefix}:missing_template_id")
        status = str(template.get("status") or "")
        if status not in _ALLOWED_STATUSES:
            errors.append(f"{prefix}:invalid_status:{status}")
        for field in (
            "temperature_c",
            "time_h",
            "concentration_m",
            "pressure_bar",
        ):
            value = template.get(field)
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, (int, float))
            ):
                errors.append(f"{prefix}:{field}_must_be_number_or_null")
        for field in ("time_h", "concentration_m", "pressure_bar"):
            value = template.get(field)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                if value < 0:
                    errors.append(f"{prefix}:{field}_must_be_nonnegative")
        slots = template.get("slots")
        if not isinstance(slots, list) or not slots:
            errors.append(f"{prefix}:slots_must_be_nonempty_array")
            continue
        slot_ids = []
        for slot_index, slot in enumerate(slots):
            slot_prefix = f"{prefix}.slots[{slot_index}]"
            if not isinstance(slot, Mapping):
                errors.append(f"{slot_prefix}:must_be_object")
                continue
            for key in _unknown_keys(slot, _SLOT_KEYS):
                errors.append(f"{slot_prefix}:unknown_key:{key}")
            slot_id = str(slot.get("slot_id") or "")
            slot_ids.append(slot_id)
            if not slot_id:
                errors.append(f"{slot_prefix}:missing_slot_id")
            role_id = str(slot.get("role_id") or "")
            if role_id not in roles:
                errors.append(f"{slot_prefix}:unknown_role:{role_id}")
            policy = str(
                slot.get("selection_policy") or "present_alternatives"
            )
            if policy not in _ALLOWED_SELECTION_POLICIES:
                errors.append(f"{slot_prefix}:invalid_selection_policy:{policy}")
            alternatives = slot.get("alternatives")
            if not isinstance(alternatives, list) or not alternatives:
                errors.append(f"{slot_prefix}:alternatives_must_be_nonempty_array")
                continue
            for option_index, option in enumerate(alternatives):
                option_prefix = f"{slot_prefix}.alternatives[{option_index}]"
                if not isinstance(option, Mapping):
                    errors.append(f"{option_prefix}:must_be_object")
                    continue
                for key in _unknown_keys(option, _OPTION_KEYS):
                    errors.append(f"{option_prefix}:unknown_key:{key}")
                substance_id = str(option.get("substance_id") or "")
                if not substance_id:
                    errors.append(f"{option_prefix}:missing_substance_id")
                    continue
                result = registry.resolve_id(substance_id)
                if result.status != "resolved" or result.substance is None:
                    errors.append(
                        f"{option_prefix}:unknown_substance_id:{substance_id}"
                    )
                elif role_id not in {
                    assignment.role_id for assignment in result.substance.roles
                }:
                    errors.append(
                        f"{option_prefix}:substance_role_mismatch:"
                        f"{substance_id}:{role_id}"
                    )
                preference = option.get("preference", 0)
                if isinstance(preference, bool) or not isinstance(preference, int):
                    errors.append(f"{option_prefix}:preference_must_be_integer")
        if len(slot_ids) != len(set(slot_ids)):
            errors.append(f"{prefix}:duplicate_slot_ids")
        for substance_id in template.get("forbidden_substance_ids") or ():
            if registry.resolve_id(str(substance_id)).status != "resolved":
                errors.append(
                    f"{prefix}:unknown_forbidden_substance_id:{substance_id}"
                )
    if len(template_ids) != len(set(template_ids)):
        errors.append("duplicate_template_ids")
    return tuple(sorted(errors))


@lru_cache(maxsize=1)
def load_recipe_template_set() -> ConditionRecipeTemplateSet:
    """Load the checked package definition or fail without partial admission."""
    with RECIPE_TEMPLATES_PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    errors = validate_recipe_template_payload(payload)
    if errors:
        raise ValueError("Invalid recipe templates: " + "; ".join(errors))
    schema_version = str(payload["schema_version"])
    return ConditionRecipeTemplateSet(
        definition_id=str(payload["definition_id"]),
        templates=tuple(
            _template(value, schema_version) for value in payload["templates"]
        ),
        schema_version=schema_version,
    )


def get_recipe_template(template_id: str) -> ConditionRecipeTemplate | None:
    """Return one template by stable ID."""
    return next(
        (
            template
            for template in load_recipe_template_set().templates
            if template.template_id == template_id
        ),
        None,
    )


__all__ = [
    "RECIPE_TEMPLATES_PATH",
    "get_recipe_template",
    "load_recipe_template_set",
    "validate_recipe_template_payload",
]
