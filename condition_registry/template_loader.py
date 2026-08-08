"""Load and validate package-owned condition-recipe templates."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Tuple

from .loader import load_role_definitions
from .resolver import ConditionRegistry
from .template_models import (
    ConditionRecipeTemplate,
    ConditionRecipeTemplateSet,
    RecipeTemplateOption,
    RecipeTemplatePartnerAmount,
    RecipeTemplateSelection,
    RecipeTemplateSlot,
    RecipeTemplateVariant,
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
    "variants",
    "temperature_c",
    "time_h",
    "concentration_m",
    "pressure_bar",
    "atmosphere",
    "partner_amounts",
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
_VARIANT_KEYS = {"variant_id", "status", "selections", "priority", "notes"}
_SELECTION_KEYS = {"slot_id", "substance_id", "amount", "amount_unit"}
_PARTNER_AMOUNT_KEYS = {"role", "minimum", "maximum", "unit"}
_ALLOWED_AMOUNT_UNITS = {"equivalent", "mol_percent"}


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


def _variant(payload: Mapping[str, Any]) -> RecipeTemplateVariant:
    selections = payload.get("selections") or ()
    return RecipeTemplateVariant(
        variant_id=str(payload.get("variant_id") or ""),
        status=str(payload.get("status") or ""),  # type: ignore[arg-type]
        selections=tuple(
            RecipeTemplateSelection(
                slot_id=str(value.get("slot_id") or ""),
                substance_id=str(value.get("substance_id") or ""),
                amount=float(value["amount"])
                if value.get("amount") is not None
                else None,
                amount_unit=str(value["amount_unit"])
                if value.get("amount_unit")
                else None,
            )
            for value in selections
        ),
        priority=int(payload.get("priority") or 0),
        notes=_tuple_strings(payload.get("notes") or ()),
    )


def _partner_amount(payload: Mapping[str, Any]) -> RecipeTemplatePartnerAmount:
    return RecipeTemplatePartnerAmount(
        role=str(payload.get("role") or ""),
        minimum=float(payload.get("minimum") or 0.0),
        maximum=float(payload.get("maximum") or 0.0),
        unit=str(payload.get("unit") or "equivalent"),  # type: ignore[arg-type]
    )


def _template(payload: Mapping[str, Any], schema_version: str) -> ConditionRecipeTemplate:
    return ConditionRecipeTemplate(
        template_id=str(payload.get("template_id") or ""),
        status=str(payload.get("status") or ""),  # type: ignore[arg-type]
        slots=tuple(_slot(value) for value in payload.get("slots") or ()),
        variants=tuple(_variant(value) for value in payload.get("variants") or ()),
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
        partner_amounts=tuple(
            _partner_amount(value)
            for value in payload.get("partner_amounts") or ()
        ),
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
    if schema_version != "1.2":
        errors.append(f"unsupported_schema_version:{schema_version}")
    if not str(payload.get("definition_id") or ""):
        errors.append("missing_definition_id")
    templates = payload.get("templates")
    if not isinstance(templates, list):
        return tuple(sorted((*errors, "templates_must_be_array")))

    registry = ConditionRegistry()
    definitions = load_role_definitions()
    roles = {
        str(value.get("id") or "")
        for value in definitions.get("roles") or ()
    }
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
        partner_amounts = template.get("partner_amounts") or []
        if not isinstance(partner_amounts, list):
            errors.append(f"{prefix}:partner_amounts_must_be_array")
            partner_amounts = ()
        partner_roles = []
        for amount_index, amount in enumerate(partner_amounts):
            amount_prefix = f"{prefix}.partner_amounts[{amount_index}]"
            if not isinstance(amount, Mapping):
                errors.append(f"{amount_prefix}:must_be_object")
                continue
            for key in _unknown_keys(amount, _PARTNER_AMOUNT_KEYS):
                errors.append(f"{amount_prefix}:unknown_key:{key}")
            role = str(amount.get("role") or "")
            partner_roles.append(role)
            if not role:
                errors.append(f"{amount_prefix}:missing_role")
            minimum = amount.get("minimum")
            maximum = amount.get("maximum")
            for field, value in (("minimum", minimum), ("maximum", maximum)):
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or value <= 0
                ):
                    errors.append(f"{amount_prefix}:{field}_must_be_positive_number")
            if (
                isinstance(minimum, (int, float))
                and not isinstance(minimum, bool)
                and isinstance(maximum, (int, float))
                and not isinstance(maximum, bool)
                and minimum > maximum
            ):
                errors.append(f"{amount_prefix}:invalid_range")
            if str(amount.get("unit") or "") != "equivalent":
                errors.append(f"{amount_prefix}:invalid_unit")
        if len(partner_roles) != len(set(partner_roles)):
            errors.append(f"{prefix}:duplicate_partner_amount_roles")
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
        slot_options = {
            str(slot.get("slot_id") or ""): {
                str(option.get("substance_id") or "")
                for option in slot.get("alternatives") or ()
                if isinstance(option, Mapping)
            }
            for slot in slots
            if isinstance(slot, Mapping)
        }
        required_slots = {
            str(slot.get("slot_id") or "")
            for slot in slots
            if isinstance(slot, Mapping) and bool(slot.get("required", True))
        }
        variants = template.get("variants") or ()
        if not isinstance(variants, list):
            errors.append(f"{prefix}:variants_must_be_array")
            variants = ()
        variant_ids = []
        for variant_index, variant in enumerate(variants):
            variant_prefix = f"{prefix}.variants[{variant_index}]"
            if not isinstance(variant, Mapping):
                errors.append(f"{variant_prefix}:must_be_object")
                continue
            for key in _unknown_keys(variant, _VARIANT_KEYS):
                errors.append(f"{variant_prefix}:unknown_key:{key}")
            variant_id = str(variant.get("variant_id") or "")
            variant_ids.append(variant_id)
            if not variant_id:
                errors.append(f"{variant_prefix}:missing_variant_id")
            variant_status = str(variant.get("status") or "")
            if variant_status not in _ALLOWED_STATUSES:
                errors.append(f"{variant_prefix}:invalid_status:{variant_status}")
            if variant_status == "active" and status != "active":
                errors.append(f"{variant_prefix}:active_variant_in_inactive_template")
            priority = variant.get("priority", 0)
            if isinstance(priority, bool) or not isinstance(priority, int):
                errors.append(f"{variant_prefix}:priority_must_be_integer")
            selections = variant.get("selections")
            if not isinstance(selections, list):
                errors.append(f"{variant_prefix}:selections_must_be_array")
                continue
            normalized_selections = []
            for selection_index, selection in enumerate(selections):
                selection_prefix = (
                    f"{variant_prefix}.selections[{selection_index}]"
                )
                if not isinstance(selection, Mapping):
                    errors.append(f"{selection_prefix}:must_be_object")
                    continue
                for key in _unknown_keys(selection, _SELECTION_KEYS):
                    errors.append(f"{selection_prefix}:unknown_key:{key}")
                slot_id = str(selection.get("slot_id") or "")
                substance_id = str(selection.get("substance_id") or "")
                if not slot_id:
                    errors.append(f"{selection_prefix}:missing_slot_id")
                if not substance_id:
                    errors.append(f"{selection_prefix}:missing_substance_id")
                amount = selection.get("amount")
                amount_unit = selection.get("amount_unit")
                if amount is not None and (
                    isinstance(amount, bool)
                    or not isinstance(amount, (int, float))
                    or amount <= 0
                ):
                    errors.append(f"{selection_prefix}:amount_must_be_positive_number")
                if (amount is None) != (amount_unit is None):
                    errors.append(
                        f"{selection_prefix}:amount_and_unit_must_appear_together"
                    )
                if amount_unit is not None and str(amount_unit) not in _ALLOWED_AMOUNT_UNITS:
                    errors.append(f"{selection_prefix}:invalid_amount_unit:{amount_unit}")
                normalized_selections.append((slot_id, substance_id, amount, amount_unit))
            selection_slots = {value[0] for value in normalized_selections}
            if len(selection_slots) != len(normalized_selections):
                errors.append(f"{variant_prefix}:duplicate_selection_slots")
            for missing in sorted(required_slots - selection_slots):
                errors.append(f"{variant_prefix}:missing_required_slot:{missing}")
            for unknown in sorted(selection_slots - set(slot_options)):
                errors.append(f"{variant_prefix}:unknown_slot:{unknown}")
            for normalized_slot, normalized_substance, amount, amount_unit in normalized_selections:
                if (
                    normalized_slot in slot_options
                    and normalized_substance not in slot_options[normalized_slot]
                ):
                    errors.append(
                        f"{variant_prefix}:selection_not_in_slot_alternatives:"
                        f"{normalized_slot}:{normalized_substance}"
                    )
                role_id = next(
                    (
                        str(slot.get("role_id") or "")
                        for slot in slots
                        if isinstance(slot, Mapping)
                        and str(slot.get("slot_id") or "") == normalized_slot
                    ),
                    "",
                )
                if role_id == "metal_catalyst" and amount_unit not in (None, "mol_percent"):
                    errors.append(
                        f"{variant_prefix}:catalyst_amount_must_be_mol_percent:"
                        f"{normalized_slot}"
                    )
                if role_id == "base" and amount_unit not in (None, "equivalent"):
                    errors.append(
                        f"{variant_prefix}:base_amount_must_be_equivalent:"
                        f"{normalized_slot}"
                    )
                if normalized_substance in {
                    str(value)
                    for value in template.get("forbidden_substance_ids") or ()
                }:
                    errors.append(
                        f"{variant_prefix}:forbidden_substance_selected:"
                        f"{normalized_substance}"
                    )
            if variant_status == "active":
                selected_roles = {
                    next(
                        (
                            str(slot.get("role_id") or "")
                            for slot in slots
                            if isinstance(slot, Mapping)
                            and str(slot.get("slot_id") or "") == slot_id
                        ),
                        "",
                    ): (amount, amount_unit)
                    for slot_id, _, amount, amount_unit in normalized_selections
                }
                for role_id in ("metal_catalyst", "base"):
                    amount, amount_unit = selected_roles.get(role_id, (None, None))
                    if amount is None or amount_unit is None:
                        errors.append(
                            f"{variant_prefix}:active_variant_missing_quantity:"
                            f"{role_id}"
                        )
        if len(variant_ids) != len(set(variant_ids)):
            errors.append(f"{prefix}:duplicate_variant_ids")
        if status == "active":
            for field in ("temperature_c", "time_h", "concentration_m", "atmosphere"):
                if template.get(field) is None:
                    errors.append(f"{prefix}:active_template_missing:{field}")
            if not partner_amounts:
                errors.append(f"{prefix}:active_template_missing:partner_amounts")
            provenance = template.get("provenance")
            sources = provenance.get("sources") if isinstance(provenance, Mapping) else None
            if not isinstance(sources, list) or not sources:
                errors.append(f"{prefix}:active_template_missing:provenance_sources")
            else:
                for source_index, source in enumerate(sources):
                    source_prefix = f"{prefix}.provenance.sources[{source_index}]"
                    if not isinstance(source, Mapping) or not all(
                        str(source.get(field) or "")
                        for field in ("doi", "url", "procedure_locator")
                    ):
                        errors.append(f"{source_prefix}:incomplete_primary_source")
            if isinstance(provenance, Mapping) and provenance.get("review_required") is not False:
                errors.append(f"{prefix}:active_template_must_be_reviewed")
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
