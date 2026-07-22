"""Strict loading and cross-package validation of condition rules."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Tuple

from condition_registry import load_recipe_template_set
from reactive_taxonomy import load_handle_patterns
from reactive_taxonomy.reaction_grammars import load_reaction_grammars

from .models import (
    ConditionRule,
    ConditionRuleSet,
    PartnerConstraint,
    RuleRecommendation,
    RuleScope,
    RuleSelection,
)


RULE_SET_PATH = (
    Path(__file__).parents[1]
    / "definitions"
    / "rule_sets"
    / "pd_sp2_cn.v1.json"
)
_ROOT_KEYS = {"schema_version", "definition_id", "selection_policy", "rules"}
_RULE_KEYS = {
    "rule_id",
    "status",
    "scope",
    "match",
    "selection",
    "recommendations",
    "rationale",
    "cautions",
    "provenance",
}
_SCOPE_KEYS = {"transformation_classes_any", "event_scopes_any"}
_MATCH_KEYS = {"partner_constraints"}
_PARTNER_KEYS = {
    "role",
    "site_type",
    "anchor_contexts_any",
    "handle_tokens_any",
    "center_tokens_any",
    "derived_families_any",
    "h_count_min",
    "h_count_max",
    "retained_contexts_all",
    "retained_contexts_allowed",
    "availability_any",
}
_SELECTION_KEYS = {"group", "tier", "priority"}
_RECOMMENDATION_KEYS = {"recipe_template_id"}
_ALLOWED_STATUSES = {"draft", "active", "retired"}
_ALLOWED_TIERS = {"specific", "fallback"}
_ALLOWED_EVENT_SCOPES = {"single_event", "multi_event", "unresolved"}


def _tuple_strings(values: Iterable[Any]) -> Tuple[str, ...]:
    return tuple(str(value) for value in values)


def _unknown_keys(value: Mapping[str, Any], allowed: set[str]) -> Tuple[str, ...]:
    return tuple(sorted(str(key) for key in set(value) - allowed))


def _partner(payload: Mapping[str, Any]) -> PartnerConstraint:
    return PartnerConstraint(
        role=str(payload.get("role") or ""),
        site_type=str(payload.get("site_type") or ""),
        anchor_contexts_any=_tuple_strings(
            payload.get("anchor_contexts_any") or ()
        ),
        handle_tokens_any=_tuple_strings(payload.get("handle_tokens_any") or ()),
        center_tokens_any=_tuple_strings(payload.get("center_tokens_any") or ()),
        derived_families_any=_tuple_strings(
            payload.get("derived_families_any") or ()
        ),
        h_count_min=int(payload["h_count_min"])
        if payload.get("h_count_min") is not None
        else None,
        h_count_max=int(payload["h_count_max"])
        if payload.get("h_count_max") is not None
        else None,
        retained_contexts_all=_tuple_strings(
            payload.get("retained_contexts_all") or ()
        ),
        retained_contexts_allowed=_tuple_strings(
            payload.get("retained_contexts_allowed") or ()
        ),
        availability_any=_tuple_strings(payload.get("availability_any") or ()),
    )


def _rule(payload: Mapping[str, Any], schema_version: str) -> ConditionRule:
    scope = payload["scope"]
    match = payload["match"]
    selection = payload["selection"]
    return ConditionRule(
        rule_id=str(payload["rule_id"]),
        status=str(payload["status"]),  # type: ignore[arg-type]
        scope=RuleScope(
            transformation_classes_any=_tuple_strings(
                scope["transformation_classes_any"]
            ),
            event_scopes_any=_tuple_strings(
                scope.get("event_scopes_any") or ("single_event",)
            ),
        ),
        partner_constraints=tuple(
            _partner(value) for value in match["partner_constraints"]
        ),
        selection=RuleSelection(
            group=str(selection["group"]),
            tier=str(selection["tier"]),  # type: ignore[arg-type]
            priority=int(selection["priority"]),
        ),
        recommendations=tuple(
            RuleRecommendation(recipe_template_id=str(value["recipe_template_id"]))
            for value in payload["recommendations"]
        ),
        rationale=str(payload["rationale"]),
        cautions=_tuple_strings(payload.get("cautions") or ()),
        provenance=dict(payload.get("provenance") or {}),
        schema_version=schema_version,
    )


def validate_condition_rule_payload(payload: Mapping[str, Any]) -> Tuple[str, ...]:
    """Validate shape, taxonomy vocabulary, and recipe-template references."""
    errors = []
    for key in _unknown_keys(payload, _ROOT_KEYS):
        errors.append(f"unknown_root_key:{key}")
    schema_version = str(payload.get("schema_version") or "")
    if schema_version != "1.1":
        errors.append(f"unsupported_schema_version:{schema_version}")
    if not str(payload.get("definition_id") or ""):
        errors.append("missing_definition_id")
    policy = payload.get("selection_policy")
    if not isinstance(policy, Mapping):
        errors.append("selection_policy_must_be_object")
    else:
        if set(policy) - {"mode", "tier_order"}:
            errors.append("selection_policy_has_unknown_keys")
        if policy.get("mode") != "first_nonempty_tier":
            errors.append("unsupported_selection_mode")
        if tuple(policy.get("tier_order") or ()) != ("specific", "fallback"):
            errors.append("tier_order_must_be_specific_then_fallback")

    grammars = load_reaction_grammars()
    grammar_roles = {}
    for grammar in grammars:
        transformation = str(grammar.get("transformation_class") or "")
        grammar_roles.setdefault(transformation, {}).update(
            {
                str(role): str(constraint.get("site_type") or "")
                for role, constraint in (grammar.get("roles") or {}).items()
            }
        )
    handle_tokens_by_site = {}
    for pattern in load_handle_patterns():
        handle_tokens_by_site.setdefault(
            str(pattern.get("site_type") or ""), set()
        ).update(str(value) for value in pattern.get("tokens") or ())
    templates = {
        template.template_id: template for template in load_recipe_template_set().templates
    }
    rules = payload.get("rules")
    if not isinstance(rules, list):
        return tuple(sorted((*errors, "rules_must_be_array")))
    rule_ids = []
    for rule_index, rule in enumerate(rules):
        prefix = f"rules[{rule_index}]"
        if not isinstance(rule, Mapping):
            errors.append(f"{prefix}:must_be_object")
            continue
        for key in _unknown_keys(rule, _RULE_KEYS):
            errors.append(f"{prefix}:unknown_key:{key}")
        rule_id = str(rule.get("rule_id") or "")
        rule_ids.append(rule_id)
        if not rule_id:
            errors.append(f"{prefix}:missing_rule_id")
        status = str(rule.get("status") or "")
        if status not in _ALLOWED_STATUSES:
            errors.append(f"{prefix}:invalid_status:{status}")
        scope = rule.get("scope")
        transformations: Tuple[str, ...] = ()
        if not isinstance(scope, Mapping):
            errors.append(f"{prefix}:scope_must_be_object")
        else:
            for key in _unknown_keys(scope, _SCOPE_KEYS):
                errors.append(f"{prefix}.scope:unknown_key:{key}")
            values = scope.get("transformation_classes_any")
            if not isinstance(values, list) or not values:
                errors.append(
                    f"{prefix}.scope:transformation_classes_must_be_nonempty_array"
                )
            else:
                transformations = _tuple_strings(values)
                for value in transformations:
                    if value not in grammar_roles:
                        errors.append(
                            f"{prefix}.scope:unknown_transformation_class:{value}"
                        )
            event_scopes = scope.get("event_scopes_any") or ("single_event",)
            if not isinstance(event_scopes, list) or not event_scopes:
                errors.append(f"{prefix}.scope:event_scopes_must_be_nonempty_array")
            else:
                for value in event_scopes:
                    if str(value) not in _ALLOWED_EVENT_SCOPES:
                        errors.append(f"{prefix}.scope:unknown_event_scope:{value}")
        match = rule.get("match")
        if not isinstance(match, Mapping):
            errors.append(f"{prefix}:match_must_be_object")
        else:
            for key in _unknown_keys(match, _MATCH_KEYS):
                errors.append(f"{prefix}.match:unknown_key:{key}")
            constraints = match.get("partner_constraints")
            if not isinstance(constraints, list) or not constraints:
                errors.append(
                    f"{prefix}.match:partner_constraints_must_be_nonempty_array"
                )
            else:
                constrained_roles = []
                for partner_index, partner in enumerate(constraints):
                    partner_prefix = (
                        f"{prefix}.match.partner_constraints[{partner_index}]"
                    )
                    if not isinstance(partner, Mapping):
                        errors.append(f"{partner_prefix}:must_be_object")
                        continue
                    for key in _unknown_keys(partner, _PARTNER_KEYS):
                        errors.append(f"{partner_prefix}:unknown_key:{key}")
                    role = str(partner.get("role") or "")
                    site_type = str(partner.get("site_type") or "")
                    constrained_roles.append(role)
                    if not role or not site_type:
                        errors.append(f"{partner_prefix}:missing_role_or_site_type")
                    for transformation in transformations:
                        expected = grammar_roles.get(transformation, {}).get(role)
                        if expected is None:
                            errors.append(
                                f"{partner_prefix}:unknown_role_for_transformation:"
                                f"{transformation}:{role}"
                            )
                        elif expected != site_type:
                            errors.append(
                                f"{partner_prefix}:site_type_mismatch:"
                                f"{transformation}:{role}:{site_type}"
                            )
                    for token in partner.get("handle_tokens_any") or ():
                        if str(token) not in handle_tokens_by_site.get(
                            site_type, set()
                        ):
                            errors.append(
                                f"{partner_prefix}:unknown_handle_token:"
                                f"{site_type}:{token}"
                            )
                    minimum = partner.get("h_count_min")
                    maximum = partner.get("h_count_max")
                    minimum_valid = minimum is None or (
                        isinstance(minimum, int) and not isinstance(minimum, bool)
                    )
                    maximum_valid = maximum is None or (
                        isinstance(maximum, int) and not isinstance(maximum, bool)
                    )
                    if not minimum_valid:
                        errors.append(f"{partner_prefix}:h_count_min_must_be_integer")
                    elif minimum is not None and minimum < 0:
                        errors.append(f"{partner_prefix}:negative_h_count_min")
                    if not maximum_valid:
                        errors.append(f"{partner_prefix}:h_count_max_must_be_integer")
                    elif maximum is not None and maximum < 0:
                        errors.append(f"{partner_prefix}:negative_h_count_max")
                    if (
                        minimum_valid
                        and maximum_valid
                        and minimum is not None
                        and maximum is not None
                        and minimum > maximum
                    ):
                        errors.append(f"{partner_prefix}:invalid_h_count_range")
                if len(constrained_roles) != len(set(constrained_roles)):
                    errors.append(f"{prefix}.match:duplicate_partner_roles")
        selection = rule.get("selection")
        if not isinstance(selection, Mapping):
            errors.append(f"{prefix}:selection_must_be_object")
        else:
            for key in _unknown_keys(selection, _SELECTION_KEYS):
                errors.append(f"{prefix}.selection:unknown_key:{key}")
            if not str(selection.get("group") or ""):
                errors.append(f"{prefix}.selection:missing_group")
            tier = str(selection.get("tier") or "")
            if tier not in _ALLOWED_TIERS:
                errors.append(f"{prefix}.selection:invalid_tier:{tier}")
            priority = selection.get("priority")
            if isinstance(priority, bool) or not isinstance(priority, int):
                errors.append(f"{prefix}.selection:priority_must_be_integer")
        recommendations = rule.get("recommendations")
        if not isinstance(recommendations, list) or not recommendations:
            errors.append(f"{prefix}:recommendations_must_be_nonempty_array")
        else:
            for recommendation_index, recommendation in enumerate(recommendations):
                recommendation_prefix = (
                    f"{prefix}.recommendations[{recommendation_index}]"
                )
                if not isinstance(recommendation, Mapping):
                    errors.append(f"{recommendation_prefix}:must_be_object")
                    continue
                for key in _unknown_keys(recommendation, _RECOMMENDATION_KEYS):
                    errors.append(f"{recommendation_prefix}:unknown_key:{key}")
                template_id = str(recommendation.get("recipe_template_id") or "")
                template = templates.get(template_id)
                if template is None:
                    errors.append(
                        f"{recommendation_prefix}:unknown_recipe_template:{template_id}"
                    )
                elif status == "active" and template.status != "active":
                    errors.append(
                        f"{recommendation_prefix}:active_rule_references_"
                        f"{template.status}_template:{template_id}"
                    )
        if not str(rule.get("rationale") or ""):
            errors.append(f"{prefix}:missing_rationale")
        if status == "active":
            provenance = rule.get("provenance")
            if not isinstance(provenance, Mapping):
                errors.append(f"{prefix}:active_rule_missing_provenance")
            else:
                if provenance.get("review_required") is not False:
                    errors.append(f"{prefix}:active_rule_must_be_reviewed")
                for field in ("doi", "procedure_locator"):
                    if not str(provenance.get(field) or ""):
                        errors.append(
                            f"{prefix}:active_rule_missing_provenance:{field}"
                        )
    if len(rule_ids) != len(set(rule_ids)):
        errors.append("duplicate_rule_ids")
    return tuple(sorted(errors))


@lru_cache(maxsize=1)
def load_condition_rule_set() -> ConditionRuleSet:
    """Load the checked rule set or fail without partial admission."""
    with RULE_SET_PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    errors = validate_condition_rule_payload(payload)
    if errors:
        raise ValueError("Invalid condition rules: " + "; ".join(errors))
    policy = payload["selection_policy"]
    schema_version = str(payload["schema_version"])
    return ConditionRuleSet(
        definition_id=str(payload["definition_id"]),
        selection_mode=str(policy["mode"]),  # type: ignore[arg-type]
        tier_order=tuple(policy["tier_order"]),  # type: ignore[arg-type]
        rules=tuple(_rule(value, schema_version) for value in payload["rules"]),
        schema_version=schema_version,
    )


__all__ = [
    "RULE_SET_PATH",
    "load_condition_rule_set",
    "validate_condition_rule_payload",
]
