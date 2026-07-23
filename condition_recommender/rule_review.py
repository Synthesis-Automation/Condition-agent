"""Generate a flat, chemist-readable review view of expert condition rules.

The CSV produced here is a derived review artifact. Runtime recommendation
continues to load the validated rule and recipe-template definitions.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

from condition_registry import (
    ConditionRecipeTemplate,
    ResolvedConditionComponent,
    ResolvedConditionRecipe,
    load_recipe_template_set,
    materialize_recipe_variant,
)

from .rules import ConditionRule, PartnerConstraint, load_condition_rule_set


RULE_REVIEW_COLUMNS: Tuple[str, ...] = (
    "rule_id",
    "rule_status",
    "rule_kind",
    "match_summary",
    "recipe_template_id",
    "template_status",
    "recipe_variant_id",
    "variant_status",
    "variant_priority",
    "catalyst",
    "ligand",
    "base",
    "solvent",
    "additive",
    "acid",
    "condensation_agent",
    "oxidant",
    "reductant",
    "other_component",
    "electrophile_equiv_min",
    "electrophile_equiv_max",
    "nucleophile_equiv_min",
    "nucleophile_equiv_max",
    "temperature_c",
    "time_h",
    "concentration_m",
    "pressure_bar",
    "atmosphere",
    "forbidden_substance_ids",
    "selection_group",
    "selection_tier",
    "selection_priority",
    "transformation_classes",
    "event_scopes",
    "evidence_qualities",
    "reaction_scopes",
    "electrophile_site_type",
    "electrophile_anchor_contexts",
    "electrophile_handle_tokens",
    "electrophile_center_tokens",
    "electrophile_families",
    "electrophile_h_count_min",
    "electrophile_h_count_max",
    "electrophile_retained_contexts_any",
    "electrophile_retained_contexts_all",
    "electrophile_retained_contexts_allowed",
    "electrophile_availability",
    "electrophile_steric_classes",
    "electrophile_electronic_classes",
    "electrophile_ortho_substituent_count_min",
    "electrophile_ortho_substituent_count_max",
    "electrophile_alpha_branched_group_count_min",
    "electrophile_alpha_branched_group_count_max",
    "nucleophile_site_type",
    "nucleophile_anchor_contexts",
    "nucleophile_handle_tokens",
    "nucleophile_center_tokens",
    "nucleophile_families",
    "nucleophile_h_count_min",
    "nucleophile_h_count_max",
    "nucleophile_retained_contexts_any",
    "nucleophile_retained_contexts_all",
    "nucleophile_retained_contexts_allowed",
    "nucleophile_availability",
    "nucleophile_steric_classes",
    "nucleophile_electronic_classes",
    "nucleophile_ortho_substituent_count_min",
    "nucleophile_ortho_substituent_count_max",
    "nucleophile_alpha_branched_group_count_min",
    "nucleophile_alpha_branched_group_count_max",
    "rule_rationale",
    "rule_cautions",
    "template_notes",
    "variant_notes",
    "legacy_source_rule_ids",
    "doi",
    "procedure_locator",
    "rule_provenance_json",
    "template_provenance_json",
    "recipe_id",
    "rule_definition_id",
    "rule_schema_version",
    "template_definition_id",
    "template_schema_version",
    "review_decision",
    "reviewer",
    "review_date",
    "review_notes",
    "proposed_rule_change",
    "proposed_recipe_change",
)

_COMPONENT_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("catalyst", "catalysts"),
    ("ligand", "ligands"),
    ("base", "bases"),
    ("solvent", "solvents"),
    ("additive", "additives"),
    ("acid", "acids"),
    ("condensation_agent", "condensation_agents"),
    ("oxidant", "oxidants"),
    ("reductant", "reductants"),
    ("other_component", "other_components"),
)


def _join(values: Iterable[object]) -> str:
    return " | ".join(str(value) for value in values)


def _value(value: object) -> str:
    return "" if value is None else str(value)


def _json(value: Mapping[str, object]) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _provenance_text(
    rule: ConditionRule,
    template: ConditionRecipeTemplate,
    field: str,
) -> str:
    direct = rule.provenance.get(field)
    if direct:
        return str(direct)
    template_direct = template.provenance.get(field)
    if template_direct:
        return str(template_direct)
    sources = template.provenance.get("sources")
    if not isinstance(sources, list):
        return ""
    values = []
    for source in sources:
        if isinstance(source, Mapping) and source.get(field):
            values.append(str(source[field]))
    return _join(dict.fromkeys(values))


def _component(component: ResolvedConditionComponent) -> str:
    identity = component.canonical_name or component.raw_identifier
    if component.substance_id:
        identity = f"{identity} [{component.substance_id}]"
    if component.amount is not None:
        amount = f"{component.amount:g}"
        if component.amount_unit:
            amount = f"{amount} {component.amount_unit}"
        identity = f"{identity} ({amount})"
    return identity


def _components(
    recipe: ResolvedConditionRecipe, attribute: str
) -> str:
    values: Sequence[ResolvedConditionComponent] = getattr(recipe, attribute)
    return _join(_component(value) for value in values)


def _constraint_columns(
    prefix: str, constraint: Optional[PartnerConstraint]
) -> Dict[str, str]:
    fields = {
        "site_type": "",
        "anchor_contexts": "",
        "handle_tokens": "",
        "center_tokens": "",
        "families": "",
        "h_count_min": "",
        "h_count_max": "",
        "retained_contexts_any": "",
        "retained_contexts_all": "",
        "retained_contexts_allowed": "",
        "availability": "",
        "steric_classes": "",
        "electronic_classes": "",
        "ortho_substituent_count_min": "",
        "ortho_substituent_count_max": "",
        "alpha_branched_group_count_min": "",
        "alpha_branched_group_count_max": "",
    }
    if constraint is not None:
        fields.update(
            {
                "site_type": constraint.site_type,
                "anchor_contexts": _join(constraint.anchor_contexts_any),
                "handle_tokens": _join(constraint.handle_tokens_any),
                "center_tokens": _join(constraint.center_tokens_any),
                "families": _join(constraint.derived_families_any),
                "h_count_min": _value(constraint.h_count_min),
                "h_count_max": _value(constraint.h_count_max),
                "retained_contexts_any": _join(
                    constraint.retained_contexts_any
                ),
                "retained_contexts_all": _join(
                    constraint.retained_contexts_all
                ),
                "retained_contexts_allowed": _join(
                    constraint.retained_contexts_allowed
                ),
                "availability": _join(constraint.availability_any),
                "steric_classes": _join(constraint.steric_classes_any),
                "electronic_classes": _join(
                    constraint.electronic_classes_any
                ),
                "ortho_substituent_count_min": _value(
                    constraint.ortho_substituent_count_min
                ),
                "ortho_substituent_count_max": _value(
                    constraint.ortho_substituent_count_max
                ),
                "alpha_branched_group_count_min": _value(
                    constraint.alpha_branched_group_count_min
                ),
                "alpha_branched_group_count_max": _value(
                    constraint.alpha_branched_group_count_max
                ),
            }
        )
    return {f"{prefix}_{key}": value for key, value in fields.items()}


def _match_summary(rule: ConditionRule) -> str:
    constraints = {value.role: value for value in rule.partner_constraints}
    electrophile = constraints.get("electrophile")
    nucleophile = constraints.get("nucleophile")

    electrophile_parts = []
    if electrophile is not None:
        electrophile_parts.extend(electrophile.anchor_contexts_any)
        electrophile_parts.extend(electrophile.handle_tokens_any)
        electrophile_parts.extend(electrophile.steric_classes_any)

    nucleophile_parts = []
    if nucleophile is not None:
        nucleophile_parts.extend(nucleophile.derived_families_any)
        nucleophile_parts.extend(nucleophile.retained_contexts_all)
        nucleophile_parts.extend(nucleophile.retained_contexts_any)
        if nucleophile.alpha_branched_group_count_min is not None:
            nucleophile_parts.append(
                "alpha_branched>="
                f"{nucleophile.alpha_branched_group_count_min}"
            )

    left = "+".join(electrophile_parts) or "electrophile"
    right = "+".join(nucleophile_parts) or "nucleophile"
    return f"{left} + {right}"


def _partner_amount_columns(
    template: ConditionRecipeTemplate,
) -> Dict[str, str]:
    row = {
        "electrophile_equiv_min": "",
        "electrophile_equiv_max": "",
        "nucleophile_equiv_min": "",
        "nucleophile_equiv_max": "",
    }
    for amount in template.partner_amounts:
        if amount.role not in {"electrophile", "nucleophile"}:
            continue
        row[f"{amount.role}_equiv_min"] = _value(amount.minimum)
        row[f"{amount.role}_equiv_max"] = _value(amount.maximum)
    return row


def _base_row(
    rule: ConditionRule,
    template: ConditionRecipeTemplate,
    *,
    rule_definition_id: str,
    rule_schema_version: str,
    template_definition_id: str,
    template_schema_version: str,
) -> Dict[str, str]:
    constraints = {value.role: value for value in rule.partner_constraints}
    row = {column: "" for column in RULE_REVIEW_COLUMNS}
    row.update(
        {
            "rule_id": rule.rule_id,
            "rule_status": rule.status,
            "rule_kind": rule.rule_kind,
            "match_summary": _match_summary(rule),
            "recipe_template_id": template.template_id,
            "template_status": template.status,
            "temperature_c": _value(template.temperature_c),
            "time_h": _value(template.time_h),
            "concentration_m": _value(template.concentration_m),
            "pressure_bar": _value(template.pressure_bar),
            "atmosphere": template.atmosphere or "",
            "forbidden_substance_ids": _join(
                template.forbidden_substance_ids
            ),
            "selection_group": rule.selection.group,
            "selection_tier": rule.selection.tier,
            "selection_priority": str(rule.selection.priority),
            "transformation_classes": _join(
                rule.scope.transformation_classes_any
            ),
            "event_scopes": _join(rule.scope.event_scopes_any),
            "evidence_qualities": _join(
                rule.scope.evidence_qualities_any
            ),
            "reaction_scopes": _join(rule.scope.reaction_scopes_any),
            "rule_rationale": rule.rationale,
            "rule_cautions": _join(rule.cautions),
            "template_notes": _join(template.notes),
            "legacy_source_rule_ids": _join(
                rule.provenance.get("source_rule_ids", ())
            ),
            "doi": _provenance_text(rule, template, "doi"),
            "procedure_locator": _provenance_text(
                rule, template, "procedure_locator"
            ),
            "rule_provenance_json": _json(rule.provenance),
            "template_provenance_json": _json(template.provenance),
            "rule_definition_id": rule_definition_id,
            "rule_schema_version": rule_schema_version,
            "template_definition_id": template_definition_id,
            "template_schema_version": template_schema_version,
        }
    )
    row.update(_partner_amount_columns(template))
    row.update(
        _constraint_columns("electrophile", constraints.get("electrophile"))
    )
    row.update(
        _constraint_columns("nucleophile", constraints.get("nucleophile"))
    )
    return row


def build_rule_review_rows() -> Tuple[Dict[str, str], ...]:
    """Build deterministic rule-by-explicit-variant rows for human review."""
    rule_set = load_condition_rule_set()
    template_set = load_recipe_template_set()
    templates = {
        template.template_id: template for template in template_set.templates
    }
    tier_rank = {
        tier: index for index, tier in enumerate(rule_set.tier_order)
    }
    rows = []
    ordered_rules = sorted(
        rule_set.rules,
        key=lambda rule: (
            tier_rank[rule.selection.tier],
            -rule.selection.priority,
            rule.rule_id,
        ),
    )
    for rule in ordered_rules:
        for recommendation in rule.recommendations:
            template = templates[recommendation.recipe_template_id]
            base = _base_row(
                rule,
                template,
                rule_definition_id=rule_set.definition_id,
                rule_schema_version=rule_set.schema_version,
                template_definition_id=template_set.definition_id,
                template_schema_version=template_set.schema_version,
            )
            variants = sorted(
                template.variants,
                key=lambda variant: (-variant.priority, variant.variant_id),
            )
            if not variants:
                rows.append(base)
                continue
            for variant in variants:
                recipe = materialize_recipe_variant(
                    template,
                    variant.variant_id,
                    transformation_class=(
                        rule.scope.transformation_classes_any[0]
                    ),
                    include_draft=True,
                )
                row = dict(base)
                row.update(
                    {
                        "recipe_variant_id": variant.variant_id,
                        "variant_status": variant.status,
                        "variant_priority": str(variant.priority),
                        "variant_notes": _join(variant.notes),
                        "recipe_id": recipe.recipe_id,
                    }
                )
                for column, attribute in _COMPONENT_COLUMNS:
                    row[column] = _components(recipe, attribute)
                rows.append(row)
    return tuple(rows)


def export_rule_review_csv(path: Path | str) -> Tuple[Dict[str, str], ...]:
    """Write the current validated expert rules as an Excel-friendly CSV."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = build_rule_review_rows()
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RULE_REVIEW_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return rows


__all__ = [
    "RULE_REVIEW_COLUMNS",
    "build_rule_review_rows",
    "export_rule_review_csv",
]
