"""Generate a concise chemist-review CSV from canonical condition rules.

The CSV is a derived review artifact. Runtime recommendation continues to load
the validated rule and recipe-template definitions.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

from condition_registry import (
    ConditionRecipeTemplate,
    RecipeTemplateVariant,
    ResolvedConditionComponent,
    ResolvedConditionRecipe,
    load_recipe_template_set,
    materialize_recipe_variant,
)

from .rules import ConditionRule, PartnerConstraint, load_condition_rule_set


RULE_REVIEW_COLUMNS: Tuple[str, ...] = (
    "rule_id",
    "status",
    "rule_kind",
    "match_summary",
    "recipe_variant_id",
    "catalyst",
    "ligand",
    "base",
    "solvent",
    "other_conditions",
    "stoichiometry",
    "process_conditions",
    "rationale",
    "definition_notes",
    "source",
    "review_decision",
    "reviewer",
    "review_date",
    "review_notes",
    "proposed_change",
)

_PRIMARY_COMPONENT_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("catalyst", "catalysts"),
    ("ligand", "ligands"),
    ("base", "bases"),
    ("solvent", "solvents"),
)

_OTHER_COMPONENT_GROUPS: Tuple[Tuple[str, str], ...] = (
    ("Additive", "additives"),
    ("Acid", "acids"),
    ("Condensation agent", "condensation_agents"),
    ("Oxidant", "oxidants"),
    ("Reductant", "reductants"),
    ("Other", "other_components"),
)

_AMOUNT_UNIT_LABELS = {
    "equivalent": "equiv",
    "mol_percent": "mol%",
}


def _join(values: Iterable[object]) -> str:
    return " | ".join(str(value) for value in values if value)


def _provenance_text(
    rule: ConditionRule,
    template: ConditionRecipeTemplate,
    field: str,
) -> str:
    """Read a citation field without exposing the full provenance payload."""
    direct = rule.provenance.get(field) or template.provenance.get(field)
    if direct:
        return str(direct)
    sources = template.provenance.get("sources")
    if not isinstance(sources, list):
        return ""
    values = (
        source[field]
        for source in sources
        if isinstance(source, Mapping) and source.get(field)
    )
    return _join(dict.fromkeys(values))


def _component(component: ResolvedConditionComponent) -> str:
    identity = component.canonical_name or component.raw_identifier
    if component.amount is None:
        return identity
    unit = _AMOUNT_UNIT_LABELS.get(
        str(component.amount_unit),
        str(component.amount_unit or ""),
    )
    return f"{identity} ({component.amount:g} {unit})".rstrip()


def _components(recipe: ResolvedConditionRecipe, attribute: str) -> str:
    values: Sequence[ResolvedConditionComponent] = getattr(recipe, attribute)
    return _join(_component(value) for value in values)


def _other_conditions(recipe: ResolvedConditionRecipe) -> str:
    return _join(
        f"{label}: {components}"
        for label, attribute in _OTHER_COMPONENT_GROUPS
        if (components := _components(recipe, attribute))
    )


def _constraint_summary(constraint: PartnerConstraint) -> str:
    details = [constraint.site_type]
    details.extend(constraint.anchor_contexts_any)
    details.extend(constraint.handle_tokens_any)
    details.extend(constraint.center_tokens_any)
    details.extend(constraint.derived_families_any)
    details.extend(constraint.retained_contexts_all)
    details.extend(constraint.retained_contexts_any)
    details.extend(constraint.availability_any)
    details.extend(constraint.context_kinds_any)
    details.extend(constraint.ring_families_any)
    details.extend(constraint.steric_accessibility_any)
    details.extend(constraint.ortho_burden_classes_any)
    details.extend(constraint.electronic_axes_any)
    details.extend(constraint.alkyl_substitutions_any)
    details.extend(constraint.beta_hydrogen_statuses_any)
    details.extend(constraint.lone_pair_availability_any)
    details.extend(constraint.reactivity_modifiers_any)
    if constraint.ortho_occupancy_min is not None:
        details.append(f"ortho occupancy >= {constraint.ortho_occupancy_min}")
    if constraint.ortho_occupancy_max is not None:
        details.append(f"ortho occupancy <= {constraint.ortho_occupancy_max}")
    if constraint.alpha_branched_group_count_min is not None:
        details.append(
            f"alpha branches >= {constraint.alpha_branched_group_count_min}"
        )
    if constraint.h_count_min is not None:
        details.append(f"H >= {constraint.h_count_min}")
    if constraint.h_count_max is not None:
        details.append(f"H <= {constraint.h_count_max}")
    summary = ", ".join(dict.fromkeys(details)) or constraint.site_type
    return f"{constraint.role}: {summary}"


def _match_summary(rule: ConditionRule) -> str:
    transformation = ", ".join(rule.scope.transformation_classes_any)
    partners = "; ".join(
        _constraint_summary(constraint)
        for constraint in rule.partner_constraints
    )
    return f"{transformation}; {partners}" if partners else transformation


def _equivalent_range(
    template: ConditionRecipeTemplate,
    role: str,
) -> str:
    amount = next(
        (value for value in template.partner_amounts if value.role == role),
        None,
    )
    if amount is None:
        return ""
    minimum = f"{amount.minimum:g}"
    maximum = f"{amount.maximum:g}"
    value = minimum if minimum == maximum else f"{minimum}-{maximum}"
    return f"{value} equiv"


def _stoichiometry(template: ConditionRecipeTemplate) -> str:
    return _join(
        f"{role} {amount}"
        for role in ("electrophile", "nucleophile")
        if (amount := _equivalent_range(template, role))
    )


def _process_conditions(template: ConditionRecipeTemplate) -> str:
    return _join(
        (
            (
                f"{template.temperature_c:g} C"
                if template.temperature_c is not None
                else ""
            ),
            f"{template.time_h:g} h" if template.time_h is not None else "",
            (
                f"{template.concentration_m:g} M"
                if template.concentration_m is not None
                else ""
            ),
            (
                f"{template.pressure_bar:g} bar"
                if template.pressure_bar is not None
                else ""
            ),
            template.atmosphere or "",
        )
    )


def _source(
    rule: ConditionRule,
    template: ConditionRecipeTemplate,
) -> str:
    doi = _provenance_text(rule, template, "doi")
    procedure = _provenance_text(rule, template, "procedure_locator")
    return _join(
        (
            f"DOI: {doi}" if doi else "",
            f"Procedure: {procedure}" if procedure else "",
        )
    )


def _combined_status(
    rule: ConditionRule,
    template: ConditionRecipeTemplate,
    variant: Optional[RecipeTemplateVariant],
) -> str:
    statuses = tuple(
        dict.fromkeys(
            value
            for value in (
                rule.status,
                template.status,
                variant.status if variant is not None else "",
            )
            if value
        )
    )
    if len(statuses) == 1:
        return statuses[0]
    return _join(
        (
            f"rule={rule.status}",
            f"template={template.status}",
            f"variant={variant.status}" if variant is not None else "",
        )
    )


def _definition_notes(
    rule: ConditionRule,
    template: ConditionRecipeTemplate,
    variant: Optional[RecipeTemplateVariant],
) -> str:
    return _join(
        (
            *(f"Caution: {value}" for value in rule.cautions),
            *(f"Template: {value}" for value in template.notes),
            *(
                f"Variant: {value}"
                for value in (variant.notes if variant is not None else ())
            ),
        )
    )


def _review_row(
    rule: ConditionRule,
    template: ConditionRecipeTemplate,
    variant: Optional[RecipeTemplateVariant],
    recipe: Optional[ResolvedConditionRecipe],
) -> Dict[str, str]:
    row = {column: "" for column in RULE_REVIEW_COLUMNS}
    row.update(
        {
            "rule_id": rule.rule_id,
            "status": _combined_status(rule, template, variant),
            "rule_kind": rule.rule_kind.replace("_", " "),
            "match_summary": _match_summary(rule),
            "recipe_variant_id": variant.variant_id if variant else "",
            "stoichiometry": _stoichiometry(template),
            "process_conditions": _process_conditions(template),
            "rationale": rule.rationale,
            "definition_notes": _definition_notes(rule, template, variant),
            "source": _source(rule, template),
        }
    )
    if recipe is not None:
        for column, attribute in _PRIMARY_COMPONENT_COLUMNS:
            row[column] = _components(recipe, attribute)
        row["other_conditions"] = _other_conditions(recipe)
    return row


def build_rule_review_rows() -> Tuple[Dict[str, str], ...]:
    """Build deterministic, concise rule-by-variant rows for chemist review."""
    rule_set = load_condition_rule_set()
    template_set = load_recipe_template_set()
    templates = {
        template.template_id: template for template in template_set.templates
    }
    tier_rank = {
        tier: index for index, tier in enumerate(rule_set.tier_order)
    }
    rows = []
    for rule in sorted(
        rule_set.rules,
        key=lambda value: (
            tier_rank[value.selection.tier],
            -value.selection.priority,
            value.rule_id,
        ),
    ):
        for recommendation in rule.recommendations:
            template = templates[recommendation.recipe_template_id]
            variants = sorted(
                template.variants,
                key=lambda value: (-value.priority, value.variant_id),
            )
            if not variants:
                rows.append(_review_row(rule, template, None, None))
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
                rows.append(_review_row(rule, template, variant, recipe))
    return tuple(rows)


def export_rule_review_csv(path: Path | str) -> Tuple[Dict[str, str], ...]:
    """Write the concise validated expert-rule review as an Excel-ready CSV."""
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
