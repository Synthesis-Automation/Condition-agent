"""CLI for structural expert-rule condition recommendations."""

from __future__ import annotations

import argparse
import json
from typing import Any, Dict, List, Optional

from .rules import recommend_rule_conditions
from .rules.models import RuleRecommendationResult


_COMPONENT_GROUPS = (
    ("catalysts", "Catalyst"),
    ("ligands", "Ligand"),
    ("bases", "Base"),
    ("acids", "Acid"),
    ("condensation_agents", "Condensation agent"),
    ("oxidants", "Oxidant"),
    ("reductants", "Reductant"),
    ("additives", "Additive"),
    ("solvents", "Solvent"),
    ("other_components", "Other"),
)

_AMOUNT_UNIT_LABELS = {
    "equivalent": "equiv",
    "mol_percent": "mol%",
}


def _format_number(value: object) -> str:
    """Render numeric recipe values without insignificant trailing zeros."""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return f"{value:g}"
    return str(value)


def _format_amount(amount: object, unit: object) -> Optional[str]:
    if amount is None:
        return None
    unit_label = _AMOUNT_UNIT_LABELS.get(str(unit), str(unit or ""))
    return f"{_format_number(amount)} {unit_label}".rstrip()


def _format_component(label: str, component: Dict[str, Any]) -> str:
    name = component.get("canonical_name") or component.get("raw_identifier")
    amount = _format_amount(component.get("amount"), component.get("amount_unit"))
    suffix = f", {amount}" if amount else ""
    return f"    - {label}: {name}{suffix}"


def _format_conditions(recipe: Dict[str, Any]) -> Optional[str]:
    conditions = []
    if recipe.get("temperature_c") is not None:
        conditions.append(f"{_format_number(recipe['temperature_c'])} C")
    if recipe.get("time_h") is not None:
        conditions.append(f"{_format_number(recipe['time_h'])} h")
    if recipe.get("concentration_m") is not None:
        conditions.append(f"{_format_number(recipe['concentration_m'])} M")
    if recipe.get("atmosphere"):
        conditions.append(str(recipe["atmosphere"]))
    if not conditions:
        return None
    return "; ".join(conditions)


def _format_partner_amounts(template: Dict[str, Any]) -> Optional[str]:
    amounts = []
    for partner in template.get("partner_amounts", []):
        minimum = partner.get("minimum")
        maximum = partner.get("maximum")
        if minimum is None:
            continue
        value = _format_number(minimum)
        if maximum is not None and maximum != minimum:
            value = f"{value}-{_format_number(maximum)}"
        unit = _AMOUNT_UNIT_LABELS.get(
            str(partner.get("unit")), str(partner.get("unit") or "")
        )
        amounts.append(f"{partner.get('role', 'partner')} {value} {unit}".rstrip())
    if not amounts:
        return None
    return "; ".join(amounts)


def _format_warning(warning: str) -> str:
    code, separator, detail = warning.partition(":")
    label = code.replace("_", " ").lower().capitalize()
    return f"{label}: {detail}" if separator else label


def format_result(result: RuleRecommendationResult) -> str:
    """Return the concise, human-readable review view for a rule result."""
    lines: List[str] = [f"Reaction: {result.query_reaction_smiles}"]
    if not result.valid:
        lines.append(f"Status: invalid ({result.error or 'unknown error'})")
        if result.warnings:
            lines.append(
                "Warnings: " + "; ".join(_format_warning(w) for w in result.warnings)
            )
        return "\n".join(lines)

    lines.append(f"Transformation: {result.transformation_class or 'unknown'}")
    if not result.recommendations:
        lines.append("Recommendations: none")
    else:
        lines.append(f"Recommendations: {len(result.recommendations)}")

    for recommendation in result.recommendations:
        rule_kind = recommendation.rule_kind.replace("_", " ")
        lines.extend(
            (
                "",
                (
                    f"{recommendation.rank}. {recommendation.rule_id} "
                    f"[{recommendation.rule_status}, {rule_kind}]"
                ),
                f"   Why: {recommendation.rationale}",
            )
        )
        for variant in recommendation.compatible_variants:
            lines.append(
                f"   Recipe {variant.rank}: {variant.variant_id} "
                f"[{variant.variant_status}; compatibility "
                f"{variant.compatibility_score:.2f}]"
            )
            recipe = variant.resolved_recipe
            component_lines = []
            for field, label in _COMPONENT_GROUPS:
                component_lines.extend(
                    _format_component(label, component)
                    for component in recipe.get(field, [])
                )
            if component_lines:
                lines.append("   Components:")
                lines.extend(component_lines)
            conditions = _format_conditions(recipe)
            if conditions:
                lines.append(f"   Conditions: {conditions}")
            partner_amounts = _format_partner_amounts(
                recommendation.recipe_template
            )
            if partner_amounts:
                lines.append(f"   Stoichiometry: {partner_amounts}")

        cautions = tuple(
            dict.fromkeys(
                (
                    *recommendation.cautions,
                    *(
                        caution
                        for variant in recommendation.compatible_variants
                        for caution in variant.cautions
                    ),
                )
            )
        )
        for caution in cautions:
            lines.append(f"   Caution: {caution}")

    if result.excluded_variants:
        lines.append(
            f"\nExcluded incompatible recipe variants: {len(result.excluded_variants)}"
        )
    if result.warnings:
        lines.append(
            "Flags: " + "; ".join(_format_warning(w) for w in result.warnings)
        )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recommend expert condition templates from taxonomy facts"
    )
    parser.add_argument("reaction_smiles")
    parser.add_argument(
        "--include-draft",
        action="store_true",
        help="Include explicitly marked draft rules and templates for review",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the complete audit payload as JSON",
    )
    args = parser.parse_args()
    result = recommend_rule_conditions(
        args.reaction_smiles,
        include_draft=args.include_draft,
    )
    if args.json:
        print(json.dumps(result.to_dict(), indent=2, ensure_ascii=False))
    else:
        print(format_result(result))


if __name__ == "__main__":
    main()
