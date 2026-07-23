"""Public expert-rule recommendation API."""

from __future__ import annotations

from dataclasses import asdict, replace

from condition_registry import (
    get_recipe_template,
    load_recipe_template_set,
    materialize_recipe_variant,
)
from reactive_taxonomy import featurize_reaction

from ..compatibility import assess_recipe_compatibility
from .facts import build_rule_query_facts
from .loader import load_condition_rule_set
from .matcher import select_condition_rules
from .models import (
    RuleConditionRecommendation,
    RuleRecipeVariantAssessment,
    RuleRecommendationResult,
)


def recommend_rule_conditions(
    reaction_smiles: str,
    *,
    include_draft: bool = False,
) -> RuleRecommendationResult:
    """Recommend expert recipe templates for one verified reaction event."""
    analysis = featurize_reaction(reaction_smiles)
    facts, error = build_rule_query_facts(analysis)
    rule_set = load_condition_rule_set()
    template_set = load_recipe_template_set()
    if facts is None:
        return RuleRecommendationResult(
            query_reaction_smiles=reaction_smiles,
            valid=False,
            warnings=tuple(analysis.warnings),
            error=error,
            rule_definition_id=rule_set.definition_id,
            rule_definition_schema_version=rule_set.schema_version,
            recipe_template_definition_id=template_set.definition_id,
            recipe_template_schema_version=template_set.schema_version,
        )

    selected, traces, selected_tiers = select_condition_rules(
        rule_set,
        facts,
        include_draft=include_draft,
    )
    trace_by_id = {trace.rule_id: trace for trace in traces}
    recommendations = []
    excluded_variants = []
    seen_templates = set()
    warnings = list(analysis.warnings)
    signature_payload = asdict(analysis.reaction_signature)
    for rule in selected:
        trace = trace_by_id[rule.rule_id]
        for recommendation in rule.recommendations:
            template_id = recommendation.recipe_template_id
            if template_id in seen_templates:
                continue
            template = get_recipe_template(template_id)
            if template is None:
                continue
            seen_templates.add(template_id)
            compatible_variants = []
            eligible_variants = tuple(
                sorted(
                    (
                        variant
                        for variant in template.variants
                        if variant.status == "active"
                        or (include_draft and variant.status == "draft")
                    ),
                    key=lambda variant: (-variant.priority, variant.variant_id),
                )
            )
            if not eligible_variants:
                warnings.append(
                    f"NO_ELIGIBLE_RECIPE_VARIANTS:{template.template_id}"
                )
            for variant in eligible_variants:
                recipe = materialize_recipe_variant(
                    template,
                    variant.variant_id,
                    transformation_class=facts.transformation_class,
                    named_family=analysis.named_family,
                    include_draft=include_draft,
                )
                assessment = assess_recipe_compatibility(
                    signature_payload,
                    recipe.to_dict(),
                )
                variant_cautions = list(assessment.evidence)
                if variant.status == "draft":
                    variant_cautions.append(
                        "Draft recipe variant; quantities and operating details "
                        "require review"
                    )
                variant_result = RuleRecipeVariantAssessment(
                    rule_id=rule.rule_id,
                    recipe_template_id=template.template_id,
                    variant_id=variant.variant_id,
                    variant_status=variant.status,
                    rank=0,
                    definition_priority=variant.priority,
                    recipe_id=recipe.recipe_id,
                    resolved_recipe=recipe.to_dict(),
                    compatible=assessment.compatible,
                    compatibility_score=assessment.score,
                    hard_conflicts=assessment.hard_conflicts,
                    penalty_ids=assessment.penalty_ids,
                    compatibility_evidence=assessment.evidence,
                    cautions=tuple(variant_cautions),
                    compatibility_definition_id=assessment.definition_id,
                )
                if assessment.compatible:
                    compatible_variants.append(variant_result)
                else:
                    excluded_variants.append(variant_result)
            compatible_variants.sort(
                key=lambda result: (
                    -result.compatibility_score,
                    -result.definition_priority,
                    result.variant_id,
                )
            )
            compatible_variants = [
                replace(result, rank=rank)
                for rank, result in enumerate(compatible_variants, start=1)
            ]
            if not compatible_variants:
                continue
            cautions = list(rule.cautions)
            if rule.status == "draft" or template.status == "draft":
                cautions.append(
                    "Draft expert definition; not active for production recommendation"
                )
            recommendations.append(
                RuleConditionRecommendation(
                    rank=len(recommendations) + 1,
                    rule_id=rule.rule_id,
                    rule_status=rule.status,
                    rule_kind=rule.rule_kind,
                    recipe_template_id=template_id,
                    recipe_template=template.to_dict(),
                    compatible_variants=tuple(compatible_variants),
                    selection_group=rule.selection.group,
                    selection_tier=rule.selection.tier,
                    rationale=rule.rationale,
                    match_evidence=trace.evidence,
                    cautions=tuple(cautions),
                )
            )

    structural_matches = tuple(trace for trace in traces if trace.matched)
    if not structural_matches:
        warnings.append("NO_STRUCTURAL_RULE_MATCH")
    elif not include_draft and not recommendations and any(
        trace.rule_status == "draft" for trace in structural_matches
    ):
        warnings.append("DRAFT_RULE_MATCHES_EXCLUDED")
    if include_draft and any(
        recommendation.rule_status == "draft"
        for recommendation in recommendations
    ):
        warnings.append("DRAFT_RULES_INCLUDED")
    if any(
        recommendation.rule_kind == "structural_override"
        for recommendation in recommendations
    ):
        warnings.append("STRUCTURAL_OVERRIDE_APPLIED")
    elif any(
        recommendation.rule_kind == "default"
        for recommendation in recommendations
    ):
        warnings.append("DEFAULT_RULE_APPLIED")
    elif any(
        recommendation.rule_kind == "fallback"
        for recommendation in recommendations
    ):
        warnings.append("FALLBACK_RULE_APPLIED")
    if excluded_variants:
        warnings.append(
            f"INCOMPATIBLE_RULE_VARIANTS_EXCLUDED:{len(excluded_variants)}"
        )
    return RuleRecommendationResult(
        query_reaction_smiles=reaction_smiles,
        valid=True,
        query_signature_id=facts.signature_id,
        reaction_signature_schema_version=facts.reaction_signature_schema_version,
        transformation_class=facts.transformation_class,
        taxonomy_definition_versions=facts.taxonomy_definition_versions,
        selected_tiers=selected_tiers,
        recommendations=tuple(recommendations),
        excluded_variants=tuple(excluded_variants),
        match_traces=traces,
        warnings=tuple(sorted(set(warnings))),
        rule_definition_id=rule_set.definition_id,
        rule_definition_schema_version=rule_set.schema_version,
        recipe_template_definition_id=template_set.definition_id,
        recipe_template_schema_version=template_set.schema_version,
    )


__all__ = ["recommend_rule_conditions"]
