"""Deterministic matching and specificity-tier selection for expert rules."""

from __future__ import annotations

from typing import Dict, Tuple

from .models import (
    ConditionRule,
    ConditionRuleSet,
    PartnerConstraint,
    PartnerRuleFacts,
    RuleMatchTrace,
    RuleQueryFacts,
)


def _partner_summary(facts: PartnerRuleFacts) -> str:
    values = [
        f"role={facts.role}",
        f"site_type={facts.site_type}",
        f"availability={facts.availability}",
    ]
    if facts.anchor_context:
        values.append(f"anchor_context={facts.anchor_context}")
    if facts.handle_token:
        values.append(f"handle_token={facts.handle_token}")
    if facts.center_token:
        values.append(f"center_token={facts.center_token}")
    if facts.derived_family:
        values.append(f"derived_family={facts.derived_family}")
    if facts.h_count is not None:
        values.append(f"h_count={facts.h_count}")
    if facts.retained_contexts:
        values.append("retained_contexts=" + ",".join(facts.retained_contexts))
    return ";".join(values)


def _match_partner(
    constraint: PartnerConstraint,
    facts: PartnerRuleFacts,
) -> Tuple[str, ...]:
    failures = []
    if facts.site_type != constraint.site_type:
        failures.append(
            f"{constraint.role}:site_type:{facts.site_type}!={constraint.site_type}"
        )
    if (
        constraint.anchor_contexts_any
        and facts.anchor_context not in constraint.anchor_contexts_any
    ):
        failures.append(
            f"{constraint.role}:anchor_context:{facts.anchor_context or 'missing'}"
        )
    if (
        constraint.handle_tokens_any
        and facts.handle_token not in constraint.handle_tokens_any
    ):
        failures.append(
            f"{constraint.role}:handle_token:{facts.handle_token or 'missing'}"
        )
    if (
        constraint.center_tokens_any
        and facts.center_token not in constraint.center_tokens_any
    ):
        failures.append(
            f"{constraint.role}:center_token:{facts.center_token or 'missing'}"
        )
    if (
        constraint.derived_families_any
        and facts.derived_family not in constraint.derived_families_any
    ):
        failures.append(
            f"{constraint.role}:derived_family:{facts.derived_family or 'missing'}"
        )
    if constraint.h_count_min is not None and (
        facts.h_count is None or facts.h_count < constraint.h_count_min
    ):
        failures.append(
            f"{constraint.role}:h_count_below_minimum:"
            f"{facts.h_count if facts.h_count is not None else 'missing'}"
        )
    if constraint.h_count_max is not None and (
        facts.h_count is None or facts.h_count > constraint.h_count_max
    ):
        failures.append(
            f"{constraint.role}:h_count_above_maximum:"
            f"{facts.h_count if facts.h_count is not None else 'missing'}"
        )
    contexts = set(facts.retained_contexts)
    required_contexts = set(constraint.retained_contexts_all)
    if not required_contexts.issubset(contexts):
        failures.append(
            f"{constraint.role}:missing_retained_contexts:"
            + ",".join(sorted(required_contexts - contexts))
        )
    allowed_contexts = set(constraint.retained_contexts_allowed)
    if allowed_contexts and not contexts.issubset(allowed_contexts):
        failures.append(
            f"{constraint.role}:disallowed_retained_contexts:"
            + ",".join(sorted(contexts - allowed_contexts))
        )
    if (
        constraint.availability_any
        and facts.availability not in constraint.availability_any
    ):
        failures.append(
            f"{constraint.role}:availability:{facts.availability or 'missing'}"
        )
    return tuple(failures)


def match_condition_rule(
    rule: ConditionRule,
    facts: RuleQueryFacts,
    *,
    include_draft: bool = False,
) -> RuleMatchTrace:
    """Match one rule without fuzzy inference or missing-value substitution."""
    evidence = []
    failures = []
    if facts.transformation_class not in rule.scope.transformation_classes_any:
        failures.append(
            f"transformation_class:{facts.transformation_class or 'missing'}"
        )
    else:
        evidence.append(f"transformation_class={facts.transformation_class}")
    if facts.event_scope not in rule.scope.event_scopes_any:
        failures.append(f"event_scope:{facts.event_scope or 'missing'}")
    else:
        evidence.append(f"event_scope={facts.event_scope}")

    partners: Dict[str, PartnerRuleFacts] = {
        partner.role: partner for partner in facts.partners
    }
    for constraint in rule.partner_constraints:
        partner = partners.get(constraint.role)
        if partner is None:
            failures.append(f"missing_partner:{constraint.role}")
            continue
        partner_failures = _match_partner(constraint, partner)
        if partner_failures:
            failures.extend(partner_failures)
        else:
            evidence.append(_partner_summary(partner))
    matched = not failures
    eligible = matched and (
        rule.status == "active" or (include_draft and rule.status == "draft")
    )
    return RuleMatchTrace(
        rule_id=rule.rule_id,
        rule_status=rule.status,
        matched=matched,
        eligible=eligible,
        selection_group=rule.selection.group,
        selection_tier=rule.selection.tier,
        priority=rule.selection.priority,
        evidence=tuple(evidence),
        failures=tuple(failures),
    )


def select_condition_rules(
    rule_set: ConditionRuleSet,
    facts: RuleQueryFacts,
    *,
    include_draft: bool = False,
) -> tuple[
    Tuple[ConditionRule, ...],
    Tuple[RuleMatchTrace, ...],
    Tuple[Tuple[str, str], ...],
]:
    """Select the first nonempty tier independently within each rule group."""
    traces = tuple(
        sorted(
            (
                match_condition_rule(rule, facts, include_draft=include_draft)
                for rule in rule_set.rules
            ),
            key=lambda trace: trace.rule_id,
        )
    )
    trace_by_id = {trace.rule_id: trace for trace in traces}
    eligible_by_group: Dict[str, list[ConditionRule]] = {}
    for rule in rule_set.rules:
        if trace_by_id[rule.rule_id].eligible:
            eligible_by_group.setdefault(rule.selection.group, []).append(rule)

    selected = []
    selected_tiers = []
    tier_positions = {tier: index for index, tier in enumerate(rule_set.tier_order)}
    for group in sorted(eligible_by_group):
        values = eligible_by_group[group]
        selected_tier = next(
            (
                tier
                for tier in rule_set.tier_order
                if any(rule.selection.tier == tier for rule in values)
            ),
            None,
        )
        if selected_tier is None:
            continue
        selected_tiers.append((group, selected_tier))
        selected.extend(
            rule for rule in values if rule.selection.tier == selected_tier
        )
    selected.sort(
        key=lambda rule: (
            tier_positions[rule.selection.tier],
            rule.selection.group,
            -rule.selection.priority,
            rule.rule_id,
        )
    )
    return tuple(selected), traces, tuple(selected_tiers)


__all__ = ["match_condition_rule", "select_condition_rules"]
