"""Declarative compatibility checks between query spectators and recipes."""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from collections import Counter
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Literal, Mapping, Tuple

from condition_registry import (
    ConditionConstraintSet,
    condition_constraint_conflicts,
    CONDITION_RECIPE_COMPONENT_BUCKETS,
    load_condition_vocabulary,
)
from reactive_taxonomy import load_molecular_motif_definitions

_RULES_PATH = Path(__file__).with_name("definitions") / "compatibility.v1.json"
_MATCH_KEYS = {
    "query_tags_any",
    "recipe_buckets_any",
    "recipe_role_ids_any",
    "recipe_substance_ids_any",
    "recipe_atmospheres_any",
    "minimum_temperature_c",
    "maximum_temperature_c",
}


@dataclass(frozen=True)
class CompatibilityCoverage:
    """Scope of executed rules, separate from admission and penalty scores.

    A supported requirement establishes only the capability named by that rule.
    It never establishes reaction success or complete experimental coverage.
    """

    schema_version: str = "compatibility_coverage.v1"
    assessment_status: Literal["not_assessed", "assessed"] = "not_assessed"
    capability_status: Literal[
        "not_assessed", "not_covered", "supported", "unresolved"
    ] = "not_assessed"
    condition_identity_status: Literal[
        "not_assessed", "no_components", "resolved", "unresolved"
    ] = "not_assessed"
    evaluated_hard_conflict_rule_ids: Tuple[str, ...] = ()
    evaluated_soft_penalty_rule_ids: Tuple[str, ...] = ()
    evaluated_regime_requirement_ids: Tuple[str, ...] = ()
    unresolved_components: Tuple[str, ...] = ()
    limitations: Tuple[str, ...] = ()
    score_meaning: str = "absence_of_known_conflicts_not_success_probability"


@dataclass(frozen=True)
class CompatibilityAssessment:
    """Hard admission decision plus auditable soft compatibility penalties."""

    compatible: bool
    score: float
    hard_conflicts: Tuple[str, ...] = ()
    penalty_ids: Tuple[str, ...] = ()
    evidence: Tuple[str, ...] = ()
    definition_id: str = "compatibility.v1"
    definition_version: str = "1.3"
    status: str = "no_known_conflict"
    checked_requirements: Tuple[str, ...] = ()
    unresolved_requirements: Tuple[str, ...] = ()
    coverage: CompatibilityCoverage = field(default_factory=CompatibilityCoverage)


@lru_cache(maxsize=1)
def load_compatibility_rules() -> Dict[str, Any]:
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    validate_compatibility_rules(rules)
    return rules


def _validate_match_vocabulary(
    rule: Mapping[str, Any],
    *,
    query_tags: set[str],
) -> None:
    unknown_tags = set(rule.get("query_tags_any") or ()) - query_tags
    if unknown_tags:
        raise ValueError(f"unknown compatibility query tags: {sorted(unknown_tags)}")
    unknown_buckets = set(rule.get("recipe_buckets_any") or ()) - set(
        CONDITION_RECIPE_COMPONENT_BUCKETS
    )
    if unknown_buckets:
        raise ValueError(
            f"unknown compatibility recipe buckets: {sorted(unknown_buckets)}"
        )
    unknown_roles = set(rule.get("recipe_role_ids_any") or ()) - set(
        load_condition_vocabulary().role_ids
    )
    if unknown_roles:
        raise ValueError(f"unknown compatibility recipe roles: {sorted(unknown_roles)}")
    for substance_id in rule.get("recipe_substance_ids_any") or ():
        from condition_registry import resolve_substance_id

        if resolve_substance_id(str(substance_id)).status != "resolved":
            raise ValueError(f"unknown compatibility recipe substance: {substance_id}")


def validate_compatibility_rules(rules: Mapping[str, Any]) -> None:
    """Validate rule structure against taxonomy and registry vocabularies."""
    if str(rules.get("schema_version") or "") != "1.3":
        raise ValueError("unsupported compatibility definition schema")
    policy = rules.get("coverage_policy") or {}
    if (
        policy.get("unknown_score") != 0.0
        or policy.get("unknown_does_not_imply_incompatible") is not True
    ):
        raise ValueError("Invalid compatibility uncertainty policy")
    for requirement in rules.get("structural_capability_requirements") or ():
        if (
            not requirement.get("id")
            or not requirement.get("order_changes")
            or not requirement.get("hydrogen_changes")
            or requirement.get("missing_policy") != "unknown_review_required"
            or set(requirement.get("recipe_buckets_any") or ())
            - set(CONDITION_RECIPE_COMPONENT_BUCKETS)
        ):
            raise ValueError("Invalid structural condition capability requirement")
    query_tags = {
        str(tag)
        for definition in load_molecular_motif_definitions()
        for tag in definition.get("tags") or ()
    }
    seen = set()
    for section in ("hard_conflicts", "soft_penalties", "regime_requirements"):
        values = rules.get(section)
        if not isinstance(values, list):
            raise ValueError(f"compatibility {section} must be a list")
        for rule in values:
            rule_id = str(rule.get("id") or "")
            if not rule_id or rule_id in seen:
                raise ValueError("compatibility rule IDs must be present and unique")
            seen.add(rule_id)
            if not str(rule.get("message") or "").strip():
                raise ValueError(f"compatibility rule {rule_id} requires a message")
            if section != "regime_requirements":
                unknown_keys = set(rule) - (
                    _MATCH_KEYS
                    | {
                        "id",
                        "message",
                        "penalty",
                        "penalty_group",
                    }
                )
                if unknown_keys:
                    raise ValueError(
                        f"unsupported keys for compatibility rule {rule_id}: "
                        f"{sorted(unknown_keys)}"
                    )
                _validate_match_vocabulary(
                    rule,
                    query_tags=query_tags,
                )
                minimum = rule.get("minimum_temperature_c")
                maximum = rule.get("maximum_temperature_c")
                if minimum is not None:
                    float(minimum)
                if maximum is not None:
                    float(maximum)
                if (
                    minimum is not None
                    and maximum is not None
                    and float(minimum) > float(maximum)
                ):
                    raise ValueError(f"invalid temperature range in {rule_id}")
            if section == "soft_penalties":
                penalty = float(rule.get("penalty") or 0.0)
                if not 0.0 < penalty < 1.0:
                    raise ValueError(
                        f"compatibility penalty must be in (0, 1): {rule_id}"
                    )
            if section == "regime_requirements":
                allowed = {
                    "id",
                    "message",
                    "named_families_any",
                    "minimum_family_confidence",
                    "required_all_buckets",
                    "unresolved_penalty",
                }
                unknown_keys = set(rule) - allowed
                if unknown_keys:
                    raise ValueError(
                        f"unsupported keys for compatibility rule {rule_id}: "
                        f"{sorted(unknown_keys)}"
                    )
                if not tuple(rule.get("named_families_any") or ()):
                    raise ValueError(
                        f"compatibility regime {rule_id} requires family evidence"
                    )
                required = set(rule.get("required_all_buckets") or ())
                if not required <= set(CONDITION_RECIPE_COMPONENT_BUCKETS):
                    raise ValueError(f"unknown required recipe bucket in {rule_id}")
                confidence = float(rule.get("minimum_family_confidence") or 0.0)
                if not 0.0 <= confidence <= 1.0:
                    raise ValueError(f"invalid family confidence in {rule_id}")
                unresolved_penalty = float(rule.get("unresolved_penalty") or 0.0)
                if not 0.0 < unresolved_penalty < 1.0:
                    raise ValueError(f"invalid unresolved penalty in {rule_id}")


def _query_tags(signature: Mapping[str, Any]) -> set[str]:
    return {
        str(tag)
        for group in signature.get("spectator_groups") or ()
        for tag in group.get("tags") or ()
    }


def _query_family_evidence(
    signature: Mapping[str, Any],
) -> tuple[tuple[str, float], ...]:
    values = []
    family = str(signature.get("named_family") or "")
    if family:
        values.append((family, float(signature.get("family_confidence") or 0.0)))
    for event in signature.get("events") or ():
        if not isinstance(event, Mapping):
            continue
        event_family = str(event.get("named_family") or "")
        if event_family:
            values.append(
                (
                    event_family,
                    float(event.get("family_confidence") or 0.0),
                )
            )
    return tuple(sorted(set(values)))


def _recipe_facts(
    recipe: Mapping[str, Any],
) -> tuple[set[str], set[str], set[str], bool, str, float | None]:
    buckets = set()
    substance_ids = set()
    role_ids = set()
    has_unresolved = False
    for bucket in CONDITION_RECIPE_COMPONENT_BUCKETS:
        components = recipe.get(bucket) or ()
        if components:
            buckets.add(bucket)
        for component in components:
            has_unresolved |= component.get("identity_status") != "resolved"
            substance_id = component.get("substance_id")
            if substance_id:
                substance_ids.add(str(substance_id))
            if component.get("role_status") == "assigned" and component.get(
                "primary_role"
            ):
                role_ids.add(str(component["primary_role"]))
    temperature_value = recipe.get("temperature_c")
    temperature = float(temperature_value) if temperature_value is not None else None
    return (
        buckets,
        substance_ids,
        role_ids,
        has_unresolved,
        str(recipe.get("atmosphere") or ""),
        temperature,
    )


def _matches(
    rule: Mapping[str, Any],
    *,
    query_tags: set[str],
    recipe_buckets: set[str],
    recipe_substance_ids: set[str],
    recipe_role_ids: set[str],
    recipe_atmosphere: str,
    recipe_temperature_c: float | None,
) -> bool:
    required_tags = set(rule.get("query_tags_any") or ())
    if required_tags and not required_tags.intersection(query_tags):
        return False
    required_buckets = set(rule.get("recipe_buckets_any") or ())
    if required_buckets and not required_buckets.intersection(recipe_buckets):
        return False
    required_substances = set(rule.get("recipe_substance_ids_any") or ())
    if required_substances and not required_substances.intersection(
        recipe_substance_ids
    ):
        return False
    required_roles = set(rule.get("recipe_role_ids_any") or ())
    if required_roles and not required_roles.intersection(recipe_role_ids):
        return False
    atmosphere_tokens = tuple(
        str(value).casefold() for value in rule.get("recipe_atmospheres_any") or ()
    )
    normalized_atmosphere = recipe_atmosphere.casefold()
    if atmosphere_tokens:
        if any(f"{token}-free" in normalized_atmosphere for token in atmosphere_tokens):
            return False
        words = set(
            normalized_atmosphere.replace(",", " ")
            .replace("/", " ")
            .replace("(", " ")
            .replace(")", " ")
            .replace("-", " ")
            .split()
        )
        if not set(atmosphere_tokens).intersection(words):
            return False
    minimum_temperature = rule.get("minimum_temperature_c")
    if minimum_temperature is not None and (
        recipe_temperature_c is None
        or recipe_temperature_c < float(minimum_temperature)
    ):
        return False
    maximum_temperature = rule.get("maximum_temperature_c")
    if maximum_temperature is not None and (
        recipe_temperature_c is None
        or recipe_temperature_c > float(maximum_temperature)
    ):
        return False
    return True


def assess_recipe_compatibility(
    signature: Mapping[str, Any], recipe: Mapping[str, Any]
) -> CompatibilityAssessment:
    """Evaluate a resolved recipe using unchanged query-group evidence."""
    rules = load_compatibility_rules()
    tags = _query_tags(signature)
    (
        buckets,
        substance_ids,
        role_ids,
        has_unresolved,
        atmosphere,
        temperature,
    ) = _recipe_facts(recipe)
    checked_requirements = []
    unresolved_requirements = []
    capability_evidence = []
    for requirement in rules.get("structural_capability_requirements") or ():
        if (
            signature.get("formed_bond_types")
            or signature.get("broken_bond_types")
            or Counter(signature.get("order_changes") or ())
            != Counter(requirement["order_changes"])
            or Counter(signature.get("hydrogen_changes") or ())
            != Counter(requirement["hydrogen_changes"])
        ):
            continue
        checked_requirements.append(requirement["id"])
        atmosphere_tokens = set(re.findall(r"[a-z0-9]+", atmosphere.casefold()))
        supported = bool(set(requirement["recipe_buckets_any"]) & buckets) or (
            bool(set(requirement["recipe_atmospheres_any"]) & atmosphere_tokens)
            and "free" not in atmosphere_tokens
        )
        if not supported:
            unresolved_requirements.append(requirement["id"])
        capability_evidence.append(
            requirement["message"]
            + (
                "; capability reported"
                if supported
                else "; supporting condition evidence is missing"
            )
        )
    hard = []
    penalties = []
    evidence = list(capability_evidence)
    evaluated_hard = []
    evaluated_soft = []
    evaluated_regime = []
    missing_regime = []
    for rule in rules.get("hard_conflicts") or ():
        evaluated_hard.append(str(rule["id"]))
        if _matches(
            rule,
            query_tags=tags,
            recipe_buckets=buckets,
            recipe_substance_ids=substance_ids,
            recipe_role_ids=role_ids,
            recipe_atmosphere=atmosphere,
            recipe_temperature_c=temperature,
        ):
            hard.append(str(rule["id"]))
            evidence.append(str(rule["message"]))
    requirement_penalties = []
    family_evidence = _query_family_evidence(signature)
    for rule in rules.get("regime_requirements") or ():
        families = set(rule.get("named_families_any") or ())
        minimum_confidence = float(rule.get("minimum_family_confidence") or 0.0)
        if families and not any(
            family in families and confidence >= minimum_confidence
            for family, confidence in family_evidence
        ):
            continue
        evaluated_regime.append(str(rule["id"]))
        required = set(rule.get("required_all_buckets") or ())
        missing = sorted(required - buckets)
        if not missing:
            continue
        missing_regime.append(str(rule["id"]))
        message = f"{rule['message']}: missing {', '.join(missing)}"
        if has_unresolved:
            requirement_penalties.append(
                (
                    rule,
                    message + "; unresolved identities prevent a hard decision",
                )
            )
        else:
            hard.append(str(rule["id"]))
            evidence.append(message)
    if hard:
        return CompatibilityAssessment(
            compatible=False,
            score=0.0,
            status="conflict",
            checked_requirements=tuple(checked_requirements),
            unresolved_requirements=tuple(unresolved_requirements),
            hard_conflicts=tuple(hard),
            evidence=tuple(evidence),
            definition_id=str(rules["definition_id"]),
            definition_version=str(rules["schema_version"]),
            coverage=_compatibility_coverage(
                recipe, evaluated_hard, evaluated_soft, evaluated_regime,
                checked_requirements, [*unresolved_requirements, *missing_regime],
            ),
        )
    matched_penalties: Dict[str, Mapping[str, Any]] = {}
    for rule in rules.get("soft_penalties") or ():
        evaluated_soft.append(str(rule["id"]))
        if _matches(
            rule,
            query_tags=tags,
            recipe_buckets=buckets,
            recipe_substance_ids=substance_ids,
            recipe_role_ids=role_ids,
            recipe_atmosphere=atmosphere,
            recipe_temperature_c=temperature,
        ):
            group = str(rule.get("penalty_group") or rule["id"])
            current = matched_penalties.get(group)
            if current is None or float(rule["penalty"]) > float(current["penalty"]):
                matched_penalties[group] = rule
    for rule, message in requirement_penalties:
        matched_penalties[f"requirement:{rule['id']}"] = {
            "id": rule["id"],
            "penalty": rule["unresolved_penalty"],
            "message": message,
        }
    total_penalty = 0.0
    for group in sorted(matched_penalties):
        rule = matched_penalties[group]
        penalties.append(str(rule["id"]))
        evidence.append(str(rule["message"]))
        total_penalty += float(rule["penalty"])
    unknown = not buckets or has_unresolved or bool(unresolved_requirements)
    if unknown:
        evidence.append(
            "Condition compatibility coverage is incomplete; no efficacy claim is supported"
        )
    return CompatibilityAssessment(
        compatible=True,
        score=float(rules["coverage_policy"]["unknown_score"])
        if unknown
        else round(max(0.0, 1.0 - total_penalty), 6),
        status="unknown" if unknown else "no_known_conflict",
        checked_requirements=tuple(checked_requirements),
        unresolved_requirements=tuple(unresolved_requirements),
        penalty_ids=tuple(penalties),
        evidence=tuple(evidence),
        definition_id=str(rules["definition_id"]),
        definition_version=str(rules["schema_version"]),
        coverage=_compatibility_coverage(
            recipe, evaluated_hard, evaluated_soft, evaluated_regime,
            checked_requirements, [*unresolved_requirements, *missing_regime],
        ),
    )


def _compatibility_coverage(
    recipe: Mapping[str, Any], hard_rules: list[str], soft_rules: list[str],
    regime_rules: list[str], capability_rules: list[str], missing_rules: list[str],
) -> CompatibilityCoverage:
    """Describe only checks that actually ran; retain the existing ranking score."""
    components = [
        (f"{bucket}[{index}]", component)
        for bucket in CONDITION_RECIPE_COMPONENT_BUCKETS
        for index, component in enumerate(recipe.get(bucket) or ())
    ]
    unresolved = tuple(
        f"{path}: {component.get('raw_identifier') or 'identity unresolved'}"
        for path, component in components if component.get("identity_status") != "resolved"
    )
    capability_status = (
        "unresolved" if missing_rules else
        "supported" if capability_rules or regime_rules else "not_covered"
    )
    limitations = [
        "Rule checks do not establish conversion, selectivity, yield "
        "or complete experimental coverage."
    ]
    if capability_status == "not_covered":
        limitations.append("No applicable reaction capability requirement was checked.")
    if unresolved:
        limitations.append("Unresolved condition identities limit the rule assessment.")
    if not components:
        limitations.append(
            "No condition components were supplied; operating settings alone "
            "do not establish a complete recipe."
        )
    return CompatibilityCoverage(
        assessment_status="assessed", capability_status=capability_status,
        condition_identity_status=(
            "unresolved" if unresolved else "resolved" if components else "no_components"
        ),
        evaluated_hard_conflict_rule_ids=tuple(hard_rules),
        evaluated_soft_penalty_rule_ids=tuple(soft_rules),
        evaluated_regime_requirement_ids=tuple(regime_rules),
        unresolved_components=unresolved, limitations=tuple(limitations),
        score_meaning=str(load_compatibility_rules()["coverage_policy"]["score_meaning"]),
    )


def filter_compatible_precedents(
    signature: Mapping[str, Any],
    precedents: Iterable[Any],
    *,
    condition_constraints: ConditionConstraintSet | None = None,
) -> tuple[
    tuple[tuple[Any, CompatibilityAssessment], ...],
    tuple[tuple[Any, CompatibilityAssessment], ...],
]:
    """Partition indexed precedents without discarding conflict evidence."""
    accepted = []
    excluded = []
    for precedent in precedents:
        assessment = assess_recipe_compatibility(signature, precedent.resolved_recipe)
        conflicts = (
            condition_constraint_conflicts(
                precedent.resolved_recipe, condition_constraints
            )
            if condition_constraints is not None
            else ()
        )
        if conflicts:
            assessment = replace(
                assessment,
                status="conflict",
                compatible=False,
                score=0.0,
                hard_conflicts=assessment.hard_conflicts + conflicts,
                evidence=assessment.evidence + conflicts,
            )
        target = accepted if assessment.compatible else excluded
        target.append((precedent, assessment))
    return tuple(accepted), tuple(excluded)


__all__ = [
    "CompatibilityAssessment",
    "CompatibilityCoverage",
    "assess_recipe_compatibility",
    "filter_compatible_precedents",
    "load_compatibility_rules",
    "validate_compatibility_rules",
]
