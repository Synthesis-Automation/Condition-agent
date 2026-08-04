"""Tiered eligibility for reaction-core queries and precedents.

This module deliberately separates a drawable/queryable reaction core from a
trusted condition precedent.  Review-core precedents are available only to an
explicit review/unrestricted index and never enter the default trusted index.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Mapping, Tuple

from reactive_taxonomy import (
    REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
    REACTION_CORE_PROJECTION_SCHEMA_VERSION,
)

from .models import CORE_ELIGIBILITY_DEFINITION_VERSION, CoreEligibility

_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "core_eligibility.v1.json"
)


@lru_cache(maxsize=1)
def load_core_eligibility_rules() -> dict[str, Any]:
    """Load and validate the declarative core-eligibility policy."""
    with _DEFINITION_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if rules.get("schema_version") != "1.0":
        raise ValueError("unsupported core-eligibility schema")
    if rules.get("definition_id") != "core_eligibility.v1":
        raise ValueError("unexpected core-eligibility definition ID")
    for field in (
        "allowed_query_evidence_statuses",
        "indexable_condition_statuses",
        "indexable_stage_statuses",
        "review_mapping_statuses",
        "blocking_warnings",
    ):
        if not tuple(rules.get(field) or ()):
            raise ValueError(f"core-eligibility policy requires {field}")
    return rules


@dataclass(frozen=True)
class CoreEligibilityAssessment:
    """One auditable core-eligibility decision."""

    tier: CoreEligibility
    reasons: Tuple[str, ...]
    query_eligible: bool
    precedent_eligible: bool
    requires_expert_review: bool
    definition_version: str = CORE_ELIGIBILITY_DEFINITION_VERSION
    schema_version: str = "1.0"


def _text(value: Any) -> str:
    return str(getattr(value, "value", value) or "")


def assess_core_eligibility(
    record: Mapping[str, Any],
) -> CoreEligibilityAssessment:
    """Classify one converted record without changing admission authority."""
    rules = load_core_eligibility_rules()
    core = record.get("reaction_core")
    if not isinstance(core, Mapping) or not core:
        return CoreEligibilityAssessment(
            CoreEligibility.UNAVAILABLE,
            ("reaction_core_unavailable",),
            False,
            False,
            False,
        )

    reasons = []
    if _text(core.get("schema_version")) != REACTION_CORE_PROJECTION_SCHEMA_VERSION:
        reasons.append("incompatible_reaction_core_schema")
    if (
        _text(core.get("algorithm_version"))
        != REACTION_CORE_PROJECTION_ALGORITHM_VERSION
    ):
        reasons.append("incompatible_reaction_core_algorithm")
    quality = core.get("quality")
    quality_status = (
        _text(quality.get("status")) if isinstance(quality, Mapping) else ""
    )
    if quality_status == "blocked":
        reasons.append("reaction_core_quality_blocked")
    warnings = {
        _text(value)
        for value in tuple(core.get("warnings") or ())
        + tuple(record.get("warnings") or ())
    }
    if warnings.intersection(set(rules["blocking_warnings"])):
        reasons.append("reaction_core_has_blocking_warning")
    evidence_status = _text(core.get("evidence_status"))
    if evidence_status not in set(rules["allowed_query_evidence_statuses"]):
        reasons.append("reaction_core_evidence_not_query_eligible")
    if int(core.get("event_count") or 0) < 1:
        reasons.append("reaction_core_event_missing")
    for field in ("exact_core_key", "typed_core_key", "shape_core_key"):
        if not _text(core.get(field)):
            reasons.append(f"{field}_missing")
    if reasons:
        return CoreEligibilityAssessment(
            CoreEligibility.BLOCKED,
            tuple(sorted(set(reasons))),
            False,
            False,
            True,
        )

    index_eligibility = _text(record.get("index_eligibility"))
    chemistry_status = _text(record.get("chemistry_status"))
    completeness = record.get("reaction_completeness")
    completeness_status = (
        _text(completeness.get("status"))
        if isinstance(completeness, Mapping)
        else ""
    )
    condition_status = _text(record.get("condition_status"))
    stage_status = _text(record.get("condition_stage_status")) or "single_stage"
    external_mapping = record.get("external_atom_mapping")
    mapping_status = (
        _text(external_mapping.get("status"))
        if isinstance(external_mapping, Mapping)
        else ""
    )
    mapping_product_coverage = (
        float(external_mapping.get("product_mapping_coverage") or 0.0)
        if isinstance(external_mapping, Mapping)
        else 0.0
    )

    external_mapping_used = bool(
        mapping_status and not mapping_status.startswith("not_requested_")
    )
    if index_eligibility == "eligible" and not external_mapping_used:
        return CoreEligibilityAssessment(
            CoreEligibility.TRUSTED_CORE, (), True, True, False
        )

    review_reasons = []
    unresolved_policy = dict(
        rules.get("mapper_supported_unresolved_requirements") or {}
    )
    mapper_supported_unresolved = bool(
        completeness_status == "unresolved"
        and mapping_status == unresolved_policy.get("mapping_status")
        and isinstance(completeness, Mapping)
        and max(
            float(completeness.get("product_mapping_coverage") or 0.0),
            mapping_product_coverage,
        )
        >= float(unresolved_policy.get("minimum_product_mapping_coverage") or 1.0)
        and (
            not unresolved_policy.get("require_no_product_element_excess")
            or not dict(completeness.get("product_element_excess") or {})
        )
        and (
            not unresolved_policy.get("require_no_suspected_missing_reactant")
            or not bool(completeness.get("suspected_missing_reactant"))
        )
        and not bool(
            unresolved_policy.get("require_no_suspected_multiplicity_problem")
            and completeness.get("suspected_insufficient_reactant_multiplicity")
        )
    )
    if completeness_status != "verified" and not mapper_supported_unresolved:
        review_reasons.append("reaction_completeness_not_verified")
    if condition_status not in set(rules["indexable_condition_statuses"]):
        review_reasons.append("conditions_not_precedent_eligible")
    if stage_status not in set(rules["indexable_stage_statuses"]):
        review_reasons.append("condition_stages_not_precedent_eligible")
    if record.get("partial_product_transformation"):
        review_reasons.append("partial_product_transformation")
    if mapping_status and mapping_status not in set(rules["review_mapping_statuses"]):
        review_reasons.append("external_mapping_lacks_internal_consensus")
    if index_eligibility not in {"eligible", "review_only"}:
        review_reasons.append("record_not_review_index_eligible")
    if chemistry_status not in {"verified", "review"}:
        review_reasons.append("chemistry_not_review_qualified")

    if not review_reasons:
        return CoreEligibilityAssessment(
            CoreEligibility.REVIEW_CORE,
            ("expert_review_required",),
            True,
            True,
            True,
        )
    return CoreEligibilityAssessment(
        CoreEligibility.QUERY_ONLY,
        tuple(sorted(set(review_reasons))),
        True,
        False,
        True,
    )


__all__ = [
    "CoreEligibilityAssessment",
    "assess_core_eligibility",
    "load_core_eligibility_rules",
]
