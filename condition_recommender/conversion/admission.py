"""Chemistry-first common admission for mixed recommendation datasets."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Optional, Tuple

from condition_registry import ResolvedConditionRecipe

from .identities import CanonicalReactionIdentity
from .input_schema import RawReactionRecord
from ..models import (
    AdmissionTier,
    ChemistryStatus,
    ConditionIdentity,
    ConditionStageStatus,
    ConditionStatus,
    IndexEligibility,
    OutcomeStatus,
)


@dataclass(frozen=True)
class AdmissionDecision:
    tier: AdmissionTier
    reasons: Tuple[str, ...]
    chemistry_status: ChemistryStatus
    condition_status: ConditionStatus
    condition_stage_status: ConditionStageStatus
    outcome_status: OutcomeStatus
    index_eligibility: IndexEligibility
    policy_version: str = "generic_admission.v1.6"


@lru_cache(maxsize=1)
def load_admission_policy() -> dict[str, Any]:
    path = Path(__file__).parents[1] / "definitions" / "generic_admission.v1.json"
    with path.open("r", encoding="utf-8") as handle:
        return dict(json.load(handle))


def decide_admission(
    *,
    record: RawReactionRecord,
    analysis: Any,
    canonical_identity: Optional[CanonicalReactionIdentity],
    conditions: ConditionIdentity,
    resolved_recipe: ResolvedConditionRecipe,
) -> AdmissionDecision:
    """Assess independent evidence dimensions and derive the legacy tier."""
    policy = load_admission_policy()
    minimum_yield, maximum_yield = policy["valid_yield_range"]

    chemistry_status = ChemistryStatus.VERIFIED
    chemistry_reasons: list[str] = []
    completeness = analysis.reaction_completeness
    fallback_descriptor = analysis.fallback_descriptor
    if not analysis.valid or canonical_identity is None:
        chemistry_status = ChemistryStatus.REJECTED
        chemistry_reasons.append("invalid_reaction_or_product")
    elif completeness is not None and completeness.status == "incomplete":
        if fallback_descriptor is not None and fallback_descriptor.retrieval_eligible:
            chemistry_status = ChemistryStatus.REVIEW
            chemistry_reasons.append("exploratory_partial_product_correspondence")
        else:
            chemistry_status = ChemistryStatus.REJECTED
            chemistry_reasons.append("unaccounted_product_heavy_atoms")
            if completeness.suspected_insufficient_reactant_multiplicity:
                chemistry_reasons.append("insufficient_reactant_multiplicity")
            elif completeness.suspected_missing_reactant:
                chemistry_reasons.append("suspected_missing_reactant")
    elif analysis.reaction_signature is None:
        if analysis.candidates or analysis.evidence_quality in {
            "reactant_grammar_only",
            "ambiguous",
        }:
            chemistry_status = ChemistryStatus.REVIEW
            chemistry_reasons.append("missing_verified_reaction_signature")
        else:
            chemistry_status = ChemistryStatus.REJECTED
            chemistry_reasons.append("no_usable_transformation_evidence")
        if completeness is not None and completeness.status == "unresolved":
            chemistry_reasons.append("reaction_completeness_unresolved")
    else:
        if analysis.evidence_quality in {
            "conflicting_edit_evidence",
            "conflicting_stereochemical_evidence",
        }:
            chemistry_reasons.append(analysis.evidence_quality)
        if any(warning in analysis.warnings for warning in policy["review_warnings"]):
            chemistry_reasons.append("multiple_products")
        if analysis.evidence_quality not in set(policy["verified_evidence"]):
            chemistry_reasons.append("insufficient_edit_evidence")
        if completeness is not None and completeness.status == "unresolved":
            chemistry_reasons.append("reaction_completeness_unresolved")
        completeness_warnings = set(
            completeness.warnings if completeness is not None else ()
        )
        if "PRODUCT_MAPS_MISSING_FROM_REACTANTS" in completeness_warnings:
            chemistry_reasons.append("inconsistent_product_atom_mapping")
        if (
            "PARTIAL_ATOM_MAPPING" in completeness_warnings
            and analysis.evidence_quality == "validated_atom_mapping"
        ):
            chemistry_reasons.append("partial_atom_mapping")
        if chemistry_reasons:
            chemistry_status = ChemistryStatus.REVIEW

    if record.yield_pct is None:
        outcome_status = OutcomeStatus.MISSING
        outcome_reasons = ["missing_yield"]
    elif minimum_yield <= record.yield_pct <= maximum_yield:
        outcome_status = OutcomeStatus.USABLE
        outcome_reasons = []
    else:
        outcome_status = OutcomeStatus.INVALID
        outcome_reasons = ["invalid_yield"]

    raw_conditions_present = bool(
        record.catalyst_cas or record.reagent_cas or record.solvent_cas
    )
    condition_reasons: list[str] = []
    condition_stage_status = _condition_stage_status(record.stages)
    if not raw_conditions_present:
        condition_status = ConditionStatus.UNUSABLE
        condition_reasons.append("no_condition_identifiers")
    elif not _has_non_solvent_component(resolved_recipe):
        condition_status = ConditionStatus.UNUSABLE
        condition_reasons.append("solvent_only_recipe")
    else:
        statuses = {
            component.identity_status for component in resolved_recipe.components
        }
        if "invalid_identifier" in statuses:
            condition_status = ConditionStatus.RESOLVED_PARTIAL
            condition_reasons.extend(
                ("condition_identifier_uncertainty", "unresolved_condition_identifiers")
            )
        elif any(status != "resolved" for status in statuses):
            condition_status = ConditionStatus.UNRESOLVED_RETAINED
            condition_reasons.append("condition_identifier_uncertainty")
        else:
            condition_status = ConditionStatus.RESOLVED_COMPLETE
    if condition_stage_status == ConditionStageStatus.UNASSIGNED_MULTISTAGE:
        condition_reasons.append("multistage_conditions_not_structured")

    raw_identifier_count = sum(
        len(values)
        for values in (
            record.catalyst_cas,
            record.reagent_cas,
            record.solvent_cas,
        )
    )
    normalized_identifier_count = sum(
        len(values)
        for values in (
            conditions.catalyst_cas,
            conditions.reagent_cas,
            conditions.solvent_cas,
        )
    )
    if (
        raw_conditions_present
        and normalized_identifier_count < raw_identifier_count
        and "condition_identifier_uncertainty" not in condition_reasons
    ):
        condition_reasons.append("condition_identifier_uncertainty")
        if condition_status == ConditionStatus.RESOLVED_COMPLETE:
            condition_status = ConditionStatus.RESOLVED_PARTIAL

    condition_is_indexable = condition_status in {
        ConditionStatus.RESOLVED_COMPLETE,
        ConditionStatus.UNRESOLVED_RETAINED,
    }
    stage_is_indexable = condition_stage_status.value in set(
        policy.get("indexable_stage_statuses") or ()
    ) or (
        condition_stage_status == ConditionStageStatus.UNASSIGNED_MULTISTAGE
        and condition_status.value
        in set(policy.get("indexable_unassigned_multistage_condition_statuses") or ())
    )
    chemistry_is_indexable = (
        chemistry_status == ChemistryStatus.VERIFIED
        or (
            chemistry_status == ChemistryStatus.REVIEW
            and analysis.reaction_signature is not None
            and completeness is not None
            and completeness.status == "verified"
            and analysis.evidence_quality
            in set(policy.get("indexable_review_evidence") or ())
        )
        or (
            chemistry_status == ChemistryStatus.REVIEW
            and fallback_descriptor is not None
            and fallback_descriptor.retrieval_eligible
        )
    )
    if (
        chemistry_status == ChemistryStatus.REJECTED
        or condition_status == ConditionStatus.UNUSABLE
    ):
        index_eligibility = IndexEligibility.INELIGIBLE
    elif chemistry_is_indexable and condition_is_indexable and stage_is_indexable:
        index_eligibility = IndexEligibility.ELIGIBLE
    else:
        index_eligibility = IndexEligibility.REVIEW_ONLY

    if index_eligibility == IndexEligibility.INELIGIBLE:
        tier = AdmissionTier.REJECTED
    elif (
        chemistry_status == ChemistryStatus.VERIFIED
        and condition_status == ConditionStatus.RESOLVED_COMPLETE
        and condition_stage_status
        in {
            ConditionStageStatus.SINGLE_STAGE,
            ConditionStageStatus.STRUCTURED_MULTISTAGE,
        }
        and outcome_status == OutcomeStatus.USABLE
        and index_eligibility == IndexEligibility.ELIGIBLE
    ):
        tier = AdmissionTier.VERIFIED
    else:
        tier = AdmissionTier.REVIEW

    reasons = sorted(set(chemistry_reasons + condition_reasons + outcome_reasons))
    if tier == AdmissionTier.VERIFIED:
        reasons = ["verified_reaction_signature_with_conditions"]
    return AdmissionDecision(
        tier=tier,
        reasons=tuple(reasons),
        chemistry_status=chemistry_status,
        condition_status=condition_status,
        condition_stage_status=condition_stage_status,
        outcome_status=outcome_status,
        index_eligibility=index_eligibility,
    )


def _condition_stage_status(raw_stages: str) -> ConditionStageStatus:
    value = str(raw_stages or "").strip()
    if not value:
        return ConditionStageStatus.SINGLE_STAGE
    try:
        stage_count = int(value)
    except ValueError:
        return ConditionStageStatus.UNASSIGNED_MULTISTAGE
    if stage_count <= 1:
        return ConditionStageStatus.SINGLE_STAGE
    return ConditionStageStatus.UNASSIGNED_MULTISTAGE


def _has_non_solvent_component(recipe: ResolvedConditionRecipe) -> bool:
    return any(component.primary_role != "solvent" for component in recipe.components)


__all__ = ["AdmissionDecision", "decide_admission", "load_admission_policy"]
