"""Chemistry-first common admission for mixed recommendation datasets."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Optional, Tuple

from .identities import CanonicalReactionIdentity
from .input_schema import RawReactionRecord
from ..models import AdmissionTier, ConditionIdentity


@dataclass(frozen=True)
class AdmissionDecision:
    tier: AdmissionTier
    reasons: Tuple[str, ...]
    policy_version: str = "generic_admission.v1"


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
) -> AdmissionDecision:
    """Classify a row without requiring a named reaction family."""
    policy = load_admission_policy()
    minimum_yield, maximum_yield = policy["valid_yield_range"]
    if not analysis.valid or canonical_identity is None:
        return AdmissionDecision(
            AdmissionTier.REJECTED, ("invalid_reaction_or_product",)
        )
    if record.yield_pct is None or not minimum_yield <= record.yield_pct <= maximum_yield:
        return AdmissionDecision(
            AdmissionTier.REJECTED, ("missing_or_invalid_yield",)
        )
    raw_conditions_present = bool(
        record.catalyst_cas or record.reagent_cas or record.solvent_cas
    )
    if not raw_conditions_present:
        return AdmissionDecision(
            AdmissionTier.REJECTED, ("no_condition_identifiers",)
        )
    signature = analysis.reaction_signature
    if signature is None:
        if analysis.candidates or analysis.evidence_quality in {
            "reactant_grammar_only",
            "ambiguous",
        }:
            return AdmissionDecision(
                AdmissionTier.REVIEW,
                ("missing_verified_reaction_signature",),
            )
        return AdmissionDecision(
            AdmissionTier.REJECTED,
            ("no_usable_transformation_evidence",),
        )
    reasons = []
    if analysis.evidence_quality == "conflicting_edit_evidence":
        reasons.append("conflicting_edit_evidence")
    if any(
        warning in analysis.warnings for warning in policy["review_warnings"]
    ):
        reasons.append("multiple_products")
    if analysis.evidence_quality not in set(policy["verified_evidence"]):
        reasons.append("insufficient_edit_evidence")
    normalized_conditions_present = bool(
        conditions.catalyst_cas
        or conditions.reagent_cas
        or conditions.solvent_cas
    )
    if not normalized_conditions_present:
        reasons.append("unresolved_condition_identifiers")
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
    if normalized_identifier_count < raw_identifier_count:
        reasons.append("condition_identifier_uncertainty")
    if reasons:
        return AdmissionDecision(
            AdmissionTier.REVIEW, tuple(sorted(set(reasons)))
        )
    return AdmissionDecision(
        AdmissionTier.VERIFIED,
        ("verified_reaction_signature_with_conditions",),
    )


__all__ = ["AdmissionDecision", "decide_admission", "load_admission_policy"]
