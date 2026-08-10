"""Condition-evidence enrichment for verified retrosynthesis candidates.

This application-layer bridge deliberately keeps graph validation authoritative.
Only candidates with a verified forward signature are sent to the condition
recommender.  Missing condition evidence is retained as uncertainty rather than
being interpreted as a chemically invalid disconnection.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Iterable, Literal, Protocol, Tuple

from .generic_models import GenericDisconnectionCandidate


ConditionSupportStatus = Literal[
    "recommended_direct",
    "recommended_fallback",
    "insufficient_evidence",
]


class ConditionRecommenderProtocol(Protocol):
    """Narrow condition-recommender interface consumed by the bridge."""

    def recommend(
        self,
        reaction_smiles: str,
        *,
        top_k: int = 5,
        minimum_pool_size: int | None = None,
        unrestricted_fallback: bool = False,
    ) -> Any: ...


def _serialized(value: Any) -> Dict[str, Any]:
    if hasattr(value, "to_dict"):
        return dict(value.to_dict())
    if hasattr(value, "__dataclass_fields__"):
        return asdict(value)
    if isinstance(value, dict):
        return dict(value)
    raise TypeError("condition recommendation must be serializable")


@dataclass(frozen=True)
class RetrosynthesisConditionEvidence:
    """Condition support attached to one structurally verified candidate."""

    status: ConditionSupportStatus
    query_reaction_smiles: str
    recommender_valid: bool
    recommendation_mode: str
    retrieval_level: str | None
    uses_type_agnostic_fallback: bool
    candidate_count: int
    independent_candidate_count: int
    compatible_candidate_count: int
    independent_compatible_candidate_count: int
    excluded_candidate_count: int
    best_recipe_score: float | None
    best_recipe_compatibility_score: float | None
    best_recipe_reference_support: int
    recommendations: Tuple[Dict[str, Any], ...]
    warnings: Tuple[str, ...]
    error: str | None

    def to_dict(self) -> Dict[str, Any]:
        """Return an auditable JSON-compatible representation."""

        return asdict(self)


@dataclass(frozen=True)
class ConditionRankedRetrosynthesisCandidate:
    """A verified retrosynthesis proposal with condition evidence and ranks."""

    candidate: GenericDisconnectionCandidate
    retrosynthesis_rank: int
    condition_informed_rank: int
    condition_evidence: RetrosynthesisConditionEvidence

    def to_dict(self) -> Dict[str, Any]:
        """Preserve the candidate schema and append condition-ranking fields."""

        return {
            **self.candidate.to_dict(),
            "retrosynthesis_rank": self.retrosynthesis_rank,
            "condition_informed_rank": self.condition_informed_rank,
            "condition_evidence": self.condition_evidence.to_dict(),
        }


def _condition_evidence(result: Any) -> RetrosynthesisConditionEvidence:
    recommendations = tuple(
        _serialized(value) for value in (result.recommendations or ())
    )
    warnings = tuple(str(value) for value in (result.warnings or ()))
    fallback = (
        str(result.recommendation_mode or "") != "verified_signature"
        or "TYPE_AGNOSTIC_FALLBACK_USED" in warnings
        or "REACTION_CORE_RETRIEVAL_USED" in warnings
        or "REACTION_TOPOLOGY_FALLBACK_USED" in warnings
    )
    if result.valid and recommendations:
        status: ConditionSupportStatus = (
            "recommended_fallback" if fallback else "recommended_direct"
        )
    else:
        status = "insufficient_evidence"
    best = recommendations[0] if recommendations else {}
    return RetrosynthesisConditionEvidence(
        status=status,
        query_reaction_smiles=str(result.query_reaction_smiles or ""),
        recommender_valid=bool(result.valid),
        recommendation_mode=str(result.recommendation_mode or ""),
        retrieval_level=(
            str(result.retrieval_level) if result.retrieval_level else None
        ),
        uses_type_agnostic_fallback=fallback,
        candidate_count=int(result.candidate_count or 0),
        independent_candidate_count=int(
            result.independent_candidate_count or 0
        ),
        compatible_candidate_count=int(result.compatible_candidate_count or 0),
        independent_compatible_candidate_count=int(
            result.independent_compatible_candidate_count or 0
        ),
        excluded_candidate_count=int(result.excluded_candidate_count or 0),
        best_recipe_score=(
            float(best["score"]) if best.get("score") is not None else None
        ),
        best_recipe_compatibility_score=(
            float(best["compatibility_score"])
            if best.get("compatibility_score") is not None
            else None
        ),
        best_recipe_reference_support=int(best.get("reference_support") or 0),
        recommendations=recommendations,
        warnings=warnings,
        error=str(result.error) if result.error else None,
    )


def _ranking_key(
    value: tuple[
        int,
        GenericDisconnectionCandidate,
        RetrosynthesisConditionEvidence,
    ],
) -> tuple[int, int, int, float, int]:
    original_rank, _, evidence = value
    status_rank = {
        "recommended_direct": 0,
        "recommended_fallback": 1,
        "insufficient_evidence": 2,
    }[evidence.status]
    return (
        status_rank,
        -evidence.independent_compatible_candidate_count,
        -evidence.best_recipe_reference_support,
        -(evidence.best_recipe_score or 0.0),
        original_rank,
    )


def rank_retrosynthesis_candidates_with_conditions(
    candidates: Iterable[GenericDisconnectionCandidate],
    recommender: ConditionRecommenderProtocol,
    *,
    condition_top_k: int = 3,
    minimum_pool_size: int | None = None,
    unrestricted_fallback: bool = False,
    rerank: bool = True,
) -> Tuple[ConditionRankedRetrosynthesisCandidate, ...]:
    """Enrich and optionally rerank structurally verified candidates.

    Ranking is lexicographic and auditable: direct condition support precedes
    fallback support, which precedes insufficient evidence.  Within a tier,
    independent compatible support, best-recipe reference support, and the
    recommender's recipe score are considered in that order.  The original
    retrosynthesis rank is the deterministic final tie-breaker.
    """

    if condition_top_k < 1:
        raise ValueError("condition_top_k must be positive")
    if minimum_pool_size is not None and minimum_pool_size < 1:
        raise ValueError("minimum_pool_size must be positive")
    assessed = []
    for original_rank, candidate in enumerate(candidates, start=1):
        if candidate.forward_validation_status != "verified_signature":
            continue
        result = recommender.recommend(
            candidate.proposed_reaction_smiles,
            top_k=condition_top_k,
            minimum_pool_size=minimum_pool_size,
            unrestricted_fallback=unrestricted_fallback,
        )
        assessed.append(
            (original_rank, candidate, _condition_evidence(result))
        )
    ordered = sorted(assessed, key=_ranking_key) if rerank else assessed
    return tuple(
        ConditionRankedRetrosynthesisCandidate(
            candidate=candidate,
            retrosynthesis_rank=original_rank,
            condition_informed_rank=condition_rank,
            condition_evidence=evidence,
        )
        for condition_rank, (original_rank, candidate, evidence) in enumerate(
            ordered,
            start=1,
        )
    )


__all__ = [
    "ConditionRankedRetrosynthesisCandidate",
    "ConditionRecommenderProtocol",
    "ConditionSupportStatus",
    "RetrosynthesisConditionEvidence",
    "rank_retrosynthesis_candidates_with_conditions",
]
