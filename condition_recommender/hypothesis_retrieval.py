"""Conservative query-only retrieval across ambiguous edit hypotheses."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Tuple

from .compatibility import filter_compatible_precedents
from .edit_prototypes import (
    AnonymousEditPrototype,
    anonymous_edit_prototype,
    anonymous_edit_prototype_from_hypothesis,
    anonymous_edit_similarity,
)
from .fallback_similarity import (
    assess_fallback_similarity,
    compatibility_signature_from_fallback,
)
from .generic_indexing import GenericIndexedReaction, GenericReactionIndex
from .generic_retrieval import CompatibleRetrievalResult
from .models import RetrievalLevelTrace
from .similarity import SimilarityAssessment
from .support import summarize_evidence_support


_RULES_PATH = (
    Path(__file__).with_name("definitions")
    / "edit_hypothesis_retrieval.v1.json"
)


@lru_cache(maxsize=1)
def load_edit_hypothesis_retrieval_rules() -> dict[str, Any]:
    """Load and validate the query-only ambiguous-edit policy."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "1.0":
        raise ValueError("unsupported edit-hypothesis retrieval schema")
    if str(rules.get("definition_id") or "") != (
        "edit_hypothesis_retrieval.v1"
    ):
        raise ValueError("unexpected edit-hypothesis retrieval definition")
    if not str(rules.get("calibration_status") or ""):
        raise ValueError("edit-hypothesis retrieval requires calibration status")
    if int(rules.get("minimum_hypothesis_count") or 0) < 2:
        raise ValueError("edit-hypothesis retrieval requires ambiguity")
    threshold = float(rules.get("minimum_prototype_similarity") or 0.0)
    if not 0.0 < threshold <= 1.0:
        raise ValueError("invalid edit-hypothesis similarity threshold")
    if int(rules.get("minimum_independent_support") or 0) < 1:
        raise ValueError("edit-hypothesis support must be positive")
    if int(rules.get("candidate_limit") or 0) < 1:
        raise ValueError("edit-hypothesis candidate limit must be positive")
    if rules.get("consensus_policy") != "intersection_across_all_hypotheses":
        raise ValueError("unsupported edit-hypothesis consensus policy")
    if rules.get("precedent_policy") != "verified_signatures_only":
        raise ValueError("unsupported edit-hypothesis precedent policy")
    weights = rules.get("similarity_weights")
    if not isinstance(weights, Mapping) or set(weights) != {
        "edit_hypothesis_consensus",
        "structure_context",
    }:
        raise ValueError("invalid edit-hypothesis similarity weights")
    normalized = {str(key): float(value) for key, value in weights.items()}
    if any(value < 0.0 for value in normalized.values()):
        raise ValueError("edit-hypothesis weights must be non-negative")
    if abs(sum(normalized.values()) - 1.0) > 1e-9:
        raise ValueError("edit-hypothesis weights must sum to one")
    rules["similarity_weights"] = normalized
    return rules


@dataclass(frozen=True)
class HypothesisRetrievalQuery:
    """Serialized ambiguous query context used for robust ranking."""

    hypotheses: Tuple[Mapping[str, Any], ...]
    prototypes: Tuple[AnonymousEditPrototype, ...]
    fallback_descriptor: Mapping[str, Any]

    def to_mapping(self) -> dict[str, Any]:
        return {
            "hypotheses": tuple(dict(value) for value in self.hypotheses),
            "prototypes": tuple(asdict(value) for value in self.prototypes),
            "fallback_descriptor": dict(self.fallback_descriptor),
        }


def build_hypothesis_retrieval_query(
    hypotheses: Tuple[Mapping[str, Any], ...],
    fallback_descriptor: Mapping[str, Any],
) -> HypothesisRetrievalQuery | None:
    """Build a query only when every retained hypothesis has an EP1 view."""
    rules = load_edit_hypothesis_retrieval_rules()
    if len(hypotheses) < int(rules["minimum_hypothesis_count"]):
        return None
    prototypes = tuple(
        prototype
        for hypothesis in hypotheses
        for prototype in (
            anonymous_edit_prototype_from_hypothesis(hypothesis),
        )
        if prototype is not None
    )
    if len(prototypes) != len(hypotheses):
        return None
    return HypothesisRetrievalQuery(
        hypotheses=hypotheses,
        prototypes=prototypes,
        fallback_descriptor=fallback_descriptor,
    )


def _candidate_positions(
    prototype: AnonymousEditPrototype,
    index: GenericReactionIndex,
    threshold: float,
) -> set[int]:
    positions = set()
    for position, row in enumerate(index.rows):
        if not row.signature:
            continue
        precedent = anonymous_edit_prototype(row.signature)
        if precedent is None:
            continue
        if anonymous_edit_similarity(prototype, precedent) >= threshold:
            positions.add(position)
    return positions


def retrieve_hypothesis_consensus_pool_with_trace(
    query: HypothesisRetrievalQuery,
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
) -> CompatibleRetrievalResult:
    """Return only verified precedents compatible with every edit hypothesis."""
    rules = load_edit_hypothesis_retrieval_rules()
    minimum = max(
        int(rules["minimum_independent_support"]),
        int(minimum_pool_size)
        if minimum_pool_size is not None
        else int(rules["minimum_independent_support"]),
    )
    threshold = float(rules["minimum_prototype_similarity"])
    traces = []
    position_sets = []
    for hypothesis, prototype in zip(query.hypotheses, query.prototypes):
        positions = _candidate_positions(prototype, index, threshold)
        position_sets.append(positions)
        rows = index.select(sorted(positions))
        support = summarize_evidence_support(rows)
        traces.append(
            RetrievalLevelTrace(
                level=(
                    "edit_hypothesis:"
                    f"{hypothesis.get('hypothesis_id') or prototype.prototype_id}"
                ),
                candidate_count=len(rows),
                independent_candidate_count=support.independent_count,
                compatible_candidate_count=len(rows),
                independent_compatible_candidate_count=support.independent_count,
                excluded_candidate_count=0,
                minimum_independent_support=minimum,
                status="scored" if rows else "empty",
            )
        )
    robust_positions = (
        set.intersection(*position_sets) if position_sets else set()
    )
    scored = []
    for position in robust_positions:
        row = index.rows[position]
        precedent = anonymous_edit_prototype(row.signature)
        if precedent is None:
            continue
        worst_case = min(
            anonymous_edit_similarity(prototype, precedent)
            for prototype in query.prototypes
        )
        scored.append(
            (
                worst_case,
                row.canonical_reaction_id,
                row.reaction_id,
                row.observation_id,
                position,
            )
        )
    scored.sort(key=lambda item: (-item[0],) + item[1:])
    limit = int(rules["candidate_limit"])
    robust_positions = {item[-1] for item in scored[:limit]}
    rows = index.select(sorted(robust_positions))
    raw_support = summarize_evidence_support(rows)
    compatibility_signature = compatibility_signature_from_fallback(
        query.fallback_descriptor
    )
    accepted, excluded = filter_compatible_precedents(
        compatibility_signature,
        rows,
    )
    accepted_rows = tuple(row for row, _ in accepted)
    accepted_support = summarize_evidence_support(accepted_rows)
    if not rows:
        status = "no_consensus_candidate"
    elif not accepted:
        status = "no_compatible_recipe"
    elif accepted_support.independent_count >= minimum:
        status = "selected"
    else:
        status = "insufficient_support"
    consensus_trace = RetrievalLevelTrace(
        level="edit_hypothesis_consensus",
        candidate_count=len(rows),
        independent_candidate_count=raw_support.independent_count,
        compatible_candidate_count=len(accepted),
        independent_compatible_candidate_count=accepted_support.independent_count,
        excluded_candidate_count=len(excluded),
        minimum_independent_support=minimum,
        status=status,
    )
    traces.append(consensus_trace)
    if not accepted:
        return CompatibleRetrievalResult(
            level=(
                "no_compatible_condition_precedent"
                if rows
                else "no_edit_hypothesis_consensus_precedent"
            ),
            pool=(),
            candidate_count=len(rows),
            independent_candidate_count=raw_support.independent_count,
            excluded_candidate_count=len(excluded),
            independent_compatible_candidate_count=0,
            trace=tuple(traces),
        )
    if accepted_support.independent_count < minimum:
        return CompatibleRetrievalResult(
            level="insufficient_edit_hypothesis_consensus_support",
            pool=(),
            candidate_count=len(rows),
            independent_candidate_count=raw_support.independent_count,
            excluded_candidate_count=len(excluded),
            independent_compatible_candidate_count=(
                accepted_support.independent_count
            ),
            trace=tuple(traces),
        )
    return CompatibleRetrievalResult(
        level="edit_hypothesis_consensus",
        pool=accepted,
        candidate_count=len(rows),
        independent_candidate_count=raw_support.independent_count,
        excluded_candidate_count=len(excluded),
        independent_compatible_candidate_count=accepted_support.independent_count,
        trace=tuple(traces),
    )


def _prototype_from_mapping(
    value: Mapping[str, Any],
) -> AnonymousEditPrototype:
    return AnonymousEditPrototype(
        prototype_id=str(value["prototype_id"]),
        formed_element_pairs=tuple(value.get("formed_element_pairs") or ()),
        broken_element_pairs=tuple(value.get("broken_element_pairs") or ()),
        order_changed_element_pairs=tuple(
            value.get("order_changed_element_pairs") or ()
        ),
        ring_count_delta=int(value.get("ring_count_delta") or 0),
        formed_ring_sizes=tuple(
            int(item) for item in value.get("formed_ring_sizes") or ()
        ),
        event_count=int(value.get("event_count") or 0),
        version=str(value.get("version") or ""),
    )


def assess_hypothesis_consensus_similarity(
    query: Mapping[str, Any],
    row: GenericIndexedReaction,
) -> SimilarityAssessment:
    """Score a precedent by worst-case edit agreement plus observed context."""
    rules = load_edit_hypothesis_retrieval_rules()
    prototypes = tuple(
        _prototype_from_mapping(value)
        for value in query.get("prototypes") or ()
        if isinstance(value, Mapping)
    )
    precedent = anonymous_edit_prototype(row.signature)
    edit_score = (
        min(
            anonymous_edit_similarity(prototype, precedent)
            for prototype in prototypes
        )
        if prototypes and precedent is not None
        else 0.0
    )
    descriptor = query.get("fallback_descriptor")
    descriptor_value = descriptor if isinstance(descriptor, Mapping) else {}
    context_score = (
        assess_fallback_similarity(
            descriptor_value,
            row.fallback_descriptor,
        ).score
        if descriptor_value and row.fallback_descriptor
        else 0.0
    )
    weights = rules["similarity_weights"]
    components = {
        "edit_hypothesis_consensus": edit_score,
        "structure_context": context_score,
    }
    available = {"edit_hypothesis_consensus"}
    if descriptor_value and row.fallback_descriptor:
        available.add("structure_context")
    denominator = sum(float(weights[name]) for name in available)
    contributions = {
        name: (
            float(weights[name]) / denominator * components[name]
            if denominator and name in available
            else 0.0
        )
        for name in components
    }
    return SimilarityAssessment(
        score=round(sum(contributions.values()), 6),
        components={
            name: round(value, 6) for name, value in components.items()
        },
        contributions={
            name: round(value, 6) for name, value in contributions.items()
        },
        definition_id=str(rules["definition_id"]),
        definition_version=str(rules["schema_version"]),
    )


__all__ = [
    "HypothesisRetrievalQuery",
    "assess_hypothesis_consensus_similarity",
    "build_hypothesis_retrieval_query",
    "load_edit_hypothesis_retrieval_rules",
    "retrieve_hypothesis_consensus_pool_with_trace",
]
