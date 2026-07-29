"""Conservative precedent retrieval for parsed reactions without signatures."""

from __future__ import annotations

from typing import Any, Mapping

from .compatibility import filter_compatible_precedents
from .fallback_similarity import (
    assess_fallback_similarity,
    compatibility_signature_from_fallback,
    fallback_index_tokens,
    load_fallback_retrieval_rules,
    shared_high_signal_count,
)
from .generic_indexing import GenericReactionIndex
from .generic_retrieval import CompatibleRetrievalResult
from .models import RetrievalLevelTrace
from .support import summarize_evidence_support


def retrieve_fallback_pool_with_trace(
    descriptor: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
) -> CompatibleRetrievalResult:
    """Retrieve high-similarity precedents without asserting query bond edits."""
    rules = load_fallback_retrieval_rules()
    configured_minimum = int(rules["minimum_independent_support"])
    minimum = max(
        configured_minimum,
        int(minimum_pool_size) if minimum_pool_size is not None else configured_minimum,
    )
    query_tokens = fallback_index_tokens(descriptor)
    positions = set()
    for token in query_tokens:
        positions.update(index.fallback_features.get(token, ()))
    raw_rows = tuple(
        index.rows[position]
        for position in sorted(positions)
        if index.rows[position].fallback_descriptor
    )
    raw_support = summarize_evidence_support(raw_rows)
    candidate_trace = RetrievalLevelTrace(
        level="fallback_descriptor_candidates",
        candidate_count=len(raw_rows),
        independent_candidate_count=raw_support.independent_count,
        compatible_candidate_count=len(raw_rows),
        independent_compatible_candidate_count=raw_support.independent_count,
        excluded_candidate_count=0,
        minimum_independent_support=minimum,
        status="scored" if raw_rows else "empty",
    )
    if not raw_rows:
        return CompatibleRetrievalResult(
            level="no_fallback_descriptor_candidate",
            pool=(),
            candidate_count=0,
            independent_candidate_count=0,
            excluded_candidate_count=0,
            independent_compatible_candidate_count=0,
            trace=(candidate_trace,),
        )

    minimum_similarity = float(rules["minimum_similarity"])
    minimum_shared = int(rules["minimum_shared_high_signal_features"])
    scored = []
    for row in raw_rows:
        precedent = row.fallback_descriptor
        shared = shared_high_signal_count(descriptor, precedent)
        if shared < minimum_shared:
            continue
        assessment = assess_fallback_similarity(descriptor, precedent)
        if assessment.score >= minimum_similarity:
            scored.append((assessment.score, shared, row))
    scored.sort(
        key=lambda item: (
            -item[0],
            -item[1],
            item[2].canonical_reaction_id,
            item[2].reaction_id,
            item[2].observation_id,
        )
    )
    matched_rows = tuple(item[2] for item in scored[: int(rules["candidate_limit"])])
    compatibility_signature = compatibility_signature_from_fallback(descriptor)
    accepted, compatibility_excluded = filter_compatible_precedents(
        compatibility_signature,
        matched_rows,
    )
    excluded_count = len(compatibility_excluded)
    accepted_rows = tuple(row for row, _ in accepted)
    accepted_support = summarize_evidence_support(accepted_rows)
    matched_support = summarize_evidence_support(matched_rows)
    if not accepted:
        status = "no_compatible_recipe" if matched_rows else "below_similarity_gate"
        trace = RetrievalLevelTrace(
            level="unverified_structure_similarity",
            candidate_count=len(matched_rows),
            independent_candidate_count=matched_support.independent_count,
            compatible_candidate_count=0,
            independent_compatible_candidate_count=0,
            excluded_candidate_count=excluded_count,
            minimum_independent_support=minimum,
            status=status,
        )
        return CompatibleRetrievalResult(
            level=(
                "no_compatible_condition_precedent"
                if matched_rows
                else "no_safe_fallback_similarity_match"
            ),
            pool=(),
            candidate_count=len(matched_rows),
            independent_candidate_count=matched_support.independent_count,
            excluded_candidate_count=excluded_count,
            independent_compatible_candidate_count=0,
            trace=(candidate_trace, trace),
        )

    selected_level = "unverified_structure_fallback"
    status = "selected"
    if accepted_support.independent_count < minimum:
        best_score = max(
            assess_fallback_similarity(
                descriptor,
                row.fallback_descriptor,
            ).score
            for row in accepted_rows
        )
        if best_score < float(rules["limited_support_minimum_similarity"]):
            trace = RetrievalLevelTrace(
                level="unverified_structure_similarity",
                candidate_count=len(matched_rows),
                independent_candidate_count=matched_support.independent_count,
                compatible_candidate_count=len(accepted),
                independent_compatible_candidate_count=(
                    accepted_support.independent_count
                ),
                excluded_candidate_count=excluded_count,
                minimum_independent_support=minimum,
                status="insufficient_independent_support",
            )
            return CompatibleRetrievalResult(
                level="insufficient_safe_fallback_support",
                pool=(),
                candidate_count=len(matched_rows),
                independent_candidate_count=matched_support.independent_count,
                excluded_candidate_count=excluded_count,
                independent_compatible_candidate_count=(
                    accepted_support.independent_count
                ),
                trace=(candidate_trace, trace),
            )
        selected_level += "_limited_support"
        status = "selected_limited_support"
    trace = RetrievalLevelTrace(
        level="unverified_structure_similarity",
        candidate_count=len(matched_rows),
        independent_candidate_count=matched_support.independent_count,
        compatible_candidate_count=len(accepted),
        independent_compatible_candidate_count=accepted_support.independent_count,
        excluded_candidate_count=excluded_count,
        minimum_independent_support=minimum,
        status=status,
    )
    return CompatibleRetrievalResult(
        level=selected_level,
        pool=accepted,
        candidate_count=len(matched_rows),
        independent_candidate_count=matched_support.independent_count,
        excluded_candidate_count=excluded_count,
        independent_compatible_candidate_count=accepted_support.independent_count,
        trace=(candidate_trace, trace),
    )


__all__ = ["retrieve_fallback_pool_with_trace"]
