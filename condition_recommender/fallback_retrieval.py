"""Conservative precedent retrieval for parsed reactions without signatures."""

from __future__ import annotations

from collections import Counter
from typing import Any, Mapping

from .compatibility import CompatibilityAssessment, filter_compatible_precedents
from .fallback_similarity import (
    assess_fallback_similarity,
    compatibility_signature_from_fallback,
    fallback_edit_index_tokens,
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
    unrestricted: bool = False,
    required_source_capability_ids: Mapping[str, str] | None = None,
    required_source_substance_ids: Mapping[str, str] | None = None,
) -> CompatibleRetrievalResult:
    """Retrieve precedents without asserting query bond edits.

    The unrestricted expert mode bypasses chemistry eligibility, similarity,
    support, and recipe-compatibility gates. Feature-index overlap and the
    configured candidate limit remain as computational bounds.
    """
    rules = load_fallback_retrieval_rules()
    configured_minimum = int(rules["minimum_independent_support"])
    minimum = max(
        configured_minimum,
        int(minimum_pool_size) if minimum_pool_size is not None else configured_minimum,
    )
    edit_tokens = fallback_edit_index_tokens(
        descriptor,
        require_eligible=not unrestricted,
    )
    partial_transformation_key = str(
        descriptor.get("partial_transformation_key") or ""
    )
    edit_position_sets = [
        set(index.fallback_features.get(token, ())) for token in edit_tokens
    ]
    if partial_transformation_key and not unrestricted:
        positions = set(
            index.partial_transformations.get(partial_transformation_key, ())
        )
        candidate_level = "partial_transformation_exact"
    elif edit_position_sets and not unrestricted:
        positions = set.intersection(*edit_position_sets)
        candidate_level = "fallback_edit_candidates"
    else:
        query_tokens = fallback_index_tokens(
            descriptor,
            require_eligible=not unrestricted,
        )
        positions = set()
        for token in query_tokens:
            positions.update(index.fallback_features.get(token, ()))
        candidate_level = (
            "unrestricted_fallback_descriptor_candidates"
            if unrestricted
            else "fallback_descriptor_candidates"
        )
    raw_rows = tuple(
        index.rows[position]
        for position in sorted(positions)
        if index.rows[position].fallback_descriptor
    )
    required_source_ids = {
        str(value.get("requirement_id") or "")
        for value in descriptor.get("source_requirements") or ()
        if isinstance(value, Mapping)
    }
    capability_constraints = dict(required_source_capability_ids or {})
    substance_constraints = dict(required_source_substance_ids or {})

    def supports_selected_sources(row: Any) -> bool:
        supports = {
            str(value.get("requirement_id") or ""): value
            for value in row.fragment_source_support
            if str(value.get("status") or "") == "supported"
        }
        if required_source_ids and not unrestricted and not (
            required_source_ids <= set(supports)
        ):
            return False
        for requirement_id, capability_id in capability_constraints.items():
            support = supports.get(requirement_id, {})
            if capability_id not in {
                str(value) for value in support.get("capability_ids") or ()
            }:
                return False
        for requirement_id, substance_id in substance_constraints.items():
            support = supports.get(requirement_id, {})
            if substance_id not in {
                str(value)
                for value in support.get("component_substance_ids") or ()
            }:
                return False
        return True

    source_excluded_count = 0
    if (
        (required_source_ids and not unrestricted)
        or capability_constraints
        or substance_constraints
    ):
        source_supported_rows = tuple(
            row
            for row in raw_rows
            if supports_selected_sources(row)
        )
        source_excluded_count = len(raw_rows) - len(source_supported_rows)
    else:
        source_supported_rows = raw_rows
    query_edits = Counter(
        str(value) for value in descriptor.get("verified_edit_tokens") or ()
    )
    if query_edits and not unrestricted:
        source_supported_rows = tuple(
            row
            for row in source_supported_rows
            if Counter(
                str(value)
                for value in row.fallback_descriptor.get("verified_edit_tokens") or ()
            )
            == query_edits
        )
    raw_support = summarize_evidence_support(raw_rows)
    source_supported_support = summarize_evidence_support(source_supported_rows)
    candidate_trace = RetrievalLevelTrace(
        level=candidate_level,
        candidate_count=len(raw_rows),
        independent_candidate_count=raw_support.independent_count,
        compatible_candidate_count=len(source_supported_rows),
        independent_compatible_candidate_count=(
            source_supported_support.independent_count
        ),
        excluded_candidate_count=source_excluded_count,
        minimum_independent_support=minimum,
        status=(
            "source_supported"
            if source_supported_rows
            else "missing_condition_fragment_source"
            if raw_rows
            else "empty"
        ),
    )
    if not source_supported_rows:
        return CompatibleRetrievalResult(
            level=(
                "no_condition_fragment_source_precedent"
                if raw_rows
                else "no_fallback_descriptor_candidate"
            ),
            pool=(),
            candidate_count=len(raw_rows),
            independent_candidate_count=raw_support.independent_count,
            excluded_candidate_count=source_excluded_count,
            independent_compatible_candidate_count=0,
            trace=(candidate_trace,),
        )
    raw_rows = source_supported_rows

    if unrestricted:
        candidate_limit = int(rules["candidate_limit"])
        scored = [
            (
                assess_fallback_similarity(
                    descriptor,
                    row.fallback_descriptor,
                ).score,
                shared_high_signal_count(
                    descriptor,
                    row.fallback_descriptor,
                ),
                row,
            )
            for row in raw_rows
        ]
        scored.sort(
            key=lambda item: (
                -item[0],
                -item[1],
                item[2].canonical_reaction_id,
                item[2].reaction_id,
                item[2].observation_id,
            )
        )
        matched_rows = tuple(item[2] for item in scored[:candidate_limit])
        neutral_compatibility = CompatibilityAssessment(
            compatible=True,
            score=1.0,
            evidence=(
                "Condition compatibility gates were explicitly bypassed",
            ),
        )
        accepted = tuple(
            (row, neutral_compatibility) for row in matched_rows
        )
        support = summarize_evidence_support(matched_rows)
        trace = RetrievalLevelTrace(
            level="unrestricted_unverified_structure_similarity",
            candidate_count=len(matched_rows),
            independent_candidate_count=support.independent_count,
            compatible_candidate_count=len(matched_rows),
            independent_compatible_candidate_count=support.independent_count,
            excluded_candidate_count=0,
            minimum_independent_support=0,
            status="selected_without_gates",
        )
        return CompatibleRetrievalResult(
            level="unrestricted_unverified_structure_fallback",
            pool=accepted,
            candidate_count=len(matched_rows),
            independent_candidate_count=support.independent_count,
            excluded_candidate_count=0,
            independent_compatible_candidate_count=support.independent_count,
            trace=(candidate_trace, trace),
        )

    exploratory_policy = rules["exploratory_partial_correspondence"]
    exploratory = bool(exploratory_policy.get("enabled")) and (
        descriptor.get("evidence_mode") == "partial_product_correspondence"
        and bool(edit_tokens)
    )
    minimum_similarity = (
        float(exploratory_policy["minimum_similarity"])
        if exploratory
        else float(rules["minimum_similarity"])
    )
    minimum_shared = (
        1 if exploratory else int(rules["minimum_shared_high_signal_features"])
    )
    candidate_limit = (
        int(exploratory_policy["candidate_limit"])
        if exploratory
        else int(rules["candidate_limit"])
    )
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
    matched_rows = tuple(item[2] for item in scored[:candidate_limit])
    compatibility_signature = compatibility_signature_from_fallback(descriptor)
    accepted, compatibility_excluded = filter_compatible_precedents(
        compatibility_signature,
        matched_rows,
    )
    excluded_count = len(compatibility_excluded) + source_excluded_count
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

    selected_level = (
        "source_supported_partial_transformation"
        if candidate_level == "partial_transformation_exact"
        else "unverified_structure_fallback_exploratory"
        if exploratory
        else "unverified_structure_fallback"
    )
    status = (
        "selected_source_supported_partial"
        if candidate_level == "partial_transformation_exact"
        else "selected_exploratory"
        if exploratory
        else "selected"
    )
    if accepted_support.independent_count < minimum:
        if exploratory and bool(
            exploratory_policy.get("support_affects_confidence_only")
        ):
            selected_level += "_limited_support"
            status = (
                "selected_source_supported_partial_limited_support"
                if candidate_level == "partial_transformation_exact"
                else "selected_exploratory_limited_support"
            )
        else:
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
