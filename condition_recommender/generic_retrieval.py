"""Type-agnostic hierarchical retrieval over generic reaction signatures."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Literal, Mapping, Tuple

from .compatibility import CompatibilityAssessment, filter_compatible_precedents
from .core_retrieval import reaction_core_query_eligible
from .edit_prototypes import (
    anonymous_edit_prototype,
    anonymous_edit_similarity,
)
from .generic_indexing import GenericIndexedReaction, GenericReactionIndex
from .models import RetrievalLevelTrace
from .signature_features import (
    environment_profile_similarity,
    environment_tokens,
)
from .similarity import generic_signature_similarity, reaction_scope
from .support import summarize_evidence_support

_RULES_PATH = Path(__file__).with_name("definitions") / "generic_retrieval.v1.json"
_SUPPORTED_RETRIEVAL_LEVELS = {
    "exact_signature",
    "handle_signature",
    "named_family",
    "transformation_signature",
    "environment_neighbors",
    "bond_edit_signature",
    "reaction_core_shape",
    "edit_graph_neighbors",
}
RetrievalStrategy = Literal[
    "hybrid",
    "family_only",
    "generic_only",
    "transformation_prior",
    "legacy_pilot",
]
_EVALUATION_STRATEGIES = {
    "family_only",
    "generic_only",
    "transformation_prior",
    "legacy_pilot",
}


@lru_cache(maxsize=1)
def load_generic_retrieval_rules() -> Dict[str, Any]:
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "1.7":
        raise ValueError("unsupported generic retrieval definition schema")
    if str(rules.get("definition_id") or "") != "generic_retrieval.v1":
        raise ValueError("unexpected generic retrieval definition ID")
    if not str(rules.get("calibration_status") or "").strip():
        raise ValueError("generic retrieval definition requires calibration status")
    ladder = tuple(str(value) for value in rules.get("retrieval_ladder") or ())
    if (
        not ladder
        or len(set(ladder)) != len(ladder)
        or set(ladder) != _SUPPORTED_RETRIEVAL_LEVELS
    ):
        raise ValueError("generic retrieval ladder is incomplete or invalid")
    evaluation_ladders = rules.get("evaluation_ladders")
    if not isinstance(evaluation_ladders, Mapping) or set(
        evaluation_ladders
    ) != _EVALUATION_STRATEGIES:
        raise ValueError("generic retrieval evaluation ladders are incomplete")
    for strategy, values in evaluation_ladders.items():
        strategy_ladder = tuple(str(value) for value in values or ())
        if (
            not strategy_ladder
            or len(set(strategy_ladder)) != len(strategy_ladder)
            or not set(strategy_ladder) <= _SUPPORTED_RETRIEVAL_LEVELS
        ):
            raise ValueError(
                f"generic retrieval evaluation ladder is invalid: {strategy}"
            )
    if int(rules["environment_neighbor_limit"]) < 1:
        raise ValueError("environment_neighbor_limit must be positive")
    environment_threshold = float(rules["environment_neighbor_min_similarity"])
    if not 0.0 < environment_threshold <= 1.0:
        raise ValueError(
            "environment_neighbor_min_similarity must be in (0, 1]"
        )
    if int(rules["edit_graph_neighbor_limit"]) < 1:
        raise ValueError("edit_graph_neighbor_limit must be positive")
    edit_threshold = float(rules["edit_graph_neighbor_min_similarity"])
    if not 0.0 < edit_threshold <= 1.0:
        raise ValueError(
            "edit_graph_neighbor_min_similarity must be in (0, 1]"
        )
    return rules


def _positions(mapping: Mapping[str, Tuple[int, ...]], key: Any) -> set[int]:
    return set(mapping.get(str(key or ""), ()))


def _compatible_edit_positions(
    signature: Mapping[str, Any], index: GenericReactionIndex
) -> set[int]:
    """Enforce the net bond-edit compatibility gate before similarity."""
    return _positions(index.bond_edits, signature.get("bond_edit_signature_key"))


def _candidate_levels(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
) -> list[tuple[str, set[int]]]:
    compatible = _compatible_edit_positions(signature, index)
    edit_graph_neighbors = _edit_graph_neighbor_positions(
        signature,
        index,
        exclude=compatible,
    )
    core_shape_positions = _core_shape_positions(reaction_core, index)
    if not compatible and not edit_graph_neighbors and not core_shape_positions:
        return []
    rules = load_generic_retrieval_rules()
    family = str(signature.get("named_family") or "")
    family_confidence = float(signature.get("family_confidence") or 0.0)
    threshold = float(rules["high_confidence_family_threshold"])
    family_positions = (
        _positions(index.families, family) & compatible
        if family and family_confidence >= threshold
        else set()
    )
    candidates = {
        "exact_signature": (
            _positions(index.exact, signature.get("exact_signature_key"))
            & compatible
        ),
        "handle_signature": (
            _positions(index.handles, signature.get("handle_signature_key"))
            & compatible
        ),
        "named_family": family_positions,
        "transformation_signature": (
            _positions(
                index.transformations,
                signature.get("transformation_signature_key"),
            )
            & compatible
        ),
        "environment_neighbors": _environment_neighbor_positions(
            signature, index, compatible
        ),
        "bond_edit_signature": compatible,
        "reaction_core_shape": core_shape_positions,
        "edit_graph_neighbors": edit_graph_neighbors,
    }
    ladder = (
        rules["retrieval_ladder"]
        if strategy == "hybrid"
        else (rules["evaluation_ladders"] or {}).get(strategy)
    )
    if not ladder:
        raise ValueError(f"Unsupported generic retrieval strategy: {strategy}")
    return [(level, candidates[level]) for level in ladder]


def _core_shape_positions(
    reaction_core: Mapping[str, Any] | None,
    index: GenericReactionIndex,
) -> set[int]:
    """Return verified precedents sharing the mapping-robust core shape."""
    if not reaction_core:
        return set()
    eligible, _ = reaction_core_query_eligible(reaction_core, index)
    if not eligible:
        return set()
    shape_key = str(reaction_core.get("shape_core_key") or "")
    if not shape_key:
        return set()
    event_count = int(reaction_core.get("event_count") or 0)
    return {
        position
        for position in index.core_shapes.get(shape_key, ())
        if index.rows[position].signature
        and index.rows[position].reaction_core
        and int(index.rows[position].reaction_core.get("event_count") or 0)
        == event_count
    }


def _edit_graph_neighbor_positions(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    exclude: set[int],
) -> set[int]:
    """Find chemistry-gated approximate edit graphs without family routing."""
    query = anonymous_edit_prototype(signature)
    if query is None:
        return set()
    rules = load_generic_retrieval_rules()
    threshold = float(rules["edit_graph_neighbor_min_similarity"])
    limit = int(rules["edit_graph_neighbor_limit"])
    scored = []
    for position, row in enumerate(index.rows):
        if position in exclude or not row.signature:
            continue
        candidate = anonymous_edit_prototype(row.signature)
        if candidate is None:
            continue
        score = anonymous_edit_similarity(query, candidate)
        if score >= threshold:
            scored.append(
                (
                    score,
                    row.canonical_reaction_id,
                    row.reaction_id,
                    row.observation_id,
                    position,
                )
            )
    scored.sort(key=lambda item: (-item[0],) + item[1:])
    return {item[-1] for item in scored[:limit]}


def _environment_neighbor_positions(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    compatible: set[int],
) -> set[int]:
    """Narrow edit-compatible rows using interpretable local environments."""
    query_tokens = environment_tokens(signature)
    if not query_tokens:
        return set()
    rules = load_generic_retrieval_rules()
    threshold = float(rules["environment_neighbor_min_similarity"])
    limit = int(rules["environment_neighbor_limit"])
    token_positions = set()
    for token in query_tokens:
        token_positions.update(index.environment_features.get(token, ()))
    scored = []
    for position in compatible & token_positions:
        score = environment_profile_similarity(
            signature,
            index.rows[position].signature,
        )
        if score >= threshold:
            scored.append((score, position))
    scored.sort(
        key=lambda item: (
            -item[0],
            index.rows[item[1]].canonical_reaction_id,
            index.rows[item[1]].reaction_id,
            index.rows[item[1]].observation_id,
        )
    )
    return {position for _, position in scored[:limit]}


def _minimum_support(minimum_pool_size: int | None) -> int:
    rules = load_generic_retrieval_rules()
    minimum = (
        int(minimum_pool_size)
        if minimum_pool_size is not None
        else int(rules["minimum_pool_size"])
    )
    if minimum < 1:
        raise ValueError("minimum_pool_size must be positive")
    return minimum


def retrieve_generic_pool_with_trace(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
) -> tuple[
    str,
    Tuple[GenericIndexedReaction, ...],
    Tuple[RetrievalLevelTrace, ...],
]:
    """Select by independent support and retain every attempted tier."""
    minimum = _minimum_support(minimum_pool_size)
    levels = _candidate_levels(
        signature,
        index,
        reaction_core=reaction_core,
        strategy=strategy,
    )
    if not levels:
        trace = RetrievalLevelTrace(
            level="bond_edit_gate",
            candidate_count=0,
            independent_candidate_count=0,
            compatible_candidate_count=0,
            independent_compatible_candidate_count=0,
            excluded_candidate_count=0,
            minimum_independent_support=minimum,
            status="no_compatible_bond_edit",
        )
        return "no_compatible_bond_edit", (), (trace,)
    fallback: tuple[str, set[int], int] | None = None
    traces = []
    for level, positions in levels:
        if not positions:
            traces.append(
                RetrievalLevelTrace(
                    level=level,
                    candidate_count=0,
                    independent_candidate_count=0,
                    compatible_candidate_count=0,
                    independent_compatible_candidate_count=0,
                    excluded_candidate_count=0,
                    minimum_independent_support=minimum,
                    status="empty",
                )
            )
            continue
        rows = index.select(sorted(positions))
        support = summarize_evidence_support(rows)
        trace = RetrievalLevelTrace(
            level=level,
            candidate_count=len(rows),
            independent_candidate_count=support.independent_count,
            compatible_candidate_count=len(rows),
            independent_compatible_candidate_count=support.independent_count,
            excluded_candidate_count=0,
            minimum_independent_support=minimum,
            status="insufficient_independent_support",
        )
        traces.append(trace)
        if fallback is None or (
            support.independent_count,
            len(rows),
        ) > (
            summarize_evidence_support(
                index.select(sorted(fallback[1]))
            ).independent_count,
            len(fallback[1]),
        ):
            fallback = (level, positions, len(traces) - 1)
        if support.independent_count >= minimum:
            traces[-1] = replace(trace, status="selected")
            return level, rows, tuple(traces)
    if fallback is None:
        return "no_compatible_precedent", (), tuple(traces)
    level, positions, trace_index = fallback
    traces[trace_index] = replace(
        traces[trace_index], status="selected_limited_support"
    )
    return (
        f"{level}_limited_support",
        index.select(sorted(positions)),
        tuple(traces),
    )


def retrieve_generic_pool(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
) -> Tuple[str, Tuple[GenericIndexedReaction, ...]]:
    """Compatibility wrapper returning the historical two-value result."""
    level, rows, _ = retrieve_generic_pool_with_trace(
        signature,
        index,
        minimum_pool_size=minimum_pool_size,
        reaction_core=reaction_core,
        strategy=strategy,
    )
    return level, rows


@dataclass(frozen=True)
class CompatibleRetrievalResult:
    """Selected compatibility-filtered pool and its complete retrieval trace."""

    level: str
    pool: tuple[tuple[GenericIndexedReaction, CompatibilityAssessment], ...]
    candidate_count: int
    independent_candidate_count: int
    excluded_candidate_count: int
    independent_compatible_candidate_count: int
    trace: Tuple[RetrievalLevelTrace, ...]


def retrieve_compatible_generic_pool_with_trace(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
) -> CompatibleRetrievalResult:
    """Apply compatibility before independent-support checks at every tier."""
    minimum = _minimum_support(minimum_pool_size)
    levels = _candidate_levels(
        signature,
        index,
        reaction_core=reaction_core,
        strategy=strategy,
    )
    if not levels:
        trace = RetrievalLevelTrace(
            level="bond_edit_gate",
            candidate_count=0,
            independent_candidate_count=0,
            compatible_candidate_count=0,
            independent_compatible_candidate_count=0,
            excluded_candidate_count=0,
            minimum_independent_support=minimum,
            status="no_compatible_bond_edit",
        )
        return CompatibleRetrievalResult(
            "no_compatible_bond_edit", (), 0, 0, 0, 0, (trace,)
        )
    fallback = None
    traces = []
    for level, positions in levels:
        if not positions:
            traces.append(
                RetrievalLevelTrace(
                    level=level,
                    candidate_count=0,
                    independent_candidate_count=0,
                    compatible_candidate_count=0,
                    independent_compatible_candidate_count=0,
                    excluded_candidate_count=0,
                    minimum_independent_support=minimum,
                    status="empty",
                )
            )
            continue
        rows = index.select(sorted(positions))
        accepted, excluded = filter_compatible_precedents(signature, rows)
        raw_support = summarize_evidence_support(rows)
        accepted_rows = tuple(row for row, _ in accepted)
        accepted_support = summarize_evidence_support(accepted_rows)
        status = (
            "insufficient_independent_support"
            if accepted
            else "no_compatible_recipe"
        )
        trace = RetrievalLevelTrace(
            level=level,
            candidate_count=len(rows),
            independent_candidate_count=raw_support.independent_count,
            compatible_candidate_count=len(accepted),
            independent_compatible_candidate_count=(
                accepted_support.independent_count
            ),
            excluded_candidate_count=len(excluded),
            minimum_independent_support=minimum,
            status=status,
        )
        traces.append(trace)
        if not accepted:
            continue
        candidate = (
            level,
            accepted,
            len(rows),
            raw_support.independent_count,
            len(excluded),
            accepted_support.independent_count,
            len(traces) - 1,
        )
        if fallback is None or (
            accepted_support.independent_count,
            len(accepted),
        ) > (
            fallback[5],
            len(fallback[1]),
        ):
            fallback = candidate
        if accepted_support.independent_count >= minimum:
            traces[-1] = replace(trace, status="selected")
            return CompatibleRetrievalResult(
                level=level,
                pool=accepted,
                candidate_count=len(rows),
                independent_candidate_count=raw_support.independent_count,
                excluded_candidate_count=len(excluded),
                independent_compatible_candidate_count=(
                    accepted_support.independent_count
                ),
                trace=tuple(traces),
            )
    if fallback is None:
        bond_rows = index.select(sorted(_compatible_edit_positions(signature, index)))
        _, excluded = filter_compatible_precedents(signature, bond_rows)
        support = summarize_evidence_support(bond_rows)
        return CompatibleRetrievalResult(
            "no_compatible_condition_precedent",
            (),
            len(bond_rows),
            support.independent_count,
            len(excluded),
            0,
            tuple(traces),
        )
    (
        level,
        accepted,
        raw_count,
        independent_raw_count,
        excluded_count,
        independent_accepted_count,
        trace_index,
    ) = fallback
    traces[trace_index] = replace(
        traces[trace_index], status="selected_limited_support"
    )
    return CompatibleRetrievalResult(
        f"{level}_limited_support",
        accepted,
        raw_count,
        independent_raw_count,
        excluded_count,
        independent_accepted_count,
        tuple(traces),
    )


def retrieve_compatible_generic_pool(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
    reaction_core: Mapping[str, Any] | None = None,
    strategy: RetrievalStrategy = "hybrid",
) -> tuple[
    str,
    tuple[tuple[GenericIndexedReaction, CompatibilityAssessment], ...],
    int,
    int,
]:
    """Compatibility wrapper returning the historical four-value result."""
    result = retrieve_compatible_generic_pool_with_trace(
        signature,
        index,
        minimum_pool_size=minimum_pool_size,
        reaction_core=reaction_core,
        strategy=strategy,
    )
    return (
        result.level,
        result.pool,
        result.candidate_count,
        result.excluded_candidate_count,
    )


__all__ = [
    "generic_signature_similarity",
    "CompatibleRetrievalResult",
    "RetrievalStrategy",
    "load_generic_retrieval_rules",
    "reaction_scope",
    "retrieve_compatible_generic_pool",
    "retrieve_compatible_generic_pool_with_trace",
    "retrieve_generic_pool",
    "retrieve_generic_pool_with_trace",
]
