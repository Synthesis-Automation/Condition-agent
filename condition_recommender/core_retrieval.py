"""Conservative query-only retrieval by minimized reaction-core identity."""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Mapping, Tuple

from reactive_taxonomy import (
    REACTION_CORE_PROJECTION_ALGORITHM_VERSION,
    REACTION_CORE_PROJECTION_SCHEMA_VERSION,
)

from .compatibility import CompatibilityAssessment, filter_compatible_precedents
from .generic_indexing import GenericIndexedReaction, GenericReactionIndex
from .models import RetrievalLevelTrace
from .support import summarize_evidence_support


_RULES_PATH = (
    Path(__file__).with_name("definitions")
    / "reaction_core_retrieval.v3.json"
)


@dataclass(frozen=True)
class CoreRetrievalResult:
    """Compatibility-filtered tier-qualified precedents for a core query."""

    level: str
    pool: Tuple[Tuple[GenericIndexedReaction, CompatibilityAssessment], ...]
    candidate_count: int
    independent_candidate_count: int
    excluded_candidate_count: int
    independent_compatible_candidate_count: int
    trace: Tuple[RetrievalLevelTrace, ...]


@lru_cache(maxsize=1)
def load_reaction_core_retrieval_rules() -> dict[str, Any]:
    """Load and validate the core-only retrieval trust policy."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "1.0":
        raise ValueError("unsupported reaction-core retrieval schema")
    if str(rules.get("definition_id") or "") != "reaction_core_retrieval.v3":
        raise ValueError("unexpected reaction-core retrieval definition ID")
    ladder = tuple(rules.get("retrieval_ladder") or ())
    expected = (
        {"level": "reaction_core_exact", "key_field": "exact_core_key", "index_map": "core_exact"},
        {"level": "reaction_core_local", "key_field": "typed_core_key", "index_map": "core_typed"},
        {"level": "reaction_core_context", "key_field": "shape_core_key", "index_map": "core_shapes"},
    )
    if ladder != expected:
        raise ValueError("reaction-core retrieval ladder is invalid")
    if not str(rules.get("calibration_status") or "").strip():
        raise ValueError("reaction-core retrieval requires calibration status")
    if int(rules.get("minimum_pool_size") or 0) < 1:
        raise ValueError("reaction-core retrieval minimum support must be positive")
    if not tuple(rules.get("allowed_query_evidence_statuses") or ()):
        raise ValueError("reaction-core retrieval requires allowed evidence")
    if tuple(rules.get("allowed_precedent_tiers") or ()) != (
        "trusted",
        "review_core",
    ):
        raise ValueError("reaction-core retrieval precedent tiers are invalid")
    if rules.get("review_core_requires_expert_mode") is not True:
        raise ValueError("review-core retrieval must require expert mode")
    return rules


def reaction_core_query_eligible(
    core: Mapping[str, Any],
    index: GenericReactionIndex,
) -> tuple[bool, str]:
    """Return whether a minimized core may query the verified core index."""
    rules = load_reaction_core_retrieval_rules()
    if str(core.get("schema_version") or "") != REACTION_CORE_PROJECTION_SCHEMA_VERSION:
        return False, "incompatible_reaction_core_schema"
    if str(core.get("algorithm_version") or "") != REACTION_CORE_PROJECTION_ALGORITHM_VERSION:
        return False, "incompatible_reaction_core_algorithm"
    if str(core.get("schema_version") or "") != index.reaction_core_schema_version:
        return False, "incompatible_reaction_core_index_schema"
    if str(core.get("algorithm_version") or "") != index.reaction_core_algorithm_version:
        return False, "incompatible_reaction_core_index_algorithm"
    if str(core.get("evidence_status") or "") not in set(rules["allowed_query_evidence_statuses"]):
        return False, "reaction_core_evidence_not_query_eligible"
    quality = core.get("quality")
    quality_status = (
        str(quality.get("status") or "")
        if isinstance(quality, Mapping)
        else ""
    )
    if quality_status == "blocked":
        return False, "reaction_core_quality_blocked"
    if set(core.get("warnings") or ()).intersection(rules["blocking_query_warnings"]):
        return False, "reaction_core_has_blocking_warning"
    for field in ("exact_core_key", "typed_core_key", "shape_core_key"):
        if not str(core.get(field) or ""):
            return False, f"{field}_missing"
    if int(core.get("event_count") or 0) < 1:
        return False, "reaction_core_event_missing"
    return True, ""


def _level_positions(
    core: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    key_field: str,
    index_map: str,
) -> set[int]:
    rules = load_reaction_core_retrieval_rules()
    allowed_tiers = set(rules["allowed_precedent_tiers"])
    key = str(core.get(key_field) or "")
    event_count = int(core.get("event_count") or 0)
    mapping = getattr(index, index_map)
    return {
        position
        for position in mapping.get(key, ())
        if index.rows[position].precedent_tier.value in allowed_tiers
        and index.rows[position].reaction_core
        and (
            index.rows[position].signature
            or index.rows[position].fallback_descriptor
        )
        and int(index.rows[position].reaction_core.get("event_count") or 0)
        == event_count
    }


def retrieve_core_pool_with_trace(
    core: Mapping[str, Any],
    compatibility_signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
) -> CoreRetrievalResult:
    """Retrieve the narrowest adequately supported exact-to-context core pool."""
    rules = load_reaction_core_retrieval_rules()
    minimum = (
        int(minimum_pool_size)
        if minimum_pool_size is not None
        else int(rules["minimum_pool_size"])
    )
    if minimum < 1:
        raise ValueError("minimum_pool_size must be positive")
    eligible, reason = reaction_core_query_eligible(core, index)
    if not eligible:
        trace = RetrievalLevelTrace(
            level="reaction_core_ladder",
            candidate_count=0,
            independent_candidate_count=0,
            compatible_candidate_count=0,
            independent_compatible_candidate_count=0,
            excluded_candidate_count=0,
            minimum_independent_support=minimum,
            status=reason,
        )
        return CoreRetrievalResult(reason, (), 0, 0, 0, 0, (trace,))

    trace = []
    limited: CoreRetrievalResult | None = None
    for definition in rules["retrieval_ladder"]:
        level = str(definition["level"])
        positions = _level_positions(
            core,
            index,
            key_field=str(definition["key_field"]),
            index_map=str(definition["index_map"]),
        )
        rows = index.select(sorted(positions))
        raw_support = summarize_evidence_support(rows)
        accepted, excluded = filter_compatible_precedents(compatibility_signature, rows)
        accepted_support = summarize_evidence_support(tuple(row for row, _ in accepted))
        if not rows:
            status = "empty"
        elif not accepted:
            status = "no_compatible_recipe"
        elif accepted_support.independent_count >= minimum:
            status = "selected"
        else:
            status = "selected_limited_support"
        trace.append(
            RetrievalLevelTrace(
                level=level,
                candidate_count=len(rows),
                independent_candidate_count=raw_support.independent_count,
                compatible_candidate_count=len(accepted),
                independent_compatible_candidate_count=accepted_support.independent_count,
                excluded_candidate_count=len(excluded),
                minimum_independent_support=minimum,
                status=status,
            )
        )
        result = CoreRetrievalResult(
            level=level if status == "selected" else f"{level}_limited_support",
            pool=accepted,
            candidate_count=len(rows),
            independent_candidate_count=raw_support.independent_count,
            excluded_candidate_count=len(excluded),
            independent_compatible_candidate_count=accepted_support.independent_count,
            trace=tuple(trace),
        )
        if status == "selected":
            return result
        if accepted and limited is None:
            limited = result

    if limited is not None:
        return replace(limited, trace=tuple(trace))
    return CoreRetrievalResult(
        "no_reaction_core_precedent",
        (),
        0,
        0,
        0,
        0,
        tuple(trace),
    )


__all__ = [
    "CoreRetrievalResult",
    "load_reaction_core_retrieval_rules",
    "reaction_core_query_eligible",
    "retrieve_core_pool_with_trace",
]
