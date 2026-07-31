"""Conservative query-only retrieval by minimized reaction-core shape."""

from __future__ import annotations

from dataclasses import dataclass
import json
from functools import lru_cache
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
    / "reaction_core_retrieval.v1.json"
)


@dataclass(frozen=True)
class CoreRetrievalResult:
    """Compatibility-filtered verified precedents for one core-only query."""

    level: str
    pool: Tuple[
        Tuple[GenericIndexedReaction, CompatibilityAssessment],
        ...,
    ]
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
    if str(rules.get("definition_id") or "") != "reaction_core_retrieval.v1":
        raise ValueError("unexpected reaction-core retrieval definition ID")
    if not str(rules.get("calibration_status") or "").strip():
        raise ValueError("reaction-core retrieval requires calibration status")
    if int(rules.get("minimum_pool_size") or 0) < 1:
        raise ValueError("reaction-core retrieval minimum support must be positive")
    if not tuple(rules.get("allowed_query_evidence_statuses") or ()):
        raise ValueError("reaction-core retrieval requires allowed evidence")
    if rules.get("verified_precedents_only") is not True:
        raise ValueError("reaction-core retrieval must use verified precedents")
    return rules


def reaction_core_query_eligible(
    core: Mapping[str, Any],
    index: GenericReactionIndex,
) -> tuple[bool, str]:
    """Return whether a minimized core may query the verified core index."""
    rules = load_reaction_core_retrieval_rules()
    if str(core.get("schema_version") or "") != (
        REACTION_CORE_PROJECTION_SCHEMA_VERSION
    ):
        return False, "incompatible_reaction_core_schema"
    if str(core.get("algorithm_version") or "") != (
        REACTION_CORE_PROJECTION_ALGORITHM_VERSION
    ):
        return False, "incompatible_reaction_core_algorithm"
    if str(core.get("schema_version") or "") != index.reaction_core_schema_version:
        return False, "incompatible_reaction_core_index_schema"
    if str(core.get("algorithm_version") or "") != (
        index.reaction_core_algorithm_version
    ):
        return False, "incompatible_reaction_core_index_algorithm"
    if str(core.get("evidence_status") or "") not in set(
        rules["allowed_query_evidence_statuses"]
    ):
        return False, "reaction_core_evidence_not_query_eligible"
    if set(core.get("warnings") or ()).intersection(
        rules["blocking_query_warnings"]
    ):
        return False, "reaction_core_has_blocking_warning"
    if not str(core.get("shape_core_key") or ""):
        return False, "reaction_core_shape_key_missing"
    if int(core.get("event_count") or 0) < 1:
        return False, "reaction_core_event_missing"
    return True, ""


def retrieve_core_shape_pool_with_trace(
    core: Mapping[str, Any],
    compatibility_signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
) -> CoreRetrievalResult:
    """Retrieve only verified precedents with the exact robust core-shape key."""
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
            level="reaction_core_shape",
            candidate_count=0,
            independent_candidate_count=0,
            compatible_candidate_count=0,
            independent_compatible_candidate_count=0,
            excluded_candidate_count=0,
            minimum_independent_support=minimum,
            status=reason,
        )
        return CoreRetrievalResult(
            reason,
            (),
            0,
            0,
            0,
            0,
            (trace,),
        )
    shape_key = str(core["shape_core_key"])
    event_count = int(core["event_count"])
    positions = {
        position
        for position in index.core_shapes.get(shape_key, ())
        if index.rows[position].signature
        and index.rows[position].reaction_core
        and int(index.rows[position].reaction_core.get("event_count") or 0)
        == event_count
    }
    rows = index.select(sorted(positions))
    raw_support = summarize_evidence_support(rows)
    accepted, excluded = filter_compatible_precedents(
        compatibility_signature,
        rows,
    )
    accepted_rows = tuple(row for row, _ in accepted)
    accepted_support = summarize_evidence_support(accepted_rows)
    if not rows:
        status = "empty"
        level = "no_reaction_core_shape_precedent"
    elif not accepted:
        status = "no_compatible_recipe"
        level = "no_compatible_condition_precedent"
    elif accepted_support.independent_count >= minimum:
        status = "selected"
        level = "reaction_core_shape"
    else:
        status = "selected_limited_support"
        level = "reaction_core_shape_limited_support"
    trace = RetrievalLevelTrace(
        level="reaction_core_shape",
        candidate_count=len(rows),
        independent_candidate_count=raw_support.independent_count,
        compatible_candidate_count=len(accepted),
        independent_compatible_candidate_count=accepted_support.independent_count,
        excluded_candidate_count=len(excluded),
        minimum_independent_support=minimum,
        status=status,
    )
    return CoreRetrievalResult(
        level=level,
        pool=accepted,
        candidate_count=len(rows),
        independent_candidate_count=raw_support.independent_count,
        excluded_candidate_count=len(excluded),
        independent_compatible_candidate_count=accepted_support.independent_count,
        trace=(trace,),
    )


__all__ = [
    "CoreRetrievalResult",
    "load_reaction_core_retrieval_rules",
    "reaction_core_query_eligible",
    "retrieve_core_shape_pool_with_trace",
]
