"""Validated, versioned policies for general retrosynthesis ranking."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Tuple

from .generic_models import GenericDisconnectionCandidate


RANKING_POLICY_PATH = (
    Path(__file__).with_name("definitions")
    / "retrosynthesis_ranking.v2.json"
)
_ALLOWED_GROUP_FIELDS = frozenset(
    {
        "operator_id",
        "disconnection_site_key",
        "synthon_signature",
        "realization_id",
    }
)
_CONDITION_STATUSES = (
    "recommended_direct",
    "recommended_fallback",
    "insufficient_evidence",
)


@dataclass(frozen=True)
class RetrosynthesisRankingPolicy:
    """One immutable structural-diversity and condition-ranking policy."""

    definition_id: str
    schema_version: str
    diversity_group_fields: Tuple[str, ...]
    candidate_pool_multiplier: int
    preserve_abstraction_level_order: bool
    diversity_score_band_width: float
    condition_score_band_width: float
    abstraction_level_order: Tuple[str, ...]
    condition_status_order: Tuple[str, ...]
    precursor_realism_band_penalties: Tuple[Tuple[float, int], ...]

    def level_rank(self, level: str) -> int:
        """Return the configured specificity rank for one abstraction level."""

        try:
            return self.abstraction_level_order.index(level)
        except ValueError:
            return len(self.abstraction_level_order)

    def condition_status_rank(self, status: str) -> int:
        """Return the configured evidence rank for one condition status."""

        try:
            return self.condition_status_order.index(status)
        except ValueError:
            return len(self.condition_status_order)

    def precursor_realism_band_penalty(self, score: float) -> int:
        """Return the configured demotion for one precursor-realism score."""

        if not math.isfinite(score) or not 0.0 <= score <= 1.0:
            raise ValueError("precursor realism score must be between zero and one")
        for minimum_score, penalty in self.precursor_realism_band_penalties:
            if score >= minimum_score:
                return penalty
        raise ValueError("precursor realism penalty policy does not cover score")


def _positive_float(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must be a positive number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized <= 0.0:
        raise ValueError(f"{field} must be a positive number")
    return normalized


@lru_cache(maxsize=1)
def load_retrosynthesis_ranking_policy() -> RetrosynthesisRankingPolicy:
    """Load and validate the canonical general ranking policy."""

    value = json.loads(RANKING_POLICY_PATH.read_text(encoding="utf-8"))
    if value.get("definition_id") != "retrosynthesis_ranking.v2":
        raise ValueError("unexpected retrosynthesis ranking definition ID")
    if value.get("schema_version") != "2.0":
        raise ValueError("unsupported retrosynthesis ranking schema")
    diversity = value.get("candidate_diversity")
    condition = value.get("condition_reranking")
    realism = value.get("precursor_realism_reranking")
    if not all(
        isinstance(section, dict)
        for section in (diversity, condition, realism)
    ):
        raise ValueError(
            "ranking policy requires diversity, condition, and realism rules"
        )
    group_fields = tuple(diversity.get("group_fields") or ())
    if (
        not group_fields
        or len(set(group_fields)) != len(group_fields)
        or any(field not in _ALLOWED_GROUP_FIELDS for field in group_fields)
    ):
        raise ValueError("ranking policy has invalid diversity group fields")
    pool_multiplier = diversity.get("candidate_pool_multiplier")
    if (
        isinstance(pool_multiplier, bool)
        or not isinstance(pool_multiplier, int)
        or pool_multiplier < 1
    ):
        raise ValueError("candidate pool multiplier must be a positive integer")
    preserve_levels = diversity.get("preserve_abstraction_level_order")
    if not isinstance(preserve_levels, bool):
        raise ValueError("abstraction-level preservation must be boolean")
    level_order = tuple(value.get("abstraction_level_order") or ())
    if not level_order or len(set(level_order)) != len(level_order):
        raise ValueError("ranking policy requires unique abstraction levels")
    status_order = tuple(condition.get("status_order") or ())
    if set(status_order) != set(_CONDITION_STATUSES):
        raise ValueError("ranking policy condition statuses are incomplete")
    raw_penalties = realism.get("band_penalties")
    if not isinstance(raw_penalties, list) or not raw_penalties:
        raise ValueError("ranking policy requires precursor-realism penalties")
    penalties = []
    for item in raw_penalties:
        if not isinstance(item, dict):
            raise ValueError("precursor-realism penalties must be objects")
        minimum_score = item.get("minimum_score")
        band_penalty = item.get("band_penalty")
        if (
            isinstance(minimum_score, bool)
            or not isinstance(minimum_score, (int, float))
            or not math.isfinite(float(minimum_score))
            or not 0.0 <= float(minimum_score) <= 1.0
        ):
            raise ValueError("realism minimum scores must be between zero and one")
        if (
            isinstance(band_penalty, bool)
            or not isinstance(band_penalty, int)
            or band_penalty < 0
        ):
            raise ValueError("realism band penalties must be nonnegative integers")
        penalties.append((float(minimum_score), band_penalty))
    if penalties != sorted(penalties, reverse=True):
        raise ValueError("realism minimum scores must be strictly descending")
    if len({minimum for minimum, _ in penalties}) != len(penalties):
        raise ValueError("realism minimum scores must be unique")
    if penalties[-1][0] != 0.0:
        raise ValueError("realism penalty policy must cover a zero score")
    if any(
        left_penalty > right_penalty
        for (_, left_penalty), (_, right_penalty) in zip(penalties, penalties[1:])
    ):
        raise ValueError("lower realism must not receive a smaller band penalty")
    return RetrosynthesisRankingPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        diversity_group_fields=group_fields,
        candidate_pool_multiplier=pool_multiplier,
        preserve_abstraction_level_order=preserve_levels,
        diversity_score_band_width=_positive_float(
            diversity.get("structural_score_band_width"),
            "candidate diversity structural score band",
        ),
        condition_score_band_width=_positive_float(
            condition.get("structural_score_band_width"),
            "condition reranking structural score band",
        ),
        abstraction_level_order=level_order,
        condition_status_order=status_order,
        precursor_realism_band_penalties=tuple(penalties),
    )


def diversity_group_key(
    candidate: GenericDisconnectionCandidate,
    policy: RetrosynthesisRankingPolicy,
) -> Tuple[str, ...]:
    """Return the chemistry-derived diversity identity for one candidate."""

    values = tuple(
        str(getattr(candidate, field, "") or "")
        for field in policy.diversity_group_fields
    )
    return values if any(values) else (candidate.precursor_smiles,)


def structural_score_bands(
    candidates: Iterable[GenericDisconnectionCandidate],
    *,
    width: float,
    separate_by_level: bool = True,
) -> dict[int, int]:
    """Assign deterministic score bands separately within every level."""

    values = tuple(candidates)
    best_by_level: dict[str, float] = {}
    for candidate in values:
        level = candidate.abstraction_level if separate_by_level else "*"
        best_by_level[level] = max(
            candidate.score,
            best_by_level.get(level, candidate.score),
        )
    bands = {}
    for candidate in values:
        level = candidate.abstraction_level if separate_by_level else "*"
        bands[id(candidate)] = int(
            math.floor(
                max(0.0, best_by_level[level] - candidate.score) / width
                + 1e-12
            )
        )
    return bands


__all__ = [
    "RANKING_POLICY_PATH",
    "RetrosynthesisRankingPolicy",
    "diversity_group_key",
    "load_retrosynthesis_ranking_policy",
    "structural_score_bands",
]
