"""Validated declarative ranking policy for deterministic route search."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


MULTISTEP_RANKING_POLICY_PATH = (
    Path(__file__).with_name("definitions") / "multistep_ranking.v3.json"
)


@dataclass(frozen=True)
class MultistepRankingPolicy:
    """Immutable route-search breadth and step-cost configuration."""

    definition_id: str
    schema_version: str
    maximum_paths_per_state: int
    minimum_candidates_per_level: int
    base_step_cost: float
    candidate_score_deficit_weight: float
    abstraction_level_penalties: tuple[tuple[str, float], ...]
    selectivity_warning_penalty: float
    heuristic_terminal_penalty: float
    precursor_compatibility_band_penalty_weight: float
    precursor_realism_band_penalty_weight: float
    strategic_progress_minimum_score: float
    strategic_progress_deficit_penalty_weight: float
    complexity_increasing_step_penalty: float
    candidate_rank_tiebreak: float
    condition_status_penalties: tuple[tuple[str, float], ...]

    def abstraction_penalty(self, level: str) -> float:
        """Return the configured fallback penalty for one abstraction level."""

        return dict(self.abstraction_level_penalties).get(level, 0.0)

    def condition_status_penalty(self, status: str) -> float:
        """Return the configured route-cost penalty for condition support."""

        return dict(self.condition_status_penalties).get(
            status,
            dict(self.condition_status_penalties)["insufficient_evidence"],
        )


def _nonnegative_float(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must be a nonnegative number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized < 0.0:
        raise ValueError(f"{field} must be a nonnegative number")
    return normalized


def _positive_integer(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{field} must be a positive integer")
    return value


def _nonnegative_integer(value: object, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{field} must be a nonnegative integer")
    return value


@lru_cache(maxsize=1)
def load_multistep_ranking_policy() -> MultistepRankingPolicy:
    """Load and validate the versioned multistep ranking definition."""

    value = json.loads(MULTISTEP_RANKING_POLICY_PATH.read_text(encoding="utf-8"))
    if value.get("definition_id") != "multistep_ranking.v3":
        raise ValueError("unexpected multistep ranking definition ID")
    if value.get("schema_version") != "3.0":
        raise ValueError("unsupported multistep ranking schema")
    search = value.get("search")
    costs = value.get("step_cost")
    condition_postprocessing = value.get("condition_postprocessing")
    if (
        not isinstance(search, dict)
        or not isinstance(costs, dict)
        or not isinstance(condition_postprocessing, dict)
    ):
        raise ValueError(
            "multistep ranking requires search, step-cost, and condition rules"
        )
    raw_penalties = costs.get("abstraction_level_penalties")
    if not isinstance(raw_penalties, dict) or not raw_penalties:
        raise ValueError("multistep ranking requires abstraction penalties")
    penalties = tuple(
        sorted(
            (
                str(level),
                _nonnegative_float(penalty, f"{level} abstraction penalty"),
            )
            for level, penalty in raw_penalties.items()
        )
    )
    raw_condition_penalties = condition_postprocessing.get("status_penalties")
    required_condition_statuses = {
        "recommended_direct",
        "recommended_fallback",
        "insufficient_evidence",
    }
    if (
        not isinstance(raw_condition_penalties, dict)
        or set(raw_condition_penalties) != required_condition_statuses
    ):
        raise ValueError("multistep ranking requires all condition statuses")
    condition_penalties = tuple(
        sorted(
            (
                str(status),
                _nonnegative_float(
                    penalty,
                    f"{status} condition-status penalty",
                ),
            )
            for status, penalty in raw_condition_penalties.items()
        )
    )
    return MultistepRankingPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        maximum_paths_per_state=_positive_integer(
            search.get("maximum_paths_per_state"),
            "maximum paths per state",
        ),
        minimum_candidates_per_level=_nonnegative_integer(
            search.get("minimum_candidates_per_level"),
            "minimum candidates per level",
        ),
        base_step_cost=_nonnegative_float(costs.get("base"), "base step cost"),
        candidate_score_deficit_weight=_nonnegative_float(
            costs.get("candidate_score_deficit_weight"),
            "candidate score-deficit weight",
        ),
        abstraction_level_penalties=penalties,
        selectivity_warning_penalty=_nonnegative_float(
            costs.get("selectivity_warning_penalty"),
            "selectivity-warning penalty",
        ),
        heuristic_terminal_penalty=_nonnegative_float(
            costs.get("heuristic_terminal_penalty"),
            "heuristic-terminal penalty",
        ),
        precursor_compatibility_band_penalty_weight=_nonnegative_float(
            costs.get("precursor_compatibility_band_penalty_weight"),
            "precursor-compatibility band-penalty weight",
        ),
        precursor_realism_band_penalty_weight=_nonnegative_float(
            costs.get("precursor_realism_band_penalty_weight"),
            "precursor-realism band-penalty weight",
        ),
        strategic_progress_minimum_score=_nonnegative_float(
            costs.get("strategic_progress_minimum_score"),
            "strategic-progress minimum score",
        ),
        strategic_progress_deficit_penalty_weight=_nonnegative_float(
            costs.get("strategic_progress_deficit_penalty_weight"),
            "strategic-progress deficit penalty weight",
        ),
        complexity_increasing_step_penalty=_nonnegative_float(
            costs.get("complexity_increasing_step_penalty"),
            "complexity-increasing step penalty",
        ),
        candidate_rank_tiebreak=_nonnegative_float(
            costs.get("candidate_rank_tiebreak"),
            "candidate-rank tiebreak",
        ),
        condition_status_penalties=condition_penalties,
    )


__all__ = [
    "MULTISTEP_RANKING_POLICY_PATH",
    "MultistepRankingPolicy",
    "load_multistep_ranking_policy",
]
