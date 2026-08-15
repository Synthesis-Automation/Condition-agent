"""Validated deterministic ranking policy for forward product candidates."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Mapping, Tuple


_PATH = Path(__file__).with_name("definitions") / "forward_ranking.v1.json"


@dataclass(frozen=True)
class ForwardRankingPolicy:
    """Versioned weights and ordering guards for forward ranking."""

    definition_id: str
    schema_version: str
    calibration_status: str
    weights: Dict[str, float]
    abstraction_priority: Tuple[str, ...]
    score_band_width: float

    def level_rank(self, level: str) -> int:
        try:
            return self.abstraction_priority.index(level)
        except ValueError:
            return len(self.abstraction_priority)


@lru_cache(maxsize=1)
def load_forward_ranking_policy() -> ForwardRankingPolicy:
    """Load and strictly validate the forward ranking definition."""

    value = json.loads(_PATH.read_text(encoding="utf-8"))
    if value.get("definition_id") != "forward_ranking.v1":
        raise ValueError("unexpected forward ranking definition ID")
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported forward ranking definition schema")
    weights = {
        str(key): float(weight) for key, weight in (value.get("weights") or {}).items()
    }
    expected = {
        "precursor_similarity",
        "template_specificity",
        "independent_support",
        "operator_edit_agreement",
        "recipe_compatibility",
    }
    if set(weights) != expected or any(weight < 0.0 for weight in weights.values()):
        raise ValueError("forward ranking weights are incomplete or invalid")
    if abs(sum(weights.values()) - 1.0) > 1e-9:
        raise ValueError("forward ranking weights must sum to one")
    priority = tuple(str(item) for item in value.get("abstraction_priority") or ())
    if not priority or len(priority) != len(set(priority)):
        raise ValueError("forward abstraction priority must be unique and nonempty")
    width = float(value.get("score_band_width") or 0.0)
    if not 0.0 < width <= 1.0:
        raise ValueError("forward score-band width must be in (0, 1]")
    return ForwardRankingPolicy(
        definition_id="forward_ranking.v1",
        schema_version="1.0",
        calibration_status=str(value.get("calibration_status") or ""),
        weights=weights,
        abstraction_priority=priority,
        score_band_width=width,
    )


def weighted_forward_score(
    components: Mapping[str, float | None],
    policy: ForwardRankingPolicy,
) -> tuple[float, Dict[str, float]]:
    """Score available components, renormalizing genuinely missing evidence."""

    available = {
        key: float(value) for key, value in components.items() if value is not None
    }
    denominator = sum(policy.weights[key] for key in available)
    if denominator <= 0.0:
        return 0.0, {key: 0.0 for key in policy.weights}
    contributions = {
        key: (
            policy.weights[key] * available[key] / denominator
            if key in available
            else 0.0
        )
        for key in policy.weights
    }
    return (
        round(sum(contributions.values()), 8),
        {key: round(value, 8) for key, value in contributions.items()},
    )


__all__ = [
    "ForwardRankingPolicy",
    "load_forward_ranking_policy",
    "weighted_forward_score",
]
