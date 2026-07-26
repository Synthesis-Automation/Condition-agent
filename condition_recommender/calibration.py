"""Development-only ranking calibration with a separate validation gate."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from statistics import mean
from typing import Any, Dict, Mapping, Tuple

from .evaluation import evaluate_generic_index
from .generic_retrieval import load_generic_retrieval_rules
from .recipe_ranking import load_generic_ranking_rules

_RULES_PATH = Path(__file__).with_name("definitions") / "generic_calibration.v1.json"
_COMPONENTS = {
    "similarity",
    "yield",
    "independent_support",
    "reaction_breadth",
    "dataset_diversity",
    "compatibility",
    "condition_certainty",
}


def _validated_weights(weights: Any) -> Dict[str, float]:
    if not isinstance(weights, Mapping) or set(weights) != _COMPONENTS:
        raise ValueError("calibration weights do not match ranking components")
    values = {str(key): float(value) for key, value in weights.items()}
    if any(value < 0.0 for value in values.values()):
        raise ValueError("calibration weights must be non-negative")
    if abs(sum(values.values()) - 1.0) > 1e-9:
        raise ValueError("calibration weights must sum to one")
    return values


def validate_generic_calibration_rules(rules: Mapping[str, Any]) -> None:
    """Validate deterministic calibration candidates and promotion gates."""
    if str(rules.get("schema_version") or "") != "1.0":
        raise ValueError("unsupported generic calibration schema")
    if str(rules.get("definition_id") or "") != "generic_calibration.v1":
        raise ValueError("unexpected generic calibration definition ID")
    seeds = tuple(int(seed) for seed in rules.get("development_seeds") or ())
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError("development calibration seeds must be unique")
    if not 0.0 < float(rules["test_fraction"]) < 1.0:
        raise ValueError("calibration test_fraction must be in (0, 1)")
    if int(rules["top_k"]) < 1:
        raise ValueError("calibration top_k must be positive")
    if int(rules["minimum_seen_recipe_queries"]) < 1:
        raise ValueError("minimum seen-recipe query count must be positive")
    if str(rules.get("split_mode") or "") not in {
        "grouped_random",
        "scaffold_disjoint",
    }:
        raise ValueError("unsupported calibration split mode")
    candidates = rules.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise ValueError("calibration candidates must be a non-empty list")
    identifiers = [str(candidate.get("candidate_id") or "") for candidate in candidates]
    if any(not identifier for identifier in identifiers) or len(set(identifiers)) != len(
        identifiers
    ):
        raise ValueError("calibration candidate IDs must be present and unique")
    for candidate in candidates:
        if int(candidate["minimum_pool_size"]) < 1:
            raise ValueError("candidate minimum_pool_size must be positive")
        _validated_weights(candidate.get("weights"))


@lru_cache(maxsize=1)
def load_generic_calibration_rules() -> Dict[str, Any]:
    """Load the versioned calibration search and promotion policy."""
    rules = json.loads(_RULES_PATH.read_text(encoding="utf-8"))
    validate_generic_calibration_rules(rules)
    return dict(rules)


def _candidate_summary(
    reports: Tuple[Mapping[str, Any], ...],
) -> Dict[str, Any]:
    metrics = tuple(report["metrics"] for report in reports)
    yield_values = [
        float(value)
        for value in (metric["yield_mae_pct"] for metric in metrics)
        if value is not None
    ]
    return {
        "run_count": len(reports),
        "query_count": sum(int(metric["query_count"]) for metric in metrics),
        "seen_recipe_query_count": sum(
            int(metric["seen_recipe_query_count"]) for metric in metrics
        ),
        "coverage_rate_mean": round(
            mean(float(metric["coverage_rate"]) for metric in metrics), 6
        ),
        "seen_top1_rate_mean": round(
            mean(
                float(metric["seen_recipe_top1_recovery_rate"])
                for metric in metrics
            ),
            6,
        ),
        "seen_topk_rate_mean": round(
            mean(
                float(metric["seen_recipe_topk_recovery_rate"])
                for metric in metrics
            ),
            6,
        ),
        "top1_rate_mean": round(
            mean(float(metric["top1_recipe_recovery_rate"]) for metric in metrics),
            6,
        ),
        "yield_mae_pct_mean": (
            round(mean(yield_values), 6) if yield_values else None
        ),
        "hard_incompatible_recommendation_count": sum(
            int(metric["hard_incompatible_recommendation_count"])
            for metric in metrics
        ),
    }


def _selection_key(summary: Mapping[str, Any], candidate_id: str) -> tuple:
    yield_mae = summary["yield_mae_pct_mean"]
    return (
        -int(summary["hard_incompatible_recommendation_count"]),
        float(summary["seen_top1_rate_mean"]),
        float(summary["seen_topk_rate_mean"]),
        float(summary["coverage_rate_mean"]),
        -(float(yield_mae) if yield_mae is not None else 1_000.0),
        candidate_id,
    )


def calibrate_generic_ranking(
    development_index_path: str | Path,
    validation_index_path: str | Path,
    output_dir: str | Path,
) -> Dict[str, Any]:
    """Select on development splits and promote only after validation."""
    rules = load_generic_calibration_rules()
    destination = Path(output_dir)
    development_reports: Dict[str, Tuple[Mapping[str, Any], ...]] = {}
    summaries = {}
    for candidate in rules["candidates"]:
        candidate_id = str(candidate["candidate_id"])
        weights = _validated_weights(candidate["weights"])
        reports = tuple(
            evaluate_generic_index(
                development_index_path,
                destination / "development" / candidate_id / f"seed-{seed}",
                test_fraction=float(rules["test_fraction"]),
                seed=int(seed),
                top_k=int(rules["top_k"]),
                minimum_pool_size=int(candidate["minimum_pool_size"]),
                split_mode=str(rules["split_mode"]),
                retrieval_strategy="hybrid",
                ranking_weights=weights,
            )
            for seed in rules["development_seeds"]
        )
        development_reports[candidate_id] = reports
        summaries[candidate_id] = _candidate_summary(reports)
    eligible = [
        candidate
        for candidate in rules["candidates"]
        if summaries[str(candidate["candidate_id"])]["seen_recipe_query_count"]
        >= int(rules["minimum_seen_recipe_queries"])
        and summaries[str(candidate["candidate_id"])][
            "hard_incompatible_recommendation_count"
        ]
        == 0
    ]
    selected = (
        max(
            eligible,
            key=lambda candidate: _selection_key(
                summaries[str(candidate["candidate_id"])],
                str(candidate["candidate_id"]),
            ),
        )
        if eligible
        else None
    )
    validation_reports = {}
    promoted = False
    promotion_reasons = []
    if selected is not None:
        selected_id = str(selected["candidate_id"])
        current = next(
            candidate
            for candidate in rules["candidates"]
            if candidate["candidate_id"] == "current_min3"
        )
        for label, candidate in (("selected", selected), ("current", current)):
            validation_reports[label] = evaluate_generic_index(
                validation_index_path,
                destination / "validation" / label,
                test_fraction=float(rules["test_fraction"]),
                seed=int(rules["validation_seed"]),
                top_k=int(rules["top_k"]),
                minimum_pool_size=int(candidate["minimum_pool_size"]),
                split_mode=str(rules["split_mode"]),
                retrieval_strategy="hybrid",
                ranking_weights=_validated_weights(candidate["weights"]),
            )
        selected_metrics = validation_reports["selected"]["metrics"]
        current_metrics = validation_reports["current"]["metrics"]
        if int(selected_metrics["hard_incompatible_recommendation_count"]):
            promotion_reasons.append("validation_hard_incompatibility")
        if int(selected_metrics["seen_recipe_query_count"]) < int(
            rules["minimum_seen_recipe_queries"]
        ):
            promotion_reasons.append("insufficient_validation_seen_recipes")
        if float(selected_metrics["coverage_rate"]) < (
            float(current_metrics["coverage_rate"])
            - float(rules["maximum_coverage_drop"])
        ):
            promotion_reasons.append("validation_coverage_regression")
        if float(selected_metrics["seen_recipe_top1_recovery_rate"]) < (
            float(current_metrics["seen_recipe_top1_recovery_rate"])
            - float(rules["maximum_seen_top1_drop"])
        ):
            promotion_reasons.append("validation_seen_top1_regression")
        promoted = not promotion_reasons
    else:
        selected_id = None
        promotion_reasons.append("insufficient_development_seen_recipes")

    report: Dict[str, Any] = {
        "schema_version": "1.0",
        "artifact_type": "generic_ranking_calibration",
        "calibration_definition_id": rules["definition_id"],
        "development_index_path": str(Path(development_index_path)),
        "validation_index_path": str(Path(validation_index_path)),
        "development_summaries": summaries,
        "selected_candidate_id": selected_id,
        "validation_metrics": {
            label: value["metrics"] for label, value in validation_reports.items()
        },
        "promoted": promoted,
        "promotion_reasons": promotion_reasons,
    }
    if promoted and selected is not None:
        ranking = dict(load_generic_ranking_rules())
        ranking["calibration_status"] = "calibrated_sample_v1"
        ranking["weights"] = _validated_weights(selected["weights"])
        ranking["calibrated_minimum_pool_size"] = int(
            selected["minimum_pool_size"]
        )
        ranking["calibration_provenance"] = {
            "definition_id": rules["definition_id"],
            "selected_candidate_id": selected_id,
        }
        (destination / "recommended_ranking_definition.json").write_text(
            json.dumps(ranking, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        retrieval = dict(load_generic_retrieval_rules())
        retrieval["minimum_pool_size"] = int(selected["minimum_pool_size"])
        retrieval["calibration_status"] = "calibrated_sample_v1"
        retrieval["calibration_provenance"] = {
            "definition_id": rules["definition_id"],
            "selected_candidate_id": selected_id,
        }
        (destination / "recommended_retrieval_definition.json").write_text(
            json.dumps(retrieval, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "calibration_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


__all__ = [
    "calibrate_generic_ranking",
    "load_generic_calibration_rules",
    "validate_generic_calibration_rules",
]
